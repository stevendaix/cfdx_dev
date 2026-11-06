// Reproducible linear-system replay, independent of the CFD flow solver.
// Usage: cfdx_linear_replay capture.cfdxlin [jacobi|ilu0|native|directional|energy|auto ...]
// All chosen candidates must pass their ORIGINAL physical residual gate.
#include "cfdx/core/linalg/linear_system_capture.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/anisotropy_analysis.h"
#include "cfdx/core/linalg/matrix_diagnostics.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;

namespace {
PreconditionerModel choose(const std::string& name) {
    if (name == "auto") return PreconditionerModel::Auto;
    if (name == "jacobi") return PreconditionerModel::Jacobi;
    if (name == "ilu0") return PreconditionerModel::ILU0;
    if (name == "native") return PreconditionerModel::NativeAMG;
    if (name == "directional") return PreconditionerModel::DirectionalSemiCoarseningAMG;
    if (name == "energy") return PreconditionerModel::EnergyMinimizingAMG;
    throw std::invalid_argument("unknown preconditioner: " + name);
}
}
int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            std::cerr << "Usage: cfdx_linear_replay file.cfdxlin "
                "[auto|jacobi|ilu0|native|directional|energy ...]\n";
            return 2;
        }
        const auto snapshot = read_linear_system_capture(argv[1]);
        const auto diagnosis = diagnose_matrix(snapshot.matrix);
        const auto anisotropy = analyze_anisotropy(snapshot.matrix);
        std::cerr << anisotropy.json(to_string(snapshot.problem)) << '\n';
        std::cerr << "CSR rows=" << diagnosis.rows << " nnz=" << diagnosis.nnz
                  << " missing_diagonal=" << diagnosis.missing_diagonal
                  << " near_zero_diagonal=" << diagnosis.near_zero_diagonal
                  << " checksum=" << std::hex << snapshot.checksum << std::dec << '\n';
        std::cout << "preconditioner,krylov,status,iterations,total_ms,true_abs_l2,true_relative\n";
        bool passed = true;
        for (int i = 2; i < (argc > 2 ? argc : 3); ++i) {
            const std::string name = argc > 2 ? argv[i] : "auto";
            auto request = LinearSolverRequest{};
            request.preconditioner = choose(name);
            request.krylov = KrylovModel::Auto;
            request.null_space = snapshot.null_space;
            request.scaling = snapshot.scaling;
            auto x = snapshot.initial_solution;
            const auto start = std::chrono::steady_clock::now();
            try {
                const auto result = solve_linear_system(snapshot.matrix, snapshot.rhs, x,
                    snapshot.problem, request, snapshot.max_iterations, snapshot.tolerance);
                const double elapsed = std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - start).count();
                const double absolute = physical_residual_norm(snapshot.matrix, x, snapshot.rhs);
                const double rhs = snapshot.rhs.norm2();
                const double initial = physical_residual_norm(snapshot.matrix,
                    snapshot.initial_solution, snapshot.rhs);
                const double reference = rhs > 0.0 ? rhs : initial;
                const double relative = reference > 0.0 ? absolute / reference :
                    (absolute == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
                const bool accepted = result.result.status == SolverStatus::CONVERGED &&
                    std::isfinite(relative) && relative <= snapshot.tolerance;
                passed = passed && accepted;
                if (!accepted && !result.preconditioner_error.empty())
                    std::cerr << "PRECONDITIONER_SETUP_ERROR " << name << ": "
                              << result.preconditioner_error << '\n';
                std::cout << name << ',' << to_string(result.plan.krylov) << ','
                    << to_string(result.result.status) << ',' << result.result.iterations << ','
                    << std::setprecision(17) << elapsed << ',' << absolute << ','
                    << relative << '\n';
            } catch (const std::exception& e) {
                passed = false;
                std::cerr << "REPLAY_SOLVE_FAILED " << name << ": " << e.what() << '\n';
                std::cout << name << ",none,exception,0,0,nan,nan\n";
            }
        }
        if (!passed) std::cerr << "CFDX_LINEAR_REPLAY_FAIL: at least one candidate failed\n";
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "CFDX_LINEAR_REPLAY_FAIL: " << e.what() << '\n';
        return 1;
    }
}
