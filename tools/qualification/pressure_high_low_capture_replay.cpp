// Offline replay of the EXACT original pressure CSR/RHS with CPU-native
// high/low iterative refinement. Never marks a CFD result physically qualified.
#include "cfdx/core/linalg/linear_system_capture.h"
#include "cfdx/core/linalg/pressure_compensated_refinement.h"
#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>

using namespace cfdx::core;
int main(int argc, char** argv) {
    if (argc < 2 || argc > 4) {
        std::cerr << "usage: pressure_high_low_capture_replay BFS_CAPTURE.cfdxlin [amg|amg_cg|ilu0|none] [OUTPUT.split]\n";
        return 2;
    }
    try {
        const auto capture = read_linear_system_capture(argv[1]);
        if (capture.problem != LinearProblemKind::PressurePoisson &&
            capture.problem != LinearProblemKind::General)
            throw std::invalid_argument("capture is not a pressure/elliptic problem");
        const std::string method = argc >= 3 ? argv[2] : "amg";
        if (method != "amg" && method != "amg_cg" && method != "ilu0" && method != "none")
            throw std::invalid_argument("unsupported preconditioner");
        ILU0Preconditioner ilu;
        NativeAMGPreconditioner amg(NativeAMGMethod::BoomerStyle);
        Preconditioner* pc = (method == "amg" || method == "amg_cg") ? static_cast<Preconditioner*>(&amg)
            : method == "ilu0" ? static_cast<Preconditioner*>(&ilu) : nullptr;
        if (pc && !pc->setup(capture.matrix))
            throw std::runtime_error("preconditioner setup failed");
        Vector high = capture.initial_solution;
        Vector low(high.size(), 0.0);
        KrylovControls policy;
        policy.adaptive_restart = false;
        policy.restart_min = 160;
        policy.restart_max = 160;
        const auto begin = std::chrono::steady_clock::now();
        const auto initial = method == "amg_cg"
            ? solve_cg_prepared(capture.matrix, capture.rhs, high, *pc, 350, capture.tolerance)
            : solve_fgmres(capture.matrix, capture.rhs, high,
                           160, 350, capture.tolerance, pc, policy, {}, false);
        const auto after_first = std::chrono::steady_clock::now();
        // The original solver's return code is not accepted as an assurance:
        // independently recompute against the captured unscaled matrix.
        const auto report = refine_pressure_high_low(
            capture.matrix, capture.rhs, high, low, capture.tolerance, 4,
            [&](const SparseMatrix& a, const Vector& residual, Vector& delta) {
                const auto next = method == "amg_cg"
                    ? solve_cg_prepared(a, residual, delta, *pc, 350, 1e-7)
                    : solve_fgmres(a, residual, delta,
                                   160, 350, 1e-7, pc, policy, {}, false);
                // Even a finite best effort must REDUCE the *original*
                // high/low true residual. Never trust nested projected values.
                return next.status != SolverStatus::NOT_APPLICABLE &&
                       next.status != SolverStatus::DIVERGED;
            });
        if (argc == 4) {
            std::ofstream out(argv[3]);
            if (!out) throw std::runtime_error("cannot open pressure split output");
            out << "CFDX_PRESSURE_SPLIT_V1\n"
                << high.size() << " " << std::hex << std::setfill('0')
                << std::setw(16) << capture.checksum << std::dec << "\n";
            for (std::size_t i = 0; i < high.size(); ++i)
                out << std::hexfloat << high(i) << " " << low(i) << "\n";
            if (!out) throw std::runtime_error("failed to save pressure split");
        }
        const auto done = std::chrono::steady_clock::now();
        std::cout << std::scientific << std::setprecision(17)
            << "PRESSURE_HILO_CAPTURE preconditioner=" << method
            << " rows=" << capture.matrix.n_rows()
            << " nnz=" << capture.matrix.nnz()
            << " checksum=" << capture.checksum
            << " requested_tolerance=" << capture.tolerance
            << " initial_solver_status=" << static_cast<int>(initial.status)
            << " initial_solver_iters=" << initial.iterations
            << " original_double_solver_relative=" << initial.residual_relative
            << " after_initial_high_low_relative=" << report.initial_relative_residual
            << " final_high_low_relative=" << report.final_relative_residual
            << " final_single_double_relative=" << report.final_double_only_residual
            << " passes_accepted=" << report.passes_accepted
            << " first_solver_seconds=" << std::chrono::duration<double>(after_first-begin).count()
            << " total_seconds=" << std::chrono::duration<double>(done-begin).count()
            << " algebraic_converged=" << report.algebraic_converged
            << " physical_fv_qualified=" << report.physical_face_flux_certified << '\n';
        // A pressure *capture* is a diagnostic even if algebraic_converged.
        // Never return 'physical success' to a CFD qualification orchestrator.
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "PRESSURE_HILO_CAPTURE_ERROR " << e.what() << '\n';
        return 2;
    }
}
