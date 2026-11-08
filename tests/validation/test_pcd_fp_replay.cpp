// Replay the EXACT nested PCD Fp operator/RHS/initial state from production.
// Compare one and two 100-iteration FGMRES passes on the same original matrix;
// no synthetic residual replaces the independent FP64/long-double CSR residual.
#include "cfdx/physics/pcd_fp_capture.h"
#include "cfdx/physics/pcd_fp_equilibration.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/amg_preconditioner.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using Clock = std::chrono::steady_clock;
static double elapsed(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct BenchmarkResult {
    std::string method;
    std::size_t rows = 0, nnz = 0;
    bool setup_ok = false;
    double setup_seconds = 0, first_seconds = 0, second_seconds = 0;
    SolverResult first{}, second{};
    PcdFpTrueResidual first_residual{}, second_residual{};
    bool preconditioner_probe_ok = false;
    double probe_correction_l2 = 0, probe_A_correction_l2 = 0, probe_cosine = 0;
    bool second_improved = false;
    bool second_attempted = false;
    std::string error;
};

static BenchmarkResult compare_one(const PcdFpCapture& input,
    const std::string& method, std::size_t budget, double tolerance) {
    BenchmarkResult report;
    report.method = method;
    report.rows = input.matrix.n_rows();
    report.nnz = input.matrix.nnz();
    const bool equil = method == "ilu0_equilibrated" || method == "amg_equilibrated";
    const bool ilu = method == "ilu0" || method == "ilu0_equilibrated";
    const bool amg = method == "amg" || method == "amg_equilibrated";
    const bool sa = method == "amg_sa";
    const bool lines = method == "amg_lines";
    const bool semicoarsen = method == "amg_semicoarsen";
    if (!ilu && !amg && !sa && !lines && !semicoarsen && method != "jacobi")
        throw std::invalid_argument("unsupported PCD Fp replay method");
    try {
        // Match production's prepared ILU0 / AMG exactly; no CG smoother gate.
        const auto original = &input.matrix;
        std::unique_ptr<PcdFpEquilibration> scaled;
        if (equil)
            scaled = std::make_unique<PcdFpEquilibration>(equilibrate_pcd_fp_rows(*original));
        const auto& A = scaled ? scaled->matrix : *original;
        Vector b = scaled ? pcd_fp_scaled_rhs(input.rhs, scaled->row_multipliers) : input.rhs;
        const auto preconditioner_model = ilu ? PreconditionerModel::ILU0 :
            (amg ? PreconditionerModel::NativeAMG :
             (sa ? PreconditionerModel::SmoothedAggregationAMG :
             (lines ? PreconditionerModel::AnisotropicBoomerAMG :
             (semicoarsen ? PreconditionerModel::DirectionalSemiCoarseningAMG
                          : PreconditionerModel::Jacobi))));
        auto preconditioner = make_scalar_preconditioner(
            preconditioner_model, false, false);
        const auto setup_begin = Clock::now();
        report.setup_ok = preconditioner && preconditioner->setup(A);
        report.setup_seconds = elapsed(setup_begin);
        if (!report.setup_ok) {
            report.error = "PRECONDITIONER_SETUP_FAILED";
            return report;
        }
        Vector x = input.initial;
        const auto initial = input.initial;
        // The same physical initial residual is probed through each prepared
        // preconditioner. A grossly amplified correction exposes a broken
        // V-cycle even before checking Krylov's status.
        Vector probe_residual(A.n_rows(), 0.0), probe(A.n_rows(), 0.0);
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            double ax0 = 0;
            const auto* ro = A.row_offsets_data();
            const auto* co = A.columns_data();
            const auto* va = A.values_data();
            for (std::size_t k = ro[i]; k < ro[i + 1]; ++k)
                ax0 += va[k] * initial(co[k]);
            probe_residual(i) = b(i) - ax0;
        }
        report.preconditioner_probe_ok = preconditioner->apply(probe_residual, probe);
        if (report.preconditioner_probe_ok && probe.is_valid()) {
            report.probe_correction_l2 = probe.norm2();
            Vector Ap(A.n_rows(), 0.0);
            const auto* ro = A.row_offsets_data();
            const auto* co = A.columns_data();
            const auto* va = A.values_data();
            double dot = 0;
            for (std::size_t i = 0; i < A.n_rows(); ++i) {
                for (std::size_t k = ro[i]; k < ro[i+1]; ++k)
                    Ap(i) += va[k] * probe(co[k]);
                dot += probe_residual(i) * Ap(i);
            }
            report.probe_A_correction_l2 = Ap.norm2();
            const double denom = probe_residual.norm2() * Ap.norm2();
            report.probe_cosine = denom > 0 ? dot / denom : 0.0;
        }

        const auto run = [&](SolverResult& status, double& seconds,
                             PcdFpTrueResidual& residual) {
            const auto start = Clock::now();
            status = solve_fgmres(A, b, x,
                std::min<int>(512, static_cast<int>(A.n_rows())),
                budget, tolerance, preconditioner.get(), {}, {}, false);
            seconds = elapsed(start);
            residual = pcd_fp_true_residual(*original, input.rhs, initial, x);
        };
        run(report.first, report.first_seconds, report.first_residual);
        // Same criterion and extra 100-step Krylov cycle as the production opt-in.
        // Always run the second cycle, but accept it only if true residual improves.
        // Do NOT allow a worse iterate to masquerade as better convergence.
        if (x.is_valid() && report.first_residual.final_l2 >
                tolerance * std::max(1.0, input.rhs.norm2())) {
            report.second_attempted = true;
            const auto best = x;
            run(report.second, report.second_seconds, report.second_residual);
            report.second_improved = x.is_valid() &&
                std::isfinite(report.second_residual.final_l2) &&
                report.second_residual.final_l2 < report.first_residual.final_l2;
            if (!report.second_improved)
                report.second_residual = report.first_residual;
        } else {
            report.second = report.first;
            report.second_residual = report.first_residual;
        }
    } catch (const std::exception& exc) {
        report.error = exc.what();
    }
    return report;
}

static void print_json(const BenchmarkResult& r, const std::string& capture) {
    auto safe = [](double x) {return std::isfinite(x) ? x : -1.0;};
    std::cout << std::setprecision(17)
        << "{\"capture\":\"" << capture << "\",\"method\":\"" << r.method
        << "\",\"rows\":" << r.rows << ",\"nnz\":" << r.nnz
        << ",\"setup_ok\":" << (r.setup_ok ? "true" : "false")
        << ",\"setup_seconds\":" << safe(r.setup_seconds)
        << ",\"first_seconds\":" << safe(r.first_seconds)
        << ",\"second_seconds\":" << safe(r.second_seconds)
        << ",\"preconditioner_probe_ok\":" << (r.preconditioner_probe_ok ? "true" : "false")
        << ",\"probe_correction_l2\":" << safe(r.probe_correction_l2)
        << ",\"probe_A_correction_l2\":" << safe(r.probe_A_correction_l2)
        << ",\"probe_cosine\":" << safe(r.probe_cosine)
        << ",\"first_status\":" << static_cast<int>(r.first.status)
        << ",\"first_iterations\":" << r.first.iterations
        << ",\"first_true_relative_rhs\":" << safe(r.first_residual.relative_rhs)
        << ",\"first_true_relative_initial\":" << safe(r.first_residual.relative_initial)
        << ",\"first_true_l2\":" << safe(r.first_residual.final_l2)
        << ",\"second_attempted\":" << (r.second_attempted ? "true" : "false")
        << ",\"second_improved\":" << (r.second_improved ? "true" : "false")
        << ",\"second_status\":" << static_cast<int>(r.second.status)
        << ",\"second_iterations\":" << r.second.iterations
        << ",\"best_true_relative_rhs\":" << safe(r.second_residual.relative_rhs)
        << ",\"best_true_l2\":" << safe(r.second_residual.final_l2)
        << ",\"total_seconds_two_passes\":"
        << safe(r.setup_seconds + r.first_seconds + r.second_seconds)
        << ",\"error\":\"";
    // Errors should be printed as one JSON-safe string.
    for (char c : r.error) {
        if (c == '"' || c == '\\') std::cout << '\\';
        if (c >= ' ' && c != '\n' && c != '\r') std::cout << c;
    }
    std::cout << "\"}\n";
}

static int self_test() {
    SparseMatrix A(3, 3);
    A.push_back(0, 0, 4.0); A.push_back(0, 1, -1.0);
    A.push_back(1, 0, -1.0); A.push_back(1, 1, 4.0);
    A.push_back(1, 2, -1.0); A.push_back(2, 1, -1.0);
    A.push_back(2, 2, 4.0); A.finalize();
    Vector b(3), initial(3, 0.0);
    b(0) = 1; b(1) = 2; b(2) = 3;
    const auto file = std::filesystem::temp_directory_path() /
        "cfdx_pcd_fp_self_test.cfp";
    write_pcd_fp_capture(file, A, b, initial);
    auto capture = read_pcd_fp_capture(file);
    std::filesystem::remove(file);
    if (capture.matrix.nnz() != A.nnz() || capture.rhs(1) != 2)
        return 2;
    for (auto method : {"ilu0", "ilu0_equilibrated", "amg"}) {
        auto result = compare_one(capture, method, 100, 1e-10);
        if (!result.setup_ok || result.first_residual.relative_rhs > 1e-8)
            return 3;
    }
    std::cout << "PCD_FP_CAPTURE_SELF_TEST_PASS\n";
    return 0;
}

int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--self-test")
            return self_test();
        if (argc < 2 || argc > 4)
            throw std::invalid_argument(
                "usage: test_pcd_fp_replay capture.cfp [budget=100] [tolerance=1e-8]");
        const std::string file = argv[1];
        const auto budget = argc > 2 ? static_cast<std::size_t>(std::stoul(argv[2])) : 100u;
        const double tolerance = argc > 3 ? std::stod(argv[3]) : 1e-8;
        if (budget == 0 || budget > 5000 || !(tolerance >= 1e-14 && tolerance <= 1e-4))
            throw std::invalid_argument("invalid nested benchmark solver budget/tolerance");
        auto capture = read_pcd_fp_capture(file);
        for (auto method : {"ilu0", "ilu0_equilibrated", "amg", "amg_equilibrated", "amg_sa", "amg_lines", "amg_semicoarsen", "jacobi"})
            print_json(compare_one(capture, method, budget, tolerance), file);
        return 0; // a negative solver result is data, never a fake qualification PASS
    } catch (const std::exception& e) {
        std::cerr << "PCD_FP_REPLAY_ERROR: " << e.what() << '\n';
        return 2;
    }
}
