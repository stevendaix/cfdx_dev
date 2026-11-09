#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/linear_system_capture.h"
#include "common/test_harness.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

struct Measurement {
    SolverResult solve;
    double true_relative = 0.0;
    double elapsed_seconds = 0.0;
};

Measurement measure(const SparseMatrix& matrix, const Vector& rhs,
                    const Vector& initial, bool periodic, int restart,
                    int period, std::size_t max_iter, double tolerance,
                    bool use_ilu0) {
    Vector x = initial;
    KrylovControls controls;
    controls.adaptive_restart = false;
    controls.restart_min = restart;
    controls.restart_max = restart;
    controls.residual_replacement = periodic;
    controls.residual_replacement_period = static_cast<std::size_t>(period);
    ILU0Preconditioner ilu;
    JacobiPreconditioner jacobi;
    Preconditioner* pc = use_ilu0 ? static_cast<Preconditioner*>(&ilu)
                                 : static_cast<Preconditioner*>(&jacobi);
    const auto before = std::chrono::steady_clock::now();
    const auto result = solve_fgmres(matrix, rhs, x, restart, max_iter,
                                     tolerance, pc, controls);
    const auto after = std::chrono::steady_clock::now();
    const auto ax = matrix.matvec(x);
    long double sq = 0.0L, rhs_sq = 0.0L;
    for (std::size_t i = 0; i < rhs.size(); ++i) {
        const long double ri = static_cast<long double>(rhs(i)) - ax[i];
        sq += ri * ri;
        rhs_sq += static_cast<long double>(rhs(i)) * rhs(i);
    }
    const double relative = rhs_sq > 0.0L
        ? static_cast<double>(std::sqrt(sq / rhs_sq))
        : static_cast<double>(std::sqrt(sq));
    return {result, relative,
            std::chrono::duration<double>(after - before).count()};
}

SparseMatrix sample_matrix(std::size_t n) {
    SparseMatrix matrix(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        matrix.push_back(i, i, 4.5);
        if (i > 0) matrix.push_back(i, i - 1, -1.2);
        if (i + 1 < n) matrix.push_back(i, i + 1, -0.7);
    }
    matrix.finalize();
    return matrix;
}

void print(const char* label, const Measurement& item) {
    std::cout << std::setprecision(17)
              << "KRYLOV_REPLAY mode=" << label
              << " status=" << to_string(item.solve.status)
              << " iterations=" << item.solve.iterations
              << " replacements=" << item.solve.residual_replacements
              << " reported_relative=" << item.solve.residual_relative
              << " independently_computed_relative=" << item.true_relative
              << " elapsed_seconds=" << item.elapsed_seconds << '\n';
}

} // namespace

int main(int argc, char** argv) {
    if (argc == 1) {
        run_case("replay_policy_preserves_nonsymmetric_solution", [] {
            const auto a = sample_matrix(36);
            Vector b(36, 0.0), initial(36, 0.0);
            for (std::size_t i = 0; i < b.size(); ++i)
                b(i) = std::sin(0.4 * static_cast<double>(i + 1));
            const auto usual = measure(a, b, initial, false, 12, 3, 400, 1e-11, true);
            const auto reliable = measure(a, b, initial, true, 12, 3, 400, 1e-11, true);
            EXPECT_TRUE(usual.solve.status == SolverStatus::CONVERGED);
            EXPECT_TRUE(reliable.solve.status == SolverStatus::CONVERGED);
            EXPECT_TRUE(reliable.solve.residual_replacements == 0);
            EXPECT_TRUE(usual.solve.residual_replacements == 0);
            EXPECT_TRUE(reliable.true_relative <= 1e-11);
            EXPECT_TRUE(usual.true_relative <= 1e-11);
            EXPECT_NEAR(reliable.true_relative, reliable.solve.residual_relative, 1e-12);
            print("standard", usual);
            print("periodic", reliable);
        });
        return run_all();
    }
    if (argc < 2 || argc > 7) {
        std::cerr << "usage: test_gmres_reliable_replay CAPTURE [restart] "
                     "[max_iters] [tolerance] [period] [ilu0|jacobi]\n";
        return 2;
    }
    try {
        const auto capture = read_linear_system_capture(argv[1]);
        const int restart = argc > 2 ? std::stoi(argv[2]) : 50;
        const auto max_iters = argc > 3 ? std::stoul(argv[3]) : 120UL;
        const double tolerance = argc > 4 ? std::stod(argv[4]) : capture.tolerance;
        const int period = argc > 5 ? std::stoi(argv[5]) : 20;
        const std::string pc = argc > 6 ? argv[6] : "ilu0";
        if (restart < 1 || restart > 512 || max_iters == 0 || period < 1 ||
            !std::isfinite(tolerance) || tolerance <= 0.0 ||
            (pc != "ilu0" && pc != "jacobi"))
            throw std::invalid_argument("invalid replay controls");
        std::cout << "KRYLOV_REPLAY_CAPTURE n=" << capture.matrix.n_rows()
                  << " nnz=" << capture.matrix.nnz()
                  << " checksum=" << capture.checksum << "\n";
        const auto usual = measure(capture.matrix, capture.rhs, capture.initial_solution,
                                   false, restart, period, max_iters, tolerance, pc == "ilu0");
        const auto reliable = measure(capture.matrix, capture.rhs, capture.initial_solution,
                                      true, restart, period, max_iters, tolerance, pc == "ilu0");
        print("standard", usual);
        print("periodic", reliable);
        // These are diagnostic experiments; a linear failure is NOT marked
        // physically qualified, and a successful replay is not a full CFD run.
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "KRYLOV_REPLAY_INPUT_ERROR " << e.what() << '\n';
        return 2;
    }
}
