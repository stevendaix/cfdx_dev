#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

SparseMatrix make_representative_coupled(std::size_t n) {
    SparseMatrix A(4 * n, 4 * n);
    for (std::size_t c = 0; c < n; ++c) {
        const double block[3][3] = {
            {4.0, 0.5, 0.0},
            {0.5, 5.0, 0.25},
            {0.0, 0.25, 3.0},
        };
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t q = 0; q < 3; ++q)
                if (block[r][q] != 0.0)
                    A.push_back(r * n + c, q * n + c, block[r][q]);

        A.push_back(3 * n + c, 3 * n + c, 1.0);
        A.push_back(3 * n + c, c, 1.0);
        A.push_back(c, 3 * n + c, 1.0);

        if (c + 1 < n) {
            A.push_back(c, 3 * n + c + 1, -1.0);
            A.push_back(3 * n + c + 1, c, -1.0);
        }
    }
    A.finalize();
    return A;
}

Vector matvec(const SparseMatrix& A, const Vector& x) {
    const auto values = A.matvec(x);
    Vector y(values.size());
    for (std::size_t i = 0; i < values.size(); ++i)
        y(i) = values[i];
    return y;
}

double true_relative_residual(const SparseMatrix& A,
                              const Vector& x,
                              const Vector& b) {
    const auto ax = A.matvec(x);
    double rr2 = 0.0;
    double bb2 = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) {
        const double r = b(i) - ax[i];
        rr2 += r * r;
        bb2 += b(i) * b(i);
    }
    return std::sqrt(rr2 / std::max(bb2, 1e-300));
}

std::size_t estimated_csr_bytes(const SparseMatrix& A) {
    // Diagnostic storage estimate only: CSR row offsets + column indices +
    // values. This is deliberately not reported as process RSS.
    return (A.n_rows() + 1) * sizeof(std::size_t) +
           A.nnz() * (sizeof(std::size_t) + sizeof(double));
}

const char* factorization_name(CoupledSchurFactorization mode) {
    switch (mode) {
        case CoupledSchurFactorization::Diagonal: return "diagonal";
        case CoupledSchurFactorization::Lower: return "lower";
        case CoupledSchurFactorization::Upper: return "upper";
        case CoupledSchurFactorization::Full: return "full";
    }
    return "unknown";
}

} // namespace

int main() {
    run_case("n8_schur_production_benchmark", [] {
        // This is a representative assembled 4N CFD block matrix. It is a
        // reproducible algebraic production-path benchmark, not a physical
        // Couette/Poiseuille/Ghia qualification oracle.
        for (const std::size_t n : {16u, 32u, 64u, 128u}) {
            const auto A = make_representative_coupled(n);

            Vector exact(4 * n);
            for (std::size_t i = 0; i < exact.size(); ++i)
                exact(i) = std::sin(0.017 * static_cast<double>(i + 1));
            const auto b = matvec(A, exact);
            Vector x(4 * n, 0.0);

            CoupledBlockSchurOptions options;
            options.factorization = CoupledSchurFactorization::Full;
            options.velocity_approximation =
                CoupledSchurVelocityApproximation::Block;
            CoupledBlockSchurAMGPreconditioner pc(n, options);

            const auto setup_begin = std::chrono::steady_clock::now();
            EXPECT_TRUE(pc.setup(A));
            const auto setup_end = std::chrono::steady_clock::now();

            const auto solve_begin = std::chrono::steady_clock::now();
            const auto result = solve_gmres(
                A, b, x, 20, 400, 1e-10, &pc, {}, false);
            const auto solve_end = std::chrono::steady_clock::now();

            const double rr = true_relative_residual(A, x, b);
            const auto setup_us = std::chrono::duration_cast<
                std::chrono::microseconds>(setup_end - setup_begin).count();
            const auto solve_us = std::chrono::duration_cast<
                std::chrono::microseconds>(solve_end - solve_begin).count();

            EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
            EXPECT_TRUE(std::isfinite(rr));
            EXPECT_TRUE(rr < 1e-9);
            EXPECT_TRUE(pc.is_ready());
            EXPECT_TRUE(pc.pressure_hierarchy_builds() == 1);
            EXPECT_TRUE(pc.pressure_numeric_updates() == 0);
            EXPECT_TRUE(pc.schur_matrix().nnz() > 0);
            EXPECT_TRUE(pc.pressure_coarse_size() > 0);

            std::cout
                << "n8_schur_benchmark"
                << " cells=" << n
                << " unknowns=" << A.n_rows()
                << " nnz=" << A.nnz()
                << " schur_nnz=" << pc.schur_matrix().nnz()
                << " csr_bytes_est=" << estimated_csr_bytes(A)
                << " schur_csr_bytes_est="
                << estimated_csr_bytes(pc.schur_matrix())
                << " factorization=" << factorization_name(pc.factorization())
                << " velocity_approximation=block"
                << " krylov=FGMRES"
                << " status=" << static_cast<int>(result.status)
                << " iterations=" << result.iterations
                << " reported_residual=" << result.residual
                << " true_residual=" << rr
                << " setup_us=" << setup_us
                << " solve_us=" << solve_us
                << " pressure_coarse_size=" << pc.pressure_coarse_size()
                << " hierarchy_builds=" << pc.pressure_hierarchy_builds()
                << " numeric_updates=" << pc.pressure_numeric_updates()
                << '\n';
        }

        // Numeric refresh must preserve the symbolic pressure hierarchy. This
        // is part of the benchmark lifecycle contract rather than a performance
        // claim.
        const std::size_t n = 64;
        const auto A = make_representative_coupled(n);
        auto A_updated = make_representative_coupled(n);

        // Preserve the exact CSR graph while changing a velocity-block
        // coefficient. Direct value access is intentional here: the test
        // exercises the values-only refresh contract, not matrix assembly.
        const std::size_t target_row = 0;
        const std::size_t target_col = 0;
        bool changed = false;
        const auto* columns = A_updated.columns_data();
        auto* values = A_updated.values_data();
        for (std::size_t k = A_updated.row_offsets_data()[target_row];
             k < A_updated.row_offsets_data()[target_row + 1]; ++k) {
            if (columns[k] == target_col) {
                values[k] *= 1.10;
                changed = true;
                break;
            }
        }
        EXPECT_TRUE(changed);
        EXPECT_TRUE(A_updated(target_row, target_col) !=
                    A(target_row, target_col));

        CoupledBlockSchurAMGPreconditioner pc(n);
        EXPECT_TRUE(pc.setup(A));
        const auto hierarchy_builds_before = pc.pressure_hierarchy_builds();
        EXPECT_TRUE(pc.update_values(A_updated));
        EXPECT_TRUE(pc.pressure_hierarchy_builds() == hierarchy_builds_before);
        EXPECT_TRUE(pc.pressure_numeric_updates() == 1);

        // The refreshed preconditioner must remain usable with the updated
        // operator; convergence is checked with an independently recomputed
        // residual, so this is not merely a counter assertion.
        Vector exact_updated(4 * n);
        for (std::size_t i = 0; i < exact_updated.size(); ++i)
            exact_updated(i) =
                std::sin(0.031 * static_cast<double>(i + 1));
        const auto b_updated = matvec(A_updated, exact_updated);
        Vector x_updated(4 * n, 0.0);
        const auto updated_result = solve_gmres(
            A_updated, b_updated, x_updated, 20, 400, 1e-10, &pc);
        const double updated_rr =
            true_relative_residual(A_updated, x_updated, b_updated);
        EXPECT_TRUE(updated_result.status == SolverStatus::CONVERGED);
        EXPECT_TRUE(std::isfinite(updated_rr));
        EXPECT_TRUE(updated_rr < 1e-9);

        std::cout << "n8_schur_benchmark_lifecycle"
                  << " cells=" << n
                  << " coefficient_changed=true"
                  << " hierarchy_builds_before=" << hierarchy_builds_before
                  << " hierarchy_builds_after="
                  << pc.pressure_hierarchy_builds()
                  << " numeric_updates=" << pc.pressure_numeric_updates()
                  << " updated_iterations=" << updated_result.iterations
                  << " updated_true_residual=" << updated_rr
                  << " graph_change_rebuild=explicit_setup_required"
                  << '\n';
    });

    return run_all();
}
