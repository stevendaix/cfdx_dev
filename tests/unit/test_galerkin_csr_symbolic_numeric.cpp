#include "cfdx/core/linalg/galerkin_csr_product.h"
#include "cfdx/core/linalg/amg_preconditioner.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/projected_gmres.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
using Rows = GalerkinCSRProduct::TransferRows;
SparseMatrix fine(bool changed_pattern = false, double scale = 1.0) {
    SparseMatrix A(5, 5);
    for (std::size_t i = 0; i < 5; ++i) {
        A.push_back(i, i, (3.4 + 0.2 * static_cast<double>(i)) * scale);
        if (i > 0) A.push_back(i, i - 1, -0.7 * scale);
        if (i + 1 < 5) A.push_back(i, i + 1, -0.6 * scale);
    }
    if (changed_pattern) A.push_back(0, 4, -0.125 * scale);
    A.finalize();
    return A;
}
Rows transfers() {
    return {{{0, 1.0}}, {{0, 0.5}, {1, 0.5}}, {{1, 1.0}},
            {{1, 0.2}, {2, 0.8}}, {{2, 1.0}}};
}
SparseMatrix legacy_reference(const SparseMatrix& A, const Rows& P,
                              std::size_t nc) {
    std::vector<std::map<std::size_t, double>> rows(nc);
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i];
             k < A.row_offsets_data()[i + 1]; ++k) {
            const auto j = A.columns_data()[k];
            for (const auto& [ci, wi] : P[i])
                for (const auto& [cj, wj] : P[j])
                    rows[ci][cj] += wi * A.values_data()[k] * wj;
        }
    SparseMatrix C(nc, nc);
    for (std::size_t i = 0; i < nc; ++i)
        for (const auto& [j, value] : rows[i])
            if (std::isfinite(value) && std::abs(value) > 1e-30)
                C.push_back(i, j, value);
    C.finalize();
    return C;
}
void compare(const SparseMatrix& lhs, const SparseMatrix& rhs) {
    EXPECT_TRUE(lhs.n_rows() == rhs.n_rows());
    EXPECT_TRUE(lhs.n_cols() == rhs.n_cols());
    for (std::size_t i = 0; i < lhs.n_rows(); ++i)
        for (std::size_t j = 0; j < lhs.n_cols(); ++j) {
            const double expected = rhs(i, j);
            EXPECT_TRUE(std::abs(lhs(i, j) - expected) <=
                        2e-13 * std::max(1.0, std::abs(expected)));
        }
}
SparseMatrix rotated_anisotropic(std::size_t nside, double main, double cross) {
    // Nine-point symmetric M-matrix: high diffusion along the 45-degree
    // diagonals and weak coupling in the perpendicular direction.
    SparseMatrix A(nside * nside, nside * nside);
    for (std::size_t y = 0; y < nside; ++y)
        for (std::size_t x = 0; x < nside; ++x) {
            const std::size_t i = y * nside + x;
            double diagonal = 1.0;
            for (const int dy : {-1, 0, 1})
                for (const int dx : {-1, 0, 1}) {
                    if (dx == 0 && dy == 0) continue;
                    const auto nx = static_cast<long>(x) + dx;
                    const auto ny = static_cast<long>(y) + dy;
                    if (nx < 0 || ny < 0 || nx >= static_cast<long>(nside) ||
                        ny >= static_cast<long>(nside)) continue;
                    const double weight = (dx == dy && dx != 0) ? main : cross;
                    const std::size_t j = static_cast<std::size_t>(ny) * nside +
                                          static_cast<std::size_t>(nx);
                    A.push_back(i, j, -weight);
                    diagonal += weight;
                }
            A.push_back(i, i, diagonal);
        }
    A.finalize();
    return A;
}
SparseMatrix anisotropic(std::size_t nside, double ax, double ay) {
    SparseMatrix A(nside * nside, nside * nside);
    for (std::size_t y = 0; y < nside; ++y)
        for (std::size_t x = 0; x < nside; ++x) {
            const auto i = y * nside + x;
            double diagonal = 1.0;
            if (x > 0) { A.push_back(i, i - 1, -ax); diagonal += ax; }
            if (x + 1 < nside) { A.push_back(i, i + 1, -ax); diagonal += ax; }
            if (y > 0) { A.push_back(i, i - nside, -ay); diagonal += ay; }
            if (y + 1 < nside) { A.push_back(i, i + nside, -ay); diagonal += ay; }
            A.push_back(i, i, diagonal);
        }
    A.finalize();
    return A;
}
}

int main() {
    run_case("csr_factory_rejects_unsorted_and_duplicate_columns", [] {
        bool invalid = false;
        try {
            (void)SparseMatrix::from_sorted_csr(2, 2, {0, 2, 2}, {1, 0}, {1, 2});
        } catch (const std::invalid_argument&) { invalid = true; }
        EXPECT_TRUE(invalid);
        invalid = false;
        try {
            (void)SparseMatrix::from_sorted_csr(2, 2, {0, 2, 2}, {0, 0}, {1, 2});
        } catch (const std::invalid_argument&) { invalid = true; }
        EXPECT_TRUE(invalid);
        invalid = false;
        try {
            (void)SparseMatrix::from_sorted_csr(1, 1, {0, 1}, {0},
                {std::numeric_limits<double>::infinity()});
        } catch (const std::invalid_argument&) { invalid = true; }
        EXPECT_TRUE(invalid);
    });
    run_case("symbolic_numeric_matches_legacy_ptap_and_reuses_pattern", [] {
        const auto P = transfers();
        auto A = fine();
        GalerkinCSRProduct plan;
        EXPECT_TRUE(plan.symbolic(A, P, 3));
        EXPECT_TRUE(plan.symbolic_ap_nnz() > 0);
        EXPECT_TRUE(plan.symbolic_coarse_nnz() > 0);
        SparseMatrix C;
        EXPECT_TRUE(plan.numeric(A, P, C));
        compare(C, legacy_reference(A, P, 3));
        auto changed = fine(false, 2.5);
        SparseMatrix refreshed;
        EXPECT_TRUE(plan.numeric(changed, P, refreshed));
        compare(refreshed, legacy_reference(changed, P, 3));
        EXPECT_TRUE(refreshed.nnz() == C.nnz());
    });
    run_case("numeric_values_only_refresh_preserves_csr_pattern", [] {
        const auto P = transfers();
        const auto original = fine();
        const auto updated = fine(false, 2.5);
        GalerkinCSRProduct plan;
        EXPECT_TRUE(plan.symbolic(original, P, 3));
        SparseMatrix coarse;
        EXPECT_TRUE(plan.numeric(original, P, coarse));
        const SparseMatrix old_coarse = coarse;
        EXPECT_TRUE(plan.numeric_into_existing(updated, P, coarse));
        compare(coarse, legacy_reference(updated, P, 3));
        EXPECT_TRUE(coarse.nnz() == old_coarse.nnz());
        EXPECT_TRUE(std::equal(coarse.row_offsets_data(),
                               coarse.row_offsets_data() + coarse.n_rows() + 1,
                               old_coarse.row_offsets_data()));
        EXPECT_TRUE(std::equal(coarse.columns_data(),
                               coarse.columns_data() + coarse.nnz(),
                               old_coarse.columns_data()));
    });
    run_case("numeric_values_only_failure_is_transactional", [] {
        const Rows identity{{{0, 1.0}}, {{1, 1.0}}};
        const auto original = SparseMatrix::from_sorted_csr(
            2, 2, {0, 2, 4}, {0, 1, 0, 1}, {2.0, -1.0, -1.0, 2.0});
        const auto pruned = SparseMatrix::from_sorted_csr(
            2, 2, {0, 2, 4}, {0, 1, 0, 1}, {2.0, 0.0, 0.0, 2.0});
        GalerkinCSRProduct plan;
        EXPECT_TRUE(plan.symbolic(original, identity, 2));
        SparseMatrix coarse;
        EXPECT_TRUE(plan.numeric(original, identity, coarse));
        const SparseMatrix retained = coarse;
        EXPECT_TRUE(!plan.numeric_into_existing(pruned, identity, coarse));
        compare(coarse, retained);
        EXPECT_TRUE(coarse.nnz() == retained.nnz());
        EXPECT_TRUE(!plan.numeric_into_existing(coarse, identity, coarse));
        compare(coarse, retained);
        const Rows reversed{{{1, 1.0}}, {{0, 1.0}}};
        EXPECT_TRUE(!plan.numeric_into_existing(original, reversed, coarse));
        compare(coarse, retained);
    });
    run_case("failed_reuse_does_not_mutate_existing_coarse_matrix", [] {
        const auto P = transfers();
        auto A = fine();
        GalerkinCSRProduct plan;
        EXPECT_TRUE(plan.symbolic(A, P, 3));
        SparseMatrix C;
        EXPECT_TRUE(plan.numeric(A, P, C));
        const SparseMatrix original = C;
        EXPECT_TRUE(!plan.numeric(fine(true), P, C));
        compare(C, original);
        auto wrong = P;
        std::reverse(wrong[1].begin(), wrong[1].end());
        EXPECT_TRUE(!plan.numeric(A, wrong, C));
        compare(C, original);
    });
    run_case("changed_coarse_sparsity_requires_new_symbolic_plan", [] {
        const Rows identity{{{0, 1.0}}, {{1, 1.0}}};
        const auto original = SparseMatrix::from_sorted_csr(
            2, 2, {0, 2, 4}, {0, 1, 0, 1}, {2.0, -1.0, -1.0, 2.0});
        const auto changed = SparseMatrix::from_sorted_csr(
            2, 2, {0, 2, 4}, {0, 1, 0, 1}, {2.0, 0.0, 0.0, 2.0});
        GalerkinCSRProduct first;
        EXPECT_TRUE(first.symbolic(original, identity, 2));
        SparseMatrix old_coarse;
        EXPECT_TRUE(first.numeric(original, identity, old_coarse));
        GalerkinCSRProduct second;
        EXPECT_TRUE(second.symbolic(old_coarse, identity, 2));
        SparseMatrix old_second;
        EXPECT_TRUE(second.numeric(old_coarse, identity, old_second));
        SparseMatrix new_coarse;
        EXPECT_TRUE(first.numeric(changed, identity, new_coarse));
        EXPECT_TRUE(new_coarse.nnz() < old_coarse.nnz());
        const SparseMatrix retained = old_second;
        EXPECT_TRUE(!second.numeric(new_coarse, identity, old_second));
        compare(old_second, retained);
    });
    run_case("constant_near_nullspace_transfer_matches_reference", [] {
        const auto P = transfers();
        auto A = fine();
        GalerkinCSRProduct plan;
        EXPECT_TRUE(plan.symbolic(A, P, 3));
        SparseMatrix C;
        EXPECT_TRUE(plan.numeric(A, P, C));
        const auto ones = std::vector<double>(3, 1.0);
        std::vector<double> prolonged(A.n_rows(), 0.0);
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            for (const auto& [col, weight] : P[i]) prolonged[i] += weight * ones[col];
            EXPECT_TRUE(std::abs(prolonged[i] - 1.0) < 1e-14);
        }
        const auto ap = A.matvec(prolonged);
        std::vector<double> ptap(3, 0.0);
        for (std::size_t i = 0; i < A.n_rows(); ++i)
            for (const auto& [col, weight] : P[i]) ptap[col] += weight * ap[i];
        const auto csr_product = C.matvec(ones);
        for (std::size_t i = 0; i < 3; ++i)
            EXPECT_TRUE(std::abs(csr_product[i] - ptap[i]) < 2e-13);
    });
    run_case("cf_cap_preserves_constant_mode_and_valid_gmres", [] {
        const auto A = anisotropic(24, 1000.0, 1.0);
        const FunctionalLinearOperator op(A.n_rows(), [&A](const Vector& x, Vector& y) {
            A.matvec_into(x, y);
        });
        MatrixFreeVcyclePreconditioner default_amg(
            op, 0.7, 3, 3, 0.25, 25, AMGInterpolationPolicy::DirectCF,
            false, false, AMGSmootherPolicy::AlgebraicLines, 0.0, 0);
        MatrixFreeVcyclePreconditioner capped_amg(
            op, 0.7, 3, 3, 0.25, 25, AMGInterpolationPolicy::DirectCF,
            false, false, AMGSmootherPolicy::AlgebraicLines, 0.0, 4);
        EXPECT_TRUE(default_amg.setup(A));
        EXPECT_TRUE(capped_amg.setup(A));
        EXPECT_TRUE(capped_amg.operator_complexity() <=
                    default_amg.operator_complexity() + 1e-12);
        EXPECT_TRUE(std::abs(capped_amg.prolongation_row_sum_min() - 1.0) < 1e-12);
        EXPECT_TRUE(std::abs(capped_amg.prolongation_row_sum_max() - 1.0) < 1e-12);
        Vector x_ref(A.n_rows());
        for (std::size_t i = 0; i < x_ref.size(); ++i)
            x_ref(i) = std::sin(0.013 * static_cast<double>(i));
        Vector rhs(A.n_rows()), x(A.n_rows(), 0.0);
        A.matvec_into(x_ref, rhs);
        const auto result = solve_fgmres(A, rhs, x, 30, 150, 1e-10, &capped_amg,
                                          KrylovControls{}, PrecisionPolicy{}, false);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        Vector ax(A.n_rows());
        A.matvec_into(x, ax);
        double norm = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i)
            norm = std::hypot(norm, rhs(i) - ax(i));
        EXPECT_TRUE(norm <= 1e-10 * rhs.norm2());
    });
    run_case("cf_cap_anisotropic_x_y_ratio_million_and_rotated_45", [] {
        const std::vector<SparseMatrix> matrices{
            anisotropic(20, 1e6, 1.0),
            anisotropic(20, 1.0, 1e6),
            rotated_anisotropic(20, 1000.0, 1.0)};
        for (const auto& A : matrices) {
            const FunctionalLinearOperator op(A.n_rows(), [&A](const Vector& x, Vector& y) {
                A.matvec_into(x, y);
            });
            MatrixFreeVcyclePreconditioner amg(
                op, 0.7, 3, 3, 0.25, 25, AMGInterpolationPolicy::DirectCF,
                false, false, AMGSmootherPolicy::AlgebraicLines, 0.0, 4);
            EXPECT_TRUE(amg.setup(A));
            EXPECT_TRUE(std::abs(amg.prolongation_row_sum_min() - 1.0) < 1e-12);
            EXPECT_TRUE(std::abs(amg.prolongation_row_sum_max() - 1.0) < 1e-12);
            Vector x_ref(A.n_rows());
            for (std::size_t i = 0; i < x_ref.size(); ++i)
                x_ref(i) = std::cos(0.037 * static_cast<double>(i));
            Vector rhs(A.n_rows()), x(A.n_rows(), 0.0);
            A.matvec_into(x_ref, rhs);
            const auto result = solve_fgmres(A, rhs, x, 30, 300, 1e-10, &amg,
                                             KrylovControls{}, PrecisionPolicy{}, false);
            EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
            Vector ax(A.n_rows());
            A.matvec_into(x, ax);
            double true_residual = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i)
                true_residual = std::hypot(true_residual, rhs(i) - ax(i));
            EXPECT_TRUE(true_residual <= rhs.norm2() * 1e-10);
        }
    });
    run_case("cf_cap_projected_neumann_anisotropy_pressure", [] {
        // Remove only the reaction from the manufactured diffusion matrix.
        // The resulting operator has a true constant pressure nullspace.
        const SparseMatrix reaction = anisotropic(16, 1000.0, 1.0);
        std::vector<SparseMatrix::Index> offsets(
            reaction.row_offsets_data(), reaction.row_offsets_data() + reaction.n_rows() + 1);
        std::vector<SparseMatrix::Index> columns(
            reaction.columns_data(), reaction.columns_data() + reaction.nnz());
        std::vector<double> values(reaction.values_data(),
                                   reaction.values_data() + reaction.nnz());
        for (std::size_t i = 0; i < reaction.n_rows(); ++i)
            for (std::size_t k = offsets[i]; k < offsets[i + 1]; ++k)
                if (columns[k] == i) values[k] -= 1.0;
        const SparseMatrix A = SparseMatrix::from_sorted_csr(
            reaction.n_rows(), reaction.n_cols(),
            std::move(offsets), std::move(columns), std::move(values));
        const auto space = NullSpaceProjector::constant(A.n_rows());
        EXPECT_TRUE(is_bilateral_null_space(A, space));
        const FunctionalLinearOperator op(A.n_rows(), [&A](const Vector& x, Vector& y) {
            A.matvec_into(x, y);
        });
        MatrixFreeVcyclePreconditioner pc(
            op, 0.7, 3, 3, 0.25, 25, AMGInterpolationPolicy::DirectCF,
            true, false, AMGSmootherPolicy::AlgebraicLines, 0.0, 4);
        Vector exact(A.n_rows());
        for (std::size_t i = 0; i < exact.size(); ++i)
            exact(i) = std::sin(0.027 * static_cast<double>(i));
        space.remove(exact);
        Vector rhs(A.n_rows()), solution(A.n_rows(), 3.0);
        A.matvec_into(exact, rhs);
        const auto result = solve_projected_gmres(
            A, rhs, solution, space, 30, 350, 1e-10, pc);
        EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
        Vector ax(A.n_rows());
        A.matvec_into(solution, ax);
        double true_residual = 0.0;
        for (std::size_t i = 0; i < solution.size(); ++i)
            true_residual = std::hypot(true_residual, rhs(i) - ax(i));
        EXPECT_TRUE(true_residual <= rhs.norm2() * 1e-10);
        EXPECT_TRUE(space.component_norm(solution) < 1e-10);
    });
    return run_all();
}
