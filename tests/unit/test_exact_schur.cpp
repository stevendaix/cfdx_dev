#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

using Dense = std::vector<std::vector<double>>;

static SparseMatrix make_sparse(std::size_t rows,
                                std::size_t cols,
                                const std::initializer_list<std::tuple<std::size_t,
                                                                         std::size_t,
                                                                         double>>& entries) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

static Dense to_dense(const SparseMatrix& A) {
    Dense out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i];
             k < A.row_offsets_data()[i + 1]; ++k)
            out[i][A.columns_data()[k]] += A.values_data()[k];
    return out;
}

static Dense dense_multiply(const Dense& A, const Dense& B) {
    if (A.empty() || B.empty() || A[0].size() != B.size())
        throw std::invalid_argument("dense_multiply: dimension mismatch");
    Dense C(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t k = 0; k < B.size(); ++k)
            for (std::size_t j = 0; j < B[0].size(); ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

static Dense dense_subtract(const Dense& A, const Dense& B) {
    Dense C = A;
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[0].size(); ++j) C[i][j] -= B[i][j];
    return C;
}

static Dense dense_identity(std::size_t n) {
    Dense I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;
    return I;
}

static Dense dense_inverse(Dense A) {
    const std::size_t n = A.size();
    if (n == 0 || A[0].size() != n)
        throw std::invalid_argument("dense_inverse: non-square matrix");
    Dense I = dense_identity(n);
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        if (std::abs(A[pivot][k]) <= 1e-14)
            throw std::invalid_argument("dense_inverse: singular matrix");
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);
        const double d = A[k][k];
        for (std::size_t j = 0; j < n; ++j) { A[k][j] /= d; I[k][j] /= d; }
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double factor = A[i][k];
            for (std::size_t j = 0; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
                I[i][j] -= factor * I[k][j];
            }
        }
    }
    return I;
}

static Vector dense_matvec(const Dense& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

int main() {
    // Reference case: symmetric positive-definite Auu, non-zero C, SPD Schur.
    // Chosen so that Auu != block_diag(Auu) (i.e., there are non-zero
    // off-diagonal Auu entries) to distinguish the exact Schur from the
    // block-local approximation.
    const auto Auu = make_sparse(2, 2, {
        {0, 0, 4.0}, {0, 1, 1.0},
        {1, 0, 1.0}, {1, 1, 3.0}
    });
    const auto G = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.5},
        {1, 0, 0.25}, {1, 1, 1.0}
    });
    // Take D = G^T so the assembled saddle-point matrix is symmetric and S
    // is SPD.
    const auto D = make_sparse(2, 2, {
        {0, 0, 1.0},  {0, 1, 0.25},
        {1, 0, 0.5},  {1, 1, 1.0}
    });
    const auto C = make_sparse(2, 2, {
        {0, 0, 2.0}, {0, 1, 0.0},
        {1, 0, 0.0}, {1, 1, 3.0}
    });

    const BlockOperator blocks(Auu, G, D, C);
    blocks.validate();

    // Independent dense oracle for the exact Schur complement.
    const Dense Auud = to_dense(Auu);
    const Dense Gd = to_dense(G);
    const Dense Dd = to_dense(D);
    const Dense Cd = to_dense(C);
    const Dense Auu_inv = dense_inverse(Auud);
    const Dense S_dense = dense_subtract(
        Cd, dense_multiply(dense_multiply(Dd, Auu_inv), Gd));
    const Dense S_dense_inv = dense_inverse(S_dense);

    // Auu solve callback: dense inverse. This is deliberately not
    // production code but it isolates the ExactSchur logic from the
    // choice of inner solver.
    auto auu_solve = [&Auu_inv](const Vector& rhs, Vector& y) {
        if (y.size() != rhs.size()) y = Vector(rhs.size(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < rhs.size(); ++j)
                s += Auu_inv[i][j] * rhs(j);
            y(i) = s;
        }
        return true;
    };

    run_case("exact_schur_setup_rejects_invalid_state", [&] {
        ExactSchurApproximation missing_solve({});
        EXPECT_TRUE(!missing_solve.setup(blocks));

        // Mismatched dimensions must fail is_valid() and be rejected.
        const auto Bad = make_sparse(3, 2, {{0, 0, 1.0}});
        const BlockOperator bad_blocks(Auu, Bad, D, C);
        EXPECT_TRUE(!bad_blocks.is_valid());
        ExactSchurApproximation ok_solve(auu_solve);
        EXPECT_TRUE(!ok_solve.setup(bad_blocks));
    });

    run_case("exact_schur_apply_matches_dense_oracle", [&] {
        ExactSchurApproximation schur(auu_solve);
        EXPECT_TRUE(schur.setup(blocks));

        // Verify S*p matches the dense oracle for a non-trivial p.
        Vector p_test(2, 0.0);
        p_test(0) = 1.3;
        p_test(1) = -0.7;
        Vector Sp(2, 0.0);
        EXPECT_TRUE(schur.apply_schur(p_test, Sp));
        const Vector Sp_ref = dense_matvec(S_dense, p_test);
        EXPECT_NEAR(Sp(0), Sp_ref(0), 1e-12);
        EXPECT_NEAR(Sp(1), Sp_ref(1), 1e-12);
    });

    run_case("exact_schur_solve_matches_dense_inverse", [&] {
        ExactSchurApproximation schur(auu_solve);
        EXPECT_TRUE(schur.setup(blocks));

        Vector rhs_p(2, 0.0);
        rhs_p(0) = 1.0;
        rhs_p(1) = 2.0;

        Vector p(2, 0.0);
        EXPECT_TRUE(schur.apply(rhs_p, p));

        // Reference: p_ref = S_dense^{-1} * rhs_p.
        const Vector p_ref = dense_matvec(S_dense_inv, rhs_p);
        EXPECT_NEAR(p(0), p_ref(0), 1e-9);
        EXPECT_NEAR(p(1), p_ref(1), 1e-9);

        // True residual check: S * p should equal rhs_p to CG tolerance.
        Vector residual(2, 0.0);
        EXPECT_TRUE(schur.apply_schur(p, residual));
        for (std::size_t i = 0; i < 2; ++i) residual(i) = rhs_p(i) - residual(i);
        EXPECT_TRUE(residual.norm2() < 1e-9 * rhs_p.norm2());
    });

    run_case("exact_schur_zero_rhs_returns_zero", [&] {
        ExactSchurApproximation schur(auu_solve);
        EXPECT_TRUE(schur.setup(blocks));
        Vector rhs_p(2, 0.0);
        Vector p(2, 42.0);
        EXPECT_TRUE(schur.apply(rhs_p, p));
        EXPECT_NEAR(p(0), 0.0, 1e-14);
        EXPECT_NEAR(p(1), 0.0, 1e-14);
    });

    run_case("exact_schur_solution_differs_from_block_local_when_Auu_is_coupled", [&] {
        // The exact Schur uses Auu^{-1}; the block-local Schur uses the
        // block-diagonal M_b^{-1}. When Auu has non-zero off-diagonals
        // (as here: entries (0,1)=(1,0)=1.0), the two Schur operators must
        // produce measurably different solutions for the same rhs.
        //
        // Reference block-local Schur assembled densely with M_b =
        // diag(Auu_diagonal).
        Dense Mb(2, std::vector<double>(2, 0.0));
        Mb[0][0] = Auud[0][0];
        Mb[1][1] = Auud[1][1];
        const Dense Mb_inv = dense_inverse(Mb);
        const Dense S_block = dense_subtract(
            Cd, dense_multiply(dense_multiply(Dd, Mb_inv), Gd));

        Vector rhs_p(2, 0.0);
        rhs_p(0) = 1.0;
        rhs_p(1) = 2.0;

        // Exact solve.
        ExactSchurApproximation exact(auu_solve);
        EXPECT_TRUE(exact.setup(blocks));
        Vector p_exact(2, 0.0);
        EXPECT_TRUE(exact.apply(rhs_p, p_exact));

        // Block-local solve via dense inverse.
        const Vector p_block = dense_matvec(dense_inverse(S_block), rhs_p);

        // The two solutions must differ by more than round-off, otherwise
        // this test case is not exercising the difference between S and S~.
        const double delta = std::sqrt(
            (p_exact(0) - p_block(0)) * (p_exact(0) - p_block(0)) +
            (p_exact(1) - p_block(1)) * (p_exact(1) - p_block(1)));
        EXPECT_TRUE(delta > 1e-6);
    });

    run_case("exact_schur_declares_the_inverse_action_and_behaves_as_it", [&] {
        ExactSchurApproximation exact(auu_solve);
        EXPECT_TRUE(exact.setup(blocks));

        // The exact oracle runs CG on S p = rhs_p, so it is the inverse action.
        // This is the opposite of SIMPLE/SIMPLEC, which apply the operator itself.
        EXPECT_TRUE(exact.action() == SchurAction::InverseOperator);
        EXPECT_TRUE(provides_action(exact, SchurAction::InverseOperator));
        EXPECT_TRUE(!provides_action(exact, SchurAction::Operator));

        Vector rhs_p(2, 0.0);
        rhs_p(0) = 1.0;
        rhs_p(1) = 2.0;

        Vector p_inv(2, 0.0);
        EXPECT_TRUE(exact.apply(rhs_p, p_inv));
        const Vector inverse_ref = dense_matvec(S_dense_inv, rhs_p);
        EXPECT_NEAR(p_inv(0), inverse_ref(0), 1e-9);
        EXPECT_NEAR(p_inv(1), inverse_ref(1), 1e-9);

        // The operator action on the same right-hand side must differ, otherwise
        // the assertions above could not distinguish the two actions.
        const Vector operator_ref = dense_matvec(S_dense, rhs_p);
        EXPECT_TRUE(std::abs(operator_ref(0) - inverse_ref(0)) > 1e-6 ||
                    std::abs(operator_ref(1) - inverse_ref(1)) > 1e-6);
    });

    return run_all();
}
