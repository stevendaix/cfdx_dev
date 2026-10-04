#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <stdexcept>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;
using Dense = std::vector<std::vector<double>>;

static SparseMatrix make_sparse(
    std::size_t rows, std::size_t cols,
    const std::initializer_list<std::tuple<std::size_t, std::size_t, double>>& entries) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

static Dense dense(const SparseMatrix& A) {
    Dense out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        for (std::size_t k = A.row_offsets_data()[i]; k < A.row_offsets_data()[i + 1]; ++k) {
            out[i][A.columns_data()[k]] += A.values_data()[k];
        }
    }
    return out;
}

static std::vector<double> dense_matvec(const Dense& A, const std::vector<double>& x) {
    std::vector<double> y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y[i] += A[i][j] * x[j];
    return y;
}

// Independent reference implementation used only by this verification test.
// It intentionally does not call CFDX Krylov/Schur code.
static std::vector<double> dense_solve(Dense A, std::vector<double> b) {
    const std::size_t n = A.size();
    if (b.size() != n) throw std::invalid_argument("dense_solve: dimension mismatch");

    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        }
        if (std::abs(A[pivot][k]) < 1e-13) throw std::runtime_error("dense_solve: singular");
        std::swap(A[k], A[pivot]);
        std::swap(b[k], b[pivot]);

        for (std::size_t i = k + 1; i < n; ++i) {
            const double factor = A[i][k] / A[k][k];
            for (std::size_t j = k; j < n; ++j) A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    std::vector<double> x(n, 0.0);
    for (std::size_t ii = 0; ii < n; ++ii) {
        const std::size_t i = n - 1 - ii;
        double rhs = b[i];
        for (std::size_t j = i + 1; j < n; ++j) rhs -= A[i][j] * x[j];
        x[i] = rhs / A[i][i];
    }
    return x;
}

static Dense block_dense(const Dense& M, const Dense& G, const Dense& D, const Dense& C) {
    const std::size_t nu = M.size();
    const std::size_t np = C.size();
    Dense A(nu + np, std::vector<double>(nu + np, 0.0));
    for (std::size_t i = 0; i < nu; ++i) {
        for (std::size_t j = 0; j < nu; ++j) A[i][j] = M[i][j];
        for (std::size_t j = 0; j < np; ++j) A[i][nu + j] = G[i][j];
    }
    for (std::size_t i = 0; i < np; ++i) {
        for (std::size_t j = 0; j < nu; ++j) A[nu + i][j] = D[i][j];
        for (std::size_t j = 0; j < np; ++j) A[nu + i][nu + j] = C[i][j];
    }
    return A;
}

static double norm2(const std::vector<double>& x) {
    double sum = 0.0;
    for (const double value : x) sum += value * value;
    return std::sqrt(sum);
}

int main() {
    // Small nonsymmetric saddle-point system. M is deliberately non-diagonal
    // and D/G are independent so the test exercises the actual 4-block
    // contract rather than a symmetric special case.
    const auto M = make_sparse(3, 3, {
        {0,0,5.0}, {0,1,0.7}, {0,2,-0.2},
        {1,0,0.4}, {1,1,4.0}, {1,2,0.6},
        {2,0,-0.1}, {2,1,0.3}, {2,2,3.5}});
    const auto G = make_sparse(3, 2, {
        {0,0,1.0}, {0,1,0.2},
        {1,0,-0.4}, {1,1,0.8},
        {2,0,0.3}, {2,1,1.1}});
    const auto D = make_sparse(2, 3, {
        {0,0,0.9}, {0,1,-0.2}, {0,2,0.4},
        {1,0,0.1}, {1,1,0.7}, {1,2,0.8}});
    const auto C = make_sparse(2, 2, {
        {0,0,1.8}, {0,1,0.15},
        {1,0,-0.1}, {1,1,2.2}});
    const BlockOperator blocks(M, G, D, C);
    blocks.validate();

    const Dense Md = dense(M);
    const Dense Gd = dense(G);
    const Dense Dd = dense(D);
    const Dense Cd = dense(C);
    const Dense block = block_dense(Md, Gd, Dd, Cd);

    // The same independent dense M inverse is supplied to the exact Schur
    // oracle. The direct block solve itself remains an independent dense solve.
    const auto solve_M = [Md](const Vector& rhs, Vector& y) {
        std::vector<double> b(rhs.size(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i) b[i] = rhs(i);
        const auto x = dense_solve(Md, b);
        y = Vector(x.size(), 0.0);
        for (std::size_t i = 0; i < x.size(); ++i) y(i) = x[i];
        return true;
    };

    run_case("n9_3_direct_block_solve_has_true_residual", [&] {
        const std::vector<double> rhs = {1.1, -0.3, 0.7, 0.4, -0.8};
        const auto x = dense_solve(block, rhs);
        const auto Ax = dense_matvec(block, x);

        std::vector<double> residual(rhs.size(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i) residual[i] = Ax[i] - rhs[i];

        EXPECT_TRUE(norm2(residual) < 1e-12);
    });

    run_case("n9_3_direct_block_solve_matches_exact_schur_reconstruction", [&] {
        const std::vector<double> rhs_u = {1.1, -0.3, 0.7};
        const std::vector<double> rhs_p = {0.4, -0.8};
        const std::vector<double> rhs = {rhs_u[0], rhs_u[1], rhs_u[2], rhs_p[0], rhs_p[1]};

        const auto direct = dense_solve(block, rhs);

        ExactSchurApproximation schur(solve_M);
        EXPECT_TRUE(schur.setup(blocks));

        Vector pressure_rhs(2, 0.0);
        for (std::size_t i = 0; i < 2; ++i) pressure_rhs(i) = rhs_p[i];

        // Schur RHS = b_p - D M^{-1} b_u.
        Vector bu(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i) bu(i) = rhs_u[i];
        Vector m_inv_bu(3, 0.0);
        EXPECT_TRUE(solve_M(bu, m_inv_bu));
        const auto Dm = Dd;
        const auto mvec = std::vector<double>{m_inv_bu(0), m_inv_bu(1), m_inv_bu(2)};
        const auto Dm_bu = dense_matvec(Dm, mvec);
        for (std::size_t i = 0; i < 2; ++i) pressure_rhs(i) -= Dm_bu[i];

        Vector pressure(2, 0.0);
        EXPECT_TRUE(schur.apply(pressure_rhs, pressure));

        Vector velocity_rhs(3, 0.0);
        const auto Gp = dense_matvec(Gd, {pressure(0), pressure(1)});
        for (std::size_t i = 0; i < 3; ++i) velocity_rhs(i) = rhs_u[i] - Gp[i];
        Vector velocity(3, 0.0);
        EXPECT_TRUE(solve_M(velocity_rhs, velocity));

        for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(velocity(i), direct[i], 1e-9);
        for (std::size_t i = 0; i < 2; ++i) EXPECT_NEAR(pressure(i), direct[3 + i], 1e-9);
    });

    run_case("n9_3_direct_and_schur_solutions_have_same_full_residual", [&] {
        const std::vector<double> rhs = {1.1, -0.3, 0.7, 0.4, -0.8};
        const auto direct = dense_solve(block, rhs);

        const auto Ax = dense_matvec(block, direct);
        std::vector<double> residual(rhs.size(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i) residual[i] = Ax[i] - rhs[i];
        EXPECT_TRUE(norm2(residual) < 1e-12);

        ExactSchurApproximation schur(solve_M);
        EXPECT_TRUE(schur.setup(blocks));
        Vector pressure_rhs(2, 0.0);
        pressure_rhs(0) = rhs[3];
        pressure_rhs(1) = rhs[4];

        Vector bu(3, 0.0);
        for (std::size_t i = 0; i < 3; ++i) bu(i) = rhs[i];
        Vector m_inv_bu(3, 0.0);
        EXPECT_TRUE(solve_M(bu, m_inv_bu));
        const auto tmp = dense_matvec(Dd, {m_inv_bu(0), m_inv_bu(1), m_inv_bu(2)});
        pressure_rhs(0) -= tmp[0];
        pressure_rhs(1) -= tmp[1];

        Vector pressure(2, 0.0);
        EXPECT_TRUE(schur.apply(pressure_rhs, pressure));
        Vector velocity_rhs(3, 0.0);
        const auto Gp = dense_matvec(Gd, {pressure(0), pressure(1)});
        for (std::size_t i = 0; i < 3; ++i) velocity_rhs(i) = rhs[i] - Gp[i];
        Vector velocity(3, 0.0);
        EXPECT_TRUE(solve_M(velocity_rhs, velocity));

        std::vector<double> coupled = {
            velocity(0), velocity(1), velocity(2), pressure(0), pressure(1)};
        const auto coupled_Ax = dense_matvec(block, coupled);
        std::vector<double> coupled_residual(rhs.size(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i) coupled_residual[i] = coupled_Ax[i] - rhs[i];
        EXPECT_TRUE(norm2(coupled_residual) < 1e-9);
    });

    return run_all();
}
