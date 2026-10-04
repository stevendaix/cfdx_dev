#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <initializer_list>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;
using Dense = std::vector<std::vector<double>>;

static SparseMatrix make_sparse(
    std::size_t rows, std::size_t cols,
    const std::initializer_list<std::tuple<std::size_t, std::size_t, double>>& e) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, v] : e) A.push_back(i, j, v);
    A.finalize();
    return A;
}

static Dense dense(const SparseMatrix& A) {
    Dense d(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i]; k < A.row_offsets_data()[i + 1]; ++k)
            d[i][A.columns_data()[k]] += A.values_data()[k];
    return d;
}

static Dense inverse(Dense A) {
    const std::size_t n = A.size();
    Dense I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        if (std::abs(A[pivot][k]) < 1e-14) throw std::runtime_error("singular");
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);
        const double d = A[k][k];
        for (std::size_t j = 0; j < n; ++j) { A[k][j] /= d; I[k][j] /= d; }
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double f = A[i][k];
            for (std::size_t j = 0; j < n; ++j) {
                A[i][j] -= f * A[k][j];
                I[i][j] -= f * I[k][j];
            }
        }
    }
    return I;
}

static Vector matvec(const Dense& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

static double max_abs(const Vector& a, const Vector& b) {
    double e = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) e = std::max(e, std::abs(a(i) - b(i)));
    return e;
}

int main() {
    // Deliberately non-diagonal M and non-transposed D/G: this is an
    // algebraic N9 oracle test, not a special symmetric saddle point.
    const auto M = make_sparse(2, 2, {{0,0,4.0},{0,1,1.0},{1,0,0.5},{1,1,3.0}});
    const auto G = make_sparse(2, 2, {{0,0,1.0},{0,1,0.5},{1,0,0.25},{1,1,1.0}});
    const auto D = make_sparse(2, 2, {{0,0,0.8},{0,1,0.1},{1,0,0.3},{1,1,0.9}});
    const auto C = make_sparse(2, 2, {{0,0,2.0},{0,1,0.2},{1,0,0.1},{1,1,2.5}});
    const BlockOperator blocks(M, G, D, C);
    blocks.validate();

    const Dense Md = dense(M);
    const Dense Gd = dense(G);
    const Dense Dd = dense(D);
    const Dense Cd = dense(C);
    const Dense Mi = inverse(Md);

    auto solve_M = [Mi](const Vector& rhs, Vector& y) {
        y = matvec(Mi, rhs);
        return true;
    };

    run_case("n9_schur_action_matches_independent_dense_oracle", [&] {
        ExactSchurApproximation schur(solve_M);
        EXPECT_TRUE(schur.setup(blocks));

        Vector p(2, 0.0);
        p(0) = 1.7; p(1) = -0.6;
        Vector got(2, 0.0);
        EXPECT_TRUE(schur.apply_schur(p, got));

        const Vector gp = matvec(Gd, p);
        const Vector y = matvec(Mi, gp);
        Vector expected = matvec(Cd, p);
        const Vector dy = matvec(Dd, y);
        for (std::size_t i = 0; i < expected.size(); ++i) expected(i) -= dy(i);
        EXPECT_TRUE(max_abs(got, expected) < 1e-12);
    });

    run_case("n9_schur_solution_has_independent_true_block_residual", [&] {
        ExactSchurApproximation schur(solve_M);
        EXPECT_TRUE(schur.setup(blocks));

        Vector rhs_p(2, 0.0);
        rhs_p(0) = 1.0; rhs_p(1) = 0.7;
        Vector p(2, 0.0);
        EXPECT_TRUE(schur.apply(rhs_p, p));

        // Reconstruct U from the momentum equation M U + G p = bU with bU=0.
        Vector gp = matvec(Gd, p);
        for (std::size_t i = 0; i < gp.size(); ++i) gp(i) = -gp(i);
        Vector u = matvec(Mi, gp);

        // Independent full block residual: [M G; D C][u p] - [0 rhs_p].
        Vector ru = matvec(Md, u);
        const Vector Gp = matvec(Gd, p);
        for (std::size_t i = 0; i < ru.size(); ++i) ru(i) += Gp(i);
        Vector rp = matvec(Dd, u);
        const Vector Cp = matvec(Cd, p);
        for (std::size_t i = 0; i < rp.size(); ++i) rp(i) += Cp(i) - rhs_p(i);

        EXPECT_TRUE(ru.norm2() < 1e-10);
        EXPECT_TRUE(rp.norm2() < 1e-9);
    });

    run_case("n9_schur_null_space_rejects_incompatible_rhs_and_projects_compatible", [&] {
        // C=0 and D*1=0 create the pressure constant null mode in the
        // Schur operator. The oracle must reject an incompatible RHS rather
        // than silently altering it.
        const auto G0 = make_sparse(2, 2, {{0,0,1.0},{0,1,-1.0},{1,0,-1.0},{1,1,1.0}});
        const auto C0 = make_sparse(2, 2, {});
        const auto D0 = make_sparse(2, 2, {{0,0,-1.0},{0,1,1.0},{1,0,1.0},{1,1,-1.0}});
        const BlockOperator singular(M, G0, D0, C0);
        singular.validate();

        ExactSchurApproximation schur(solve_M, {}, NullSpaceProjector::constant(2));
        EXPECT_TRUE(schur.setup(singular));

        Vector incompatible(2, 0.0);
        incompatible(0) = 1.0; incompatible(1) = 0.0;
        Vector p(2, 0.0);
        EXPECT_TRUE(!schur.apply(incompatible, p));

        Vector compatible(2, 0.0);
        compatible(0) = 1.0; compatible(1) = -1.0;
        EXPECT_TRUE(schur.apply(compatible, p));
        EXPECT_TRUE(std::abs(p(0) + p(1)) < 1e-9);
    });

    run_case("n9_schur_reference_pressure_is_explicit_reduced_system", [&] {
        // A reference-cell gauge is represented explicitly by eliminating
        // pressure cell 0. The remaining 1x1 Schur is solved by the same
        // exact algebraic definition; no hidden gauge is applied.
        const auto Gr = make_sparse(2, 1, {{0,0,0.5},{1,0,1.0}});
        const auto Dr = make_sparse(1, 2, {{0,0,0.3},{0,1,0.9}});
        const auto Cr = make_sparse(1, 1, {{0,0,2.5}});
        const BlockOperator reduced(M, Gr, Dr, Cr);
        reduced.validate();

        ExactSchurApproximation schur(solve_M);
        EXPECT_TRUE(schur.setup(reduced));
        Vector rhs(1, 0.0); rhs(0) = 1.2;
        Vector p(1, 0.0);
        EXPECT_TRUE(schur.apply(rhs, p));

        const Dense Gred = dense(Gr);
        const Dense Dred = dense(Dr);
        const Dense Cred = dense(Cr);
        const double s = Cred[0][0] - (Dred[0][0] * Mi[0][0] * Gred[0][0]
            + Dred[0][0] * Mi[0][1] * Gred[1][0]
            + Dred[0][1] * Mi[1][0] * Gred[0][0]
            + Dred[0][1] * Mi[1][1] * Gred[1][0]);
        EXPECT_NEAR(s * p(0), rhs(0), 1e-10);
    });

    return run_all();
}
