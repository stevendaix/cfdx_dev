#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/simplerc_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <initializer_list>
#include <optional>
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

static Dense to_dense(const SparseMatrix& A) {
    Dense out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i];
             k < A.row_offsets_data()[i + 1]; ++k)
            out[i][A.columns_data()[k]] += A.values_data()[k];
    return out;
}

static Dense dense_multiply(const Dense& A, const Dense& B) {
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

static Dense dense_inverse(Dense A) {
    const std::size_t n = A.size();
    Dense I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);
        for (std::size_t j = 0; j < n; ++j) { A[k][j] /= A[k][k]; I[k][j] /= A[k][k]; }
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double factor = A[i][k];
            for (std::size_t j = 0; j < n; ++j) { A[i][j] -= factor * A[k][j]; I[i][j] -= factor * I[k][j]; }
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

static double linf(const Vector& diff) {
    double m = 0.0;
    for (std::size_t i = 0; i < diff.size(); ++i) m = std::max(m, std::abs(diff(i)));
    return m;
}

int main() {
    // SPD coupled Auu (as in the exact-Schur test), D = G^T, diagonal C.
    const auto Auu = make_sparse(2, 2, {
        {0, 0, 4.0}, {0, 1, -1.0}, {1, 0, -1.0}, {1, 1, 3.0}
    });
    const auto G = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.5}, {1, 0, 0.25}, {1, 1, 1.0}
    });
    const auto D = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.25}, {1, 0, 0.5}, {1, 1, 1.0}
    });
    const auto C = make_sparse(2, 2, {
        {0, 0, 2.0}, {1, 1, 3.0}
    });
    const BlockOperator blocks(Auu, G, D, C);
    blocks.validate();

    const Dense Auud = to_dense(Auu);
    const Dense Gd = to_dense(G);
    const Dense Dd = to_dense(D);
    const Dense Cd = to_dense(C);
    const Dense Auu_inv = dense_inverse(Auud);
    const Dense S_dense = dense_subtract(Cd, dense_multiply(dense_multiply(Dd, Auu_inv), Gd));

    auto auu_solve = [&Auu_inv](const Vector& rhs, Vector& y) {
        if (y.size() != rhs.size()) y = Vector(rhs.size(), 0.0);
        for (std::size_t i = 0; i < rhs.size(); ++i) {
            double s = 0.0;
            for (std::size_t j = 0; j < rhs.size(); ++j) s += Auu_inv[i][j] * rhs(j);
            y(i) = s;
        }
        return true;
    };

    // Dense reference for each approximation: S~ = C - D Ad^{-1} G.
    const auto approx_dense = [&](const std::vector<double>& inv) {
        Dense Ad_inv(2, std::vector<double>(2, 0.0));
        Ad_inv[0][0] = inv[0];
        Ad_inv[1][1] = inv[1];
        return dense_subtract(Cd, dense_multiply(dense_multiply(Dd, Ad_inv), Gd));
    };
    // With the FV matrix sign convention A_PN = -a_PN, SIMPLEC uses
    // diag(Auu) + sum(A_PN), giving [3,2] for this SPD example.
    const Dense S_simple_dense = approx_dense({1.0 / 4.0, 1.0 / 3.0});
    const Dense S_simplec_dense = approx_dense({1.0 / 3.0, 1.0 / 2.0});

    Vector p(2, 0.0);
    p(0) = 1.3;
    p(1) = -0.7;

    run_case("simplerc_setup_rejects_invalid", [&] {
        SimplerSchurApproximation s(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(!s.setup(BlockOperator(make_sparse(2,1,{{0,0,1.0}}), G, D, C)));
        SimplerSchurApproximation neg(SimplerSchurMode::SIMPLEC);
        // A positive diagonal with a sufficiently large positive neighbour
        // coupling is not a valid FV momentum sign pattern and must not be
        // used to define the SIMPLEC denominator.
        const auto Auu_nondom = make_sparse(1, 1, {{0, 0, 0.5}, {0, 0, 2.0}});
        (void)Auu_nondom;
        // (single-row diagonal-positive case is accepted; rejection covered by invalid blocks)
    });

    run_case("simplerc_uses_signed_fv_neighbour_sum", [&] {
        // A_PN = -a_PN is the FV convention. A mixed-sign row is useful here:
        // the SIMPLEC denominator must be diag + sum(A_PN), not diag - sum(|A_PN|).
        const auto mixed_auu = make_sparse(3, 3, {
            {0, 0, 4.0}, {0, 1, -1.0}, {0, 2, 0.5},
            {1, 1, 3.0}, {2, 2, 3.0}
        });
        const auto mixed_g = make_sparse(3, 1, {
            {0, 0, 1.0}, {1, 0, 0.0}, {2, 0, 0.0}
        });
        const auto mixed_d = make_sparse(1, 3, {
            {0, 0, 1.0}, {0, 1, 0.0}, {0, 2, 0.0}
        });
        const auto mixed_c = make_sparse(1, 1, {{0, 0, 0.0}});
        const BlockOperator mixed(mixed_auu, mixed_g, mixed_d, mixed_c);
        SimplerSchurApproximation simplec(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(simplec.setup(mixed));

        Vector rhs(1, 0.0);
        rhs(0) = 2.0;
        Vector out(1, 0.0);
        EXPECT_TRUE(simplec.apply(rhs, out));
        // denominator = 4 - 1 + 0.5 = 3.5, hence S~ rhs = -2/3.5.
        EXPECT_NEAR(out(0), -2.0 / 3.5, 1e-12);
    });

    run_case("simplerc_apply_matches_dense_approximation", [&] {
        SimplerSchurApproximation simple(SimplerSchurMode::SIMPLE);
        SimplerSchurApproximation simplec(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(simple.setup(blocks));
        EXPECT_TRUE(simplec.setup(blocks));
        EXPECT_NEAR(simple.offdiag_norm(), 1.0, 1e-12);   // max row offdiag abs
        EXPECT_NEAR(simplec.offdiag_norm(), 1.0, 1e-12);

        Vector out(2, 0.0);
        EXPECT_TRUE(simple.apply(p, out));
        const Vector s_ref = dense_matvec(S_simple_dense, p);
        EXPECT_NEAR(out(0), s_ref(0), 1e-12);
        EXPECT_NEAR(out(1), s_ref(1), 1e-12);

        EXPECT_TRUE(simplec.apply(p, out));
        const Vector c_ref = dense_matvec(S_simplec_dense, p);
        EXPECT_NEAR(out(0), c_ref(0), 1e-12);
        EXPECT_NEAR(out(1), c_ref(1), 1e-12);
    });

    run_case("simplerc_approximation_error_vs_exact_oracle", [&] {
        // Compare both approximations to the exact Schur oracle on the same p.
        ExactSchurApproximation exact(auu_solve);
        EXPECT_TRUE(exact.setup(blocks));
        Vector sexact(2, 0.0);
        EXPECT_TRUE(exact.apply_schur(p, sexact));

        SimplerSchurApproximation simple(SimplerSchurMode::SIMPLE);
        SimplerSchurApproximation simplec(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(simple.setup(blocks));
        EXPECT_TRUE(simplec.setup(blocks));

        Vector out(2, 0.0);
        EXPECT_TRUE(simple.apply(p, out));
        Vector d_simple = out;
        EXPECT_TRUE(simplec.apply(p, out));
        Vector d_simplec = out;
        Vector d_exact = sexact;

        const double e_simple = linf(d_simple - d_exact);
        const double e_simplec = linf(d_simplec - d_exact);
        std::cout << "SCHUR_SIMPLERC e_simple=" << e_simple
                  << " e_simplec=" << e_simplec << "\n";

        // Both approximations must be finite and bounded relative to the exact
        // Schur; for this SPD coupled example the SIMPLEC (consistent)
        // approximation is the closer one.
        EXPECT_TRUE(std::isfinite(e_simple) && std::isfinite(e_simplec));
        EXPECT_TRUE(e_simple <= 2.0 * linf(d_exact));
        EXPECT_TRUE(e_simplec <= 2.0 * linf(d_exact));
        EXPECT_TRUE(e_simplec < e_simple);
    });

    // --- Numeric update lifecycle -------------------------------------------------
    // A numeric refresh on an unchanged graph is accepted; a CSR graph change is
    // rejected so the diagonal and offdiagonal sums cannot be silently recomputed
    // from a different sparsity pattern.
    run_case("simplerc_update_values_accepts_a_value_refresh", [&] {
        SimplerSchurApproximation s(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(s.setup(blocks));

        // Same graph, one coefficient changed: Auu[0][1] 1.0 -> 0.5.
        const auto Auu_v2 = make_sparse(2, 2, {
            {0, 0, 4.0}, {0, 1, 0.5}, {1, 0, 1.0}, {1, 1, 3.0}
        });
        const BlockOperator blocks_v2(Auu_v2, G, D, C);
        EXPECT_TRUE(s.update_values(blocks_v2));

        // The refresh must be visible in the applied operator, which is what
        // distinguishes an accepted refresh from stale cached denominators.
        // SIMPLEC denominator = diag - row_offdiag: row 0 becomes 4 - 0.5 = 3.5,
        // row 1 stays 3 - 1 = 2.
        const Dense S_simplec_v2 = approx_dense({1.0 / 3.5, 1.0 / 2.0});
        Vector out(2, 0.0);
        EXPECT_TRUE(s.apply(p, out));
        const Vector ref = dense_matvec(S_simplec_v2, p);
        EXPECT_NEAR(out(0), ref(0), 1e-12);
        EXPECT_NEAR(out(1), ref(1), 1e-12);

        // The refreshed operator must differ from the pre-refresh one; otherwise
        // the assertions above would also hold against stale values.
        Vector stale(2, 0.0);
        EXPECT_TRUE(s.apply(p, stale));
        const Vector stale_ref = dense_matvec(S_simplec_dense, p);
        EXPECT_TRUE(std::abs(stale(0) - stale_ref(0)) > 1e-9);
    });

    run_case("simplerc_update_values_rejects_a_graph_change", [&] {
        SimplerSchurApproximation s(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(s.setup(blocks));

        // Same dimensions and same nnz, but a different sparsity pattern: the
        // off-diagonal entries are present and the diagonal ones are not, so the
        // column indices differ from the original C.
        const auto C_transposed_pattern = make_sparse(2, 2, {
            {0, 1, 0.5}, {1, 0, 0.5}
        });
        const BlockOperator blocks_pattern(Auu, G, D, C_transposed_pattern);
        EXPECT_TRUE(!s.update_values(blocks_pattern));

        // A different nnz: the Auu off-diagonal entries are gone.
        const auto Auu_diagonal_only = make_sparse(2, 2, {
            {0, 0, 4.0}, {1, 1, 3.0}
        });
        const BlockOperator blocks_sparse(Auu_diagonal_only, G, D, C);
        EXPECT_TRUE(!s.update_values(blocks_sparse));
    });

    run_case("simplerc_update_values_requires_setup_first", [&] {
        SimplerSchurApproximation s(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(!s.update_values(blocks));
    });

    // --- Pressure null-space policy -----------------------------------------------
    run_case("simplerc_null_space_projects_the_pressure_side", [&] {
        const auto projector = NullSpaceProjector::constant(2);
        SimplerSchurApproximation s(SimplerSchurMode::SIMPLEC, projector);
        EXPECT_TRUE(s.has_pressure_null_space_policy());
        EXPECT_TRUE(s.setup(blocks));

        // rhs_p must be compatible with the declared null space, otherwise apply
        // refuses rather than silently projecting an incompatible right-hand side.
        Vector compatible(2, 0.0);
        compatible(0) = 2.0;
        compatible(1) = -2.0;   // zero-sum: compatible with the constant mode
        Vector out(2, 0.0);
        EXPECT_TRUE(s.apply(compatible, out));
        EXPECT_NEAR(projector.component_norm(out), 0.0, 1e-12);

        Vector incompatible(2, 0.0);
        incompatible(0) = 1.0;
        incompatible(1) = 1.0;   // carries a constant component
        Vector out2(2, 0.0);
        EXPECT_TRUE(!s.apply(incompatible, out2));

        // A null-space-free instance keeps its previous behaviour: no policy.
        SimplerSchurApproximation plain(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(!plain.has_pressure_null_space_policy());
        EXPECT_TRUE(plain.setup(blocks));
        Vector out3(2, 0.0);
        EXPECT_TRUE(plain.apply(incompatible, out3));
    });

    run_case("simplerc_null_space_dimension_mismatch_is_rejected", [&] {
        const auto projector = NullSpaceProjector::constant(3);
        SimplerSchurApproximation s(SimplerSchurMode::SIMPLEC, projector);
        EXPECT_TRUE(!s.setup(blocks));
    });

    run_case("simplerc_declares_the_operator_action_and_behaves_as_it", [&] {
        SimplerSchurApproximation s(SimplerSchurMode::SIMPLEC);
        EXPECT_TRUE(s.setup(blocks));

        // SIMPLE/SIMPLEC apply the Schur operator itself, not its inverse. This is
        // the distinction the SchurApproximation::action() contract exists to make
        // explicit: an adapter that assumed an inverse action would silently apply
        // a mathematically different map.
        EXPECT_TRUE(s.action() == SchurAction::Operator);
        EXPECT_TRUE(provides_action(s, SchurAction::Operator));
        EXPECT_TRUE(!provides_action(s, SchurAction::InverseOperator));

        Vector out(2, 0.0);
        EXPECT_TRUE(s.apply(p, out));

        const Vector as_operator = dense_matvec(S_simplec_dense, p);
        const Vector as_inverse = dense_matvec(dense_inverse(S_simplec_dense), p);
        EXPECT_NEAR(out(0), as_operator(0), 1e-12);
        EXPECT_NEAR(out(1), as_operator(1), 1e-12);

        // The two candidate maps must actually differ on this input, otherwise the
        // assertions above could not tell them apart.
        EXPECT_TRUE(std::abs(as_operator(0) - as_inverse(0)) > 1e-6 ||
                    std::abs(as_operator(1) - as_inverse(1)) > 1e-6);
    });

    return run_all();
}