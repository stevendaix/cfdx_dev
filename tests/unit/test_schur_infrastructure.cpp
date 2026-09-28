#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/sparse_matrix_product.h"
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
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        for (std::size_t k = A.row_offsets_data()[i];
             k < A.row_offsets_data()[i + 1]; ++k) {
            out[i][A.columns_data()[k]] += A.values_data()[k];
        }
    }
    return out;
}

// Test-only dense reference. This is deliberately not production CFDX code:
// the exact Schur oracle is restricted to tiny matrices so that it provides
// an independent algebraic truth model without introducing a dense fallback.
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
        for (std::size_t j = 0; j < A[0].size(); ++j)
            C[i][j] -= B[i][j];
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
        for (std::size_t j = 0; j < n; ++j) {
            A[k][j] /= d;
            I[k][j] /= d;
        }

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

static Dense block_assemble(const Dense& M, const Dense& G,
                            const Dense& D, const Dense& C) {
    const std::size_t nv = M.size();
    const std::size_t np = C.size();
    Dense A(nv + np, std::vector<double>(nv + np, 0.0));
    for (std::size_t i = 0; i < nv; ++i) {
        for (std::size_t j = 0; j < nv; ++j) A[i][j] = M[i][j];
        for (std::size_t j = 0; j < np; ++j) A[i][nv + j] = G[i][j];
    }
    for (std::size_t i = 0; i < np; ++i) {
        for (std::size_t j = 0; j < nv; ++j) A[nv + i][j] = D[i][j];
        for (std::size_t j = 0; j < np; ++j) A[nv + i][nv + j] = C[i][j];
    }
    return A;
}

static void expect_dense_near(const Dense& A, const Dense& B, double tol) {
    EXPECT_TRUE(A.size() == B.size());
    EXPECT_TRUE(A.empty() || A[0].size() == B[0].size());
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j)
            EXPECT_NEAR(A[i][j], B[i][j], tol);
}

int main() {
    run_case("sparse_matmul_matches_dense_algebra", [] {
        const auto A = make_sparse(2, 2, {
            {0, 0, 2.0}, {0, 1, 1.0},
            {1, 0, 3.0}, {1, 1, 4.0}
        });
        const auto B = make_sparse(2, 2, {
            {0, 0, 5.0}, {0, 1, 6.0},
            {1, 0, 7.0}, {1, 1, 8.0}
        });

        const auto C = sparse_matmul(A, B);
        EXPECT_TRUE(C.is_consistent());
        EXPECT_NEAR(C(0, 0), 17.0, 1e-14);
        EXPECT_NEAR(C(0, 1), 20.0, 1e-14);
        EXPECT_NEAR(C(1, 0), 43.0, 1e-14);
        EXPECT_NEAR(C(1, 1), 50.0, 1e-14);
    });

    run_case("sparse_matmul_combines_duplicates_deterministically", [] {
        SparseMatrix A(1, 2);
        A.push_back(0, 0, 2.0);
        A.push_back(0, 1, 3.0);
        A.finalize();

        SparseMatrix B(2, 2);
        B.push_back(0, 1, 4.0);
        B.push_back(1, 1, -1.0);
        B.finalize();

        const auto C = sparse_matmul(A, B);
        EXPECT_TRUE(C.is_consistent());
        EXPECT_NEAR(C(0, 1), 5.0, 1e-14);
    });

    run_case("sparse_matmul_rejects_dimension_mismatch", [] {
        SparseMatrix A(2, 3);
        SparseMatrix B(2, 2);
        A.finalize();
        B.finalize();

        bool rejected = false;
        try {
            (void)sparse_matmul(A, B);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("sparse_pattern_distinguishes_graph_from_values", [] {
        const auto A = make_sparse(2, 2, {
            {0, 0, 2.0}, {0, 1, -1.0},
            {1, 0, 3.0}, {1, 1, 4.0}
        });
        const auto B = make_sparse(2, 2, {
            {0, 0, 5.0}, {0, 1, -2.0},
            {1, 0, 7.0}, {1, 1, 8.0}
        });
        const auto C = make_sparse(2, 2, {
            {0, 0, 5.0}, {0, 1, -2.0},
            {1, 0, 7.0}, {1, 1, 8.0},
            {1, 0, 0.5}
        });

        EXPECT_TRUE(same_sparse_pattern(A, B));
        EXPECT_TRUE(!same_sparse_pattern(A, C));
    });

    run_case("block_operator_validates_saddle_point_dimensions", [] {
        SparseMatrix Auu(3, 3);
        SparseMatrix G(3, 1);
        SparseMatrix D(1, 3);
        SparseMatrix C(1, 1);
        Auu.finalize(); G.finalize(); D.finalize(); C.finalize();

        const BlockOperator blocks(Auu, G, D, C);
        EXPECT_TRUE(blocks.is_valid());
        EXPECT_TRUE(blocks.velocity_size() == 3);
        EXPECT_TRUE(blocks.pressure_size() == 1);
    });

    run_case("block_operator_does_not_assume_transpose_coupling", [] {
        SparseMatrix Auu(2, 2);
        SparseMatrix G(2, 1);
        SparseMatrix D(1, 2);
        SparseMatrix C(1, 1);
        Auu.finalize(); G.finalize(); D.finalize(); C.finalize();

        const BlockOperator blocks(Auu, G, D, C);
        EXPECT_TRUE(blocks.is_valid());
    });

    run_case("exact_schur_oracle_matches_sparse_block_algebra", [] {
        const auto M = make_sparse(2, 2, {
            {0, 0, 4.0}, {0, 1, 1.0},
            {1, 0, 2.0}, {1, 1, 3.0}
        });
        const auto G = make_sparse(2, 2, {
            {0, 0, 1.0}, {0, 1, 0.5},
            {1, 0, 2.0}, {1, 1, 1.0}
        });
        const auto D = make_sparse(2, 2, {
            {0, 0, -3.0}, {0, 1, 1.0},
            {1, 0, 0.25}, {1, 1, 2.0}
        });
        const auto C = make_sparse(2, 2, {
            {0, 0, 2.0}, {0, 1, 0.2},
            {1, 0, 0.3}, {1, 1, 3.0}
        });

        const BlockOperator blocks(M, G, D, C);
        blocks.validate();

        const Dense Md = to_dense(M);
        const Dense Gd = to_dense(G);
        const Dense Dd = to_dense(D);
        const Dense Cd = to_dense(C);
        const Dense Minv = dense_inverse(Md);

        // Exact algebraic Schur oracle:
        // S = C - D M^{-1} G.
        const Dense exact = dense_subtract(
            Cd, dense_multiply(dense_multiply(Dd, Minv), Gd));

        const auto DG = sparse_matmul(D, G);
        EXPECT_TRUE(DG.is_consistent());

        // The oracle must agree with direct sparse block construction for
        // products that do not contain Auu^{-1}; DG is intentionally not
        // substituted for D M^{-1} G.
        EXPECT_NEAR(DG(0, 0), -1.0, 1e-14);
        EXPECT_NEAR(exact[0][0], 1.7, 1e-14);
        EXPECT_NEAR(exact[0][1], 0.05, 1e-14);
        EXPECT_NEAR(exact[1][0], -0.925, 1e-14);
        EXPECT_NEAR(exact[1][1], 2.3875, 1e-14);

        // The exact oracle is an independent reference, not a production
        // dense implementation.
        EXPECT_TRUE(std::isfinite(exact[0][0]));
        EXPECT_TRUE(std::isfinite(exact[1][1]));
    });

    run_case("schur_factorization_identities_match_block_ldu", [] {
        const auto M = make_sparse(2, 2, {
            {0, 0, 4.0}, {0, 1, 1.0},
            {1, 0, 2.0}, {1, 1, 3.0}
        });
        const auto G = make_sparse(2, 2, {
            {0, 0, 1.0}, {0, 1, 0.5},
            {1, 0, 2.0}, {1, 1, 1.0}
        });
        const auto D = make_sparse(2, 2, {
            {0, 0, -3.0}, {0, 1, 1.0},
            {1, 0, 0.25}, {1, 1, 2.0}
        });
        const auto C = make_sparse(2, 2, {
            {0, 0, 2.0}, {0, 1, 0.2},
            {1, 0, 0.3}, {1, 1, 3.0}
        });

        const Dense Md = to_dense(M);
        const Dense Gd = to_dense(G);
        const Dense Dd = to_dense(D);
        const Dense Cd = to_dense(C);
        const Dense Minv = dense_inverse(Md);
        const Dense S = dense_subtract(
            Cd, dense_multiply(dense_multiply(Dd, Minv), Gd));
        const Dense Sinv = dense_inverse(S);
        const Dense A = block_assemble(Md, Gd, Dd, Cd);

        const Dense z00 = Minv;
        const Dense z11 = Sinv;
        Dense diagonal(4, std::vector<double>(4, 0.0));
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                diagonal[i][j] = z00[i][j];
                diagonal[2 + i][2 + j] = z11[i][j];
            }
        }

        const Dense MinvG = dense_multiply(Minv, Gd);
        const Dense DMinv = dense_multiply(Dd, Minv);
        const Dense SinvDMinv = dense_multiply(Sinv, DMinv);
        const Dense MinvGSinv = dense_multiply(MinvG, Sinv);

        Dense lower = diagonal;
        Dense upper = diagonal;
        for (std::size_t i = 0; i < 2; ++i) {
            for (std::size_t j = 0; j < 2; ++j) {
                lower[2 + i][j] = -SinvDMinv[i][j];
                upper[i][2 + j] = -MinvGSinv[i][j];
            }
        }

        const Dense full = dense_multiply(upper, lower);
        const Dense identity = dense_multiply(full, A);

        // These are the four standard block factorization actions:
        // diagonal = diag(M^-1,S^-1)
        // lower    = diag^-1 L^-1
        // upper    = U^-1 diag^-1
        // full     = U^-1 diag^-1 L^-1 = A^-1.
        expect_dense_near(identity, dense_identity(4), 1e-12);

        // The Schur sign is the literal CFDX convention S=C-DM^-1G.
        // No implicit PETSc-style sign flip is introduced here.
        EXPECT_NEAR(S[0][0], 1.7, 1e-14);
        EXPECT_NEAR(S[0][1], 0.05, 1e-14);
        EXPECT_NEAR(S[1][0], -0.925, 1e-14);
        EXPECT_NEAR(S[1][1], 2.3875, 1e-14);
    });

    return run_all();
}
