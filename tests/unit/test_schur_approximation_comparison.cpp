#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
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
    const std::initializer_list<std::tuple<std::size_t, std::size_t, double>>& entries) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

static Dense dense(const SparseMatrix& A) {
    Dense out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i]; k < A.row_offsets_data()[i + 1]; ++k)
            out[i][A.columns_data()[k]] += A.values_data()[k];
    return out;
}

static Dense inverse(Dense A) {
    const std::size_t n = A.size();
    Dense I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        EXPECT_TRUE(std::abs(A[pivot][k]) > 1e-14);
        std::swap(A[k], A[pivot]);
        std::swap(I[k], I[pivot]);
        const double d = A[k][k];
        for (std::size_t j = 0; j < n; ++j) {
            A[k][j] /= d;
            I[k][j] /= d;
        }
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

static Dense multiply(const Dense& A, const Dense& B) {
    Dense C(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t k = 0; k < B.size(); ++k)
            for (std::size_t j = 0; j < B[0].size(); ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

static Dense subtract(const Dense& A, const Dense& B) {
    Dense C = A;
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) C[i][j] -= B[i][j];
    return C;
}

static double inf_norm(const Dense& A) {
    double m = 0.0;
    for (const auto& row : A) {
        double s = 0.0;
        for (double v : row) s += std::abs(v);
        m = std::max(m, s);
    }
    return m;
}

static Vector matvec(const Dense& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

static double relative_error(const Vector& a, const Vector& b) {
    const double scale = std::max(1.0, b.norm2());
    return (a - b).norm2() / scale;
}

int main() {
    // This is a comparison harness, not a production qualification gate.
    // The family deliberately varies Auu coupling/conditioning while keeping
    // G and D fixed. The exact Schur is the reference; LSC/BFBt are measured
    // against it without declaring a method acceptable from one threshold.
    const auto G = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.5},
        {1, 0, 0.25}, {1, 1, 1.0}
    });
    const auto D = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.25},
        {1, 0, 0.5}, {1, 1, 1.0}
    });
    const auto C = make_sparse(2, 2, {
        {0, 0, 0.5}, {1, 1, 0.75}
    });

    const std::vector<SparseMatrix> velocity_matrices{
        make_sparse(2, 2, {{0, 0, 4.0}, {0, 1, 0.25}, {1, 0, 0.25}, {1, 1, 3.0}}),
        make_sparse(2, 2, {{0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 3.0}}),
        make_sparse(2, 2, {{0, 0, 4.0}, {0, 1, 1.75}, {1, 0, 1.75}, {1, 1, 3.0}})
    };

    const Vector rhs = [] {
        Vector b(2, 0.0);
        b(0) = 1.0;
        b(1) = -0.35;
        return b;
    }();

    run_case("exact_schur_is_reference_for_conditioning_family", [&] {
        for (const auto& Auu : velocity_matrices) {
            const Dense Ad = dense(Auu);
            const Dense Ainv = inverse(Ad);
            const double cond_inf = inf_norm(Ad) * inf_norm(Ainv);
            EXPECT_TRUE(std::isfinite(cond_inf) && cond_inf > 1.0);

            const BlockOperator blocks(Auu, G, D, C);
            blocks.validate();

            auto auu_solve = [Ainv](const Vector& r, Vector& y) {
                y = matvec(Ainv, r);
                return true;
            };
            ExactSchurApproximation exact(auu_solve);
            EXPECT_TRUE(exact.setup(blocks));

            Vector exact_action(2, 0.0);
            EXPECT_TRUE(exact.apply_schur(rhs, exact_action));

            // Independently assemble S = C - D Auu^-1 G and check the
            // implicit exact operator. This guards the comparison reference
            // itself before any approximate method is assessed.
            const Dense S = subtract(
                dense(C),
                multiply(multiply(dense(D), Ainv), dense(G)));
            const Vector assembled_action = matvec(S, rhs);
            EXPECT_TRUE(relative_error(exact_action, assembled_action) < 1e-12);
        }
    });

    run_case("lsc_and_bfbt_are_compared_to_exact_reference", [&] {
        for (const auto& Auu : velocity_matrices) {
            const Dense Ad = dense(Auu);
            const Dense Ainv = inverse(Ad);
            const BlockOperator blocks(Auu, G, D, C);
            blocks.validate();

            const auto run = [&](LscBfbtSchurApproximation::Mode mode,
                                 const std::vector<double>& qdiag) {
                std::vector<double> qinv(2, 1.0);
                if (qdiag.empty()) {
                    qinv[0] = 1.0 / Ad[0][0];
                    qinv[1] = 1.0 / Ad[1][1];
                } else {
                    qinv[0] = 1.0 / qdiag[0];
                    qinv[1] = 1.0 / qdiag[1];
                }

                Dense P = multiply(multiply(dense(D),
                                             {{qinv[0], 0.0}, {0.0, qinv[1]}}),
                                   dense(G));
                const Dense Qinv{{qinv[0], 0.0}, {0.0, qinv[1]}};
                const Dense E = multiply(
                    multiply(multiply(multiply(dense(D), Qinv), Ad), Qinv),
                    dense(G));
                const Dense H = subtract(P, dense(C));
                const Dense Hinv = inverse(H);

                auto pressure_solve = [Hinv](const Vector& r, Vector& z) {
                    z = matvec(Hinv, r);
                    return true;
                };
                LscBfbtSchurApproximation approx(mode, pressure_solve, qdiag);
                EXPECT_TRUE(approx.setup(blocks));

                Vector estimated(2, 0.0);
                EXPECT_TRUE(approx.apply(rhs, estimated));

                const Dense K = subtract(E, dense(C));
                const Dense KHinv = multiply(K, Hinv);
                const Dense HinvKHinv = multiply(Hinv, KHinv);
                Vector reference(2, 0.0);
                for (std::size_t i = 0; i < 2; ++i)
                    for (std::size_t j = 0; j < 2; ++j)
                        reference(i) += HinvKHinv[i][j] * rhs(j);
                reference *= -1.0;

                // Verify implementation against its documented -H^-1 K H^-1
                // algebra before comparing that operator to the exact Schur.
                EXPECT_TRUE(relative_error(estimated, reference) < 1e-12);

                const Dense S = subtract(
                    dense(C), multiply(multiply(dense(D), Ainv), dense(G)));
                const Dense Sinv = inverse(S);
                const Vector exact_solution = matvec(Sinv, rhs);
                const double comparison_error =
                    relative_error(estimated, exact_solution);
                EXPECT_TRUE(std::isfinite(comparison_error));

                // The value is intentionally not turned into a qualification
                // threshold here. Recordability and conditioning-awareness are
                // the purpose of this harness; production gates belong to the
                // benchmark/qualification PR.
                (void)comparison_error;
                (void)mode;
            };

            run(LscBfbtSchurApproximation::Mode::LSC, {});
            run(LscBfbtSchurApproximation::Mode::BFBT, {2.0, 3.0});
        }
    });

    run_case("comparison_does_not_assume_transpose_coupling", [&] {
        const auto nonsymmetric_D = make_sparse(2, 2, {
            {0, 0, 1.0}, {0, 1, 0.1},
            {1, 0, 0.7}, {1, 1, 1.0}
        });
        const auto Auu = velocity_matrices[1];
        const BlockOperator blocks(Auu, G, nonsymmetric_D, C);
        blocks.validate();

        const auto solve_identity = [](const Vector& r, Vector& y) {
            y = r;
            return true;
        };
        LscBfbtSchurApproximation approx(
            LscBfbtSchurApproximation::Mode::BFBT, solve_identity, {2.0, 3.0});
        EXPECT_TRUE(approx.setup(blocks));
        Vector out(2, 0.0);
        EXPECT_TRUE(approx.apply(rhs, out));
        EXPECT_TRUE(std::isfinite(out(0)) && std::isfinite(out(1)));
    });

    return run_all();
}
