#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <limits>
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
        if (!(std::abs(A[pivot][k]) > 0.0) || !std::isfinite(A[pivot][k])) return {};
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
    Dense C(A.size(), std::vector<double>(B.front().size(), 0.0));
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t k = 0; k < B.size(); ++k)
            for (std::size_t j = 0; j < B.front().size(); ++j)
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
    double result = 0.0;
    for (const auto& row : A) {
        double sum = 0.0;
        for (double value : row) sum += std::abs(value);
        result = std::max(result, sum);
    }
    return result;
}

static Vector matvec(const Dense& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

static double relative_error(const Vector& a, const Vector& b) {
    return (a - b).norm2() / std::max(1.0, b.norm2());
}

static double backward_error(const Dense& A, const Vector& x, const Vector& b) {
    const Vector residual = matvec(A, x) - b;
    const double denominator = std::max(
        std::numeric_limits<double>::min(), inf_norm(A) * x.norm2() + b.norm2());
    return residual.norm2() / denominator;
}

int main() {
    // N8 quantitative qualification harness:
    // condition numbers and measured FP64 floors are evidence; approximate
    // Schur quality remains diagnostic until representative CFD matrices have
    // established an acceptance envelope.
    const auto G = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.5}, {1, 0, 0.25}, {1, 1, 1.0}});
    const auto D = make_sparse(2, 2, {
        {0, 0, 1.0}, {0, 1, 0.25}, {1, 0, 0.5}, {1, 1, 1.0}});
    const auto C = make_sparse(2, 2, {{0, 0, 0.5}, {1, 1, 0.75}});
    const std::vector<SparseMatrix> family{
        make_sparse(2, 2, {{0, 0, 4.0}, {0, 1, 0.25}, {1, 0, 0.25}, {1, 1, 3.0}}),
        make_sparse(2, 2, {{0, 0, 4.0}, {0, 1, 1.0}, {1, 0, 1.0}, {1, 1, 3.0}}),
        make_sparse(2, 2, {{0, 0, 4.0}, {0, 1, 1.75}, {1, 0, 1.75}, {1, 1, 3.0}})};

    const Vector rhs = [] {
        Vector b(2, 0.0);
        b(0) = 1.0;
        b(1) = -0.35;
        return b;
    }();

    run_case("schur_conditioning_and_fp64_floor", [&] {
        for (std::size_t case_id = 0; case_id < family.size(); ++case_id) {
            const Dense Ad = dense(family[case_id]);
            const Dense Ainv = inverse(Ad);
            EXPECT_TRUE(!Ainv.empty());

            const double cond_inf = inf_norm(Ad) * inf_norm(Ainv);
            EXPECT_TRUE(std::isfinite(cond_inf) && cond_inf >= 1.0);

            const BlockOperator blocks(family[case_id], G, D, C);
            blocks.validate();

            auto solve_auu = [Ainv](const Vector& r, Vector& x) {
                x = matvec(Ainv, r);
                return true;
            };
            ExactSchurApproximation exact(solve_auu);
            EXPECT_TRUE(exact.setup(blocks));

            Vector exact_action(rhs.size(), 0.0);
            EXPECT_TRUE(exact.apply_schur(rhs, exact_action));

            const Dense S = subtract(
                dense(C), multiply(multiply(dense(D), Ainv), dense(G)));
            const Dense Sinv = inverse(S);
            EXPECT_TRUE(!Sinv.empty());

            const Vector exact_solution = matvec(Sinv, rhs);
            const double oracle_discrepancy =
                relative_error(exact_action, matvec(S, exact_solution));
            const double floor = backward_error(S, exact_solution, rhs);

            EXPECT_TRUE(std::isfinite(oracle_discrepancy));
            EXPECT_TRUE(std::isfinite(floor));

            std::cout << "N8_SCHUR case=" << case_id
                      << " cond_inf_Auu=" << cond_inf
                      << " exact_oracle_discrepancy=" << oracle_discrepancy
                      << " exact_solve_backward_error=" << floor
                      << " machine_epsilon="
                      << std::numeric_limits<double>::epsilon() << '\n';

            for (const auto mode : {LscBfbtSchurApproximation::Mode::LSC,
                                    LscBfbtSchurApproximation::Mode::BFBT}) {
                const std::vector<double> qdiag =
                    mode == LscBfbtSchurApproximation::Mode::LSC
                        ? std::vector<double>{}
                        : std::vector<double>{2.0, 3.0};
                const std::vector<double> qinv =
                    mode == LscBfbtSchurApproximation::Mode::LSC
                        ? std::vector<double>{1.0 / Ad[0][0], 1.0 / Ad[1][1]}
                        : std::vector<double>{1.0 / qdiag[0], 1.0 / qdiag[1]};

                Dense Qinv{{qinv[0], 0.0}, {0.0, qinv[1]}};
                const Dense P = multiply(multiply(dense(D), Qinv), dense(G));
                const Dense E = multiply(
                    multiply(multiply(multiply(dense(D), Qinv), Ad), Qinv), dense(G));
                const Dense Pinv = inverse(P);
                EXPECT_TRUE(!Pinv.empty());

                auto pressure_solve = [Pinv](const Vector& r, Vector& x) {
                    x = matvec(Pinv, r);
                    return true;
                };
                LscBfbtSchurApproximation approx(mode, pressure_solve, qdiag);
                EXPECT_TRUE(approx.setup(blocks));

                Vector estimated(rhs.size(), 0.0);
                EXPECT_TRUE(approx.apply(rhs, estimated));

                const Vector documented = matvec(Pinv, matvec(E, matvec(Pinv, rhs))) * -1.0;
                const double algebra_error = relative_error(estimated, documented);
                const double exact_action_error = relative_error(estimated, exact_action);

                EXPECT_TRUE(std::isfinite(algebra_error));
                EXPECT_TRUE(std::isfinite(exact_action_error));

                std::cout << "N8_SCHUR case=" << case_id
                          << " method="
                          << (mode == LscBfbtSchurApproximation::Mode::LSC ? "LSC" : "BFBT")
                          << " algebra_error=" << algebra_error
                          << " exact_schur_error=" << exact_action_error << '\n';
            }
        }
    });

    return run_all();
}
