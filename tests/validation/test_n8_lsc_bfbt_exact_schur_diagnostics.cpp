#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

using Dense = std::vector<std::vector<double>>;

struct Blocks {
    SparseMatrix Auu;
    SparseMatrix G;
    SparseMatrix D;
    SparseMatrix C;
};

SparseMatrix make_sparse(
    std::size_t rows, std::size_t cols,
    const std::vector<std::tuple<std::size_t, std::size_t, double>>& entries) {
    SparseMatrix A(rows, cols);
    for (const auto& [i, j, value] : entries) A.push_back(i, j, value);
    A.finalize();
    return A;
}

Blocks make_case(std::size_t n, double convection) {
    const std::size_t np = n * n;
    const std::size_t nu = 2 * np;
    std::vector<std::tuple<std::size_t, std::size_t, double>> auu, g, d, c;
    const auto cell = [n](std::size_t x, std::size_t y) { return y * n + x; };

    for (std::size_t y = 0; y < n; ++y) {
        for (std::size_t x = 0; x < n; ++x) {
            const std::size_t p = cell(x, y);
            for (std::size_t comp = 0; comp < 2; ++comp) {
                const std::size_t row = comp * np + p;
                auu.emplace_back(row, row, 5.0 + convection);
                if (x > 0) auu.emplace_back(row, comp * np + cell(x - 1, y), -1.0 - convection);
                if (x + 1 < n) auu.emplace_back(row, comp * np + cell(x + 1, y), -1.0);
                if (y > 0) auu.emplace_back(row, comp * np + cell(x, y - 1), -1.0);
                if (y + 1 < n) auu.emplace_back(row, comp * np + cell(x, y + 1), -1.0);
            }

            const std::size_t ux = p;
            const std::size_t uy = np + p;
            if (x > 0) g.emplace_back(ux, cell(x - 1, y), -0.5);
            if (x + 1 < n) g.emplace_back(ux, cell(x + 1, y), 0.5);
            if (y > 0) g.emplace_back(uy, cell(x, y - 1), -0.5);
            if (y + 1 < n) g.emplace_back(uy, cell(x, y + 1), 0.5);
            if (x == 0) g.emplace_back(ux, p, -1.0);
            else if (x + 1 == n) g.emplace_back(ux, p, 1.0);
            if (y == 0) g.emplace_back(uy, p, -1.0);
            else if (y + 1 == n) g.emplace_back(uy, p, 1.0);

            // Deliberately non-transposed D: the oracle must not rely on
            // D = -G^T.
            d.emplace_back(p, ux, -0.95);
            d.emplace_back(p, uy, -1.0);

            // Pin one pressure row so the diagnostic inverse is nonsingular.
            c.emplace_back(p, p, p == 0 ? 1.0 : 0.02);
        }
    }
    return {make_sparse(nu, nu, auu), make_sparse(nu, np, g),
            make_sparse(np, nu, d), make_sparse(np, np, c)};
}

Dense dense(const SparseMatrix& A) {
    Dense out(A.n_rows(), std::vector<double>(A.n_cols(), 0.0));
    for (std::size_t i = 0; i < A.n_rows(); ++i)
        for (std::size_t k = A.row_offsets_data()[i]; k < A.row_offsets_data()[i + 1]; ++k)
            out[i][A.columns_data()[k]] += A.values_data()[k];
    return out;
}

Dense inverse(Dense A) {
    const std::size_t n = A.size();
    Dense I(n, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) I[i][i] = 1.0;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) pivot = i;
        EXPECT_TRUE(std::abs(A[pivot][k]) > 1e-12);
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

Dense multiply(const Dense& A, const Dense& B) {
    Dense C(A.size(), std::vector<double>(B[0].size(), 0.0));
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t k = 0; k < B.size(); ++k)
            for (std::size_t j = 0; j < B[0].size(); ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

Dense subtract(const Dense& A, const Dense& B) {
    Dense C = A;
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) C[i][j] -= B[i][j];
    return C;
}

Vector matvec(const Dense& A, const Vector& x) {
    Vector y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A[i].size(); ++j) y(i) += A[i][j] * x(j);
    return y;
}

double relative_error(const Vector& a, const Vector& b) {
    return (a - b).norm2() / std::max(1.0, b.norm2());
}

double inf_norm(const Dense& A) {
    double result = 0.0;
    for (const auto& row : A) {
        double sum = 0.0;
        for (double x : row) sum += std::abs(x);
        result = std::max(result, sum);
    }
    return result;
}

double symmetry_error(const Dense& A) {
    double num = 0.0;
    double den = 1.0;
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A.size(); ++j) {
            num = std::max(num, std::abs(A[i][j] - A[j][i]));
            den = std::max(den, std::abs(A[i][j]));
        }
    return num / den;
}

Vector make_rhs(std::size_t n) {
    Vector b(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        b(i) = std::sin(0.31 * static_cast<double>(i + 1))
             + 0.2 * std::cos(0.17 * static_cast<double>(i + 1));
    b(0) = 0.0;
    return b;
}

} // namespace

int main() {
    run_case("n8_lsc_bfbt_exact_schur_diagnostics", [] {
        // Diagnostic-only campaign. No approximation-quality threshold is
        // invented here: the Exact Schur is the reference and all measured
        // discrepancies are emitted for later N8 qualification.
        for (const auto& cfg : std::vector<std::pair<std::size_t, double>>{
                 {4, 0.0}, {4, 1.0}, {4, 4.0}, {6, 0.0}, {6, 1.0}, {6, 4.0}}) {
            const auto data = make_case(cfg.first, cfg.second);
            const BlockOperator blocks(data.Auu, data.G, data.D, data.C);
            blocks.validate();

            const Dense A = dense(blocks.Auu());
            const Dense G = dense(blocks.G());
            const Dense D = dense(blocks.D());
            const Dense C = dense(blocks.C());
            const Dense Ainv = inverse(A);
            const Dense Qinv_lsc = [&] {
                Dense q(A.size(), std::vector<double>(A.size(), 0.0));
                for (std::size_t i = 0; i < A.size(); ++i) q[i][i] = 1.0 / A[i][i];
                return q;
            }();
            const Dense Iinv = [&] {
                Dense q(A.size(), std::vector<double>(A.size(), 0.0));
                for (std::size_t i = 0; i < A.size(); ++i) q[i][i] = 1.0;
                return q;
            }();

            const Dense P_lsc = multiply(multiply(D, Qinv_lsc), G);
            const Dense P_bfbt = multiply(multiply(D, Iinv), G);
            const Dense S = subtract(C, multiply(multiply(D, Ainv), G));
            const Dense P_lsc_inv = inverse(P_lsc);
            const Dense P_bfbt_inv = inverse(P_bfbt);
            const Dense S_inv = inverse(S);
            const Vector b = make_rhs(cfg.first * cfg.first);

            ExactSchurApproximation exact(
                [Ainv](const Vector& r, Vector& x) {
                    x = matvec(Ainv, r);
                    return true;
                });
            EXPECT_TRUE(exact.setup(blocks));

            Vector exact_action(b.size(), 0.0);
            EXPECT_TRUE(exact.apply_schur(b, exact_action));
            const Vector dense_exact_action = matvec(S, b);
            const double exact_action_self_error =
                relative_error(exact_action, dense_exact_action);

            auto make_pressure_solve = [](const Dense& Pinv) {
                return [Pinv](const Vector& r, Vector& x) {
                    x = matvec(Pinv, r);
                    return true;
                };
            };

            double action_error_lsc = 0.0;
            double action_error_bfbt = 0.0;
            double inverse_error_lsc = 0.0;
            double inverse_error_bfbt = 0.0;
            double pressure_residual_lsc = 0.0;
            double pressure_residual_bfbt = 0.0;

            {
                LscBfbtSchurApproximation lsc(
                    LscBfbtSchurApproximation::Mode::LSC,
                    make_pressure_solve(P_lsc_inv));
                LscBfbtSchurApproximation bfbt(
                    LscBfbtSchurApproximation::Mode::BFBT,
                    make_pressure_solve(P_bfbt_inv));
                EXPECT_TRUE(lsc.setup(blocks));
                EXPECT_TRUE(bfbt.setup(blocks));

                Vector lsc_inverse(b.size(), 0.0);
                Vector bfbt_inverse(b.size(), 0.0);
                EXPECT_TRUE(lsc.apply(b, lsc_inverse));
                EXPECT_TRUE(bfbt.apply(b, bfbt_inverse));

                const Vector exact_inverse = matvec(S_inv, b);
                inverse_error_lsc = relative_error(lsc_inverse, exact_inverse);
                inverse_error_bfbt = relative_error(bfbt_inverse, exact_inverse);
                action_error_lsc = relative_error(matvec(P_lsc, b), matvec(S, b));
                action_error_bfbt = relative_error(matvec(P_bfbt, b), matvec(S, b));

                const Vector lsc_z = matvec(P_lsc_inv, b);
                const Vector bfbt_z = matvec(P_bfbt_inv, b);
                pressure_residual_lsc =
                    (matvec(P_lsc, lsc_z) - b).norm2() / std::max(1.0, b.norm2());
                pressure_residual_bfbt =
                    (matvec(P_bfbt, bfbt_z) - b).norm2() / std::max(1.0, b.norm2());
            }

            // Sign/scaling diagnostics: report the best direct alignment of P
            // with S and -S. This distinguishes a sign convention error from
            // a genuinely poor approximation without silently accepting either.
            const double direct = action_error_lsc;
            const double opposite_sign = relative_error(
                matvec(P_lsc, b), matvec(S, b) * -1.0);
            const double lsc_bfbt_separation =
                relative_error(matvec(P_lsc, b), matvec(P_bfbt, b));

            std::cout << "N8_SCHUR_ORACLE_DIAGNOSTIC"
                      << " n=" << cfg.first
                      << " convection=" << cfg.second
                      << " nu=" << blocks.velocity_size()
                      << " np=" << blocks.pressure_size()
                      << " Auu_inf=" << inf_norm(A)
                      << " Auu_symmetry_error=" << symmetry_error(A)
                      << " P_lsc_inf=" << inf_norm(P_lsc)
                      << " P_bfbt_inf=" << inf_norm(P_bfbt)
                      << " S_inf=" << inf_norm(S)
                      << " exact_action_self_error=" << exact_action_self_error
                      << " P_lsc_vs_S_action_error=" << direct
                      << " P_lsc_vs_minusS_action_error=" << opposite_sign
                      << " P_bfbt_vs_S_action_error=" << action_error_bfbt
                      << " Sinv_vs_lsc_inverse_error=" << inverse_error_lsc
                      << " Sinv_vs_bfbt_inverse_error=" << inverse_error_bfbt
                      << " P_inverse_true_residual_lsc=" << pressure_residual_lsc
                      << " P_inverse_true_residual_bfbt=" << pressure_residual_bfbt
                      << " LSC_vs_BFBT_action_separation=" << lsc_bfbt_separation
                      << "\n";

            EXPECT_TRUE(std::isfinite(exact_action_self_error));
            EXPECT_TRUE(std::isfinite(direct));
            EXPECT_TRUE(std::isfinite(opposite_sign));
            EXPECT_TRUE(std::isfinite(action_error_bfbt));
            EXPECT_TRUE(std::isfinite(inverse_error_lsc));
            EXPECT_TRUE(std::isfinite(inverse_error_bfbt));
            EXPECT_TRUE(std::isfinite(pressure_residual_lsc));
            EXPECT_TRUE(std::isfinite(pressure_residual_bfbt));
            EXPECT_TRUE(std::isfinite(lsc_bfbt_separation));
        }
    });
    return run_all();
}
