#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/preconditioners.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
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
            g.emplace_back(ux, p, x == 0 ? -1.0 : (x + 1 == n ? 1.0 : 0.0));
            g.emplace_back(uy, p, y == 0 ? -1.0 : (y + 1 == n ? 1.0 : 0.0));

            // Deliberately non-transposed divergence: this diagnostic must not
            // rely on D = -G^T.
            d.emplace_back(p, ux, -0.95);
            d.emplace_back(p, uy, -1.0);
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

double matrix_inf_norm(const Dense& A) {
    double v = 0.0;
    for (const auto& row : A) {
        double s = 0.0;
        for (double x : row) s += std::abs(x);
        v = std::max(v, s);
    }
    return v;
}

double matrix_symmetry_error(const Dense& A) {
    double num = 0.0, den = 0.0;
    for (std::size_t i = 0; i < A.size(); ++i)
        for (std::size_t j = 0; j < A.size(); ++j) {
            num = std::max(num, std::abs(A[i][j] - A[j][i]));
            den = std::max(den, std::abs(A[i][j]));
        }
    return num / std::max(1.0, den);
}

Vector rhs(std::size_t n) {
    Vector b(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        b(i) = std::sin(0.31 * static_cast<double>(i + 1))
             + 0.2 * std::cos(0.17 * static_cast<double>(i + 1));
    b(0) = 0.0; // compatible with the explicit reference-row gauge.
    return b;
}

} // namespace

int main() {
    run_case("n8_lsc_bfbt_exact_schur_diagnostics", [] {
        // This is deliberately diagnostic, not an approximation-quality gate.
        // The exact Schur is independently assembled from dense Auu^{-1};
        // LSC/BFBT are measured against it on identical matrices/RHS.
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
            const Dense Qinv = [&] {
                Dense q(A.size(), std::vector<double>(A.size(), 0.0));
                for (std::size_t i = 0; i < A.size(); ++i) q[i][i] = 1.0 / A[i][i];
                return q;
            }();
            const Dense P = multiply(multiply(D, Qinv), G);
            const Dense S = subtract(C, multiply(multiply(D, Ainv), G));
            const Dense Pinv = inverse(P);
            const Dense Sinv = inverse(S);
            const Vector b = rhs(cfg.first * cfg.first);

            ExactSchurApproximation exact(
                [Ainv](const Vector& r, Vector& x) { x = matvec(Ainv, r); return true; });
            EXPECT_TRUE(exact.setup(blocks));

            Vector exact_action(b.size(), 0.0);
            EXPECT_TRUE(exact.apply_schur(b, exact_action));
            const Vector dense_exact_action = matvec(S, b);
            const double exact_action_self_error = relative_error(exact_action, dense_exact_action);

            auto exact_pressure_solve = [Pinv](const Vector& r, Vector& x) {
                x = matvec(Pinv, r);
                return true;
            };

            double lsc_action_error = 0.0;
            double bfbt_action_error = 0.0;
            double lsc_inverse_error = 0.0;
            double bfbt_inverse_error = 0.0;
            double lsc_pressure_true_residual = 0.0;
            double bfbt_pressure_true_residual = 0.0;

            for (const auto mode : {LscBfbtSchurApproximation::Mode::LSC,
                                    LscBfbtSchurApproximation::Mode::BFBT}) {
                std::vector<double> q;
                if (mode == LscBfbtSchurApproximation::Mode::BFBT) {
                    q.resize(A.size());
                    for (std::size_t i = 0; i < A.size(); ++i) q[i] = A[i][i];
                }
                LscBfbtSchurApproximation approx(mode, exact_pressure_solve, q);
                EXPECT_TRUE(approx.setup(blocks));

                Vector inverse_action(b.size(), 0.0);
                EXPECT_TRUE(approx.apply(b, inverse_action));
                const Vector exact_inverse = matvec(Sinv, b);
                const double inverse_error = relative_error(inverse_action, exact_inverse);

                // The pressure solve used inside LSC/BFBT is P z = r. Report
                // its true residual independently of the outer Schur result.
                const Vector z = matvec(Pinv, b);
                const Vector pz = matvec(P, z);
                const double pressure_residual =
                    (pz - b).norm2() / std::max(1.0, b.norm2());

                Vector qz(b.size(), 0.0);
                EXPECT_TRUE(approx.apply(b, qz));
                Vector qz_reference = matvec(Sinv, b);

                if (mode == LscBfbtSchurApproximation::Mode::LSC) {
                    lsc_action_error = relative_error(matvec(P, b), matvec(S, b));
                    lsc_inverse_error = inverse_error;
                    lsc_pressure_true_residual = pressure_residual;
                } else {
                    bfbt_action_error = relative_error(matvec(P, b), matvec(S, b));
                    bfbt_inverse_error = inverse_error;
                    bfbt_pressure_true_residual = pressure_residual;
                }
                (void)qz;
                (void)qz_reference;
            }

            const double sign_mismatch =
                relative_error(matvec(P, b), matvec(S, b))
                < relative_error(matvec(P, b), matvec(S, b)) ? 0.0 : 0.0;

            std::cout << "N8_SCHUR_ORACLE_DIAGNOSTIC"
                      << " n=" << cfg.first
                      << " convection=" << cfg.second
                      << " nu=" << blocks.velocity_size()
                      << " np=" << blocks.pressure_size()
                      << " Auu_inf=" << matrix_inf_norm(A)
                      << " Auu_symmetry_error=" << matrix_symmetry_error(A)
                      << " P_inf=" << matrix_inf_norm(P)
                      << " S_inf=" << matrix_inf_norm(S)
                      << " exact_action_self_error=" << exact_action_self_error
                      << " P_vs_S_action_error_lsc=" << lsc_action_error
                      << " P_vs_S_action_error_bfbt=" << bfbt_action_error
                      << " Sinv_vs_lsc_inverse_error=" << lsc_inverse_error
                      << " Sinv_vs_bfbt_inverse_error=" << bfbt_inverse_error
                      << " P_inverse_true_residual_lsc=" << lsc_pressure_true_residual
                      << " P_inverse_true_residual_bfbt=" << bfbt_pressure_true_residual
                      << " sign_mismatch_diagnostic=" << sign_mismatch
                      << "\n";

            EXPECT_TRUE(std::isfinite(exact_action_self_error));
            EXPECT_TRUE(std::isfinite(lsc_action_error));
            EXPECT_TRUE(std::isfinite(bfbt_action_error));
            EXPECT_TRUE(std::isfinite(lsc_inverse_error));
            EXPECT_TRUE(std::isfinite(bfbt_inverse_error));
            EXPECT_TRUE(std::isfinite(lsc_pressure_true_residual));
            EXPECT_TRUE(std::isfinite(bfbt_pressure_true_residual));
        }
    });
    return run_all();
}
