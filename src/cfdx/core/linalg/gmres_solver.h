#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "cg_solver.h"
#include <cstddef>
#include <cmath>
#include <vector>

namespace cfdx {
namespace core {

inline SolverResult solve_gmres(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    int restart = 30,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12)
{
    SolverResult result;

    if (b.size() != A.n_rows() || x.size() != A.n_cols()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = A.n_rows();
    const double* Av = A.values_data();
    const auto* Ac = A.columns_data();
    const auto* Ar = A.row_offsets_data();

    const double b_norm = b.norm2();
    const double tol_abs = tolerance * std::max(b_norm, 1e-15);

    // Helper: matrix-vector product
    auto matvec = [&](const std::vector<double>& v) -> std::vector<double> {
        std::vector<double> y(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) y[i] += Av[k] * v[Ac[k]];
        }
        return y;
    };

    // Compute initial residual: r = b - A*x
    std::vector<double> r(n);
    for (std::size_t i = 0; i < n; ++i) {
        r[i] = b(i);
        for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) r[i] -= Av[k] * x(Ac[k]);
    }

    double res = 0.0;
    for (std::size_t i = 0; i < n; ++i) res += r[i] * r[i];
    res = std::sqrt(res);

    if (res < tol_abs) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = res;
        result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
        return result;
    }

    int total_iter = 0;
    int m = (restart > 0) ? restart : 30;

    for (int restart_iter = 0; restart_iter < static_cast<int>(max_iter); ++restart_iter) {
        // Arnoldi: V_{m+1} and upper Hessenberg H_m
        std::vector<std::vector<double>> V(m + 1, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> H(m + 1, std::vector<double>(m + 1, 0.0));

        // v0 = r0 / ||r0||
        double beta = res;
        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i] / beta;

        // g for least squares (after Givens rotations)
        std::vector<double> g(m + 1, 0.0);
        g[0] = beta;

        int j = 0;
        bool converged = false;

        for (j = 1; j <= m; ++j) {
            // w = A * v_{j-1}
            std::vector<double> w = matvec(V[j - 1]);

            // Modified Gram-Schmidt with accumulated Givens rotations
            for (int i = 0; i < j; ++i) {
                double h_ij = 0.0;
                for (std::size_t k = 0; k < n; ++k) h_ij += V[i][k] * w[k];
                H[i][j] = h_ij;
                for (std::size_t k = 0; k < n; ++k) w[k] -= h_ij * V[i][k];
            }

            // h_{j,j} = ||w||
            double h_jj = 0.0;
            for (std::size_t k = 0; k < n; ++k) h_jj += w[k] * w[k];
            h_jj = std::sqrt(h_jj);
            H[j][j] = h_jj;

            if (h_jj < 1e-14 * beta) {
                // Happy breakdown
                break;
            }

            // v_j = w / h_{j,j}
            for (std::size_t k = 0; k < n; ++k) V[j][k] = w[k] / h_jj;

            ++total_iter;

            // Apply Givens rotation to eliminate H[j-1][j]
            double h_prev = H[j - 1][j];
            double h_curr = H[j][j];
            double norm = std::sqrt(h_prev * h_prev + h_curr * h_curr);

            if (norm > 1e-14) {
                double c = h_prev / norm;
                double s = h_curr / norm;

                // Rotate H
                H[j - 1][j] = norm;
                H[j][j] = 0.0;

                // Rotate g
                double g_prev = g[j - 1];
                g[j - 1] = c * g_prev + s * g[j];
                g[j] = -s * g_prev + c * g[j];
            }

            // Current residual estimate
            res = std::abs(g[j]);

            if (res < tol_abs) {
                converged = true;
                break;
            }
        }

        // Back substitution to solve H y = g
        int jj_end = (converged) ? j - 1 : m;
        std::vector<double> y(m + 1, 0.0);

        if (jj_end >= 0 && jj_end <= m && std::abs(H[jj_end][jj_end]) > 1e-14) {
            y[jj_end] = g[jj_end] / H[jj_end][jj_end];
            for (int jj = jj_end - 1; jj >= 0; --jj) {
                double sum = g[jj];
                for (int kk = jj + 1; kk <= jj_end; ++kk) sum -= H[jj][kk] * y[kk];
                y[jj] = sum / H[jj][jj];
            }
        }

        // x = x + V * y
        for (std::size_t i = 0; i < n; ++i) {
            for (int jj = 0; jj <= jj_end; ++jj) {
                x(i) += V[jj][i] * y[jj];
            }
        }

        // Compute new residual
        for (std::size_t i = 0; i < n; ++i) {
            r[i] = b(i);
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) r[i] -= Av[k] * x(Ac[k]);
        }

        res = 0.0;
        for (std::size_t i = 0; i < n; ++i) res += r[i] * r[i];
        res = std::sqrt(res);

        if (res < tol_abs) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = total_iter;
            result.residual = res;
            result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
            return result;
        }

        if (static_cast<std::size_t>(total_iter) >= max_iter) break;
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = total_iter;
    result.residual = res;
    result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
    return result;
}

}  // namespace core
}  // namespace cfdx