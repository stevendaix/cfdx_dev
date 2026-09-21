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
    const double tol = tolerance * std::max(b_norm, 1.0);

    std::vector<double> r(n);
    for (std::size_t i = 0; i < n; ++i) {
        r[i] = b(i);
        for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) r[i] -= Av[k] * x(Ac[k]);
    }

    double res = 0.0;
    for (std::size_t i = 0; i < n; ++i) res += r[i] * r[i];
    res = std::sqrt(res);

    if (res < tol) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = res;
        result.residual_relative = res / std::max(b_norm, 1.0);
        return result;
    }

    int iter = 0;
    int m = restart > 0 ? restart : 30;
    if (static_cast<std::size_t>(m) > n) m = static_cast<int>(n);

    for (int restart_iter = 0; iter < static_cast<int>(max_iter); ++restart_iter) {
        int restart_limit = static_cast<int>(std::min(
            static_cast<std::size_t>(m),
            static_cast<std::size_t>(max_iter - iter)));

        if (restart_limit <= 0) break;

        std::vector<std::vector<double>> V(m + 1, std::vector<double>(n, 0.0));
        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i] / res;

        std::vector<std::vector<double>> H(m + 1, std::vector<double>(m + 1, 0.0));

        for (int j = 1; j <= restart_limit; ++j) {
            std::vector<double> w(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                    w[i] += Av[k] * V[j - 1][Ac[k]];
                }
            }

            for (int i = 0; i < j; ++i) {
                double h = 0.0;
                for (std::size_t k = 0; k < n; ++k) h += V[i][k] * w[k];
                H[i][j] = h;
                for (std::size_t k = 0; k < n; ++k) w[k] -= h * V[i][k];
            }

            double h_jj = 0.0;
            for (std::size_t k = 0; k < n; ++k) h_jj += w[k] * w[k];
            h_jj = std::sqrt(h_jj);
            H[j][j] = h_jj;

            if (h_jj < 1e-14 * res) {
                break;
            }

            for (std::size_t k = 0; k < n; ++k) V[j][k] = w[k] / h_jj;

            ++iter;
        }

        std::vector<double> g(restart_limit + 2, 0.0);
        g[0] = res;

        for (int j = 1; j <= restart_limit; ++j) {
            double h1 = H[j - 1][j];
            double h2 = H[j][j];
            double norm = std::sqrt(h1 * h1 + h2 * h2);

            if (norm > 1e-14) {
                double c = h1 / norm;
                double s = h2 / norm;

                H[j - 1][j] = norm;
                H[j][j] = 0.0;

                double g1 = g[j - 1];
                double g2 = g[j];
                g[j - 1] = c * g1 + s * g2;
                g[j] = -s * g1 + c * g2;
            }
        }

        res = std::abs(g[restart_limit]);

        if (res < tol) {
            std::vector<double> y(restart_limit + 1, 0.0);
            if (std::abs(H[restart_limit][restart_limit]) > 1e-14) {
                y[restart_limit] = g[restart_limit] / H[restart_limit][restart_limit];
                for (int i = restart_limit - 1; i >= 0; --i) {
                    double sum = g[i];
                    for (int j = i + 1; j <= restart_limit; ++j) sum -= H[i][j] * y[j];
                    y[i] = sum / H[i][i];
                }
            }

            for (std::size_t i = 0; i < n; ++i) {
                for (int j = 0; j < restart_limit; ++j) {
                    x(i) += V[j][i] * y[j];
                }
            }

            result.status = SolverStatus::CONVERGED;
            result.iterations = iter;
            result.residual = res;
            result.residual_relative = res / std::max(b_norm, 1.0);
            return result;
        }

        std::vector<double> y(restart_limit + 1, 0.0);
        if (std::abs(H[restart_limit][restart_limit]) > 1e-14) {
            y[restart_limit] = g[restart_limit] / H[restart_limit][restart_limit];
            for (int i = restart_limit - 1; i >= 0; --i) {
                double sum = g[i];
                for (int j = i + 1; j <= restart_limit; ++j) sum -= H[i][j] * y[j];
                if (std::abs(H[i][i]) > 1e-14) y[i] = sum / H[i][i];
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            for (int j = 0; j < restart_limit; ++j) {
                x(i) += V[j][i] * y[j];
            }
        }

        for (std::size_t i = 0; i < n; ++i) {
            r[i] = b(i);
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) r[i] -= Av[k] * x(Ac[k]);
        }

        res = 0.0;
        for (std::size_t i = 0; i < n; ++i) res += r[i] * r[i];
        res = std::sqrt(res);

        if (res < tol) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iter;
            result.residual = res;
            result.residual_relative = res / std::max(b_norm, 1.0);
            return result;
        }
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = iter;
    result.residual = res;
    result.residual_relative = res / std::max(b_norm, 1.0);
    return result;
}

}  // namespace core
}  // namespace cfdx