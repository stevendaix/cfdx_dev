#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace cfdx::core {

inline SolverResult solve_gmres(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    int restart = 30,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    const Preconditioner* preconditioner = nullptr)
{
    SolverResult result;
    if (A.n_rows() != A.n_cols() || b.size() != A.n_rows() || x.size() != A.n_cols() ||
        restart <= 0 || max_iter == 0 || tolerance <= 0.0) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = A.n_rows();
    const auto* values = A.values_data();
    const auto* columns = A.columns_data();
    const auto* rows = A.row_offsets_data();
    const double b_norm = b.norm2();
    const double tol = tolerance * std::max(b_norm, 1.0);
    const int m = std::min<int>(restart, static_cast<int>(n));

        auto matvec = [&](const std::vector<double>& in, std::vector<double>& out) {
        out.assign(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = rows[i]; k < rows[i + 1]; ++k)
                out[i] += values[k] * in[columns[k]];
    };

    auto residual = [&](std::vector<double>& r) {
        r.assign(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            r[i] = b(i);
            for (std::size_t k = rows[i]; k < rows[i + 1]; ++k)
                r[i] -= values[k] * x(columns[k]);
        }
    };

    std::vector<double> r, w(n), z(n);
    residual(r);
    auto norm = [](const std::vector<double>& v) {
        double s = 0.0;
        for (double a : v) s += a * a;
        return std::sqrt(s);
    };

    double beta = norm(r);
    if (beta <= tol) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = beta;
        result.residual_relative = beta / std::max(b_norm, 1.0);
        return result;
    }

    std::vector<std::vector<double>> V(static_cast<std::size_t>(m + 1), std::vector<double>(n));
    std::vector<std::vector<double>> Z(static_cast<std::size_t>(m), std::vector<double>(n));
    std::vector<std::vector<double>> H(static_cast<std::size_t>(m + 1), std::vector<double>(m, 0.0));
    std::vector<double> cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0), y(m, 0.0);

    std::size_t iterations = 0;
    while (iterations < max_iter) {
        residual(r);
        beta = norm(r);
        if (beta <= tol) break;

        for (std::size_t i = 0; i < n; ++i) V[0][i] = r[i] / beta;
        for (auto& row : H) std::fill(row.begin(), row.end(), 0.0);
        std::fill(cs.begin(), cs.end(), 0.0);
        std::fill(sn.begin(), sn.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;

        int used = 0;
        double estimated_residual = beta;

        for (int j = 0; j < m && iterations < max_iter; ++j, ++iterations) {
            if (preconditioner) {
                for (std::size_t i = 0; i < n; ++i) z[i] = V[j][i];
                Vector rv(n), zv(n);
                for (std::size_t i = 0; i < n; ++i) rv(i) = z[i];
                if (!preconditioner->apply(rv, zv)) {
                    result.status = SolverStatus::NOT_APPLICABLE;
                    result.iterations = iterations;
                    return result;
                }
                for (std::size_t i = 0; i < n; ++i) Z[j][i] = zv(i);
            } else {
                Z[j] = V[j];
            }

            matvec(Z[j], w);

            for (int i = 0; i <= j; ++i) {
                double h = 0.0;
                for (std::size_t k = 0; k < n; ++k) h += V[i][k] * w[k];
                H[i][j] = h;
                for (std::size_t k = 0; k < n; ++k) w[k] -= h * V[i][k];
            }

            H[j + 1][j] = norm(w);
            if (H[j + 1][j] > 0.0) {
                for (std::size_t k = 0; k < n; ++k) V[j + 1][k] = w[k] / H[j + 1][j];
            }

            for (int i = 0; i < j; ++i) {
                const double h0 = H[i][j];
                const double h1 = H[i + 1][j];
                H[i][j] = cs[i] * h0 + sn[i] * h1;
                H[i + 1][j] = -sn[i] * h0 + cs[i] * h1;
            }

            const double a = H[j][j];
            const double b2 = H[j + 1][j];
            const double rho = std::hypot(a, b2);
            if (rho > 0.0) {
                cs[j] = a / rho;
                sn[j] = b2 / rho;
                H[j][j] = rho;
                H[j + 1][j] = 0.0;
                const double gj = g[j];
                const double gj1 = g[j + 1];
                g[j] = cs[j] * gj + sn[j] * gj1;
                g[j + 1] = -sn[j] * gj + cs[j] * gj1;
            }

            estimated_residual = std::abs(g[j + 1]);
            used = j + 1;
            if (estimated_residual <= tol || H[j + 1][j] == 0.0) break;
        }

        if (used == 0) break;

        std::fill(y.begin(), y.end(), 0.0);
        for (int i = used - 1; i >= 0; --i) {
            double sum = g[i];
            for (int j = i + 1; j < used; ++j) sum -= H[i][j] * y[j];
            if (std::abs(H[i][i]) <= 1e-30) {
                result.status = SolverStatus::DIVERGED;
                result.iterations = iterations;
                result.residual = estimated_residual;
                result.residual_relative = estimated_residual / std::max(b_norm, 1.0);
                return result;
            }
            y[i] = sum / H[i][i];
        }

        for (int j = 0; j < used; ++j)
            for (std::size_t i = 0; i < n; ++i)
                x(i) += Z[j][i] * y[j];

        residual(r);
        beta = norm(r);
        if (beta <= tol) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iterations;
            result.residual = beta;
            result.residual_relative = beta / std::max(b_norm, 1.0);
            return result;
        }
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = iterations;
    result.residual = beta;
    result.residual_relative = beta / std::max(b_norm, 1.0);
    return result;
}

} // namespace cfdx::core
