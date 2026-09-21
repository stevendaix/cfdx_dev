#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "cg_solver.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace cfdx {
namespace core {

/**
 * Restarted GMRES with modified Gram-Schmidt orthogonalization and
 * incremental Givens rotations.
 *
 * The Hessenberg matrix uses the conventional indexing H(i,j), where
 * Arnoldi step j (zero based) produces H(0..j+1,j).  Givens rotations
 * are applied immediately to the current column, so the least-squares
 * right-hand side is available without a second pass.
 */
inline SolverResult solve_gmres(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    int restart = 30,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12)
{
    SolverResult result;

    if (A.n_rows() != A.n_cols() ||
        b.size() != A.n_rows() ||
        x.size() != A.n_cols() ||
        restart <= 0 ||
        tolerance <= 0.0 ||
        max_iter == 0) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = A.n_rows();
    const std::size_t m = std::min<std::size_t>(
        static_cast<std::size_t>(restart), n);

    const double b_norm = b.norm2();
    const double abs_tol = tolerance * std::max(b_norm, 1.0);

    auto residual = [&](std::vector<double>& r) {
        r.assign(n, 0.0);
        const double* Av = A.values_data();
        const auto* Ac = A.columns_data();
        const auto* Ar = A.row_offsets_data();
        for (std::size_t i = 0; i < n; ++i) {
            double value = b(i);
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                value -= Av[k] * x(Ac[k]);
            }
            r[i] = value;
        }
    };

    std::vector<double> r;
    residual(r);
    double beta = 0.0;
    for (double v : r) beta += v * v;
    beta = std::sqrt(beta);

    if (beta <= abs_tol) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = beta;
        result.residual_relative = (b_norm > 0.0) ? beta / b_norm : 0.0;
        return result;
    }

    std::size_t total_iter = 0;

    while (total_iter < max_iter) {
        const std::size_t cycle_m = std::min(m, max_iter - total_iter);

        std::vector<std::vector<double>> V(
            cycle_m + 1, std::vector<double>(n, 0.0));
        std::vector<std::vector<double>> H(
            cycle_m + 1, std::vector<double>(cycle_m, 0.0));

        for (std::size_t i = 0; i < n; ++i) {
            V[0][i] = r[i] / beta;
        }

        std::vector<double> cs(cycle_m, 0.0);
        std::vector<double> sn(cycle_m, 0.0);
        std::vector<double> g(cycle_m + 1, 0.0);
        g[0] = beta;

        std::size_t used = 0;
        double estimated_residual = beta;

        for (std::size_t j = 0; j < cycle_m; ++j) {
            std::vector<double> w(n, 0.0);
            const double* Av = A.values_data();
            const auto* Ac = A.columns_data();
            const auto* Ar = A.row_offsets_data();

            for (std::size_t row = 0; row < n; ++row) {
                double value = 0.0;
                for (std::size_t k = Ar[row]; k < Ar[row + 1]; ++k) {
                    value += Av[k] * V[j][Ac[k]];
                }
                w[row] = value;
            }

            // Modified Gram-Schmidt.
            for (std::size_t i = 0; i <= j; ++i) {
                double hij = 0.0;
                for (std::size_t k = 0; k < n; ++k) {
                    hij += V[i][k] * w[k];
                }
                H[i][j] = hij;
                for (std::size_t k = 0; k < n; ++k) {
                    w[k] -= hij * V[i][k];
                }
            }

            double hnext = 0.0;
            for (double value : w) hnext += value * value;
            hnext = std::sqrt(hnext);
            H[j + 1][j] = hnext;

            if (hnext > std::numeric_limits<double>::epsilon() * beta) {
                for (std::size_t k = 0; k < n; ++k) {
                    V[j + 1][k] = w[k] / hnext;
                }
            }

            // Apply all previously computed Givens rotations.
            for (std::size_t i = 0; i < j; ++i) {
                const double h0 = H[i][j];
                const double h1 = H[i + 1][j];
                H[i][j] = cs[i] * h0 + sn[i] * h1;
                H[i + 1][j] = -sn[i] * h0 + cs[i] * h1;
            }

            // Generate the new rotation. This handles both normal Arnoldi
            // steps and a happy breakdown (hnext == 0).
            const double h0 = H[j][j];
            const double h1 = H[j + 1][j];
            const double rho = std::hypot(h0, h1);

            if (rho > std::numeric_limits<double>::epsilon()) {
                cs[j] = h0 / rho;
                sn[j] = h1 / rho;
            } else {
                cs[j] = 1.0;
                sn[j] = 0.0;
            }

            H[j][j] = cs[j] * h0 + sn[j] * h1;
            H[j + 1][j] = 0.0;

            const double g0 = g[j];
            const double g1 = g[j + 1];
            g[j] = cs[j] * g0 + sn[j] * g1;
            g[j + 1] = -sn[j] * g0 + cs[j] * g1;

            estimated_residual = std::abs(g[j + 1]);
            used = j + 1;
            ++total_iter;

            if (estimated_residual <= abs_tol ||
                hnext <= std::numeric_limits<double>::epsilon() * beta) {
                break;
            }

            if (total_iter >= max_iter) break;
        }

        if (used == 0) break;

        // Solve the upper triangular least-squares system H y = g.
        std::vector<double> y(used, 0.0);
        bool singular = false;
        for (std::size_t ii = used; ii-- > 0;) {
            double rhs = g[ii];
            for (std::size_t j = ii + 1; j < used; ++j) {
                rhs -= H[ii][j] * y[j];
            }
            const double diag = H[ii][ii];
            const double scale = std::max(1.0, std::abs(H[0][0]));
            if (std::abs(diag) <= std::numeric_limits<double>::epsilon() * scale) {
                singular = true;
                break;
            }
            y[ii] = rhs / diag;
        }

        if (!singular) {
            for (std::size_t j = 0; j < used; ++j) {
                for (std::size_t i = 0; i < n; ++i) {
                    x(i) += V[j][i] * y[j];
                }
            }
        }

        // Always recompute the true residual after updating x.  This avoids
        // accepting a false convergence caused by loss of orthogonality.
        residual(r);
        beta = 0.0;
        for (double value : r) beta += value * value;
        beta = std::sqrt(beta);

        if (beta <= abs_tol) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = total_iter;
            result.residual = beta;
            result.residual_relative = (b_norm > 0.0) ? beta / b_norm : 0.0;
            return result;
        }

        if (singular || used < cycle_m) {
            // A breakdown without convergence indicates that the current
            // Krylov space cannot provide another independent direction.
            break;
        }
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = total_iter;
    result.residual = beta;
    result.residual_relative = (b_norm > 0.0) ? beta / b_norm : 0.0;
    return result;
}

}  // namespace core
}  // namespace cfdx
