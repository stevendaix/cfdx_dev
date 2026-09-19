// M0.8-T05 — BiCGStab solver
//
// Spécification CFDX v0.7 §35 :
//   BiCGStab (van der Vorst) pour les systèmes non symétriques.

#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "cg_solver.h"
#include <cstddef>
#include <stdexcept>
#include <cmath>

namespace cfdx {
namespace core {

// Résout A x = b par BiCGStab.
//
// Args:
//   A          : matrice (pas nécessairement symétrique).
//   b          : vecteur second membre.
//   x          : solution (initialisé à l'entrée, modifié en sortie).
//   max_iter   : nombre maximum d'itérations.
//   tolerance  : tolérance sur le résidu relatif.
//
// Retourne un SolverResult.
inline SolverResult solve_bicgstab(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
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

    // Préconditionneur diagonal.
    std::vector<double> M(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
            if (Ac[k] == static_cast<std::uint32_t>(i)) {
                M[i] = Av[k];
                break;
            }
        }
        if (M[i] == 0.0) {
            result.status = SolverStatus::NOT_APPLICABLE;
            return result;
        }
    }

    // r0 = b - A x
    std::vector<double> r(n);
    const double* xd = x.data();
    for (std::size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
            sum += Av[k] * xd[Ac[k]];
        }
        r[i] = b(i) - sum;
    }

    // r_hat = r0
    std::vector<double> r_hat = r;

    // rho0 = alpha0 = omega0 = 1.0
    double rho = 1.0;
    double alpha = 1.0;
    double omega = 1.0;

    // v = p = 0
    std::vector<double> v(n, 0.0);
    std::vector<double> p(n, 0.0);

    const double b_norm = b.norm2();
    const double tol_abs = tolerance * std::max(b_norm, 1e-15);

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

    for (std::size_t iter = 1; iter <= max_iter; ++iter) {
        // rho_i = <r_hat, r>
        double rho_new = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            rho_new += r_hat[i] * r[i];
        }

        // beta = (rho_i / rho_{i-1}) * (alpha / omega)
        const double beta = (rho != 0.0)
            ? (rho_new / rho) * (alpha / omega)
            : 0.0;

        // p = r + beta (p - omega v)
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }

        // z = M^{-1} p
        std::vector<double> z(n);
        for (std::size_t i = 0; i < n; ++i) {
            z[i] = p[i] / M[i];
        }

        // v = A z
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                sum += Av[k] * z[Ac[k]];
            }
            v[i] = sum;
        }

        // alpha = rho_i / <r_hat, v>
        double rhat_v = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            rhat_v += r_hat[i] * v[i];
        }
        if (rhat_v == 0.0) {
            result.status = SolverStatus::NOT_APPLICABLE;
            result.iterations = iter;
            return result;
        }
        alpha = rho_new / rhat_v;

        // s = r - alpha v
        std::vector<double> s(n);
        for (std::size_t i = 0; i < n; ++i) {
            s[i] = r[i] - alpha * v[i];
        }

        // z_s = M^{-1} s  (variable distincte de z)
        std::vector<double> z_s(n);
        for (std::size_t i = 0; i < n; ++i) {
            z_s[i] = s[i] / M[i];
        }

        // t = A z_s
        std::vector<double> t(n);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                sum += Av[k] * z_s[Ac[k]];
            }
            t[i] = sum;
        }

        // omega = <t, s> / <t, t>
        double ts = 0.0;
        double tt = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            ts += t[i] * s[i];
            tt += t[i] * t[i];
        }
        if (tt == 0.0) {
            // s == 0 → convergence : x = x + alpha * z (z = M^{-1} p).
            for (std::size_t i = 0; i < n; ++i) {
                x(i) += alpha * z[i];
            }
            result.status = SolverStatus::CONVERGED;
            result.iterations = iter;
            result.residual = 0.0;
            result.residual_relative = 0.0;
            return result;
        }
        omega = ts / tt;

        // x = x + alpha z + omega z_s
        for (std::size_t i = 0; i < n; ++i) {
            x(i) += alpha * z[i] + omega * z_s[i];
        }

        // r = s - omega t
        for (std::size_t i = 0; i < n; ++i) {
            r[i] = s[i] - omega * t[i];
        }

        // Vérification de la convergence.
        res = 0.0;
        for (std::size_t i = 0; i < n; ++i) res += r[i] * r[i];
        res = std::sqrt(res);

        if (res < tol_abs) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iter;
            result.residual = res;
            result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
            return result;
        }

        rho = rho_new;
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = max_iter;
    result.residual = res;
    result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
    return result;
}

}  // namespace core
}  // namespace cfdx