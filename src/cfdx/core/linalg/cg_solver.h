// M0.8-T04 — Conjugate Gradient solver
//
// Spécification CFDX v0.7 §35 :
//   Premiers solveurs : CG, BiCGStab, GMRES.
//
// CG (Hestenes-Stiefel) :
//   Résout A x = b pour A symétrique définie positive.
//  avec préconditionneur diagonal (Jacobi) par défaut.

#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include <cstddef>
#include <stdexcept>
#include <cmath>

namespace cfdx {
namespace core {

enum class SolverStatus : std::uint8_t {
    CONVERGED = 0,
    MAX_ITER_REACHED,
    DIVERGED,
    NOT_APPLICABLE
};

inline const char* to_string(SolverStatus s) {
    switch (s) {
        case SolverStatus::CONVERGED:       return "converged";
        case SolverStatus::MAX_ITER_REACHED: return "max_iter_reached";
        case SolverStatus::DIVERGED:        return "diverged";
        case SolverStatus::NOT_APPLICABLE:  return "not_applicable";
        default:                             return "unknown";
    }
}

struct SolverResult {
    SolverStatus status = SolverStatus::NOT_APPLICABLE;
    std::size_t iterations = 0;
    double residual = 0.0;
    double residual_relative = 0.0;
};

// Résout A x = b par la méthode du gradient conjugué (CG).
//
// Args:
//   A          : matrice symétrique définie positive.
//   b          : vecteur second membre.
//   x          : solution (initialisé à l'entrée, modifié en sortie).
//   max_iter   : nombre maximum d'itérations.
//   tolerance  : tolérance sur le résidu relatif.
//
// Retourne un SolverResult.
inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12)
{
    SolverResult result;

    if (A.n_rows() != A.n_cols()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }
    if (b.size() != A.n_rows()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }
    if (x.size() != A.n_cols()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = A.n_rows();
    for (std::size_t k = 0; k < A.nnz(); ++k)
        if (!std::isfinite(A.values_data()[k])) {
            result.status = SolverStatus::DIVERGED;
            return result;
        }
    for (std::size_t i = 0; i < b.size(); ++i)
        if (!std::isfinite(b(i)) || !std::isfinite(x(i))) {
            result.status = SolverStatus::DIVERGED;
            return result;
        }
    const double* Av = A.values_data();
    const auto* Ac = A.columns_data();
    const auto* Ar = A.row_offsets_data();

    // Préconditionneur diagonal : M = diag(A).
    std::vector<double> M(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
            if (Ac[k] == static_cast<std::uint32_t>(i)) {
                M[i] = Av[k];
                break;
            }
        }
        if (!(M[i] > 0.0) || !std::isfinite(M[i])) {
            result.status = SolverStatus::NOT_APPLICABLE;
            return result;
        }
    }

    // r = b - A x
    std::vector<double> r(n);
    const double* xd = x.data();
    for (std::size_t i = 0; i < n; ++i) {
        double sum = 0.0;
        for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
            sum += Av[k] * xd[Ac[k]];
        }
        r[i] = b(i) - sum;
    }

    // z = M^{-1} r
    std::vector<double> z(n);
    for (std::size_t i = 0; i < n; ++i) {
        z[i] = r[i] / M[i];
        if (!std::isfinite(z[i])) { result.status=SolverStatus::DIVERGED; return result; }
    }

    // p = z
    std::vector<double> p = z;

    double rsold = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        rsold += r[i] * z[i];
    }

    const double b_norm = b.norm2();
    const double b_scale = std::max(b_norm, 1e-15);
    const double tol_abs = tolerance * b_scale;

    if (!std::isfinite(rsold)) { result.status=SolverStatus::DIVERGED; return result; }

    double true_residual = 0.0;
    for (std::size_t i = 0; i < n; ++i)
        true_residual += r[i] * r[i];
    true_residual = std::sqrt(true_residual);
    if (!std::isfinite(true_residual)) {
        result.status = SolverStatus::DIVERGED;
        return result;
    }
    if (true_residual <= tol_abs) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = true_residual;
        result.residual_relative = true_residual / b_scale;
        return result;
    }

    for (std::size_t iter = 1; iter <= max_iter; ++iter) {
        // Ap = A p
        std::vector<double> Ap(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                sum += Av[k] * p[Ac[k]];
            }
            Ap[i] = sum;
        }

        // alpha = rsold / (p · Ap)
        double pAp = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            pAp += p[i] * Ap[i];
        }
        if (!(pAp > 0.0) || !std::isfinite(pAp)) {
            result.status = SolverStatus::NOT_APPLICABLE;
            result.iterations = iter;
            return result;
        }
        const double alpha = rsold / pAp;
        if (!std::isfinite(alpha)) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }

        // x = x + alpha p
        // r = r - alpha Ap
        for (std::size_t i = 0; i < n; ++i) {
            x(i) += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        for (std::size_t i = 0; i < n; ++i)
            if (!std::isfinite(x(i)) || !std::isfinite(r[i])) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }

        // z = M^{-1} r
        for (std::size_t i = 0; i < n; ++i) {
            z[i] = r[i] / M[i];
        }

        // rsnew = r · z
        double rsnew = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            rsnew += r[i] * z[i];
        }

        if (!std::isfinite(rsnew)) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }
        true_residual = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            true_residual += r[i] * r[i];
        true_residual = std::sqrt(true_residual);
        if (!std::isfinite(true_residual)) {
            result.status = SolverStatus::DIVERGED;
            result.iterations = iter;
            return result;
        }
        if (true_residual <= tol_abs) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iter;
            result.residual = true_residual;
            result.residual_relative = true_residual / b_scale;
            return result;
        }

        // beta = rsnew / rsold
        const double beta = rsnew / rsold;

        // p = z + beta p
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = z[i] + beta * p[i];
        }

        rsold = rsnew;
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = max_iter;
    result.residual = true_residual;
    result.residual_relative = (b_norm > 0.0) ? true_residual / b_norm : 0.0;
    return result;
}

}  // namespace core
}  // namespace cfdx