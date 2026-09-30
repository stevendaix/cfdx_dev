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
#include "mixed_precision.h"
#include "krylov_reductions.h"
#include "null_space.h"
#include "preconditioner.h"
#include <cstddef>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <vector>

namespace cfdx {
namespace core {

enum class CGResidualReplacementPolicy : std::uint8_t {
    Adaptive,
    Disabled,
    Periodic
};

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

struct SolverDiagnosticSample {
    std::size_t iteration = 0;
    double true_residual = 0.0;
    double recursive_residual = 0.0;
    double preconditioned_dot = 0.0;
    double pAp = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
};

struct SolverResult {
    SolverStatus status = SolverStatus::NOT_APPLICABLE;
    std::size_t iterations = 0;
    double residual = 0.0;
    double residual_relative = 0.0;
    // Diagnostic evidence sampled during the iteration. This is intentionally
    // observational: it does not alter stopping criteria or Krylov updates.
    std::vector<SolverDiagnosticSample> diagnostics;
    // Residual-replacement evidence is observational and does not alter
    // convergence criteria. It is used to distinguish Krylov recurrence drift
    // from preconditioner/coarse-space limitations.
    std::size_t residual_replacements = 0;
    double max_true_recursive_gap = 0.0;
    double min_true_recursive_ratio = std::numeric_limits<double>::infinity();
    double max_true_recursive_ratio = 0.0;
    double max_true_residual = 0.0;
    double min_true_residual = std::numeric_limits<double>::infinity();
    std::size_t first_residual_replacement = 0;
    std::size_t last_residual_replacement = 0;
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
namespace detail {
inline SolverResult solve_cg_impl(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    std::size_t max_iter,
    double tolerance,
    PrecisionPolicy precision,
    KrylovReductionPolicy reduction,
    Preconditioner* preconditioner,
    const NullSpaceProjector* null_space,
    bool setup_preconditioner,
    CGResidualReplacementPolicy residual_policy = CGResidualReplacementPolicy::Adaptive,
    std::size_t periodic_replacement_interval = 0)
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

    if (null_space &&
        (null_space->dimension() != A.n_rows() ||
         !null_space->is_compatible(b) ||
         !null_space->is_null_space(A))) {
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
    if (null_space) null_space->remove(x);
    const double* Av = A.values_data();
    const bool mp32 = precision.enabled && precision.operator_precision == SolverPrecision::FP32;
    const SolverPrecision redp = precision.enabled ? precision.reduction_precision : SolverPrecision::FP64;
    const auto* Ac = A.columns_data();
    const auto* Ar = A.row_offsets_data();

    // Preserve the historical Jacobi path when no external preconditioner is
    // supplied. This keeps the native solver behavior and API unchanged.
    std::vector<double> M(n, 0.0);
    if (preconditioner) {
        if (setup_preconditioner && !preconditioner->setup(A)) {
            result.status = SolverStatus::NOT_APPLICABLE;
            return result;
        }
    } else {
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
    }

    // r = b - A x
    std::vector<double> r(n);
    Vector rv(n); mixed_precision_true_residual(A,b,x,rv);
    if (null_space) null_space->remove(rv);
    for (std::size_t i=0;i<n;++i) r[i]=rv(i);

    // z = M^{-1} r
    std::vector<double> z(n);
    Vector z_vector(n);
    Vector r_vector(n);
    const auto apply_preconditioner = [&]() {
        if (preconditioner) {
            for (std::size_t i = 0; i < n; ++i) r_vector(i) = r[i];
            if (!preconditioner->apply(r_vector, z_vector)) return false;
            for (std::size_t i = 0; i < n; ++i) z[i] = z_vector(i);
        } else {
            for (std::size_t i = 0; i < n; ++i) z[i] = r[i] / M[i];
        }
        if (null_space) null_space->remove(z);
        for (const double value : z)
            if (!std::isfinite(value)) return false;
        return true;
    };
    if (!apply_preconditioner()) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    // p = z
    std::vector<double> p = z;

    for (std::size_t i = 0; i < n; ++i) z_vector(i) = z[i];
    for (std::size_t i = 0; i < n; ++i) r_vector(i) = r[i];
    double rsold = krylov_dot(r_vector, z_vector, redp, reduction);

    const double b_norm = krylov_norm2(b, SolverPrecision::FP64, reduction);
    const double tol_abs = tolerance * std::max(b_norm, 1e-15);

    if (!std::isfinite(rsold)) { result.status=SolverStatus::DIVERGED; return result; }

    const double initial_true_residual =
        krylov_norm2(r_vector, SolverPrecision::FP64, reduction);
    if (preconditioner) {
        result.diagnostics.push_back(
            {0, initial_true_residual, initial_true_residual,
             rsold, 0.0, 0.0, 0.0});
    }
    if ((preconditioner && initial_true_residual <= tol_abs) ||
        (!preconditioner && rsold < tol_abs * tol_abs)) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = preconditioner
            ? initial_true_residual
            : krylov_norm2(rv, redp, reduction);
        result.residual_relative = (b_norm > 0.0) ? result.residual / b_norm : 0.0;
        return result;
    }

    // Residual replacement only pays for itself while b-Ax is still
    // decreasing. Once the iteration reaches the attainable accuracy of the
    // system -- about eps*cond(A) in relative terms, which no solver, direct or
    // iterative, can improve upon in FP64 -- the recursive residual keeps
    // shrinking while b-Ax does not. The gap criterion below is then satisfied
    // at every remaining iteration, and the restart that follows each
    // replacement discards the Krylov direction, degrading CG into
    // preconditioned steepest descent for the rest of the budget. Requiring
    // measurable progress since the previous replacement bounds the number of
    // replacements to O(log(r0/attainable_floor)).
    double residual_at_last_replacement = std::numeric_limits<double>::infinity();

    for (std::size_t iter = 1; iter <= max_iter; ++iter) {
        // Ap = A p
        std::vector<double> Ap(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                sum += mp32 ? static_cast<double>(static_cast<float>(Av[k])*static_cast<float>(p[Ac[k]])) : Av[k]*p[Ac[k]];
            }
            Ap[i] = sum;
        }
        if (null_space) null_space->remove(Ap);

        // alpha = rsold / (p · Ap)
        Vector p_vector(n), Ap_vector(n);
        for (std::size_t i = 0; i < n; ++i) {
            p_vector(i) = p[i];
            Ap_vector(i) = Ap[i];
        }
        const double pAp = krylov_dot(p_vector, Ap_vector, redp, reduction);
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
        if (null_space) {
            null_space->remove(x);
            null_space->remove(r);
        }

        for (std::size_t i = 0; i < n; ++i)
            if (!std::isfinite(x(i)) || !std::isfinite(r[i])) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }

        // z = M^{-1} r
        if (!apply_preconditioner()) {
            result.status = SolverStatus::NOT_APPLICABLE;
            result.iterations = iter;
            return result;
        }

        // rsnew = r · z
        for (std::size_t i = 0; i < n; ++i) {
            r_vector(i) = r[i];
            z_vector(i) = z[i];
        }
        const double rsnew = krylov_dot(r_vector, z_vector, redp, reduction);

        if (!std::isfinite(rsnew)) { result.status=SolverStatus::DIVERGED; result.iterations=iter; return result; }
        // A preconditioned norm is not a true residual norm. Re-evaluate
        // b-A*x in FP64 for the external-preconditioner path before deciding
        // convergence, just as the reported final residual is evaluated.
        const double res = preconditioner
            ? mixed_precision_true_residual(A, b, x, rv)
            : std::sqrt(std::abs(rsnew));
        if (preconditioner &&
            (iter == 1 || iter % 100 == 0 || iter == max_iter)) {
            const double recursive_residual =
                std::sqrt(std::max(0.0, krylov_dot(r_vector, r_vector,
                                                   SolverPrecision::FP64,
                                                   reduction)));
            result.diagnostics.push_back(
                {iter, res, recursive_residual, rsnew, pAp, alpha,
                 rsold != 0.0 ? rsnew / rsold : 0.0});
        }
        if (res < tol_abs) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iter;
            result.residual = res;
            result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
            return result;
        }

        // The recursive CG residual can drift away from b-Ax in finite
        // precision, especially once an effective AMG preconditioner has
        // reduced the true residual by many orders of magnitude. Continuing
        // with the stale recurrence can eventually destroy conjugacy and make
        // p^T A p appear non-positive even though A is SPD. Replacing the
        // recursive residual when the gap becomes material is a standard
        // residual-replacement safeguard; the true residual is already
        // available here, so this adds no extra fine-grid matvec.
        if (preconditioner) {
            double recursive_r2 = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                recursive_r2 += r[i] * r[i];
            const double recursive_res = std::sqrt(recursive_r2);
            const double gap = std::abs(res - recursive_res);
            const double scale = std::max(res, recursive_res);
            result.max_true_recursive_gap =
                std::max(result.max_true_recursive_gap, gap);
            if (recursive_res > 0.0 && std::isfinite(recursive_res)) {
                const double ratio = res / recursive_res;
                if (std::isfinite(ratio)) {
                    result.min_true_recursive_ratio =
                        std::min(result.min_true_recursive_ratio, ratio);
                    result.max_true_recursive_ratio =
                        std::max(result.max_true_recursive_ratio, ratio);
                }
            }
            result.max_true_residual = std::max(result.max_true_residual, res);
            result.min_true_residual = std::min(result.min_true_residual, res);
            constexpr double residual_gap_fraction = 0.25;
            const bool made_progress =
                res < 0.5 * residual_at_last_replacement;
            // A non-positive r^T z means the recursive residual has been
            // annihilated by cancellation. z = M^{-1}r is then numerically
            // zero, so the next search direction is zero and p^T A p vanishes.
            // That is a dead recurrence, not an unusable preconditioner, and it
            // must be recognized here: otherwise it surfaces one iteration
            // later as a spurious NOT_APPLICABLE from the p^T A p guard, which
            // wrongly tells the caller to abandon a perfectly valid
            // preconditioner.
            const bool recurrence_collapsed = !(rsnew > 0.0);
            const bool adaptive_replacement =
                residual_policy == CGResidualReplacementPolicy::Adaptive &&
                std::isfinite(recursive_res) && scale > 0.0 &&
                gap > residual_gap_fraction * scale && made_progress;
            const bool periodic_replacement =
                residual_policy == CGResidualReplacementPolicy::Periodic &&
                periodic_replacement_interval > 0 &&
                (iter % periodic_replacement_interval == 0);
            // Repairing a collapsed recurrence is itself a replacement, so the
            // Disabled policy must not do it, and it is pointless once b-Ax has
            // stopped improving.
            const bool repair_collapse =
                recurrence_collapsed && made_progress &&
                residual_policy != CGResidualReplacementPolicy::Disabled;
            if (recurrence_collapsed && !repair_collapse &&
                !periodic_replacement) {
                // The recurrence is dead and replacing it cannot help, so the
                // iteration has extracted all the accuracy this system allows.
                // Report the accuracy actually achieved, at the iteration where
                // it was reached, instead of spending the remaining budget on a
                // zero search direction.
                result.status = SolverStatus::MAX_ITER_REACHED;
                result.iterations = iter;
                result.residual = res;
                result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
                return result;
            }
            if (adaptive_replacement || periodic_replacement ||
                repair_collapse) {
                ++result.residual_replacements;
                residual_at_last_replacement = res;
                if (result.first_residual_replacement == 0)
                    result.first_residual_replacement = iter;
                result.last_residual_replacement = iter;
                for (std::size_t i = 0; i < n; ++i) r[i] = rv(i);
                if (null_space) null_space->remove(r);
                if (!apply_preconditioner()) {
                    result.status = SolverStatus::NOT_APPLICABLE;
                    result.iterations = iter;
                    return result;
                }
                for (std::size_t i = 0; i < n; ++i) {
                    r_vector(i) = r[i];
                    z_vector(i) = z[i];
                }
                const double replaced_rsnew =
                    krylov_dot(r_vector, z_vector, redp, reduction);
                if (!std::isfinite(replaced_rsnew) || !(replaced_rsnew > 0.0)) {
                    result.status = SolverStatus::NOT_APPLICABLE;
                    result.iterations = iter;
                    return result;
                }
                // The replacement recomputes b-Ax after the normal convergence
                // check. Re-test the true residual here so a converged iterate
                // produced exactly at a replacement is not reported as
                // MAX_ITER_REACHED merely because the replacement was the last
                // operation of the iteration.
                if (res <= tol_abs) {
                    result.status = SolverStatus::CONVERGED;
                    result.iterations = iter;
                    result.residual = res;
                    result.residual_relative = (b_norm > 0.0) ? res / b_norm : 0.0;
                    return result;
                }
                // Restart the search direction after replacement. This avoids
                // carrying a direction built from a residual that no longer
                // represents b-Ax while retaining the current iterate.
                p = z;
                rsold = replaced_rsnew;
                continue;
            }
        }

        // beta = rsnew / rsold
        const double beta = rsnew / rsold;

        // p = z + beta p
        for (std::size_t i = 0; i < n; ++i) {
            p[i] = z[i] + beta * p[i];
        }
        if (null_space) null_space->remove(p);

        rsold = rsnew;
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = max_iter;
    result.residual = mixed_precision_true_residual(A,b,x,rv);
    result.residual_relative = (b_norm > 0.0) ? result.residual / b_norm : 0.0;
    return result;
}
} // namespace detail

// Diagnostic/qualification entry point: keeps the production solve_cg API
// unchanged while allowing V&V campaigns to isolate residual-replacement
// behavior without changing the stopping criterion.
inline SolverResult solve_cg_controlled(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    std::size_t max_iter,
    double tolerance,
    CGResidualReplacementPolicy residual_policy,
    std::size_t periodic_replacement_interval = 0)
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, {}, {},
        &preconditioner, nullptr, true,
        residual_policy, periodic_replacement_interval);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        nullptr, nullptr, true);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        &preconditioner, nullptr, true);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    const NullSpaceProjector& null_space,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        nullptr, &null_space, true);
}

inline SolverResult solve_cg(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    const NullSpaceProjector& null_space,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {})
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        &preconditioner, &null_space, true);
}

// Solve with a preconditioner that has already been prepared for A. This is
// the low-level path used by reusable nonlinear/pressure-solver contexts.
inline SolverResult solve_cg_prepared(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    Preconditioner& preconditioner,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    PrecisionPolicy precision = {},
    KrylovReductionPolicy reduction = {},
    const NullSpaceProjector* null_space = nullptr)
{
    return detail::solve_cg_impl(
        A, b, x, max_iter, tolerance, precision, reduction,
        &preconditioner, null_space, false);
}

}  // namespace core
}  // namespace cfdx
