#pragma once

#include "sparse_matrix.h"
#include "cg_solver.h"
#include "vector.h"
#include "preconditioner.h"
#include "solver_workspace.h"
#include "krylov_controls.h"
#include "mixed_precision.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <cstdlib>
#include <stdexcept>

namespace cfdx::core {

struct LinearOperator {
    std::size_t size = 0;
    std::function<void(const Vector&, Vector&)> apply;
};

inline SolverResult solve_gmres(
    const LinearOperator& op,
    const Vector& b,
    Vector& x,
    int restart = 30,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    Preconditioner* preconditioner = nullptr,
    KrylovControls controls = {})
{
    SolverResult result;
    if (op.size == 0 || !op.apply || b.size() != op.size || x.size() != op.size ||
        restart <= 0 || max_iter == 0 || tolerance <= 0.0) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const std::size_t n = op.size;
    int current_restart = std::clamp(restart, controls.restart_min, controls.restart_max);
    current_restart = std::min<int>(current_restart, static_cast<int>(n));
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || max_iter == 0) {
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }

    const double b_norm = krylov_norm2(b, SolverPrecision::FP64, controls.reduction);
    if (!std::isfinite(b_norm)) {
        result.status = SolverStatus::DIVERGED;
        return result;
    }
    const double tol = tolerance * std::max(b_norm, 1.0);

    GmresWorkspace w;
    w.resize(n, static_cast<std::size_t>(current_restart));

    auto apply_operator = [&](const double* in, Vector& out) {
        for (std::size_t i = 0; i < n; ++i) {
            w.vin(i) = in[i];
        }
        op.apply(w.vin, out);
    };
    auto true_residual = [&]() {
        apply_operator(x.data(), w.ax);
        for (std::size_t i = 0; i < n; ++i) w.r(i) = b(i) - w.ax(i);
        return krylov_norm2(w.r, SolverPrecision::FP64, controls.reduction);
    };

    double beta = true_residual();
    if (beta <= tol) {
        result.status = SolverStatus::CONVERGED;
        result.iterations = 0;
        result.residual = beta;
        result.residual_relative = beta / std::max(b_norm, 1.0);
        return result;
    }

    std::size_t iterations = 0;
    double previous_cycle_residual = beta;

    while (iterations < max_iter) {
        w.resize(n, static_cast<std::size_t>(current_restart));
        beta = true_residual();
        if (beta <= tol) break;

        for (std::size_t i = 0; i < n; ++i) w.v(0)[i] = w.r(i) / beta;
        std::fill(w.h.begin(), w.h.end(), 0.0);
        std::fill(w.cs.begin(), w.cs.end(), 0.0);
        std::fill(w.sn.begin(), w.sn.end(), 0.0);
        std::fill(w.g.begin(), w.g.end(), 0.0);
        w.g[0] = beta;

        int used = 0;
        double estimated_residual = beta;
        bool arnoldi_breakdown = false;

        for (int j = 0; j < current_restart && iterations < max_iter; ++j) {
            if (preconditioner) {
                for (std::size_t i = 0; i < n; ++i) w.z(i) = w.v(j)[i];
                if (!preconditioner->apply(w.z, w.vout)) {
                    result.status = SolverStatus::NOT_APPLICABLE;
                    result.iterations = iterations;
                    return result;
                }
                for (std::size_t i = 0; i < n; ++i) w.zv(static_cast<std::size_t>(j))[i] = w.vout(i);
            } else {
                std::copy(w.v(j), w.v(j) + n, w.zv(static_cast<std::size_t>(j)));
            }

            apply_operator(w.zv(static_cast<std::size_t>(j)), w.w);
            const double arnoldi_action_norm =
                krylov_norm2(w.w, SolverPrecision::FP64, controls.reduction);
            // Modified Gram-Schmidt once is vulnerable to loss of
            // orthogonality on nonsymmetric saddle-point systems. A second
            // orthogonalization pass is inexpensive for the small acceptance
            // systems and substantially improves the Hessenberg relation.
            for (int i = 0; i <= j; ++i) {
                const double* vi = w.v(static_cast<std::size_t>(i));
                Vector vi_vector(n);
                for (std::size_t k = 0; k < n; ++k) vi_vector(k) = vi[k];
                const double h = krylov_dot(vi_vector, w.w, SolverPrecision::FP64, controls.reduction);
                w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = h;
                for (std::size_t k = 0; k < n; ++k) w.w(k) -= h * vi[k];
            }
            for (int i = 0; i <= j; ++i) {
                const double* vi = w.v(static_cast<std::size_t>(i));
                Vector vi_vector(n);
                for (std::size_t k = 0; k < n; ++k) vi_vector(k) = vi[k];
                const double correction = krylov_dot(
                    vi_vector, w.w, SolverPrecision::FP64, controls.reduction);
                w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) += correction;
                for (std::size_t k = 0; k < n; ++k) w.w(k) -= correction * vi[k];
            }

            // Scale the breakdown test with ||A z_j||, not with an arbitrary
            // absolute value of one. The operator can legitimately be very
            // small after nondimensionalisation or physical scaling; an
            // absolute epsilon floor would then misclassify a valid Krylov
            // vector as a breakdown.
            const double hnext = krylov_norm2(w.w, SolverPrecision::FP64, controls.reduction);
            w.H(static_cast<std::size_t>(j + 1), static_cast<std::size_t>(j)) = hnext;
            // Arnoldi vectors are normalized, so loss of the new direction is
            // a breakdown. Use a relative test against the current operator
            // action, but never let a tiny absolute value trigger it.
            const double arnoldi_breakdown_floor =
                128.0 * std::numeric_limits<double>::epsilon() *
                std::max(arnoldi_action_norm, 1e-300);
            const bool arnoldi_breakdown_detected =
                hnext <= arnoldi_breakdown_floor && arnoldi_action_norm > 0.0;
            if (hnext > 0.0) {
                double* vnext = w.v(static_cast<std::size_t>(j + 1));
                for (std::size_t k = 0; k < n; ++k) vnext[k] = w.w(k) / hnext;
            }

            for (int i = 0; i < j; ++i) {
                const double h0 = w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
                const double h1 = w.H(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j));
                w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) =
                    w.cs[static_cast<std::size_t>(i)] * h0 + w.sn[static_cast<std::size_t>(i)] * h1;
                w.H(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j)) =
                    -w.sn[static_cast<std::size_t>(i)] * h0 + w.cs[static_cast<std::size_t>(i)] * h1;
            }

            const double a = w.H(static_cast<std::size_t>(j), static_cast<std::size_t>(j));
            const double b2 = w.H(static_cast<std::size_t>(j + 1), static_cast<std::size_t>(j));
            const double rho = std::hypot(a, b2);
            if (rho > 0.0) {
                w.cs[static_cast<std::size_t>(j)] = a / rho;
                w.sn[static_cast<std::size_t>(j)] = b2 / rho;
                w.H(static_cast<std::size_t>(j), static_cast<std::size_t>(j)) = rho;
                w.H(static_cast<std::size_t>(j + 1), static_cast<std::size_t>(j)) = 0.0;
                const double gj = w.g[static_cast<std::size_t>(j)];
                const double gj1 = w.g[static_cast<std::size_t>(j + 1)];
                w.g[static_cast<std::size_t>(j)] =
                    w.cs[static_cast<std::size_t>(j)] * gj + w.sn[static_cast<std::size_t>(j)] * gj1;
                w.g[static_cast<std::size_t>(j + 1)] =
                    -w.sn[static_cast<std::size_t>(j)] * gj + w.cs[static_cast<std::size_t>(j)] * gj1;
            }

            ++iterations;
            ++used;
            estimated_residual = std::abs(w.g[static_cast<std::size_t>(j + 1)]);

            if (controls.residual_replacement &&
                controls.residual_replacement_period > 0 &&
                iterations % controls.residual_replacement_period == 0) {
                const double exact = true_residual();
                if (exact <= tol) {
                    result.status = SolverStatus::CONVERGED;
                    result.iterations = iterations;
                    result.residual = exact;
                    result.residual_relative = exact / std::max(b_norm, 1.0);
                    return result;
                }
            }

            if (estimated_residual <= tol)
                break;
            if (arnoldi_breakdown_detected) {
                arnoldi_breakdown = true;
                break;
            }
        }

        if (used == 0) break;

        std::fill(w.y.begin(), w.y.end(), 0.0);
        for (int i = used - 1; i >= 0; --i) {
            double sum = w.g[static_cast<std::size_t>(i)];
            for (int j = i + 1; j < used; ++j)
                sum -= w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) *
                       w.y[static_cast<std::size_t>(j)];
            const double diag = w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(i));
            if (std::abs(diag) <= 1e-30) {
                // Happy Arnoldi breakdown can produce a zero diagonal in the
                // Givens-reduced Hessenberg system. If the corresponding
                // transformed RHS is also zero, that coefficient is free and
                // must not be reported as a solver divergence. The exact
                // residual is checked immediately after the Krylov update.
                if (std::abs(sum) <= tol) {
                    w.y[static_cast<std::size_t>(i)] = 0.0;
                    continue;
                }
                if (std::getenv("CFDX_DEBUG_COUPLED"))
                    std::cerr << "GMRES_BREAKDOWN i=" << i
                              << " diag=" << diag
                              << " rhs=" << sum
                              << " estimated_residual=" << estimated_residual
                              << " iterations=" << iterations << "\n";
                // The Givens/Hessenberg estimate can be exactly zero at a
                // non-happy breakdown. It is not the physical residual.
                // Report the true ||b-Ax||_2 at the last admissible iterate.
                const double breakdown_residual = true_residual();
                result.status = SolverStatus::DIVERGED;
                result.iterations = iterations;
                result.residual = breakdown_residual;
                result.residual_relative =
                    breakdown_residual / std::max(b_norm, 1.0);
                if (std::getenv("CFDX_DEBUG_COUPLED"))
                    std::cerr << "GMRES_BREAKDOWN_TRUE_RESIDUAL value="
                              << breakdown_residual << "\n";
                return result;
            }
            w.y[static_cast<std::size_t>(i)] = sum / diag;
        }

        for (int j = 0; j < used; ++j) {
            const double alpha = w.y[static_cast<std::size_t>(j)];
            const double* zj = w.zv(static_cast<std::size_t>(j));
            for (std::size_t i = 0; i < n; ++i) x(i) += alpha * zj[i];
        }

        beta = true_residual();
        if (beta <= tol) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iterations;
            result.residual = beta;
            result.residual_relative = beta / std::max(b_norm, 1.0);
            return result;
        }

        // Arnoldi breakdown means that the Krylov space stopped growing.
        // If the updated iterate is still not converged, continuing into a
        // fresh cycle can only repeat the same subspace and misclassify a
        // singular/inconsistent system as MAX_ITER_REACHED. Report the
        // physical residual as a genuine solver failure instead.
        if (arnoldi_breakdown) {
            result.status = SolverStatus::DIVERGED;
            result.iterations = iterations;
            result.residual = beta;
            result.residual_relative = beta / std::max(b_norm, 1.0);
            return result;
        }

        if (std::getenv("CFDX_DEBUG_COUPLED"))
            std::cerr << "GMRES_CYCLE iterations=" << iterations
                      << " true_residual=" << beta
                      << " relative=" << beta / std::max(b_norm, 1.0)
                      << " restart=" << current_restart << "\n";

        const double reduction = beta / std::max(previous_cycle_residual, 1e-300);
        if (controls.adaptive_restart)
            current_restart = choose_gmres_restart(current_restart, reduction, controls);
        previous_cycle_residual = beta;
    }

    result.status = SolverStatus::MAX_ITER_REACHED;
    result.iterations = iterations;
    result.residual = beta;
    result.residual_relative = beta / std::max(b_norm, 1.0);
    return result;
}

inline SolverResult solve_gmres(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    int restart,
    std::size_t max_iter,
    double tolerance,
    Preconditioner* preconditioner,
    KrylovControls controls,
    PrecisionPolicy precision = {}) {
    LinearOperator op;
    op.size = A.n_rows();
    op.apply = [&A, precision](const Vector& in, Vector& out) {
        if (out.size() != A.n_rows()) out.resize(A.n_rows());
        if (precision.enabled && precision.operator_precision == SolverPrecision::FP32) { mixed_precision_matvec(A, in, out, SolverPrecision::FP32); return; }
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            double sum = 0.0;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                sum += val[k] * in(col[k]);
            out(i) = sum;
        }
    };
    if (preconditioner && !preconditioner->setup(A)) {
        SolverResult result;
        result.status = SolverStatus::NOT_APPLICABLE;
        return result;
    }
    return solve_gmres(op, b, x, restart, max_iter, tolerance, preconditioner, controls);
}

inline SolverResult solve_gmres(
    const SparseMatrix& A,
    const Vector& b,
    Vector& x,
    int restart = 30,
    std::size_t max_iter = 1000,
    double tolerance = 1e-12,
    Preconditioner* preconditioner = nullptr,
    PrecisionPolicy precision = {})
{
    return solve_gmres(
        A, b, x, restart, max_iter, tolerance, preconditioner,
        KrylovControls{}, precision);
}

} // namespace cfdx::core
