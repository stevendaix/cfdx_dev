#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>

namespace cfdx::core {

// Exact algebraic Schur complement approximation:
//
//     S = C - D Auu^{-1} G
//
// Unlike CoupledBlockSchurAMGPreconditioner, which uses the cell-local 3x3
// block inverse M_b^{-1} in place of Auu^{-1} and assembles S~ explicitly,
// this class applies the true Schur operator implicitly: each matvec S*p
// costs one Auu solve, one G matvec, one D matvec and one C matvec.
//
// The caller supplies the Auu solve as a callback so that any velocity
// preconditioner (AMG, dense LU on tiny systems, Krylov, ...) can be
// plugged in without this class depending on a specific solver stack.
//
// The pressure system S*p = rhs_p is solved by an internal preconditioner-
// free conjugate gradient. Two reasons:
//   1. It keeps this class header-only and dependency-light.
//   2. When Auu is SPD and A is a standard incompressible saddle-point
//      operator, S is SPD, so CG is the correct outer Krylov choice. The
//      SPD assumption is checked at runtime through p^T S p > 0.
//
// This is intentionally the reference/oracle implementation. It is not
// meant to replace SIMPLE/SIMPLEC/LSC/BFBt/PCD preconditioners in
// production: the inner Auu solve to full accuracy per outer matvec is
// deliberately expensive. It is meant to be the source of truth against
// which cheaper approximations are compared. Tracked in #481.
class ExactSchurApproximation final : public SchurApproximation {
public:
    struct Controls {
        std::size_t max_iter = 200;
        double tolerance = 1e-10;
    };

    // Solve Auu * y = rhs, filling y. Return false on failure.
    // The callback owns any preconditioner/hierarchy state; it must be
    // consistent with the Auu block passed at setup().
    using AuuSolve = std::function<bool(const Vector& rhs, Vector& y)>;

    explicit ExactSchurApproximation(AuuSolve auu_solve,
                                     const Controls& controls = Controls())
        : auu_solve_(std::move(auu_solve)), controls_(controls) {}

    const char* name() const noexcept override { return "exact_schur"; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || !auu_solve_) return false;
        blocks_ = &blocks;
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        return setup(blocks);
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (!blocks_ || !auu_solve_) return false;
        const std::size_t n = blocks_->pressure_size();
        if (rhs_p.size() != n) return false;

        if (pressure.size() != n) pressure = Vector(n, 0.0);
        else pressure.fill(0.0);

        // Preconditioner-free CG on the implicit Schur operator.
        Vector r(n), p(n), Ap(n);
        for (std::size_t i = 0; i < n; ++i) r(i) = rhs_p(i);
        for (std::size_t i = 0; i < n; ++i) p(i) = r(i);
        double r_dot_r = r.dot(r);
        const double r0 = std::sqrt(r_dot_r);
        if (!(r0 > 0.0)) return true;  // rhs is zero, pressure=0 is the solution

        const double abs_tol = controls_.tolerance * r0;

        for (std::size_t it = 0; it < controls_.max_iter; ++it) {
            if (!apply_schur(p, Ap)) return false;

            double pAp = 0.0;
            for (std::size_t i = 0; i < n; ++i) pAp += p(i) * Ap(i);
            if (!(pAp > 0.0) || !std::isfinite(pAp)) return false;

            const double alpha = r_dot_r / pAp;
            for (std::size_t i = 0; i < n; ++i) pressure(i) += alpha * p(i);
            for (std::size_t i = 0; i < n; ++i) r(i) -= alpha * Ap(i);

            const double new_r_dot_r = r.dot(r);
            if (std::sqrt(new_r_dot_r) < abs_tol) return true;

            const double beta = new_r_dot_r / r_dot_r;
            for (std::size_t i = 0; i < n; ++i) p(i) = r(i) + beta * p(i);
            r_dot_r = new_r_dot_r;
        }
        return false;  // max_iter reached without meeting tolerance
    }

    // Public so tests can compare S*p to a dense oracle without running the
    // outer CG. Not part of the SchurApproximation interface.
    bool apply_schur(const Vector& p, Vector& out) const {
        if (!blocks_) return false;
        const std::size_t np = blocks_->pressure_size();
        const std::size_t nu = blocks_->velocity_size();
        if (p.size() != np) return false;
        if (out.size() != np) out = Vector(np, 0.0);

        // Gp = G * p
        const auto Gp = blocks_->G().matvec(p);
        Vector Gp_v(nu, 0.0);
        for (std::size_t i = 0; i < nu; ++i) Gp_v(i) = Gp[i];

        // y = Auu^{-1} * Gp
        Vector y(nu, 0.0);
        if (!auu_solve_(Gp_v, y)) return false;

        // out = C*p - D*y
        const auto Cp = blocks_->C().matvec(p);
        const auto Dy = blocks_->D().matvec(y);
        for (std::size_t i = 0; i < np; ++i) out(i) = Cp[i] - Dy[i];
        return true;
    }

private:
    const BlockOperator* blocks_ = nullptr;
    AuuSolve auu_solve_;
    Controls controls_;
};

} // namespace cfdx::core
