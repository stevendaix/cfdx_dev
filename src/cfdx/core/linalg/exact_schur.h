#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <optional>

namespace cfdx::core {

struct ExactSchurControls {
    std::size_t max_iter = 200;
    double tolerance = 1e-10;
};

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
    using Controls = ExactSchurControls;

    // Solve Auu * y = rhs, filling y. Return false on failure.
    // The callback owns any preconditioner/hierarchy state; it must be
    // consistent with the Auu block passed at setup().
    using AuuSolve = std::function<bool(const Vector& rhs, Vector& y)>;

    explicit ExactSchurApproximation(AuuSolve auu_solve,
                                     const Controls& controls = Controls{},
                                     std::optional<NullSpaceProjector> pressure_null_space = std::nullopt)
        : auu_solve_(std::move(auu_solve)),
          controls_(controls),
          pressure_null_space_(std::move(pressure_null_space)) {}

    const char* name() const noexcept override { return "exact_schur"; }

    // apply() runs CG to solve S p = rhs_p, so it is the inverse action.
    SchurAction action() const noexcept override { return SchurAction::InverseOperator; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || !auu_solve_) return false;
        if (pressure_null_space_ && pressure_null_space_->dimension() != blocks.pressure_size())
            return false;
        blocks_ = &blocks;
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        if (!blocks_ || graph_signature(blocks) != graph_signature_)
            return false;
        blocks_ = &blocks;
        return true;
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (!blocks_ || !auu_solve_) return false;
        const std::size_t n = blocks_->pressure_size();
        if (rhs_p.size() != n) return false;
        if (pressure_null_space_ && !pressure_null_space_->is_compatible(rhs_p)) return false;

        if (pressure.size() != n) pressure = Vector(n, 0.0);
        else pressure.fill(0.0);

        // Preconditioner-free CG on the implicit Schur operator.
        Vector r(n), p(n), Ap(n);
        for (std::size_t i = 0; i < n; ++i) r(i) = rhs_p(i);
        if (pressure_null_space_) pressure_null_space_->remove(r);
        for (std::size_t i = 0; i < n; ++i) p(i) = r(i);
        double r_dot_r = r.dot(r);
        const double r0 = std::sqrt(r_dot_r);
        if (!(r0 > 0.0)) return true;  // rhs is zero, pressure=0 is the solution

        const double abs_tol = controls_.tolerance * r0;

        for (std::size_t it = 0; it < controls_.max_iter; ++it) {
            if (pressure_null_space_) pressure_null_space_->remove(p);
            if (!apply_schur(p, Ap)) return false;
            if (pressure_null_space_) pressure_null_space_->remove(Ap);

            double pAp = 0.0;
            for (std::size_t i = 0; i < n; ++i) pAp += p(i) * Ap(i);
            if (!(pAp > 0.0) || !std::isfinite(pAp)) return false;

            const double alpha = r_dot_r / pAp;
            for (std::size_t i = 0; i < n; ++i) pressure(i) += alpha * p(i);
            for (std::size_t i = 0; i < n; ++i) r(i) -= alpha * Ap(i);
            if (pressure_null_space_) pressure_null_space_->remove(r);

            const double new_r_dot_r = r.dot(r);
            if (std::sqrt(new_r_dot_r) < abs_tol) {
                if (pressure_null_space_) pressure_null_space_->remove(pressure);
                return true;
            }

            const double beta = new_r_dot_r / r_dot_r;
            for (std::size_t i = 0; i < n; ++i) p(i) = r(i) + beta * p(i);
            if (pressure_null_space_) pressure_null_space_->remove(p);
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
        Vector projected_p = p;
        if (pressure_null_space_) pressure_null_space_->remove(projected_p);
        if (out.size() != np) out = Vector(np, 0.0);

        // Gp = G * p
        const auto Gp = blocks_->G().matvec(projected_p);
        Vector Gp_v(nu, 0.0);
        for (std::size_t i = 0; i < nu; ++i) Gp_v(i) = Gp[i];

        // y = Auu^{-1} * Gp
        Vector y(nu, 0.0);
        if (!auu_solve_(Gp_v, y)) return false;

        // out = C*p - D*y
        const auto Cp = blocks_->C().matvec(projected_p);
        const auto Dy = blocks_->D().matvec(y);
        for (std::size_t i = 0; i < np; ++i) out(i) = Cp[i] - Dy[i];
        if (pressure_null_space_) pressure_null_space_->remove(out);
        return true;
    }

private:
    struct GraphSignature {
        std::size_t hash = 0;
        bool operator==(const GraphSignature& other) const noexcept {
            return hash == other.hash;
        }
        bool operator!=(const GraphSignature& other) const noexcept {
            return !(*this == other);
        }
    };

    static GraphSignature graph_signature(const BlockOperator& blocks) {
        std::size_t h = 1469598103934665603ULL;
        const auto mix = [&h](std::size_t value) {
            h ^= value;
            h *= 1099511628211ULL;
        };
        const auto add = [&mix](const SparseMatrix& A) {
            mix(A.n_rows());
            mix(A.n_cols());
            mix(A.nnz());
            for (std::size_t i = 0; i < A.n_rows() + 1; ++i)
                mix(A.row_offsets_data()[i]);
            for (std::size_t k = 0; k < A.nnz(); ++k)
                mix(A.columns_data()[k]);
        };
        add(blocks.Auu());
        add(blocks.G());
        add(blocks.D());
        add(blocks.C());
        return GraphSignature{h};
    }

    const BlockOperator* blocks_ = nullptr;
    AuuSolve auu_solve_;
    Controls controls_;
    std::optional<NullSpaceProjector> pressure_null_space_;
    GraphSignature graph_signature_{};
};

} // namespace cfdx::core
