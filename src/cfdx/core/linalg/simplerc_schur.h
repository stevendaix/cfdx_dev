#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <map>
#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cfdx::core {

enum class SimplerSchurMode {
    SIMPLE,   // S~ = C - D diag(Auu)^{-1} G
    SIMPLEC   // S~ = C - D (diag(Auu) - offdiag(Auu))^{-1} G  (consistent)
};

inline const char* to_string(SimplerSchurMode m) {
    return m == SimplerSchurMode::SIMPLE ? "simple" : "simplec";
}

// Algebraic SIMPLE / SIMPLEC Schur-complement approximation.
//
// The exact Schur complement S = C - D Auu^{-1} G is approximated by replacing
// Auu^{-1} with a cheap diagonal operator:
//   SIMPLE  : Ad = diag(Auu)
//   SIMPLEC : Ad = diag(Auu) - row_offdiag(Auu)      (the "consistent" form)
// S~ = C - D Ad^{-1} G is built implicitly (one D matvec, one G matvec) so the
// only cost per application is two sparse matvecs and one diagonal division.
// This is the standard pressure-preconditioner approximation used by the
// SIMPLE/SIMPLEC pressure-correction algorithms and is intended to be
// qualified against the ExactSchur oracle (#492).
class SimplerSchurApproximation final : public SchurApproximation {
public:
    // The optional pressure null-space projector is the same mean-zero policy the
    // exact and LSC/BFBt approximations accept. It is opt-in so an existing caller
    // keeps its previous algebra exactly.
    explicit SimplerSchurApproximation(
        SimplerSchurMode mode,
        std::optional<NullSpaceProjector> pressure_null_space = std::nullopt)
        : mode_(mode), pressure_null_space_(std::move(pressure_null_space)) {}

    const char* name() const noexcept override {
        return mode_ == SimplerSchurMode::SIMPLE ? "simple_schur" : "simplec_schur";
    }

    // apply() returns C rhs_p - D Ad^{-1} G rhs_p = S~ rhs_p: the approximated
    // Schur operator itself, not its inverse. This differs from the exact and
    // LSC/BFBt approximations, which apply an inverse action, and is why
    // SchurApproximation::action() is part of the contract.
    SchurAction action() const noexcept override { return SchurAction::Operator; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid()) return false;
        const SparseMatrix& auu = blocks.Auu();
        const std::size_t n = auu.n_rows();
        if (n == 0) return false;
        if (pressure_null_space_ &&
            pressure_null_space_->dimension() != blocks.pressure_size())
            return false;

        diagonal_.assign(n, 0.0);
        denominator_.assign(n, 0.0);
        offdiag_norm_ = 0.0;

        for (std::size_t i = 0; i < n; ++i) {
            double diag = 0.0;
            double row_off = 0.0;
            for (std::size_t k = auu.row_offsets_data()[i];
                 k < auu.row_offsets_data()[i + 1]; ++k) {
                const std::size_t j = auu.columns_data()[k];
                const double v = auu.values_data()[k];
                if (j == i) diag += v;
                else row_off += std::abs(v);
            }
            if (!(diag > 0.0) || !std::isfinite(diag)) return false;
            diagonal_[i] = diag;
            offdiag_norm_ = std::max(offdiag_norm_, row_off);

            const double ad = (mode_ == SimplerSchurMode::SIMPLEC)
                ? diag - row_off
                : diag;
            if (!(ad > 1e-14) || !std::isfinite(ad)) return false;
            denominator_[i] = ad;
        }
        blocks_ = &blocks;
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    // A CSR graph change invalidates the diagonal and offdiagonal sums captured
    // by setup, so it is rejected rather than silently recomputed from a different
    // sparsity pattern. This matches the exact and LSC/BFBt approximations: a
    // numeric refresh is accepted, a graph change requires an explicit setup.
    bool update_values(const BlockOperator& blocks) override {
        if (!blocks_ || graph_signature(blocks) != graph_signature_)
            return false;
        return setup(blocks);
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (!blocks_) return false;
        const std::size_t np = blocks_->pressure_size();
        const std::size_t nu = blocks_->velocity_size();
        if (rhs_p.size() != np) return false;
        if (pressure.size() != np)
            pressure = Vector(np, 0.0);

        Vector projected_rhs = rhs_p;
        if (pressure_null_space_) {
            if (!pressure_null_space_->is_compatible(projected_rhs))
                return false;
            pressure_null_space_->remove(projected_rhs);
        }

        // Gp = G * rhs_p
        const auto Gp = blocks_->G().matvec(projected_rhs);
        Vector Gp_v(nu, 0.0);
        for (std::size_t i = 0; i < nu; ++i) Gp_v(i) = Gp[i];

        // y = Ad^{-1} Gp
        Vector y(nu, 0.0);
        for (std::size_t i = 0; i < nu && i < denominator_.size(); ++i)
            y(i) = Gp_v(i) / denominator_[i];

        // out = C*rhs_p - D*y
        const auto Cp = blocks_->C().matvec(projected_rhs);
        const auto Dy = blocks_->D().matvec(y);
        for (std::size_t i = 0; i < np; ++i)
            pressure(i) = Cp[i] - Dy[i];
        if (pressure_null_space_) pressure_null_space_->remove(pressure);
        return true;
    }

    // Assemble the algebraic SIMPLE/SIMPLEC Schur operator
    // S~=C-D A_d^{-1} G for production pressure AMG. The matrix is built
    // from the same coefficients used by apply(), so the preconditioner never
    // treats the operator action as if it were an inverse action.
    SparseMatrix assembled_operator() const {
        if (!blocks_) return SparseMatrix();
        const std::size_t np = blocks_->pressure_size();
        SparseMatrix out(np, np);
        for (std::size_t row = 0; row < np; ++row) {
            std::map<std::size_t, double> entries;
            for (std::size_t k = blocks_->C().row_offsets_data()[row];
                 k < blocks_->C().row_offsets_data()[row + 1]; ++k)
                entries[blocks_->C().columns_data()[k]] += blocks_->C().values_data()[k];
            for (std::size_t dk = blocks_->D().row_offsets_data()[row];
                 dk < blocks_->D().row_offsets_data()[row + 1]; ++dk) {
                const std::size_t u = blocks_->D().columns_data()[dk];
                const double d = blocks_->D().values_data()[dk] / denominator_[u];
                for (std::size_t gk = blocks_->G().row_offsets_data()[u];
                     gk < blocks_->G().row_offsets_data()[u + 1]; ++gk)
                    entries[blocks_->G().columns_data()[gk]] -=
                        d * blocks_->G().values_data()[gk];
            }
            for (const auto& [col, value] : entries) {
                if (!std::isfinite(value)) return SparseMatrix();
                if (value != 0.0) out.push_back(row, col, value);
            }
        }
        out.finalize();
        return out;
    }

    double offdiag_norm() const noexcept { return offdiag_norm_; }

    bool has_pressure_null_space_policy() const noexcept {
        return pressure_null_space_.has_value();
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

    SimplerSchurMode mode_;
    std::optional<NullSpaceProjector> pressure_null_space_;
    const BlockOperator* blocks_ = nullptr;
    std::vector<double> diagonal_;
    std::vector<double> denominator_;
    double offdiag_norm_ = 0.0;
    GraphSignature graph_signature_{};
};

} // namespace cfdx::core