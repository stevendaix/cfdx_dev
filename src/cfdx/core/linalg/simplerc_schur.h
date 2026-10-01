#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstddef>
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
    explicit SimplerSchurApproximation(SimplerSchurMode mode) : mode_(mode) {}

    const char* name() const noexcept override {
        return mode_ == SimplerSchurMode::SIMPLE ? "simple_schur" : "simplec_schur";
    }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid()) return false;
        const SparseMatrix& auu = blocks.Auu();
        const std::size_t n = auu.n_rows();
        if (n == 0) return false;

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
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        return setup(blocks);
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (!blocks_) return false;
        const std::size_t np = blocks_->pressure_size();
        const std::size_t nu = blocks_->velocity_size();
        if (rhs_p.size() != np) return false;
        if (pressure.size() != np)
            pressure = Vector(np, 0.0);

        // Gp = G * rhs_p
        const auto Gp = blocks_->G().matvec(rhs_p);
        Vector Gp_v(nu, 0.0);
        for (std::size_t i = 0; i < nu; ++i) Gp_v(i) = Gp[i];

        // y = Ad^{-1} Gp
        Vector y(nu, 0.0);
        for (std::size_t i = 0; i < nu && i < denominator_.size(); ++i)
            y(i) = Gp_v(i) / denominator_[i];

        // out = C*rhs_p - D*y
        const auto Cp = blocks_->C().matvec(rhs_p);
        const auto Dy = blocks_->D().matvec(y);
        for (std::size_t i = 0; i < np; ++i)
            pressure(i) = Cp[i] - Dy[i];
        return true;
    }

    double offdiag_norm() const noexcept { return offdiag_norm_; }

private:
    SimplerSchurMode mode_;
    const BlockOperator* blocks_ = nullptr;
    std::vector<double> diagonal_;
    std::vector<double> denominator_;
    double offdiag_norm_ = 0.0;
};

} // namespace cfdx::core