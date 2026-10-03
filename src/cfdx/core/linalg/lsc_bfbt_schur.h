#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::core {

// Algebraic pressure-Schur inverse approximations for the coupled operator
//
//     A = [ Auu  G ]
//         [ D    C ]
//
// For incompressible problems the pressure Schur complement is commonly
// dominated by S = -D Auu^{-1} G.  LSC/BFBt approximate the inverse action
// without forming Auu^{-1}.  With a positive diagonal velocity scaling Q,
// define
//
//     P = D Q^{-1} G
//     E = D Q^{-1} Auu Q^{-1} G
//
// and use
//
//     S^{-1} r ~= - P^{-1} E P^{-1} r.
//
// This is the LSC/BFBt algebraic form documented by PETSc and the CFD
// literature.  The pressure solve P^{-1} is supplied by the caller so that
// CFDX can use AMG/Krylov without coupling the Schur approximation to a
// particular pressure solver.
//
// LSC uses a momentum-based diagonal scaling Q = diag(Auu) by default.
// BFBt uses an explicitly supplied velocity-space scaling; the identity
// scaling is the classical unscaled BFBt form.
//
// These classes deliberately do not silently assume G = -D^T and do not
// claim that C is represented: as in the standard LSC/BFBt construction,
// A11/C is outside this approximation.  Production qualification therefore
// requires a pressure-side solve and a benchmark campaign in addition to
// these algebraic operator tests.
class LscBfbtSchurApproximation final : public SchurApproximation {
public:
    enum class Mode {
        LSC,
        BFBT
    };

    // Solve P*z = rhs, where P = D Q^{-1} G for the current setup.
    using PressureSolve = std::function<bool(const Vector& rhs, Vector& z)>;

    LscBfbtSchurApproximation(Mode mode,
                              PressureSolve pressure_solve,
                              std::vector<double> q_diagonal = {})
        : mode_(mode),
          pressure_solve_(std::move(pressure_solve)),
          q_diagonal_(std::move(q_diagonal)) {}

    const char* name() const noexcept override {
        return mode_ == Mode::LSC ? "lsc_schur" : "bfbt_schur";
    }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || !pressure_solve_) return false;
        const std::size_t nu = blocks.velocity_size();
        if (nu == 0 || blocks.pressure_size() == 0) return false;

        q_inverse_.assign(nu, 1.0);
        if (q_diagonal_.empty()) {
            if (mode_ == Mode::LSC) {
                const auto& A = blocks.Auu();
                for (std::size_t i = 0; i < nu; ++i) {
                    double d = 0.0;
                    for (std::size_t k = A.row_offsets_data()[i];
                         k < A.row_offsets_data()[i + 1]; ++k) {
                        if (A.columns_data()[k] == i) d += A.values_data()[k];
                    }
                    if (!(d > 0.0) || !std::isfinite(d)) return false;
                    q_inverse_[i] = 1.0 / d;
                }
            }
        } else {
            if (q_diagonal_.size() != nu) return false;
            for (std::size_t i = 0; i < nu; ++i) {
                if (!(q_diagonal_[i] > 0.0) || !std::isfinite(q_diagonal_[i]))
                    return false;
                q_inverse_[i] = 1.0 / q_diagonal_[i];
            }
        }

        blocks_ = &blocks;
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        // A numerical refresh is valid only for the same algebraic dimensions.
        // The caller owns CSR graph-generation checks; a different graph must
        // use setup() explicitly rather than being silently rebuilt here.
        if (!blocks_ || blocks.velocity_size() != blocks_->velocity_size() ||
            blocks.pressure_size() != blocks_->pressure_size())
            return false;
        return setup(blocks);
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (!blocks_ || !pressure_solve_) return false;
        const std::size_t np = blocks_->pressure_size();
        if (rhs_p.size() != np) return false;

        Vector y(np, 0.0);
        if (!pressure_solve_(rhs_p, y)) return false;

        // E*y = D Q^-1 Auu Q^-1 G*y.
        const auto Gy = blocks_->G().matvec(y);
        Vector qg(Gy.size(), 0.0);
        for (std::size_t i = 0; i < Gy.size(); ++i) qg(i) = q_inverse_[i] * Gy[i];

        const auto Aqg = blocks_->Auu().matvec(qg);
        Vector qAqg(Aqg.size(), 0.0);
        for (std::size_t i = 0; i < Aqg.size(); ++i) qAqg(i) = q_inverse_[i] * Aqg[i];

        const auto Eq = blocks_->D().matvec(qAqg);
        Vector z(np, 0.0);
        for (std::size_t i = 0; i < np; ++i) z(i) = Eq[i];

        if (!pressure_solve_(z, pressure)) return false;
        for (std::size_t i = 0; i < np; ++i) pressure(i) = -pressure(i);
        return true;
    }

    Mode mode() const noexcept { return mode_; }

    bool uses_default_scaling() const noexcept { return q_diagonal_.empty(); }

private:
    Mode mode_;
    const BlockOperator* blocks_ = nullptr;
    PressureSolve pressure_solve_;
    std::vector<double> q_diagonal_;
    std::vector<double> q_inverse_;
};

using LeastSquaresCommutatorSchurApproximation = LscBfbtSchurApproximation;

} // namespace cfdx::core
