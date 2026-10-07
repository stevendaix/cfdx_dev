#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <map>
#include <optional>
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

    // Solve H*z = rhs, where H = D Q^{-1} G - C for the stabilized setup.
    using PressureSolve = std::function<bool(const Vector& rhs, Vector& z)>;

    LscBfbtSchurApproximation(Mode mode,
                              PressureSolve pressure_solve,
                              std::vector<double> q_diagonal = {},
                              std::optional<NullSpaceProjector> pressure_null_space = std::nullopt)
        : mode_(mode),
          pressure_solve_(std::move(pressure_solve)),
          q_diagonal_(std::move(q_diagonal)),
          pressure_null_space_(std::move(pressure_null_space)) {}

    const char* name() const noexcept override {
        return mode_ == Mode::LSC ? "lsc_schur" : "bfbt_schur";
    }

    // apply() returns -H^{-1} K H^{-1} r, which approximates S^{-1} r.
    SchurAction action() const noexcept override { return SchurAction::InverseOperator; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || !pressure_solve_) return false;
        const std::size_t nu = blocks.velocity_size();
        if (nu == 0 || blocks.pressure_size() == 0) return false;
        if (pressure_null_space_ && pressure_null_space_->dimension() != blocks.pressure_size())
            return false;

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
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        if (!blocks_ || graph_signature(blocks) != graph_signature_)
            return false;
        return setup(blocks);
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (!blocks_ || !pressure_solve_) return false;
        const std::size_t np = blocks_->pressure_size();
        if (rhs_p.size() != np) return false;

        Vector projected_rhs = rhs_p;
        if (pressure_null_space_) {
            if (!pressure_null_space_->is_compatible(projected_rhs))
                return false;
            pressure_null_space_->remove(projected_rhs);
        }

        Vector y(np, 0.0);
        if (!pressure_solve_(projected_rhs, y)) return false;
        if (pressure_null_space_) pressure_null_space_->remove(y);

        // E*y = D Q^-1 Auu Q^-1 G*y.
        const auto Gy = blocks_->G().matvec(y);
        Vector qg(Gy.size(), 0.0);
        for (std::size_t i = 0; i < Gy.size(); ++i) qg(i) = q_inverse_[i] * Gy[i];

        const auto Aqg = blocks_->Auu().matvec(qg);
        Vector qAqg(Aqg.size(), 0.0);
        for (std::size_t i = 0; i < Aqg.size(); ++i) qAqg(i) = q_inverse_[i] * Aqg[i];

        const auto Eq = blocks_->D().matvec(qAqg);
        const auto Cy = blocks_->C().matvec(y);
        Vector z(np, 0.0);
        for (std::size_t i = 0; i < np; ++i) z(i) = Eq[i] - Cy[i];
        if (pressure_null_space_) pressure_null_space_->remove(z);

        if (!pressure_solve_(z, pressure)) return false;
        if (pressure_null_space_) pressure_null_space_->remove(pressure);
        for (std::size_t i = 0; i < np; ++i) pressure(i) = -pressure(i);
        if (pressure_null_space_) pressure_null_space_->remove(pressure);
        return true;
    }

    // Assemble the classical LSC/BFBT component P = D Q^{-1} G.  It is
    // retained as a diagnostic/reference operator.
    static SparseMatrix assemble_pressure_operator(
        const BlockOperator& blocks,
        Mode mode,
        const std::vector<double>& q_diagonal = {}) {
        return assemble_scaled_divergence_gradient(blocks, mode, q_diagonal);
    }

    // Assemble the stabilized pressure-side operator
    //
    //     H = D Q^{-1} G - C.
    //
    // The physical Schur approximation is S_hat = -H.  Keeping H as the
    // pressure solve operator gives NativeAMG the same sign-normalized,
    // positive-diagonal system used by the existing CFDX pressure solvers.
    static SparseMatrix assemble_stabilized_pressure_operator(
        const BlockOperator& blocks,
        Mode mode,
        const std::vector<double>& q_diagonal = {}) {
        SparseMatrix H = assemble_scaled_divergence_gradient(blocks, mode, q_diagonal);
        const auto& C = blocks.C();
        for (std::size_t row = 0; row < C.n_rows(); ++row) {
            for (std::size_t k = C.row_offsets_data()[row];
                 k < C.row_offsets_data()[row + 1]; ++k) {
                H.push_back(row, C.columns_data()[k], -C.values_data()[k]);
            }
        }
        H.finalize();
        return H;
    }

    Mode mode() const noexcept { return mode_; }

    bool uses_default_scaling() const noexcept { return q_diagonal_.empty(); }

    bool has_pressure_null_space_policy() const noexcept {
        return pressure_null_space_.has_value();
    }

private:
    static SparseMatrix assemble_scaled_divergence_gradient(
        const BlockOperator& blocks,
        Mode mode,
        const std::vector<double>& q_diagonal = {}) {
        const auto q_inverse = make_q_inverse(blocks, mode, q_diagonal);
        const std::size_t nu = blocks.velocity_size();
        const auto& D = blocks.D();
        const auto& G = blocks.G();
        SparseMatrix P(blocks.pressure_size(), blocks.pressure_size());
        for (std::size_t row = 0; row < blocks.pressure_size(); ++row) {
            std::map<std::size_t, double> entries;
            for (std::size_t dk = D.row_offsets_data()[row];
                 dk < D.row_offsets_data()[row + 1]; ++dk) {
                const std::size_t velocity_col = D.columns_data()[dk];
                if (velocity_col >= nu)
                    throw std::invalid_argument("LSC/BFBT D column is outside velocity space");
                const double d_value = D.values_data()[dk] * q_inverse[velocity_col];
                for (std::size_t gk = G.row_offsets_data()[velocity_col];
                     gk < G.row_offsets_data()[velocity_col + 1]; ++gk) {
                    entries[G.columns_data()[gk]] += d_value * G.values_data()[gk];
                }
            }
            for (const auto& [col, value] : entries) {
                if (std::isfinite(value) && value != 0.0)
                    P.push_back(row, col, value);
            }
        }
        P.finalize();
        return P;
    }

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

    Mode mode_;
    const BlockOperator* blocks_ = nullptr;
    PressureSolve pressure_solve_;
    std::vector<double> q_diagonal_;
    std::vector<double> q_inverse_;
    std::optional<NullSpaceProjector> pressure_null_space_;
    GraphSignature graph_signature_{};
};

using LeastSquaresCommutatorSchurApproximation = LscBfbtSchurApproximation;

} // namespace cfdx::core
