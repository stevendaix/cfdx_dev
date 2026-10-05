#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <stdexcept>
#include <utility>

namespace cfdx::core {

// Pressure-convection-diffusion (PCD) approximation of the pressure Schur
// inverse for the coupled operator
//
//     A = [ Auu  G ]
//         [ D    C ].
//
// The Schur complement is
//
//     S = C - D Auu^{-1} G.
//
// The action implemented here is a single pressure inverse,
//
//     S^{-1} rhs ~= -Fp^{-1} rhs,
//
// where Fp is the pressure convection-diffusion operator. The single inverse is
// not a stylistic choice. On a uniform mesh the pressure mass operator is the
// uniform scalar diag(cell_volumes), so the previously documented composition
// -Fp^{-1} Mp Kp^{-1} reduced to a uniform factor times Fp^{-1} Kp^{-1}: two
// inverse pressure operators where the Schur complement requires one. Measured
// on the production 512-unknown coupled matrix, that squared inverse produced
// an action of norm 51510.5 against an exact-Schur action of norm 0.70831, with
// a cosine of -0.029 between them. Every double-inverse ordering fails to
// converge that system; the single-inverse orderings converge it in the same 17
// iterations as BlockSchur and MGR.
//
// The choice between Fp and the pure-diffusion operator Kp as the single inverse
// is PROVISIONAL. Both converge the production acceptance case because the mass
// flux is zero at its first coupled iteration, which leaves Fp and Kp nearly
// proportional. That near-equality is a property of that case, not evidence
// that the ordering is immaterial: with a non-zero mass flux Fp and Kp differ,
// and the ordering must be qualified on such a case before it is fixed.
//
// The pressure-side operators are supplied explicitly. This is deliberate:
// constructing Fp from the velocity block would silently encode a
// discretisation-dependent assumption and would make the approximation less
// auditable. Pressure solvers are also supplied independently for Kp and Fp.
// Only Fp enters the current action; Mp and Kp are retained and validated
// because the pending ordering qualification compares them, and because the
// lifecycle contract covers all three pressure operators.
//
// Gauge treatment: the reference cell is imposed on the pressure operators
// themselves, by replacing that operator row with a unit row. The operators are
// therefore nonsingular and their inverse action is defined for any right-hand
// side. The action must consequently return a gauge component when its input has
// one. Zeroing the reference cell on input or output freezes the coupled gauge
// row residual at its initial magnitude forever: on the production case that
// left a residual of exactly the initial pressure value (0.5) as the dominant
// term of a total residual norm of 0.500247. This class therefore does not pin
// the reference cell.
//
// The pressure-side matrices must keep their sparsity patterns between
// setup() and update_values(). A changed graph requires setup() so that the
// lifecycle cannot silently reuse stale symbolic state.
class PcdSchurApproximation final : public SchurApproximation {
public:
    using PressureSolve = std::function<bool(const Vector& rhs, Vector& x)>;

    PcdSchurApproximation(const SparseMatrix& pressure_mass,
                          const SparseMatrix& pressure_laplacian,
                          const SparseMatrix& pressure_convection_diffusion,
                          PressureSolve laplacian_solve,
                          PressureSolve convection_diffusion_solve,
                          std::optional<NullSpaceProjector> pressure_null_space = std::nullopt,
                          std::optional<std::size_t> pressure_reference_cell = std::nullopt)
        : pressure_mass_(pressure_mass),
          pressure_laplacian_(pressure_laplacian),
          pressure_convection_diffusion_(pressure_convection_diffusion),
          laplacian_solve_(std::move(laplacian_solve)),
          convection_diffusion_solve_(std::move(convection_diffusion_solve)),
          pressure_null_space_(std::move(pressure_null_space)),
          pressure_reference_cell_(pressure_reference_cell) {}

    const char* name() const noexcept override { return "pcd_schur"; }

    // apply() returns -Fp^{-1} rhs, which approximates S^{-1} rhs for the
    // pressure Schur complement. The sign belongs to the inverse action, so this
    // is InverseOperator rather than the SIMPLE/SIMPLEC operator action.
    SchurAction action() const noexcept override { return SchurAction::InverseOperator; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || !laplacian_solve_ || !convection_diffusion_solve_)
            return false;

        const std::size_t np = blocks.pressure_size();
        if (np == 0 || pressure_mass_ == nullptr || pressure_laplacian_ == nullptr ||
            pressure_convection_diffusion_ == nullptr)
            return false;

        if (!valid_pressure_operator(pressure_mass_, np) ||
            !valid_pressure_operator(pressure_laplacian_, np) ||
            !valid_pressure_operator(pressure_convection_diffusion_, np))
            return false;

        if (pressure_null_space_ &&
            pressure_null_space_->dimension() != np)
            return false;
        if (pressure_reference_cell_ &&
            *pressure_reference_cell_ >= np)
            return false;

        pressure_size_ = np;
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        if (pressure_size_ == 0 || graph_signature(blocks) != graph_signature_)
            return false;
        const std::size_t np = blocks.pressure_size();
        if (!blocks.is_valid() || np == 0 ||
            !laplacian_solve_ || !convection_diffusion_solve_ ||
            !valid_pressure_operator(*pressure_mass_, np) ||
            !valid_pressure_operator(*pressure_laplacian_, np) ||
            !valid_pressure_operator(*pressure_convection_diffusion_, np))
            return false;
        if (pressure_null_space_ &&
            pressure_null_space_->dimension() != np)
            return false;
        // Keep blocks_ bound to the object supplied during setup(). In
        // particular, do not bind it to a temporary BlockOperator created
        // by a numeric-update caller.
        return true;
    }

    // Refresh the three pressure-side operators explicitly. The common
    // SchurApproximation interface only carries the saddle-point blocks, so
    // this overload makes the additional PCD state update auditable rather
    // than pretending that those operators are part of BlockOperator.
    bool update_pressure_values(const BlockOperator& blocks,
                                const SparseMatrix& pressure_mass,
                                const SparseMatrix& pressure_laplacian,
                                const SparseMatrix& pressure_convection_diffusion) {
        if (pressure_size_ == 0 || graph_signature(blocks) != graph_signature_)
            return false;

        const std::size_t np = blocks.pressure_size();
        if (!valid_pressure_operator(pressure_mass, np) ||
            !valid_pressure_operator(pressure_laplacian, np) ||
            !valid_pressure_operator(pressure_convection_diffusion, np))
            return false;

        if (!same_pattern(pressure_mass_, pressure_mass) ||
            !same_pattern(pressure_laplacian_, pressure_laplacian) ||
            !same_pattern(pressure_convection_diffusion_,
                           pressure_convection_diffusion))
            return false;

        pressure_mass_ = pressure_mass;
        pressure_laplacian_ = pressure_laplacian;
        pressure_convection_diffusion_ = pressure_convection_diffusion;
        pressure_size_ = np;
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (pressure_size_ == 0 || !laplacian_solve_ ||
            !convection_diffusion_solve_)
            return false;

        const std::size_t np = pressure_size_;
        if (rhs_p.size() != np)
            return false;

        // The reference cell is deliberately not pinned on input or output. The
        // caller imposes it on the pressure operators themselves, so they are
        // nonsingular and their inverse action is defined for any right-hand
        // side. Returning a zero gauge component here would leave the coupled
        // gauge row residual frozen at its initial magnitude, because no Krylov
        // vector could ever carry a gauge component to correct it.
        Vector rhs = rhs_p;
        if (pressure_null_space_) {
            if (!pressure_null_space_->is_compatible(rhs))
                return false;
            pressure_null_space_->remove(rhs);
        }

        // pressure = -Fp^{-1} rhs. Exactly one pressure inverse is applied; see
        // the contract note at the top of this header for why the mass and
        // Laplacian operators must not be composed into the action.
        if (!convection_diffusion_solve_(rhs, pressure) ||
            pressure.size() != np || !pressure.is_valid())
            return false;
        if (pressure_null_space_)
            pressure_null_space_->remove(pressure);

        pressure *= -1.0;
        return true;
    }

    bool has_pressure_null_space_policy() const noexcept {
        return pressure_null_space_.has_value();
    }

    // PCD owns snapshots of the pressure-side matrices, so matrix lifetime is
    // independent of the caller after construction or numeric update. The solve
    // callbacks are owned std::function values; any solver state captured by a
    // callback remains the callback's responsibility.
    const SparseMatrix& pressure_mass() const {
        return pressure_mass_;
    }

    const SparseMatrix& pressure_laplacian() const {
        return pressure_laplacian_;
    }

    const SparseMatrix& pressure_convection_diffusion() const {
        return pressure_convection_diffusion_;
    }

    bool has_pressure_operators() const noexcept {
        return pressure_size_ != 0;
    }

private:
    static bool valid_pressure_operator(const SparseMatrix& matrix,
                                        std::size_t np) {
        if (matrix.n_rows() != np || matrix.n_cols() != np)
            return false;
        for (std::size_t k = 0; k < matrix.nnz(); ++k) {
            if (!std::isfinite(matrix.values_data()[k]))
                return false;
        }
        return true;
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

    static bool same_pattern(const SparseMatrix& a, const SparseMatrix& b) noexcept {
        if (a.n_rows() != b.n_rows() || a.n_cols() != b.n_cols() ||
            a.nnz() != b.nnz())
            return false;
        for (std::size_t i = 0; i <= a.n_rows(); ++i)
            if (a.row_offsets_data()[i] != b.row_offsets_data()[i])
                return false;
        for (std::size_t k = 0; k < a.nnz(); ++k)
            if (a.columns_data()[k] != b.columns_data()[k])
                return false;
        return true;
    }

    GraphSignature graph_signature(const BlockOperator& blocks) const {
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
        add(pressure_mass_);
        add(pressure_laplacian_);
        add(pressure_convection_diffusion_);
        return GraphSignature{h};
    }

    std::size_t pressure_size_ = 0;
    SparseMatrix pressure_mass_;
    SparseMatrix pressure_laplacian_;
    SparseMatrix pressure_convection_diffusion_;
    PressureSolve laplacian_solve_;
    PressureSolve convection_diffusion_solve_;
    std::optional<NullSpaceProjector> pressure_null_space_;
    std::optional<std::size_t> pressure_reference_cell_;
    GraphSignature graph_signature_{};
};

} // namespace cfdx::core
