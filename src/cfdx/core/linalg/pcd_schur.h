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
// For incompressible flow, the pressure Schur complement is commonly
//
//     S = C - D Auu^{-1} G ~= -Kp Mp^{-1} Fp,
//
// where Mp is the pressure mass operator, Kp the pressure Laplacian/diffusion
// operator and Fp the pressure convection-diffusion operator. Consequently
//
//     S^{-1} ~= -Fp^{-1} Mp Kp^{-1}.
//
// The three pressure-side operators are supplied explicitly. This is
// deliberate: constructing Fp from the velocity block would silently encode
// a discretisation-dependent assumption and would make the approximation less
// auditable. Pressure solvers are also supplied independently for Kp and Fp.
//
// This class implements the algebraic PCD action only. It does not claim that
// a particular pressure mass, diffusion or convection discretisation is
// universally appropriate; those choices belong to the caller/case.
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
        : pressure_mass_(&pressure_mass),
          pressure_laplacian_(&pressure_laplacian),
          pressure_convection_diffusion_(&pressure_convection_diffusion),
          laplacian_solve_(std::move(laplacian_solve)),
          convection_diffusion_solve_(std::move(convection_diffusion_solve)),
          pressure_null_space_(std::move(pressure_null_space)),
          pressure_reference_cell_(pressure_reference_cell) {}

    const char* name() const noexcept override { return "pcd_schur"; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || !laplacian_solve_ || !convection_diffusion_solve_)
            return false;

        const std::size_t np = blocks.pressure_size();
        if (np == 0 || pressure_mass_ == nullptr || pressure_laplacian_ == nullptr ||
            pressure_convection_diffusion_ == nullptr)
            return false;

        if (!valid_pressure_operator(*pressure_mass_, np) ||
            !valid_pressure_operator(*pressure_laplacian_, np) ||
            !valid_pressure_operator(*pressure_convection_diffusion_, np))
            return false;

        if (pressure_null_space_ &&
            pressure_null_space_->dimension() != np)
            return false;
        if (pressure_reference_cell_ &&
            *pressure_reference_cell_ >= np)
            return false;

        blocks_ = &blocks;
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    bool update_values(const BlockOperator& blocks) override {
        if (!blocks_ || graph_signature(blocks) != graph_signature_)
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
        if (!blocks_ || graph_signature(blocks) != graph_signature_)
            return false;

        const std::size_t np = blocks.pressure_size();
        if (!valid_pressure_operator(pressure_mass, np) ||
            !valid_pressure_operator(pressure_laplacian, np) ||
            !valid_pressure_operator(pressure_convection_diffusion, np))
            return false;

        if (!same_pattern(*pressure_mass_, pressure_mass) ||
            !same_pattern(*pressure_laplacian_, pressure_laplacian) ||
            !same_pattern(*pressure_convection_diffusion_,
                           pressure_convection_diffusion))
            return false;

        pressure_mass_ = &pressure_mass;
        pressure_laplacian_ = &pressure_laplacian;
        pressure_convection_diffusion_ = &pressure_convection_diffusion;
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    bool apply(const Vector& rhs_p, Vector& pressure) const override {
        if (!blocks_ || !pressure_mass_ || !pressure_laplacian_ ||
            !pressure_convection_diffusion_ || !laplacian_solve_ ||
            !convection_diffusion_solve_)
            return false;

        const std::size_t np = blocks_->pressure_size();
        if (rhs_p.size() != np)
            return false;

        Vector rhs = rhs_p;
        if (pressure_reference_cell_)
            rhs(*pressure_reference_cell_) = 0.0;
        if (pressure_null_space_) {
            if (!pressure_null_space_->is_compatible(rhs))
                return false;
            pressure_null_space_->remove(rhs);
        }

        // z = Kp^{-1} rhs.
        Vector z(np, 0.0);
        if (!laplacian_solve_(rhs, z) || z.size() != np || !z.is_valid())
            return false;
        if (pressure_null_space_)
            pressure_null_space_->remove(z);

        // y = Mp z.
        Vector y = pressure_mass_->matvec(z);
        if (y.size() != np || !y.is_valid())
            return false;
        if (pressure_null_space_)
            pressure_null_space_->remove(y);

        // pressure = -Fp^{-1} y.
        if (!convection_diffusion_solve_(y, pressure) ||
            pressure.size() != np || !pressure.is_valid())
            return false;
        if (pressure_null_space_)
            pressure_null_space_->remove(pressure);

        pressure *= -1.0;
        if (pressure_reference_cell_)
            pressure(*pressure_reference_cell_) = 0.0;
        if (pressure_null_space_)
            pressure_null_space_->remove(pressure);
        return true;
    }

    bool has_pressure_null_space_policy() const noexcept {
        return pressure_null_space_.has_value();
    }

    // PCD is intentionally exposed as an explicit Schur approximation object.
    // The production coupled preconditioner must own the pressure operators and
    // solver objects; this class only owns their non-owning views and therefore
    // requires those objects to outlive the approximation.
    const SparseMatrix& pressure_mass() const {
        if (!pressure_mass_) throw std::logic_error("PCD: pressure mass operator is not configured");
        return *pressure_mass_;
    }

    const SparseMatrix& pressure_laplacian() const {
        if (!pressure_laplacian_) throw std::logic_error("PCD: pressure Laplacian is not configured");
        return *pressure_laplacian_;
    }

    const SparseMatrix& pressure_convection_diffusion() const {
        if (!pressure_convection_diffusion_)
            throw std::logic_error("PCD: pressure convection-diffusion operator is not configured");
        return *pressure_convection_diffusion_;
    }

    bool has_pressure_operators() const noexcept {
        return pressure_mass_ && pressure_laplacian_ &&
               pressure_convection_diffusion_;
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
        add(*pressure_mass_);
        add(*pressure_laplacian_);
        add(*pressure_convection_diffusion_);
        return GraphSignature{h};
    }

    const BlockOperator* blocks_ = nullptr;
    const SparseMatrix* pressure_mass_ = nullptr;
    const SparseMatrix* pressure_laplacian_ = nullptr;
    const SparseMatrix* pressure_convection_diffusion_ = nullptr;
    PressureSolve laplacian_solve_;
    PressureSolve convection_diffusion_solve_;
    std::optional<NullSpaceProjector> pressure_null_space_;
    std::optional<std::size_t> pressure_reference_cell_;
    GraphSignature graph_signature_{};
};

} // namespace cfdx::core
