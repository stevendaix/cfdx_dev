#pragma once

#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/null_space.h"
#include "cfdx/core/linalg/schur_approximation.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

// Pressure Schur inverse approximation for a coupled operator
//
//     A = [ Auu  G ]  with  S = C - D Auu^{-1} G.
//         [ D    C ]
//
// The DEFAULT legacy CFDX action is -Fp^{-1} rhs. It passed an earlier
// 512-unknown production acceptance case, but must not be called the classical
// PCD commutator formula. The older -Fp^{-1} Mp Kp^{-1} composition failed
// that case; this is not a general proof against all composite PCD operators.
//
// The optional literature-inspired action is -Mp^{-1} Fp Kp^{-1} rhs:
// first SOLVE Kp, then MULTIPLY Fp, then invert the positive LUMPED mass.
// See Kay/Loghin/Wathen (2002) and Elman et al. (2006). The latter
// multiplication is crucial: it does not introduce an Fp inverse. This
// order is not interchangeable with -Kp^{-1} Fp Mp^{-1} in general.
//
// Pressure-side operators are supplied explicitly, because inferring them
// from the velocity block would hide discretization and boundary assumptions.
// The literature model is a research option until the FV pressure-space
// operators, their boundary conditions, gauge scaling and true coupled
// convergence have been qualified independently on non-zero-flow cases.
//
// Gauge treatment: pressure reference rows belong in the supplied pressure
// operators, not zeroed in the apply() input/output. A frozen reference-cell
// correction prevents the outer Krylov solver from fixing a nonzero gauge
// residual. The optional commutator keeps the same input/output treatment;
// its physical gauge scaling is NOT yet qualified for production.
//
// The pressure operators must keep their CSR patterns between setup() and
// numeric update. A changed graph requires a fresh setup(). An updated
// operator used by a nested inverse also requires updating its callback.
// The old CFDX single-inverse approximation is preserved as the default.
// The optional commutator model implements, for a lumped positive pressure
// mass, -Mp^{-1} Fp Kp^{-1}. This is the order derived in the PCD literature;
// it is NOT equivalent to -Kp^{-1} Fp Mp^{-1} on noncommuting matrices.
//
// The optional model validates its Kp inverse callback against the actual
// Kp matrix on every apply. This detects stale callbacks after numeric updates.
// It is an algebraic research option, not a claim of finite-volume or physical
// qualification for the coupled CFD solver.
enum class PcdSchurActionMode {
    LegacySingleConvectionDiffusionInverse,
    LumpedMassCommutator,
    // Research only: the full (non-lumped) pressure mass is inverted by an
    // explicit callback, whose true residual is checked on every action.
    ConsistentMassCommutator
};

class PcdSchurApproximation final : public SchurApproximation {
public:
    using PressureSolve = std::function<bool(const Vector& rhs, Vector& x)>;

    PcdSchurApproximation(const SparseMatrix& pressure_mass,
                          const SparseMatrix& pressure_laplacian,
                          const SparseMatrix& pressure_convection_diffusion,
                          PressureSolve laplacian_solve,
                          PressureSolve convection_diffusion_solve,
                          std::optional<NullSpaceProjector> pressure_null_space = std::nullopt,
                          std::optional<std::size_t> pressure_reference_cell = std::nullopt,
                          PcdSchurActionMode mode =
                              PcdSchurActionMode::LegacySingleConvectionDiffusionInverse,
                          PressureSolve pressure_mass_solve = {})
        : pressure_mass_(pressure_mass),
          pressure_laplacian_(pressure_laplacian),
          pressure_convection_diffusion_(pressure_convection_diffusion),
          laplacian_solve_(std::move(laplacian_solve)),
          convection_diffusion_solve_(std::move(convection_diffusion_solve)),
          pressure_null_space_(std::move(pressure_null_space)),
          pressure_reference_cell_(pressure_reference_cell),
          mode_(mode),
          pressure_mass_solve_(std::move(pressure_mass_solve)) {}

    const char* name() const noexcept override {
        if (mode_ == PcdSchurActionMode::LumpedMassCommutator)
            return "pcd_lumped_mass_commutator";
        if (mode_ == PcdSchurActionMode::ConsistentMassCommutator)
            return "pcd_consistent_mass_commutator";
        return "pcd_schur";
    }
    PcdSchurActionMode mode() const noexcept { return mode_; }

    // Both actions approximate S^{-1} rhs for the pressure Schur complement.
    // The sign belongs to the inverse action, so this
    // is InverseOperator rather than the SIMPLE/SIMPLEC operator action.
    SchurAction action() const noexcept override { return SchurAction::InverseOperator; }

    bool setup(const BlockOperator& blocks) override {
        if (!blocks.is_valid() || !laplacian_solve_ || !convection_diffusion_solve_)
            return false;

        const std::size_t np = blocks.pressure_size();
        if (np == 0)
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
        if (mode_ == PcdSchurActionMode::LumpedMassCommutator &&
            !positive_diagonal_lumped_mass(pressure_mass_, np))
            return false;
        if (mode_ == PcdSchurActionMode::ConsistentMassCommutator &&
            (!pressure_mass_solve_ || !certified_positive_definite_mass(pressure_mass_, np)))
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
            !valid_pressure_operator(pressure_mass_, np) ||
            !valid_pressure_operator(pressure_laplacian_, np) ||
            !valid_pressure_operator(pressure_convection_diffusion_, np))
            return false;
        if (pressure_null_space_ &&
            pressure_null_space_->dimension() != np)
            return false;
        if (mode_ == PcdSchurActionMode::LumpedMassCommutator &&
            !positive_diagonal_lumped_mass(pressure_mass_, np))
            return false;
        if (mode_ == PcdSchurActionMode::ConsistentMassCommutator &&
            (!pressure_mass_solve_ || !certified_positive_definite_mass(pressure_mass_, np)))
            return false;
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
        // Any matrix used by a nested inverse must keep the numerical state
        // of its callback. Updating only the matrix would silently continue
        // solving the old operator, even though the graph is unchanged.
        if ((mode_ == PcdSchurActionMode::LumpedMassCommutator &&
             (!positive_diagonal_lumped_mass(pressure_mass, np) ||
              !same_values(pressure_laplacian_, pressure_laplacian))) ||
            (mode_ == PcdSchurActionMode::ConsistentMassCommutator &&
             (!pressure_mass_solve_ ||
              !certified_positive_definite_mass(pressure_mass, np) ||
              !same_values(pressure_mass_, pressure_mass) ||
              !same_values(pressure_laplacian_, pressure_laplacian))) ||
            (mode_ == PcdSchurActionMode::LegacySingleConvectionDiffusionInverse &&
             !same_values(pressure_convection_diffusion_,
                          pressure_convection_diffusion)))
            return false; // Use update_pressure_values_and_solves instead.

        pressure_mass_ = pressure_mass;
        pressure_laplacian_ = pressure_laplacian;
        pressure_convection_diffusion_ = pressure_convection_diffusion;
        pressure_size_ = np;
        graph_signature_ = graph_signature(blocks);
        return true;
    }

    // A complete numeric update including the solvers whose captured state
    // must correspond to the new matrices. Validate before modifying any
    // member; pressure operators and callbacks change as one transaction.
    bool update_pressure_values_and_solves(
        const BlockOperator& blocks,
        const SparseMatrix& pressure_mass,
        const SparseMatrix& pressure_laplacian,
        const SparseMatrix& pressure_convection_diffusion,
        PressureSolve laplacian_solve,
        PressureSolve convection_diffusion_solve,
        PressureSolve pressure_mass_solve = {}) {
        if (pressure_size_ == 0 || !blocks.is_valid() ||
            graph_signature(blocks) != graph_signature_ ||
            !laplacian_solve || !convection_diffusion_solve)
            return false;
        const std::size_t np = blocks.pressure_size();
        if (!valid_pressure_operator(pressure_mass, np) ||
            !valid_pressure_operator(pressure_laplacian, np) ||
            !valid_pressure_operator(pressure_convection_diffusion, np) ||
            !same_pattern(pressure_mass_, pressure_mass) ||
            !same_pattern(pressure_laplacian_, pressure_laplacian) ||
            !same_pattern(pressure_convection_diffusion_,
                          pressure_convection_diffusion) ||
            (mode_ == PcdSchurActionMode::LumpedMassCommutator &&
             !positive_diagonal_lumped_mass(pressure_mass, np)) ||
            (mode_ == PcdSchurActionMode::ConsistentMassCommutator &&
             (!pressure_mass_solve ||
              !certified_positive_definite_mass(pressure_mass, np))))
            return false;

        SparseMatrix new_mass(pressure_mass);
        SparseMatrix new_laplacian(pressure_laplacian);
        SparseMatrix new_convection_diffusion(pressure_convection_diffusion);
        pressure_mass_ = std::move(new_mass);
        pressure_laplacian_ = std::move(new_laplacian);
        pressure_convection_diffusion_ = std::move(new_convection_diffusion);
        laplacian_solve_ = std::move(laplacian_solve);
        convection_diffusion_solve_ = std::move(convection_diffusion_solve);
        pressure_mass_solve_ = std::move(pressure_mass_solve);
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

        Vector candidate(np, 0.0);
        if (mode_ != PcdSchurActionMode::LegacySingleConvectionDiffusionInverse) {
            // Apply -Mp^{-1} Fp Kp^{-1} in the mathematically correct order.
            // This requires an explicitly verified *lumped* Mp; a general
            // mass matrix must not be silently inverted using its diagonal.
            Vector kp_inverse_rhs(np, 0.0);
            if (!laplacian_solve_(rhs, kp_inverse_rhs) ||
                kp_inverse_rhs.size() != np || !kp_inverse_rhs.is_valid() ||
                !laplacian_inverse_is_current(rhs, kp_inverse_rhs))
                return false;

            Vector fp_product(np, 0.0);
            pressure_convection_diffusion_.matvec_into(
                kp_inverse_rhs, fp_product);
            if (mode_ == PcdSchurActionMode::LumpedMassCommutator) {
                for (std::size_t row = 0; row < np; ++row) {
                    const double mass = pressure_mass_(row, row);
                    candidate(row) = -fp_product(row) / mass;
                }
            } else {
                // A consistent/full pressure mass cannot be inverted using
                // its diagonal. Fail closed on stale or inaccurate callbacks.
                if (!pressure_mass_solve_ ||
                    !pressure_mass_solve_(fp_product, candidate) ||
                    candidate.size() != np || !candidate.is_valid() ||
                    !inverse_is_current(pressure_mass_, fp_product, candidate))
                    return false;
                candidate *= -1.0;
            }
            if (!candidate.is_valid())
                return false;
        } else {
            // Historical single-inverse action retained by default. Write to
            // a temporary so even a failing callback cannot corrupt output.
            if (!convection_diffusion_solve_(rhs, candidate) ||
                candidate.size() != np || !candidate.is_valid())
                return false;
            candidate *= -1.0;
        }
        if (pressure_null_space_)
            pressure_null_space_->remove(candidate);
        pressure = std::move(candidate);
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
    static bool positive_diagonal_lumped_mass(const SparseMatrix& mass,
                                               std::size_t np) {
        if (!valid_pressure_operator(mass, np))
            return false;
        for (std::size_t row = 0; row < np; ++row) {
            const double diagonal = mass(row, row);
            if (!(diagonal > 0.0) || !std::isfinite(diagonal))
                return false;
            for (std::size_t k = mass.row_offsets_data()[row];
                 k < mass.row_offsets_data()[row + 1]; ++k) {
                if (mass.columns_data()[k] != row &&
                    mass.values_data()[k] != 0.0)
                    return false;
            }
        }
        return true;
    }

    // Independent true-residual check; hypot avoids overflow of residual^2.
    // Always verify a nested inverse against the current stored operator.
    bool inverse_is_current(const SparseMatrix& matrix, const Vector& rhs,
                            const Vector& solution) const {
        if (rhs.size() != pressure_size_ || solution.size() != pressure_size_ ||
            !rhs.is_valid() || !solution.is_valid())
            return false;
        Vector actual(pressure_size_, 0.0);
        matrix.matvec_into(solution, actual);
        if (!actual.is_valid()) return false;
        double residual_norm = 0.0;
        double rhs_norm = 0.0;
        for (std::size_t row = 0; row < pressure_size_; ++row) {
            residual_norm = std::hypot(residual_norm, actual(row) - rhs(row));
            rhs_norm = std::hypot(rhs_norm, rhs(row));
        }
        constexpr double relative_tol = 1e-9;
        constexpr double absolute_tol = 1e-12;
        return std::isfinite(residual_norm) && std::isfinite(rhs_norm) &&
            residual_norm <= absolute_tol + relative_tol * rhs_norm;
    }

    bool laplacian_inverse_is_current(const Vector& rhs,
                                      const Vector& solution) const {
        return inverse_is_current(pressure_laplacian_, rhs, solution);
    }

    // A symmetric matrix with positive diagonal can still be indefinite.
    // The consistent-mass PCD path requires an SPD pressure mass. Check
    // symmetry, then use an actual Cholesky factorization on small matrices;
    // for larger matrices accept only the sufficient (conservative) criterion
    // of strict symmetric diagonal dominance. An unproved SPD matrix is
    // rejected, never treated as valid merely because its diagonals are >0.
    // This research-only setup check is not a substitute for mesh/FE mass
    // provenance or a scalable sparse SPD factorization.
    static bool certified_positive_definite_mass(const SparseMatrix& matrix,
                                                 std::size_t np) {
        if (!valid_pressure_operator(matrix, np)) return false;
        double maximum_diagonal = 0.0;
        bool strict_diagonal_dominance = true;
        for (std::size_t row = 0; row < np; ++row) {
            const double diagonal = matrix(row, row);
            if (!(diagonal > 0.0) || !std::isfinite(diagonal)) return false;
            maximum_diagonal = std::max(maximum_diagonal, diagonal);
            double off_diagonal_sum = 0.0;
            for (std::size_t k = matrix.row_offsets_data()[row];
                 k < matrix.row_offsets_data()[row + 1]; ++k) {
                const std::size_t col = matrix.columns_data()[k];
                const double val = matrix.values_data()[k];
                if (col >= np) return false;
                const double opposite = matrix(col, row);
                if (!std::isfinite(opposite) ||
                    std::abs(val - opposite) >
                        1e-12 * std::max({1.0, std::abs(val), std::abs(opposite)}))
                    return false;
                if (col != row) off_diagonal_sum += std::abs(val);
            }
            if (!std::isfinite(off_diagonal_sum) || !(diagonal > off_diagonal_sum))
                strict_diagonal_dominance = false;
        }
        // Avoid quadratic memory on a large pressure system. Symmetric,
        // strictly diagonally dominant matrices with positive diagonal are SPD.
        constexpr std::size_t dense_limit = 96;
        if (np > dense_limit) return strict_diagonal_dominance;

        std::vector<double> lower(np * np, 0.0);
        const double pivot_floor = 1e-12 * maximum_diagonal;
        for (std::size_t i = 0; i < np; ++i) {
            for (std::size_t j = 0; j <= i; ++j) {
                double sum = matrix(i, j);
                for (std::size_t k = 0; k < j; ++k)
                    sum -= lower[i * np + k] * lower[j * np + k];
                if (!std::isfinite(sum)) return false;
                if (i == j) {
                    if (!(sum > pivot_floor)) return false;
                    lower[i * np + j] = std::sqrt(sum);
                } else {
                    lower[i * np + j] = sum / lower[j * np + j];
                    if (!std::isfinite(lower[i * np + j])) return false;
                }
            }
        }
        return true;
    }

    static bool same_values(const SparseMatrix& a, const SparseMatrix& b) {
        if (!same_pattern(a, b))
            return false;
        for (std::size_t k = 0; k < a.nnz(); ++k)
            if (a.values_data()[k] != b.values_data()[k])
                return false;
        return true;
    }

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
    PcdSchurActionMode mode_ =
        PcdSchurActionMode::LegacySingleConvectionDiffusionInverse;
    PressureSolve pressure_mass_solve_;
    GraphSignature graph_signature_{};
};

} // namespace cfdx::core
