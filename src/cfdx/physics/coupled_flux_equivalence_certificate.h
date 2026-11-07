#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

namespace cfdx::physics {

// Independent continuity row / authoritative face-flux certificate.
// The CSR matrix and RHS come from the actual Krylov assembly; face fluxes
// must come from the production FV reconstruction, NOT from A*x-b itself.
// This distinguishes an algebraically converged system from continuity of the
// face flux used by the rest of the CFD solver (Rhie-Chow, boundaries, etc.).
struct CoupledFluxEquivalenceControls {
    double relative_tolerance = 1e-8;
    double absolute_tolerance = 1e-10;
    std::optional<std::size_t> pressure_gauge_cell = std::nullopt;
    // Optional independent conservation requirement; no implicit tolerance
    // relaxation and no inference from linear residual alone.
    std::optional<double> continuity_linf_gate = std::nullopt;
};

struct CoupledFluxEquivalenceCertificate {
    bool complete = false;
    bool passed = false;
    bool continuity_gate_passed = false;
    std::size_t checked_cells = 0;
    std::size_t worst_cell = 0;
    double max_identity_error = std::numeric_limits<double>::infinity();
    double max_identity_tolerance_ratio = std::numeric_limits<double>::infinity();
    double max_conservative_flux_divergence = std::numeric_limits<double>::infinity();
    double net_boundary_flux = std::numeric_limits<double>::infinity();
    double net_cell_divergence = std::numeric_limits<double>::infinity();
};

inline CoupledFluxEquivalenceCertificate certify_coupled_face_flux_equivalence(
    const cfdx::core::SparseMatrix& A,
    const cfdx::core::Vector& rhs,
    const cfdx::core::Vector& state,
    std::size_t nc,
    const std::vector<std::size_t>& face_owner,
    const std::vector<std::int64_t>& face_neighbour,
    const std::vector<double>& authoritative_flux,
    const CoupledFluxEquivalenceControls& controls = {})
{
    CoupledFluxEquivalenceCertificate result;
    if (nc == 0 || nc > std::numeric_limits<std::size_t>::max() / 4 ||
        A.n_rows() != 4 * nc || A.n_cols() != 4 * nc ||
        rhs.size() != 4 * nc || state.size() != 4 * nc ||
        !rhs.is_valid() || !state.is_valid() ||
        face_owner.size() != face_neighbour.size() ||
        face_owner.size() != authoritative_flux.size() ||
        face_owner.empty() ||
        !std::isfinite(controls.absolute_tolerance) ||
        !std::isfinite(controls.relative_tolerance) ||
        controls.absolute_tolerance < 0.0 || controls.relative_tolerance < 0.0 ||
        (controls.absolute_tolerance == 0.0 && controls.relative_tolerance == 0.0) ||
        (controls.pressure_gauge_cell && *controls.pressure_gauge_cell >= nc) ||
        (controls.continuity_linf_gate &&
         (!std::isfinite(*controls.continuity_linf_gate) ||
          *controls.continuity_linf_gate < 0.0))) return result;

    std::vector<double> divergence(nc, 0.0);
    double boundary_net = 0.0;
    for (std::size_t face = 0; face < authoritative_flux.size(); ++face) {
        const std::size_t owner = face_owner[face];
        const std::int64_t neighbour = face_neighbour[face];
        const double phi = authoritative_flux[face];
        if (owner >= nc || !std::isfinite(phi) ||
            (neighbour >= 0 && (static_cast<std::uint64_t>(neighbour) >= nc ||
                                static_cast<std::size_t>(neighbour) == owner)) ||
            neighbour < -1) return result;
        divergence[owner] += phi;
        if (neighbour >= 0)
            divergence[static_cast<std::size_t>(neighbour)] -= phi;
        else boundary_net += phi;
    }
    double global_divergence = 0.0;
    double max_divergence = 0.0;
    for (const double value : divergence) {
        if (!std::isfinite(value)) return result;
        global_divergence += value;
        max_divergence = std::max(max_divergence, std::abs(value));
    }
    if (!std::isfinite(global_divergence) || !std::isfinite(boundary_net))
        return result;
    result.net_boundary_flux = boundary_net;
    result.net_cell_divergence = global_divergence;
    result.max_conservative_flux_divergence = max_divergence;
    result.max_identity_error = 0.0;
    result.max_identity_tolerance_ratio = 0.0;
    for (std::size_t cell = 0; cell < nc; ++cell) {
        if (controls.pressure_gauge_cell && cell == *controls.pressure_gauge_cell)
            continue; // Gauge equation is NOT a continuity equation.
        const std::size_t row = 3 * nc + cell;
        double ax = 0.0;
        for (std::size_t k = A.row_offsets_data()[row];
             k < A.row_offsets_data()[row + 1]; ++k) {
            const std::size_t col = A.columns_data()[k];
            if (col >= 4 * nc || !std::isfinite(A.values_data()[k])) return result;
            ax += A.values_data()[k] * state(col);
        }
        const double matrix_continuity = ax - rhs(row);
        const double error = std::abs(matrix_continuity - divergence[cell]);
        const double tolerance = controls.absolute_tolerance +
            controls.relative_tolerance *
                std::max(std::abs(matrix_continuity), std::abs(divergence[cell]));
        // Relative-only tolerances have a meaningful zero/zero case:
        // an exact zero residual satisfies the gate without a positive floor.
        const double ratio = tolerance > 0.0 ? error / tolerance
            : (error == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
        if (!std::isfinite(ratio)) return result;
        if (ratio > result.max_identity_tolerance_ratio) {
            result.max_identity_tolerance_ratio = ratio;
            result.worst_cell = cell;
        }
        result.max_identity_error = std::max(result.max_identity_error, error);
        ++result.checked_cells;
    }
    if (result.checked_cells != nc - (controls.pressure_gauge_cell ? 1U : 0U))
        return result;
    // This equality is a topology/orientation check: internal flux must
    // cancel between owner/neighbour. It is NOT a mass-convergence criterion.
    const double budget_tolerance = controls.absolute_tolerance +
        controls.relative_tolerance *
            std::max(std::abs(global_divergence), std::abs(boundary_net));
    const bool budget_consistent = std::abs(global_divergence - boundary_net) <=
        budget_tolerance;
    result.complete = budget_consistent;
    result.continuity_gate_passed = !controls.continuity_linf_gate ||
        max_divergence <= *controls.continuity_linf_gate;
    result.passed = result.complete && result.max_identity_tolerance_ratio <= 1.0 &&
        result.continuity_gate_passed;
    return result;
}

} // namespace cfdx::physics
