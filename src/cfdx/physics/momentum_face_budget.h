#pragma once

#include "cfdx/physics/finite_volume_transport.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Diagnostic ONLY. Every value is an integrated signed cell contribution to
// div(rho U u_d) - div(mu grad u_d) + grad_d p - rho f_d = 0.
// `collocated_difference` measures the convection difference between the
// AUTHORITATIVE mass flux and the direct interpolation of cell velocities.
// This includes both Rhie-Chow and pressure-relaxation splitting; it must
// NOT be mislabelled as the pure Rhie-Chow correction.
struct MomentumFaceCellBudget {
    double convection = 0.0;
    double viscous = 0.0;
    double pressure = 0.0;
    double body_force = 0.0;
    double bounded_divergence = 0.0;
    double transpose_viscous = 0.0;
    double mass_divergence = 0.0;
    double collocated_difference = 0.0;
    // Positive physical time derivative, zero for steady/pseudo-time solves.
    double temporal = 0.0;
    double physical_residual = 0.0;
    double bounded_operator_residual = 0.0;
    // Independently evaluated CSR row A(U)U-b(U) at the same accepted state.
    // Filled only by reconcile_momentum_budget_with_matrix().
    double algebraic_residual = std::numeric_limits<double>::quiet_NaN();
    double matrix_face_difference = std::numeric_limits<double>::quiet_NaN();
    double matrix_face_scale = std::numeric_limits<double>::quiet_NaN();
};

struct MomentumFaceBudget {
    std::vector<MomentumFaceCellBudget> cells;
    std::size_t worst_cell = 0;
    double worst_absolute_residual = 0.0;
    double worst_bounded_operator_residual = 0.0;
    double flux_decomposition_error = 0.0;
    bool matrix_closure_checked = false;
    double matrix_closure_linf = std::numeric_limits<double>::infinity();
    double matrix_closure_relative_linf = std::numeric_limits<double>::infinity();
    std::size_t matrix_closure_worst_cell = 0;
};

// Inputs represent the ACCEPTED iterate; the mass flux must be the same
// conservative pressure-corrected flux consumed by the solver. The optional
// collocated flux is an alternate observable, never a replacement for it.
inline MomentumFaceBudget diagnose_momentum_face_budget(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_force,
    const std::vector<double>& cell_viscosity,
    const ScalarBoundaryConditions& component_bcs,
    ConvectionScheme scheme,
    double high_order_fraction,
    bool bounded_convection,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>* collocated_mass_flux = nullptr,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>* transpose_stress_flux = nullptr,
    const std::vector<double>* wall_face_viscosity = nullptr,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* negative_accumulation = nullptr,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
{
    using cfdx::core::Field;
    using cfdx::core::Location;
    const auto nc = mesh.n_cells();
    const auto nf = mesh.n_faces();
    if (mass_flux.size() != nf || mass_flux.dimension() != 1 ||
        component.size() != nc || component.dimension() != 1 ||
        pressure_gradient.size() != nc || pressure_gradient.dimension() != 1 ||
        body_force.size() != nc || body_force.dimension() != 1 ||
        cell_viscosity.size() != nc || geometry.cell_volumes.size() != nc ||
        (collocated_mass_flux && (collocated_mass_flux->size() != nf ||
                                  collocated_mass_flux->dimension() != 1)) ||
        (transpose_stress_flux && (transpose_stress_flux->size() != nf ||
                                   transpose_stress_flux->dimension() != 1)) ||
        (wall_face_viscosity && wall_face_viscosity->size() != nf) ||
        (negative_accumulation && (negative_accumulation->size() != nc ||
                                   negative_accumulation->dimension() != 1)))
        throw std::invalid_argument("momentum face budget: field shape mismatch");
    if (!(high_order_fraction >= 0.0 && high_order_fraction <= 1.0) ||
        !std::isfinite(high_order_fraction))
        throw std::invalid_argument("momentum face budget: invalid order fraction");

    // Run the physical face reconstructor twice with identical convection
    // inputs: the difference is the complete orthogonal + deferred viscous
    // flux, so internal faces are accounted for once with opposite signs.
    const auto full_flux = reconstruct_scalar_transport_flux(
        mesh, geometry, mass_flux, component, 0.0, component_bcs,
        bounded_convection, scheme, face_conditions, &cell_viscosity, nullptr,
        wall_face_viscosity, high_order_fraction);
    const auto convection_flux = reconstruct_scalar_transport_flux(
        mesh, geometry, mass_flux, component, 0.0, component_bcs,
        bounded_convection, scheme, face_conditions, nullptr, nullptr,
        nullptr, high_order_fraction);
    Field<double, Location::FACE> collocated_convection;
    if (collocated_mass_flux) {
        collocated_convection = reconstruct_scalar_transport_flux(
            mesh, geometry, *collocated_mass_flux, component, 0.0,
            component_bcs, bounded_convection, scheme, face_conditions, nullptr,
            nullptr, nullptr, high_order_fraction);
    }

    MomentumFaceBudget out;
    out.cells.resize(nc);
    for (std::size_t cell = 0; cell < nc; ++cell) {
        const auto volume = geometry.cell_volumes[cell];
        if (!std::isfinite(volume) || !(volume > 0.0) ||
            !std::isfinite(pressure_gradient(cell)) ||
            !std::isfinite(body_force(cell)) || !std::isfinite(component(cell)))
            throw std::invalid_argument("momentum face budget: invalid cell state");
        out.cells[cell].pressure = pressure_gradient(cell) * volume;
        out.cells[cell].body_force = -body_force(cell) * volume;
    }
    const auto& ownership = mesh.ownership();
    for (std::size_t face = 0; face < nf; ++face) {
        const double F = mass_flux(face);
        const double adv = convection_flux(face);
        const double visc = full_flux(face) - adv;
        const double trans = transpose_stress_flux
            ? -(*transpose_stress_flux)(face) : 0.0;
        const double collocated_delta = collocated_mass_flux
            ? adv - collocated_convection(face) : 0.0;
        if (!std::isfinite(F) || !std::isfinite(adv) ||
            !std::isfinite(visc) || !std::isfinite(trans) ||
            !std::isfinite(collocated_delta))
            throw std::invalid_argument("momentum face budget: nonfinite face transport");
        const auto owner = ownership.owner(face);
        if (owner >= nc)
            throw std::invalid_argument("momentum face budget: invalid face owner");
        const auto apply = [&](std::size_t cell, double sign) {
            auto& b = out.cells[cell];
            b.convection += sign * adv;
            b.viscous += sign * visc;
            b.transpose_viscous += sign * trans;
            b.mass_divergence += sign * F;
            b.collocated_difference += sign * collocated_delta;
        };
        apply(owner, 1.0);
        const auto neighbour = ownership.neighbour(face);
        if (neighbour >= 0) {
            if (static_cast<std::size_t>(neighbour) >= nc)
                throw std::invalid_argument("momentum face budget: invalid face neighbour");
            apply(static_cast<std::size_t>(neighbour), -1.0);
        }
    }
    for (std::size_t cell = 0; cell < nc; ++cell) {
        auto& b = out.cells[cell];
        b.bounded_divergence = bounded_convection
            ? -component(cell) * b.mass_divergence : 0.0;
        const double temporal = negative_accumulation
            ? (*negative_accumulation)(cell) : 0.0;
        if (!std::isfinite(temporal))
            throw std::invalid_argument("momentum face budget: nonfinite accumulation");
        b.temporal = -temporal;
        b.physical_residual = b.convection + b.viscous + b.transpose_viscous +
                              b.pressure + b.body_force + b.temporal;
        b.bounded_operator_residual = b.physical_residual + b.bounded_divergence;
        const double magnitude = std::abs(b.physical_residual);
        if (magnitude > out.worst_absolute_residual) {
            out.worst_absolute_residual = magnitude;
            out.worst_cell = cell;
        }
        out.worst_bounded_operator_residual = std::max(
            out.worst_bounded_operator_residual,
            std::abs(b.bounded_operator_residual));
    }
    return out;
}

// This is NOT a second pressure or convection discretization. It is an
// independent arithmetic path (CSR matrix row versus signed face loop)
// which detects missing/double-counted source, TVD deferred correction,
// transpose stress, non-orthogonal diffusion and bounded-divergence terms.
// Critically, compare SIGNED values for EVERY cell, not just two maxima.
inline void reconcile_momentum_budget_with_matrix(
    MomentumFaceBudget& budget,
    const ScalarEquation& equation,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& component)
{
    const std::size_t nc = budget.cells.size();
    const auto& A = equation.matrix;
    if (nc == 0 || component.size() != nc || component.dimension() != 1 ||
        equation.rhs.size() != nc || A.n_rows() != nc || A.n_cols() != nc)
        throw std::invalid_argument("momentum matrix/face closure: field shape mismatch");

    budget.matrix_closure_checked = true;
    budget.matrix_closure_linf = 0.0;
    budget.matrix_closure_relative_linf = 0.0;
    budget.matrix_closure_worst_cell = 0;
    const auto* offsets = A.row_offsets_data();
    const auto* columns = A.columns_data();
    const auto* values = A.values_data();
    for (std::size_t row = 0; row < nc; ++row) {
        // Evaluate A*U-b by walking the CSR row independently of the
        // face-flux reconstruction. Sum absolute terms to set a round-off
        // SCALE without relaxing any actual physical convergence criterion.
        double ax = 0.0;
        double scale = std::abs(equation.rhs(row));
        for (auto k = offsets[row]; k < offsets[row+1]; ++k) {
            if (columns[k] >= nc)
                throw std::invalid_argument("momentum matrix/face closure: invalid CSR column");
            const double term = values[k] * component(columns[k]);
            if (!std::isfinite(term))
                throw std::invalid_argument("momentum matrix/face closure: nonfinite CSR term");
            ax += term;
            scale += std::abs(term);
        }
        const double algebraic = ax - equation.rhs(row);
        auto& face = budget.cells[row];
        scale = std::max(scale,
            std::abs(face.convection) + std::abs(face.viscous) +
            std::abs(face.transpose_viscous) + std::abs(face.pressure) +
            std::abs(face.body_force) + std::abs(face.temporal) +
            std::abs(face.bounded_divergence));
        const double error = algebraic - face.bounded_operator_residual;
        if (!std::isfinite(algebraic) || !std::isfinite(error) || !std::isfinite(scale))
            throw std::invalid_argument("momentum matrix/face closure: nonfinite result");
        face.algebraic_residual = algebraic;
        face.matrix_face_difference = error;
        face.matrix_face_scale = scale;
        const double relative = std::abs(error) / std::max(scale, 1e-300);
        if (std::abs(error) > budget.matrix_closure_linf) {
            budget.matrix_closure_linf = std::abs(error);
            budget.matrix_closure_worst_cell = row;
        }
        budget.matrix_closure_relative_linf = std::max(
            budget.matrix_closure_relative_linf, relative);
    }
}

inline bool momentum_matrix_face_closure_consistent(
    const MomentumFaceBudget& budget, double arithmetic_relative_tolerance = 2e-10)
{
    if (!budget.matrix_closure_checked || !(arithmetic_relative_tolerance > 0.0) ||
        !std::isfinite(arithmetic_relative_tolerance))
        return false;
    for (const auto& cell : budget.cells) {
        if (!std::isfinite(cell.matrix_face_difference) ||
            !std::isfinite(cell.matrix_face_scale) ||
            std::abs(cell.matrix_face_difference) >
                arithmetic_relative_tolerance * std::max(cell.matrix_face_scale, 1e-15))
            return false;
    }
    return true;
}

} // namespace cfdx::physics
