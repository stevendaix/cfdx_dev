#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace cfdx::physics {

// Checks A*x-b against a second evaluation of the FROZEN Picard equations.
// This is deliberately NOT a Newton-Jacobian test.  The momentum source and
// transport coefficients are those received by the actual COUPLED assembler,
// while the pressure gradient and Rhie-Chow flux are recomputed by FV faces.
// It catches inconsistencies in block coupling, boundary offsets and gauge
// treatment, but cannot validate an upstream erroneous momentum discretization.
struct CoupledFrozenFvReport {
    bool complete = false;
    bool passed = false;
    std::size_t checked_rows = 0;
    std::size_t worst_row = 0;
    double max_absolute_error = std::numeric_limits<double>::infinity();
    double max_tolerance_ratio = std::numeric_limits<double>::infinity();
    double momentum_linf = std::numeric_limits<double>::infinity();
    double continuity_linf = std::numeric_limits<double>::infinity();
};

inline CoupledFrozenFvReport certify_coupled_frozen_fv_operator(
    const CoupledAssemblySnapshot& evidence,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    double density,
    double reference_value,
    const cfdx::core::Vector& state,
    double relative_tolerance = 1e-8,
    double absolute_tolerance = 1e-9)
{
    using namespace cfdx::core;
    CoupledFrozenFvReport report;
    const std::size_t nc = mesh.n_cells();
    if (!nc || nc > std::numeric_limits<std::size_t>::max() / 4 ||
        !evidence.matrix || !evidence.rhs || !evidence.pressure_old ||
        !std::all_of(evidence.momentum_equations.begin(),
                     evidence.momentum_equations.end(),
                     [](const ScalarEquation* e) { return e != nullptr; }) ||
        evidence.matrix->n_rows() != 4*nc || evidence.matrix->n_cols() != 4*nc ||
        evidence.rhs->size() != 4*nc || state.size() != 4*nc ||
        !evidence.rhs->is_valid() || !state.is_valid() ||
        evidence.pressure_old->size() != nc ||
        evidence.reference_cell >= nc ||
        !(density > 0.0) || !std::isfinite(density) ||
        !std::isfinite(reference_value) ||
        !std::isfinite(relative_tolerance) || !std::isfinite(absolute_tolerance) ||
        relative_tolerance < 0.0 || absolute_tolerance < 0.0 ||
        (relative_tolerance == 0.0 && absolute_tolerance == 0.0) ||
        geometry.cell_volumes.size() != nc ||
        geometry.cell_centres.size() != nc ||
        geometry.face_centres.size() != mesh.n_faces() ||
        geometry.face_area_vectors.size() != mesh.n_faces() ||
        geometry.face_patch.size() != mesh.n_faces()) return report;

    for (const ScalarEquation* eq : evidence.momentum_equations)
        if (eq->matrix.n_rows() != nc || eq->matrix.n_cols() != nc ||
            eq->rhs.size() != nc || eq->diagonal.size() != nc ||
            !eq->rhs.is_valid()) return report;

    Field<double,Location::CELL> pressure(nc, "frozen_p", "Pa", 1);
    for (std::size_t c = 0; c < nc; ++c)
        pressure(c) = (!evidence.has_fixed_pressure && c == evidence.reference_cell)
            ? reference_value : state(3*nc + c);

    // Separate re-evaluation of the actual Gauss pressure gradient (including
    // nonzero prescribed normal gradients) with all Picard data held frozen.
    const auto grad_trial = gauss_gradient_with_boundary(
        pressure, mesh, geometry, pressure_bcs);
    const auto grad_old = gauss_gradient_with_boundary(
        *evidence.pressure_old, mesh, geometry, pressure_bcs);

    std::vector<double> fv(4*nc, 0.0);
    for (std::size_t component = 0; component < 3; ++component) {
        const ScalarEquation& eq = *evidence.momentum_equations[component];
        for (std::size_t c = 0; c < nc; ++c) {
            double transport = 0.0;
            bool diagonal_present = false;
            for (std::size_t k = eq.matrix.row_offsets_data()[c];
                 k < eq.matrix.row_offsets_data()[c+1]; ++k) {
                const std::size_t col = eq.matrix.columns_data()[k];
                if (col >= nc) return report;
                diagonal_present |= (col == c);
                transport += eq.matrix.values_data()[k] * state(component*nc + col);
            }
            if (!diagonal_present)
                transport += eq.diagonal[c] * state(component*nc + c);
            fv[component*nc + c] = transport - eq.rhs(c) +
                geometry.cell_volumes[c] *
                (grad_trial.component_data(component)[c] -
                 grad_old.component_data(component)[c]);
        }
    }

    // Reconstruct each face flux ONCE, oriented owner -> neighbour.  This is
    // independent of the COUPLED D/C matrix and its continuity RHS, but uses
    // the same documented Rhie-Chow face coefficient and frozen correction.
    std::vector<double> divergence(nc, 0.0);
    const auto& own = mesh.ownership();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t o = own.owner(f);
        if (o >= nc) return report;
        const auto neighbour = own.neighbour(f);
        const Vec3 Sf = geometry.face_area_vectors[f];
        double phi = 0.0;
        if (neighbour >= 0) {
            const std::size_t n = static_cast<std::size_t>(neighbour);
            if (n >= nc || n == o) return report;
            const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
            const Vec3 response{
                0.5 * (geometry.cell_volumes[o] / evidence.momentum_equations[0]->diagonal[o] +
                       geometry.cell_volumes[n] / evidence.momentum_equations[0]->diagonal[n]),
                0.5 * (geometry.cell_volumes[o] / evidence.momentum_equations[1]->diagonal[o] +
                       geometry.cell_volumes[n] / evidence.momentum_equations[1]->diagonal[n]),
                0.5 * (geometry.cell_volumes[o] / evidence.momentum_equations[2]->diagonal[o] +
                       geometry.cell_volumes[n] / evidence.momentum_equations[2]->diagonal[n])};
            const Vec3 grad{
                0.5*(grad_old.component_data(0)[o] + grad_old.component_data(0)[n]),
                0.5*(grad_old.component_data(1)[o] + grad_old.component_data(1)[n]),
                0.5*(grad_old.component_data(2)[o] + grad_old.component_data(2)[n])};
            phi = density * 0.5 * (
                (state(o) + state(n))*Sf.x +
                (state(nc+o) + state(nc+n))*Sf.y +
                (state(2*nc+o) + state(2*nc+n))*Sf.z);
            phi += density * pressure_face_orthogonal_coefficient(Sf,dvec,response) *
                (pressure(o) - pressure(n));
            phi -= density * pressure_face_nonorthogonal_flux(Sf,dvec,response,grad);
            phi += density * (response.x*grad.x*Sf.x +
                              response.y*grad.y*Sf.y + response.z*grad.z*Sf.z);
            divergence[n] -= phi;
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY)
                continue;
            const std::string name = patch < mesh.boundary().n_patches()
                ? mesh.boundary().patch(patch).name : "";
            const auto vit = velocity_bcs.find(name);
            if (vit != velocity_bcs.end() &&
                vit->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE) {
                phi = density * vit->second.value.dot(Sf);
            } else {
                phi = density * (state(o)*Sf.x + state(nc+o)*Sf.y +
                                 state(2*nc+o)*Sf.z);
            }
            const auto pit = pressure_bcs.find(name);
            if (pit != pressure_bcs.end() &&
                pit->second.type == ScalarBoundaryType::FIXED_VALUE) {
                const Vec3 dvec = geometry.face_centres[f] - geometry.cell_centres[o];
                const Vec3 response{
                    geometry.cell_volumes[o] / evidence.momentum_equations[0]->diagonal[o],
                    geometry.cell_volumes[o] / evidence.momentum_equations[1]->diagonal[o],
                    geometry.cell_volumes[o] / evidence.momentum_equations[2]->diagonal[o]};
                const Vec3 grad{
                    grad_old.component_data(0)[o], grad_old.component_data(1)[o],
                    grad_old.component_data(2)[o]};
                phi += density * pressure_face_orthogonal_coefficient(Sf,dvec,response) *
                    (pressure(o) - pit->second.value);
                phi -= density * pressure_face_nonorthogonal_flux(Sf,dvec,response,grad);
                phi += density * (response.x*grad.x*Sf.x +
                                  response.y*grad.y*Sf.y + response.z*grad.z*Sf.z);
            }
        }
        if (!std::isfinite(phi)) return report;
        divergence[o] += phi;
    }
    for (std::size_t c = 0; c < nc; ++c)
        fv[3*nc+c] = (!evidence.has_fixed_pressure && c == evidence.reference_cell)
            ? state(3*nc+c) - reference_value : divergence[c];

    report.max_absolute_error = 0.0;
    report.max_tolerance_ratio = 0.0;
    report.momentum_linf = 0.0;
    report.continuity_linf = 0.0;
    const SparseMatrix& A = *evidence.matrix;
    const Vector& rhs = *evidence.rhs;
    for (std::size_t row = 0; row < 4*nc; ++row) {
        double algebraic = -rhs(row);
        for (std::size_t k = A.row_offsets_data()[row];
             k < A.row_offsets_data()[row+1]; ++k) {
            if (A.columns_data()[k] >= 4*nc) return CoupledFrozenFvReport{};
            algebraic += A.values_data()[k] * state(A.columns_data()[k]);
        }
        const double err = std::abs(algebraic - fv[row]);
        const double scale = absolute_tolerance + relative_tolerance *
            std::max(std::abs(algebraic),std::abs(fv[row]));
        const double ratio = scale > 0.0 ? err/scale :
            (err == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
        if (!std::isfinite(ratio)) return CoupledFrozenFvReport{};
        if (ratio > report.max_tolerance_ratio) {
            report.max_tolerance_ratio = ratio;
            report.worst_row = row;
        }
        report.max_absolute_error = std::max(report.max_absolute_error,err);
        if (row < 3*nc) report.momentum_linf = std::max(report.momentum_linf,err);
        else report.continuity_linf = std::max(report.continuity_linf,err);
        ++report.checked_rows;
    }
    report.complete = report.checked_rows == 4*nc;
    report.passed = report.complete && report.max_tolerance_ratio <= 1.0;
    return report;
}

} // namespace cfdx::physics
