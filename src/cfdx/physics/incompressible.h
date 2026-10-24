#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/convection.h"
#include "cfdx/core/numerics/flux.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/numerics/laplacian.h"
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

struct IncompressibleMomentumTerms {
    cfdx::core::Field<double, cfdx::core::Location::CELL> convection;
    cfdx::core::Field<double, cfdx::core::Location::CELL> diffusion;
    cfdx::core::Field<double, cfdx::core::Location::CELL> pressure_gradient;
};

inline IncompressibleMomentumTerms evaluate_momentum_terms(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure,
    const cfdx::core::Mesh& mesh,
    double kinematic_viscosity,
    cfdx::core::InterpScheme convection_scheme = cfdx::core::InterpScheme::UPWIND)
{
    if (velocity.dimension() != 3 || pressure.dimension() != 1) {
        throw std::runtime_error("evaluate_momentum_terms: expected U(dim=3) and p(dim=1)");
    }
    if (!(kinematic_viscosity >= 0.0) || !std::isfinite(kinematic_viscosity)) {
        throw std::invalid_argument("evaluate_momentum_terms: invalid viscosity");
    }

    auto face_velocity = cfdx::core::interpolate_cell_to_face(
        velocity, mesh, cfdx::core::InterpScheme::LINEAR);
    auto face_flux = cfdx::core::compute_flux(face_velocity, mesh);

    auto convection = cfdx::core::compute_vector_convection(
        velocity, face_flux, mesh, convection_scheme);
    // F12 (audit 2026-10-24): the three returned terms must be homogeneous.
    // compute_convection returns an integrated face-flux divergence while
    // the laplacian and the gradient are already per-unit-volume operators,
    // so adding the three fields mixed m^4/s^2 with m/s^2. Divide the
    // convective divergence by the cell volume: every term below is now a
    // body term in m/s^2 for the local equation
    //   dU/dt = -convection + diffusion - pressure_gradient/rho
    // with pressure a dynamic pressure in Pa (gradient in Pa/m).
    const auto f12_geometry = cfdx::core::make_geometry_cache(mesh);
    for (std::size_t d = 0; d < 3; ++d) {
        double* data = convection.component_data(d);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const double volume = f12_geometry.cell_volumes[c];
            if (!(volume > 0.0) || !std::isfinite(volume))
                throw std::runtime_error("evaluate_momentum_terms: invalid cell volume");
            data[c] /= volume;
        }
    }

    cfdx::core::Field<double, cfdx::core::Location::CELL> diffusion(
        mesh.n_cells(), "U_diffusion", velocity.metadata().unit + "/s", 3);
    for (std::size_t d = 0; d < 3; ++d) {
        cfdx::core::Field<double, cfdx::core::Location::CELL> component(
            mesh.n_cells(), "U_component", velocity.metadata().unit, 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            component.component_data(0)[c] = velocity.component_data(d)[c];
        const auto lap = cfdx::core::compute_laplacian(component, mesh);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            diffusion.component_data(d)[c] = kinematic_viscosity * lap.component_data(0)[c];
    }

    auto pressure_gradient = cfdx::core::compute_gradient_gauss(pressure, mesh);
    return {std::move(convection), std::move(diffusion), std::move(pressure_gradient)};
}

inline cfdx::core::Field<double, cfdx::core::Location::CELL>
compute_continuity_residual(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Mesh& mesh)
{
    if (velocity.dimension() != 3 || velocity.size() != mesh.n_cells()) {
        throw std::runtime_error("compute_continuity_residual: invalid velocity");
    }
    auto face_velocity = cfdx::core::interpolate_cell_to_face(
        velocity, mesh, cfdx::core::InterpScheme::LINEAR);
    auto flux = cfdx::core::compute_flux(face_velocity, mesh);

    cfdx::core::Field<double, cfdx::core::Location::CELL> result(
        mesh.n_cells(), "continuity_residual", "m^3/s", 1);
    const auto& own = mesh.ownership();
    const auto& cells = mesh.cells();
    const auto* cf = cells.faces_data();
    const auto* co = cells.offsets_data();
    const double* f = flux.component_data(0);
    double* r = result.component_data(0);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        r[c] = 0.0;
        for (cfdx::core::Offset k = co[c]; k < co[c + 1]; ++k) {
            const std::size_t face = cf[k];
            r[c] += (own.owner(face) == c) ? f[face] : -f[face];
        }
    }
    return result;
}

}  // namespace cfdx::physics
