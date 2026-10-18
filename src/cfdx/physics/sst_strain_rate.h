#pragma once

#include "cfdx/core/numerics/gradient.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

// Return the SST strain-rate magnitude |S| = sqrt(2 S_ij S_ij) from the
// nine components of grad(U), ordered as du_i/dx_j with j varying fastest.
inline double sst_strain_rate_magnitude(const std::array<double, 9>& grad_u)
{
    for (const double value : grad_u) {
        if (!std::isfinite(value))
            throw std::invalid_argument("sst_strain_rate_magnitude: non-finite velocity gradient");
    }

    const double sxx = grad_u[0];
    const double syy = grad_u[4];
    const double szz = grad_u[8];
    const double sxy = 0.5 * (grad_u[1] + grad_u[3]);
    const double sxz = 0.5 * (grad_u[2] + grad_u[6]);
    const double syz = 0.5 * (grad_u[5] + grad_u[7]);
    const double norm_squared = 2.0 * (
        sxx * sxx + syy * syy + szz * szz +
        2.0 * (sxy * sxy + sxz * sxz + syz * syz));
    const double magnitude = std::sqrt(norm_squared);
    if (!std::isfinite(magnitude))
        throw std::invalid_argument("sst_strain_rate_magnitude: strain-rate magnitude overflow");
    return magnitude;
}

// Compute the cell-centred SST strain-rate magnitude from the production
// velocity field using the configured baseline Gauss gradient operator.
// This function deliberately does not invent boundary gradient values: it uses
// the same gradient contract as compute_gradient_gauss.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
compute_sst_strain_rate(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Mesh& mesh)
{
    using namespace cfdx::core;
    if (velocity.dimension() != 3 || velocity.size() != mesh.n_cells())
        throw std::invalid_argument("compute_sst_strain_rate: velocity must be a 3-component cell field");

    const GeometryCache geometry = make_geometry_cache(mesh);
    std::array<Field<double, Location::CELL>, 3> component_gradients;
    for (std::size_t component = 0; component < 3; ++component) {
        Field<double, Location::CELL> scalar(
            mesh.n_cells(), "velocity_component", "m/s", 1);
        const double* source = velocity.component_data(component);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            if (!std::isfinite(source[cell]))
                throw std::invalid_argument("compute_sst_strain_rate: non-finite velocity");
            scalar(cell) = source[cell];
        }
        component_gradients[component] =
            compute_gradient_gauss(scalar, mesh, geometry);
    }

    Field<double, Location::CELL> strain(
        mesh.n_cells(), "sst_strain_rate", "1/s", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        std::array<double, 9> grad{};
        for (std::size_t velocity_component = 0; velocity_component < 3; ++velocity_component) {
            for (std::size_t direction = 0; direction < 3; ++direction) {
                grad[3 * velocity_component + direction] =
                    component_gradients[velocity_component].component_data(direction)[cell];
            }
        }
        strain(cell) = sst_strain_rate_magnitude(grad);
    }
    return strain;
}

} // namespace cfdx::physics
