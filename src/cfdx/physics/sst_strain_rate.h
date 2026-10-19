#pragma once

#include "cfdx/core/numerics/gradient.h"
#include <array>
#include <cmath>
#include <stdexcept>

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
// velocity field using the baseline Gauss gradient operator. The gradient
// operator's existing boundary treatment is preserved, not replaced here.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
compute_sst_strain_rate(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Mesh& mesh)
{
    using namespace cfdx::core;
    if (velocity.dimension() != 3 || velocity.size() != mesh.n_cells())
        throw std::invalid_argument("compute_sst_strain_rate: velocity must be a 3-component cell field");

    const GeometryCache geometry = make_geometry_cache(mesh);
    Field<double, Location::CELL> u(mesh.n_cells(), "u", "m/s", 1);
    Field<double, Location::CELL> v(mesh.n_cells(), "v", "m/s", 1);
    Field<double, Location::CELL> w(mesh.n_cells(), "w", "m/s", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        for (std::size_t component = 0; component < 3; ++component) {
            if (!std::isfinite(velocity.component_data(component)[cell]))
                throw std::invalid_argument("compute_sst_strain_rate: non-finite velocity");
        }
        u(cell) = velocity.component_data(0)[cell];
        v(cell) = velocity.component_data(1)[cell];
        w(cell) = velocity.component_data(2)[cell];
    }
    const auto grad_u = compute_gradient_gauss(u, mesh, geometry);
    const auto grad_v = compute_gradient_gauss(v, mesh, geometry);
    const auto grad_w = compute_gradient_gauss(w, mesh, geometry);

    Field<double, Location::CELL> strain(
        mesh.n_cells(), "sst_strain_rate", "1/s", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const std::array<double, 9> grad{
            grad_u.component_data(0)[cell], grad_u.component_data(1)[cell],
            grad_u.component_data(2)[cell], grad_v.component_data(0)[cell],
            grad_v.component_data(1)[cell], grad_v.component_data(2)[cell],
            grad_w.component_data(0)[cell], grad_w.component_data(1)[cell],
            grad_w.component_data(2)[cell]};
        strain(cell) = sst_strain_rate_magnitude(grad);
    }
    return strain;
}

} // namespace cfdx::physics
