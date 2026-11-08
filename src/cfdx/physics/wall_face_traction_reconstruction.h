#pragma once

// Opt-in wall traction reconstruction for verification/qualification.
// Implements a linearly exact face-Dirichlet velocity correction and a
// tangential extrapolation of pressure at a homogeneous-Neumann wall.
// Both components require independent mesh-convergence validation before use
// in production force gates.

#include "cfdx/core/field/field.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

struct WallFaceTraction {
    double pressure = 0.0;
    cfdx::core::Vec3 viscous_on_fluid_normal{}; // tau * n, fluid-outward normal
};

inline WallFaceTraction reconstruct_wall_face_traction(
    const cfdx::core::Vec3& d_owner_to_face,
    const cfdx::core::Vec3& fluid_outward_unit_normal,
    double cell_pressure,
    const cfdx::core::Vec3& cell_pressure_gradient,
    const cfdx::core::Vec3& cell_velocity,
    const cfdx::core::Vec3& wall_velocity,
    const std::array<cfdx::core::Vec3, 3>& cell_velocity_gradients,
    double dynamic_viscosity,
    double prescribed_wall_pressure_normal_derivative = 0.0)
{
    using cfdx::core::Vec3;
    const Vec3& d = d_owner_to_face;
    const Vec3& n = fluid_outward_unit_normal;
    const double d2 = d.dot(d);
    const double n2 = n.dot(n);
    if (!(d2 > 0.0) || !std::isfinite(d2) ||
        !(std::abs(n2-1.0) <= 1e-10) || !std::isfinite(n2) ||
        !(dynamic_viscosity >= 0.0) || !std::isfinite(dynamic_viscosity) ||
        !std::isfinite(cell_pressure) ||
        !std::isfinite(prescribed_wall_pressure_normal_derivative))
        throw std::invalid_argument("wall face traction: invalid geometry or material");
    const auto finite_vec = [](const Vec3& v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    if (!finite_vec(d) || !finite_vec(n) || !finite_vec(cell_pressure_gradient) ||
        !finite_vec(cell_velocity) || !finite_vec(wall_velocity))
        throw std::invalid_argument("wall face traction: non-finite input");
    for (const auto& g : cell_velocity_gradients)
        if (!finite_vec(g)) throw std::invalid_argument("wall face traction: non-finite gradient");

    // Enforce the *declared* wall-normal pressure derivative rather than
    // silently extending a possibly erroneous cell gradient through the wall.
    // The legacy call (last argument omitted) explicitly assumes dp/dn=0.
    // That is a pressure-BC modelling assumption, NOT the exact NSE wall
    // identity on a curved viscous body; use the explicit derivative for MMS
    // and for A/B diagnostics when a consistent derivative is available.
    const Vec3 d_tangent = d - n * d.dot(n);
    const double p_face = cell_pressure + cell_pressure_gradient.dot(d_tangent)
                          + prescribed_wall_pressure_normal_derivative * d.dot(n);

    // Affine-exact correction: G_face * d == U_wall - U_owner. This enforces
    // the actual wall Dirichlet value without replacing the tangent gradients.
    std::array<Vec3, 3> corrected = cell_velocity_gradients;
    const double cell_u[3] = {cell_velocity.x, cell_velocity.y, cell_velocity.z};
    const double wall_u[3] = {wall_velocity.x, wall_velocity.y, wall_velocity.z};
    for (std::size_t i = 0; i < 3; ++i) {
        const double defect = wall_u[i] - cell_u[i] - corrected[i].dot(d);
        corrected[i] += d * (defect/d2);
    }
    const Vec3 t{
        dynamic_viscosity * (2.0*corrected[0].x*n.x +
                             (corrected[0].y+corrected[1].x)*n.y +
                             (corrected[0].z+corrected[2].x)*n.z),
        dynamic_viscosity * ((corrected[1].x+corrected[0].y)*n.x +
                             2.0*corrected[1].y*n.y +
                             (corrected[1].z+corrected[2].y)*n.z),
        dynamic_viscosity * ((corrected[2].x+corrected[0].z)*n.x +
                             (corrected[2].y+corrected[1].z)*n.y +
                             2.0*corrected[2].z*n.z)};
    if (!std::isfinite(p_face) || !finite_vec(t))
        throw std::runtime_error("wall face traction: non-finite result");
    return {p_face,t};
}

}  // namespace cfdx::physics
