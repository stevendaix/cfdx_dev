#pragma once

// Experimental wall-face post-processing for incompressible flow.
// Gradient consistency is verified by manufactured solutions. This helper
// does not change the momentum/pressure discretization used by the solver.

#include "cfdx/core/field/field.h"
#include <array>
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

struct WallFaceReconstruction {
    double pressure = 0.0;
    // Row i is grad(U_i). Units: 1/s for velocity measured in m/s.
    std::array<cfdx::core::Vec3, 3> velocity_gradient{};
    // Fluid -> solid normal. Force exerted by fluid on solid is
    // (p I - tau) . n, with tau = mu (grad(U) + grad(U)^T).
    cfdx::core::Vec3 force_per_area{};
};

inline WallFaceReconstruction reconstruct_wall_face_linear(
    const cfdx::core::Vec3& owner_center,
    const cfdx::core::Vec3& face_center,
    const cfdx::core::Vec3& area_vector_fluid_outward,
    double owner_pressure,
    const cfdx::core::Vec3& owner_pressure_gradient,
    const cfdx::core::Vec3& owner_velocity,
    const cfdx::core::Vec3& prescribed_wall_velocity,
    const std::array<cfdx::core::Vec3, 3>& owner_velocity_gradient,
    double dynamic_viscosity)
{
    const auto finite_vec = [](const cfdx::core::Vec3& v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    if (!finite_vec(owner_center) || !finite_vec(face_center) ||
        !finite_vec(area_vector_fluid_outward) ||
        !finite_vec(owner_pressure_gradient) || !finite_vec(owner_velocity) ||
        !finite_vec(prescribed_wall_velocity) ||
        !std::isfinite(owner_pressure) || !std::isfinite(dynamic_viscosity) ||
        !(dynamic_viscosity >= 0.0))
        throw std::invalid_argument("wall face reconstruction: invalid input");

    const cfdx::core::Vec3 displacement = face_center - owner_center;
    const double dist2 = displacement.mag2();
    const double area = area_vector_fluid_outward.mag();
    if (!(dist2 > 0.0) || !std::isfinite(dist2) ||
        !(area > 0.0) || !std::isfinite(area) ||
        !(displacement.dot(area_vector_fluid_outward) > 0.0))
        throw std::invalid_argument("wall face reconstruction: invalid wall geometry");
    for (const auto& g : owner_velocity_gradient)
        if (!finite_vec(g))
            throw std::invalid_argument("wall face reconstruction: invalid velocity gradient");

    WallFaceReconstruction result;
    result.pressure = owner_pressure + owner_pressure_gradient.dot(displacement);
    // Correct the cell gradient with a minimum-norm rank-one update so that
    // G_face (x_face - x_owner) = U_wall - U_owner. Linear exactness holds if
    // the input gradient is affine-exact. This is NOT a substitute for a
    // validated LS/WLS gradient and a high-quality wall mesh.
    const double u_owner[3] = {owner_velocity.x, owner_velocity.y, owner_velocity.z};
    const double u_wall[3] = {
        prescribed_wall_velocity.x, prescribed_wall_velocity.y, prescribed_wall_velocity.z};
    for (std::size_t i = 0; i < 3; ++i) {
        const double correction =
            (u_wall[i] - u_owner[i] - owner_velocity_gradient[i].dot(displacement)) / dist2;
        result.velocity_gradient[i] = owner_velocity_gradient[i] + displacement * correction;
    }
    const auto& g = result.velocity_gradient;
    const cfdx::core::Vec3 n = area_vector_fluid_outward * (1.0 / area);
    // tau . n, using the full symmetric Newtonian stress tensor.
    const cfdx::core::Vec3 tau_n{
        dynamic_viscosity * (2.0*g[0].x*n.x + (g[0].y + g[1].x)*n.y + (g[0].z + g[2].x)*n.z),
        dynamic_viscosity * ((g[1].x + g[0].y)*n.x + 2.0*g[1].y*n.y + (g[1].z + g[2].y)*n.z),
        dynamic_viscosity * ((g[2].x + g[0].z)*n.x + (g[2].y + g[1].z)*n.y + 2.0*g[2].z*n.z)};
    result.force_per_area = n * result.pressure - tau_n;
    if (!std::isfinite(result.pressure) || !finite_vec(result.force_per_area))
        throw std::runtime_error("wall face reconstruction: non-finite face traction");
    return result;
}

} // namespace cfdx::physics
