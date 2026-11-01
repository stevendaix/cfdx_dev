#pragma once

// Independent, opt-in wall-shear estimator for qualification only.
// The near-wall tangential velocity is reconstructed to the wall-normal line
// through the face, then a two-point quadratic with a no-slip wall value is
// differentiated at the wall. The estimate must be mesh-converged before it
// can replace a production wall traction model.

#include "cfdx/core/field/field.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

struct WallShearSample {
    cfdx::core::Vec3 cell_position;
    cfdx::core::Vec3 cell_velocity;
    // Row i is gradient of velocity component i.
    std::array<cfdx::core::Vec3, 3> velocity_gradients{};
};

struct WallShearQuadraticResult {
    double shear_on_body_tangent = 0.0; // Pa, positive in the supplied tangent
    double first_distance = 0.0;
    double second_distance = 0.0;
};

inline WallShearQuadraticResult reconstruct_wall_shear_quadratic(
    const cfdx::core::Vec3& face_position,
    const cfdx::core::Vec3& outward_fluid_normal,
    const cfdx::core::Vec3& tangent,
    const cfdx::core::Vec3& wall_velocity,
    const WallShearSample& first,
    const WallShearSample& second,
    double dynamic_viscosity)
{
    using cfdx::core::Vec3;
    const auto finite = [](const Vec3& v) {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    };
    if (!finite(face_position) || !finite(outward_fluid_normal) ||
        !finite(tangent) || !finite(wall_velocity) ||
        !finite(first.cell_position) || !finite(first.cell_velocity) ||
        !finite(second.cell_position) || !finite(second.cell_velocity) ||
        !std::isfinite(dynamic_viscosity) || !(dynamic_viscosity > 0.0))
        throw std::invalid_argument("quadratic wall shear: non-finite inputs");
    for (const auto* sample : {&first, &second})
        for (const auto& g : sample->velocity_gradients)
            if (!finite(g))
                throw std::invalid_argument("quadratic wall shear: non-finite gradient");
    if (std::abs(outward_fluid_normal.mag2() - 1.0) > 1e-10 ||
        std::abs(tangent.mag2() - 1.0) > 1e-10 ||
        std::abs(outward_fluid_normal.dot(tangent)) > 1e-10)
        throw std::invalid_argument("quadratic wall shear: nonorthonormal directions");
    const Vec3 inward = outward_fluid_normal * -1.0;
    const auto projected = [&](const WallShearSample& sample) {
        const Vec3 r = sample.cell_position - face_position;
        const double h = r.dot(inward);
        if (!std::isfinite(h) || !(h > 0.0))
            throw std::invalid_argument("quadratic wall shear: cell outside fluid");
        const Vec3 r_tangent = r - inward * h;
        // Project each velocity component along the supplied surface tangent.
        const Vec3 grad_t = sample.velocity_gradients[0] * tangent.x +
                            sample.velocity_gradients[1] * tangent.y +
                            sample.velocity_gradients[2] * tangent.z;
        const double u = (sample.cell_velocity - wall_velocity).dot(tangent) -
                         grad_t.dot(r_tangent);
        if (!std::isfinite(u))
            throw std::invalid_argument("quadratic wall shear: non-finite projected velocity");
        return std::pair<double, double>{h, u};
    };
    const auto [h1, u1] = projected(first);
    const auto [h2, u2] = projected(second);
    // This is a genuine two-distance extrapolation. Never solve an almost
    // singular fit or silently fall back to an unverified first-order value.
    // Reject near-coincident depths: the second sample must represent a distinct
    // wall-normal layer. Also avoid very remote samples on curved surfaces.
    if (!(h2 >= 1.5 * h1) || !(h2 <= 6.0 * h1))
        throw std::invalid_argument("quadratic wall shear: ill-conditioned normal stencil");
    const double du_dn = ((u1 / h1) * h2 - (u2 / h2) * h1) / (h2 - h1);
    const double shear = dynamic_viscosity * du_dn;
    if (!std::isfinite(shear))
        throw std::runtime_error("quadratic wall shear: non-finite result");
    return {shear, h1, h2};
}

} // namespace cfdx::physics
