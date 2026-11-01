#pragma once

#include "cfdx/core/field/field.h"

#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

// Experimental decomposition for pressure-response fluxes on skew meshes.
// This is a local algebraic operator, not a replacement for the coupled
// pressure / Rhie-Chow formulation. Both terms must use the SAME scheme on
// every face and be assembled with opposite signs in its two neighbour cells.
// The deferred term must reach its full value before qualification.
enum class PressureFluxSplitScheme {
    MINIMUM_CORRECTION,  // production-compatible orthogonal projection
    OVER_RELAXED         // stronger implicit two-point coupling
};

struct PressureFluxSplit {
    double implicit_coefficient = 0.0;
    double deferred_flux = 0.0;
    double exact_affine_flux = 0.0;
    double cosine_angle = 0.0;
};

inline PressureFluxSplit split_pressure_flux(
    const cfdx::core::Vec3& oriented_area,
    const cfdx::core::Vec3& owner_to_neighbour,
    const cfdx::core::Vec3& positive_directional_response,
    const cfdx::core::Vec3& face_gradient,
    PressureFluxSplitScheme scheme)
{
    using cfdx::core::Vec3;
    const Vec3& S = oriented_area;
    const Vec3& d = owner_to_neighbour;
    const Vec3& r = positive_directional_response;
    const double d2 = d.mag2();
    const double area2 = S.mag2();
    const double alignment = S.dot(d);
    if (!(d2 > 0.0) || !(area2 > 0.0) || !(alignment > 0.0) ||
        !std::isfinite(d2) || !std::isfinite(area2) || !std::isfinite(alignment) ||
        !(r.x > 0.0) || !(r.y > 0.0) || !(r.z > 0.0) ||
        !std::isfinite(r.x) || !std::isfinite(r.y) || !std::isfinite(r.z) ||
        !std::isfinite(face_gradient.x) || !std::isfinite(face_gradient.y) ||
        !std::isfinite(face_gradient.z))
        throw std::invalid_argument("split_pressure_flux: invalid oriented face or response");

    // Use e.r.e as the directional pressure response. Anisotropic response
    // generally leaves an additional *physical* deferred contribution.
    const double rnormal = (r.x*d.x*d.x + r.y*d.y*d.y + r.z*d.z*d.z)/d2;
    double coefficient = rnormal * alignment / d2;
    if (scheme == PressureFluxSplitScheme::OVER_RELAXED)
        coefficient = rnormal * area2 / alignment;
    const double full = r.x*face_gradient.x*S.x +
                        r.y*face_gradient.y*S.y +
                        r.z*face_gradient.z*S.z;
    const double deferred = full - coefficient*face_gradient.dot(d);
    const double cos_angle = alignment/std::sqrt(d2*area2);
    if (!(coefficient > 0.0) || !std::isfinite(coefficient) ||
        !std::isfinite(deferred) || !std::isfinite(full) ||
        !(cos_angle > 0.0) || cos_angle > 1.0 + 1.0e-12)
        throw std::runtime_error("split_pressure_flux: non-finite or inconsistent split");
    return {coefficient, deferred, full, cos_angle};
}

}  // namespace cfdx::physics
