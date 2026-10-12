// Issue #461 — N3 diffusion mesh-quality diagnostics.
//
// The diagnostics are deliberately independent of the Laplacian assembly.
// They expose the geometry population seen by the diffusion operator and
// identify faces on which the over-relaxed decomposition is singular or
// numerically ill-conditioned.
//
// No acceptance decision is made here: the caller owns the thresholds.

#pragma once

#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/laplacian.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx {
namespace core {

struct DiffusionQualityDiagnostics {
    std::size_t internal_faces = 0;
    std::size_t non_orthogonal_faces = 0;
    std::size_t faces_above_60_deg = 0;
    std::size_t faces_above_75_deg = 0;
    std::size_t faces_above_85_deg = 0;
    std::size_t faces_above_skew_05 = 0;
    std::size_t faces_above_skew_10 = 0;
    std::size_t singular_over_relaxed_faces = 0;

    double mean_non_orthogonality_deg = 0.0;
    double max_non_orthogonality_deg = 0.0;
    double mean_skewness = 0.0;
    double max_skewness = 0.0;
};

inline DiffusionQualityDiagnostics diagnose_diffusion_quality(
    const Mesh& mesh,
    const GeometryCache& geometry,
    double singular_relative_threshold = 1e-12)
{
    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("diagnose_diffusion_quality: invalid geometry cache");
    if (!(singular_relative_threshold > 0.0) ||
        !std::isfinite(singular_relative_threshold))
        throw std::invalid_argument(
            "diagnose_diffusion_quality: singular threshold must be finite and > 0");

    DiffusionQualityDiagnostics out;
    const FaceOwnership& ownership = mesh.ownership();

    double sum_angle = 0.0;
    double sum_skew = 0.0;

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::int64_t neighbour = ownership.neighbour(f);
        if (neighbour < 0)
            continue;

        const std::size_t owner = ownership.owner(f);
        const std::size_t nb = static_cast<std::size_t>(neighbour);
        if (owner >= mesh.n_cells() || nb >= mesh.n_cells())
            throw std::runtime_error(
                "diagnose_diffusion_quality: internal-face ownership out of range");

        ++out.internal_faces;

        const double angle = geometry.face_non_orthogonality_deg[f];
        const double skew = geometry.face_skewness[f];
        if (!std::isfinite(angle) || !std::isfinite(skew))
            throw std::runtime_error(
                "diagnose_diffusion_quality: non-finite face-quality metric");

        sum_angle += angle;
        sum_skew += skew;
        out.max_non_orthogonality_deg = std::max(
            out.max_non_orthogonality_deg, angle);
        out.max_skewness = std::max(out.max_skewness, skew);

        if (angle > 0.0)
            ++out.non_orthogonal_faces;
        if (angle > 60.0)
            ++out.faces_above_60_deg;
        if (angle > 75.0)
            ++out.faces_above_75_deg;
        if (angle > 85.0)
            ++out.faces_above_85_deg;
        if (skew > 0.5)
            ++out.faces_above_skew_05;
        if (skew > 1.0)
            ++out.faces_above_skew_10;

        const Vec3 d = geometry.cell_centres[nb] - geometry.cell_centres[owner];
        const Vec3 sf = geometry.face_Sf[f];
        const double d2 = d.dot(d);
        const double sf2 = sf.dot(sf);
        const double sd = sf.dot(d);
        if (!(d2 > 0.0) || !(sf2 > 0.0) ||
            !std::isfinite(d2) || !std::isfinite(sf2) || !std::isfinite(sd))
            throw std::runtime_error(
                "diagnose_diffusion_quality: invalid internal-face geometry");

        if (std::abs(sd) <= singular_relative_threshold * std::sqrt(sf2 * d2))
            ++out.singular_over_relaxed_faces;
    }

    if (out.internal_faces > 0) {
        const double inv = 1.0 / static_cast<double>(out.internal_faces);
        out.mean_non_orthogonality_deg = sum_angle * inv;
        out.mean_skewness = sum_skew * inv;
    }

    return out;
}

}  // namespace core
}  // namespace cfdx
