#pragma once

#include "cfdx/core/mesh/mesh.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct InitialFaceFluxContinuityReport {
    double local_l1 = 0.0;
    double face_flux_l1 = 0.0;
    double normalization_scale = 0.0;
    double normalized = 0.0;
};

// Sum of absolute cell-wise imbalances (NOT merely the global signed balance),
// normalized by prescribed boundary throughput whenever it exists (preserves
// the original open-flow certification), otherwise by actual face throughput.
// Enclosed recirculation therefore has a scale even with zero wall throughput.
// A truly motionless domain is admissible only with zero local imbalance.
inline InitialFaceFluxContinuityReport initial_face_flux_continuity_report(
    const cfdx::core::Mesh& mesh, const std::vector<double>& flux,
    double prescribed_boundary_throughput = 0.0)
{
    if (!std::isfinite(prescribed_boundary_throughput) ||
        prescribed_boundary_throughput < 0.0)
        throw std::invalid_argument("initial face flux: invalid boundary throughput");
    if (flux.size() != mesh.n_faces())
        throw std::invalid_argument("initial face flux: wrong face count");

    std::vector<long double> divergence(mesh.n_cells(), 0.0L);
    long double face_l1 = 0.0L;
    for (std::size_t face = 0; face < flux.size(); ++face) {
        const double mdot = flux[face];
        if (!std::isfinite(mdot))
            throw std::invalid_argument("initial face flux: non-finite value");
        const std::size_t owner = mesh.ownership().owner(face);
        if (owner >= mesh.n_cells())
            throw std::invalid_argument("initial face flux: invalid owner");
        const auto neighbour = mesh.ownership().neighbour(face);
        if (neighbour >= 0 &&
            static_cast<std::size_t>(neighbour) >= mesh.n_cells())
            throw std::invalid_argument("initial face flux: invalid neighbour");
        const long double value = static_cast<long double>(mdot);
        face_l1 += std::abs(value);
        divergence[owner] += value;
        if (neighbour >= 0)
            divergence[static_cast<std::size_t>(neighbour)] -= value;
    }
    long double local_l1 = 0.0L;
    for (const auto imbalance : divergence)
        local_l1 += std::abs(imbalance);

    // Do not use an arbitrary 1e-30 fallback: it makes a closed cavity's
    // normalized imbalance depend on physical units instead of its fluxes.
    const long double scale = prescribed_boundary_throughput > 0.0
        ? static_cast<long double>(prescribed_boundary_throughput)
        : face_l1;
    const long double ratio = scale > 0.0L
        ? local_l1 / scale
        : (local_l1 == 0.0L ? 0.0L
                            : std::numeric_limits<long double>::infinity());
    if (!std::isfinite(ratio) ||
        local_l1 > static_cast<long double>(std::numeric_limits<double>::max()) ||
        face_l1 > static_cast<long double>(std::numeric_limits<double>::max()) ||
        scale > static_cast<long double>(std::numeric_limits<double>::max()))
        throw std::invalid_argument("initial face flux: invalid continuity normalization");
    return {static_cast<double>(local_l1),
            static_cast<double>(face_l1), static_cast<double>(scale),
            static_cast<double>(ratio)};
}

} // namespace cfdx::physics
