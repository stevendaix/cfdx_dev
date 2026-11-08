#pragma once

#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Conservative two-state defect-correction *iteration filter*, not a change
// of the final TVD operator. The previous input is either the last committed
// raw correction (default) or the last committed effective correction
// (optional recursive mode). Both are equal to the raw flux at a fixed point.
//   D_used = (1-weight)*D(raw) + weight*D(previous).
// The accepted-state momentum residual MUST always use the unfiltered,
// second-order reconstructed face flux. The optional source-consistency gate
// below rejects an apparently converged state if this cancellation is stale.
struct TvdDeferredDefectMemoryControls {
    bool enabled = false;
    double previous_weight = 0.5;
    bool recursive = false;
    // Apply blending ONLY on faces with a verified two-step correction
    // oscillation, to leave locally monotone high-order updates untouched.
    bool oscillation_only = false;
    bool require_source_consistency = true;
};

inline void validate_tvd_deferred_defect_memory(
    const TvdDeferredDefectMemoryControls& c)
{
    if ((c.recursive && c.oscillation_only) ||
        !std::isfinite(c.previous_weight) || c.previous_weight < 0.0 ||
        c.previous_weight >= 1.0)
        throw std::invalid_argument("TVD defect memory weight must be in [0,1)");
}

struct TvdDeferredMemoryResult {
    std::vector<double> current_raw;
    std::vector<double> used;
    double max_face_difference = 0.0;
    double max_cell_source_difference = 0.0;
    double signed_source_difference = 0.0;
    std::size_t filtered_faces = 0;
};

// Apply a two-level, conservative correction to the ALREADY assembled
// high-order explicit RHS. Both owner and neighbour receive opposite
// increments, so summing over all cells gives zero even off orthogonality.
// The input current_raw is the exact FV face correction (TVD - UPWIND)
// computed from the predictor U and its frozen conservative Rhie-Chow flux.
inline TvdDeferredMemoryResult apply_tvd_deferred_defect_memory(
    const cfdx::core::Mesh& mesh,
    ScalarEquation& equation,
    const std::vector<double>& current_raw,
    const std::vector<double>* previous_raw,
    double previous_weight,
    const std::vector<double>* two_steps_ago = nullptr,
    bool oscillation_only = false)
{
    const std::size_t nc = mesh.n_cells(), nf = mesh.n_faces();
    if (equation.rhs.size() != nc || current_raw.size() != nf ||
        (previous_raw && previous_raw->size() != nf) ||
        (two_steps_ago && two_steps_ago->size() != nf) ||
        (oscillation_only && previous_raw && !two_steps_ago) ||
        !std::isfinite(previous_weight) || previous_weight < 0.0 ||
        previous_weight >= 1.0)
        throw std::invalid_argument("TVD defect memory: incompatible sizes/weight");
    TvdDeferredMemoryResult result;
    result.current_raw = current_raw;
    result.used = current_raw;
    std::vector<double> cell_difference(nc, 0.0);
    for (std::size_t face = 0; face < nf; ++face) {
        const double value = current_raw[face];
        if (!std::isfinite(value))
            throw std::invalid_argument("TVD defect memory: nonfinite current correction");
        const auto neighbour = mesh.ownership().neighbour(face);
        if (neighbour < 0) {
            // High-order deferred correction of an internal face only. No
            // extrapolation across physical walls, fixed BCs or EMPTY faces.
            if (previous_raw && !std::isfinite((*previous_raw)[face]))
                throw std::invalid_argument("TVD defect memory: nonfinite old boundary correction");
            continue;
        }
        if (mesh.ownership().owner(face) >= nc ||
            static_cast<std::size_t>(neighbour) >= nc)
            throw std::invalid_argument("TVD defect memory: invalid mesh face ownership");
        if (!previous_raw || previous_weight == 0.0) continue;
        const double old = (*previous_raw)[face];
        if (!std::isfinite(old))
            throw std::invalid_argument("TVD defect memory: nonfinite old correction");
        double effective_weight = previous_weight;
        if (oscillation_only) {
            const double prior = (*two_steps_ago)[face];
            if (!std::isfinite(prior))
                throw std::invalid_argument("TVD defect memory: nonfinite two-step correction");
            const double now_change = value - old;
            const double past_change = old - prior;
            const double flux_scale = std::max({1e-12, std::abs(value),
                                                std::abs(old), std::abs(prior)});
            if (!(now_change * past_change < 0.0 &&
                  std::abs(now_change) > 1e-10*flux_scale &&
                  std::abs(past_change) > 1e-10*flux_scale))
                effective_weight = 0.0;
        }
        const double delta = effective_weight * (old - value);
        if (effective_weight > 0.0) ++result.filtered_faces;
        if (!std::isfinite(delta))
            throw std::invalid_argument("TVD defect memory: nonfinite filtered correction");
        // Original RHS contribution is -D on owner, +D on neighbour.
        // Replace D_current with D_current + delta without touching A.
        const auto owner = mesh.ownership().owner(face);
        equation.rhs(owner) -= delta;
        equation.rhs(static_cast<std::size_t>(neighbour)) += delta;
        cell_difference[owner] -= delta;
        cell_difference[static_cast<std::size_t>(neighbour)] += delta;
        result.used[face] = value + delta;
        result.max_face_difference = std::max(result.max_face_difference, std::abs(delta));
    }
    for (const double source : cell_difference) {
        result.max_cell_source_difference =
            std::max(result.max_cell_source_difference, std::abs(source));
        result.signed_source_difference += source;
    }
    return result;
}

} // namespace cfdx::physics

