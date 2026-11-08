#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

// Opt-in Armijo-like acceptance filter for transactionally retried SIMPLE
// steps. The direction is a *recomputed Picard step* (not a Newton direction):
// the solver must not pretend that every trial is a descent direction.
// Never adjust a physical convergence threshold to make a trial pass.
struct PicardPhysicalGlobalizationControls {
    bool enabled = false;
    // Picard is not generally a descent direction during boundary-layer startup.
    // Begin safeguarding only after a warm-up of fully accepted iterates.
    std::size_t start_iteration = 2;
    // Required fractional reduction in physical merit per unit of the retry
    // relaxation scale. Zero still imposes monotonicity of the true residual.
    double sufficient_decrease = 1e-4;
    // Reject excessive deterioration of any individual conservation channel,
    // even when a different channel reduces the max-norm merit.
    double channel_growth_cap = 1.25;
    // Once every channel is within this multiple of the ORIGINAL convergence
    // gate, stop backtracking over floating-point noise, but do not relax
    // the original convergence gate itself.
    double noise_floor_gate_multiple = 2.0;
    // Optional bounded nonmonotone merit: an admissible trial may rise by a
    // SMALL factor relative to the worst of N recently accepted iterates.
    // This is NOT Newton Armijo and does not assert Picard is descent.
    // Window=1, growth=1 retain strict legacy behavior.
    std::size_t nonmonotone_window = 1;
    double maximum_window_growth = 1.0;
};

inline void validate_picard_physical_globalization(
    const PicardPhysicalGlobalizationControls& c)
{
    if (c.start_iteration < 2 ||
        !std::isfinite(c.sufficient_decrease) || c.sufficient_decrease < 0.0 ||
        c.sufficient_decrease >= 1.0 ||
        !std::isfinite(c.channel_growth_cap) || c.channel_growth_cap < 1.0 ||
        !std::isfinite(c.noise_floor_gate_multiple) ||
        c.noise_floor_gate_multiple < 1.0 ||
        c.nonmonotone_window < 1 || c.nonmonotone_window > 32 ||
        !std::isfinite(c.maximum_window_growth) ||
        c.maximum_window_growth < 1.0 || c.maximum_window_growth > 1.25)
        throw std::invalid_argument("invalid Picard physical globalization controls");
}

struct PicardPhysicalSample {
    std::array<double,3> momentum_absolute{};
    std::array<double,3> momentum_rhs_scale{};
    double continuity_l1 = 0.0;
    double conservative_face_flux_l1 = 0.0;
};

struct PicardPhysicalMerit {
    double momentum = 0.0;  // maximum of active components, scaled to old RHS
    double continuity = 0.0; // local absolute mass defect / old face throughput
    double combined = 0.0;   // maximum, each divided by its UNCHANGED gate
};

// Compare before and after under the SAME RHS and flux normalization: using
// the trial's potentially increased RHS/flux as denominator would let a poor
// trial masquerade as improved. A component with a vanishing reference RHS is
// scored on the larger old/new absolute scale, not silently dropped.
inline PicardPhysicalMerit picard_physical_merit(
    const PicardPhysicalSample& sample,
    const PicardPhysicalSample& reference,
    double momentum_gate,
    double continuity_gate)
{
    if (!(momentum_gate > 0.0 && continuity_gate > 0.0) ||
        !std::isfinite(momentum_gate) || !std::isfinite(continuity_gate))
        throw std::invalid_argument("invalid Picard globalization convergence gates");
    PicardPhysicalMerit out;
    for (std::size_t d = 0; d < 3; ++d) {
        const auto a = sample.momentum_absolute[d];
        const auto scale = reference.momentum_rhs_scale[d];
        const auto old = reference.momentum_absolute[d];
        if (!(a >= 0.0 && old >= 0.0 && scale >= 0.0) ||
            !std::isfinite(a) || !std::isfinite(old) || !std::isfinite(scale))
            throw std::invalid_argument("invalid Picard physical momentum evidence");
        if (std::max({a, old, scale}) < 1e-280) continue; // inactive 2-D component
        const double stable_scale = std::max({scale, old, 1e-30});
        out.momentum = std::max(out.momentum, a / stable_scale);
    }
    if (!(sample.continuity_l1 >= 0.0 && reference.continuity_l1 >= 0.0) ||
        !(reference.conservative_face_flux_l1 > 0.0) ||
        !std::isfinite(sample.continuity_l1) ||
        !std::isfinite(reference.continuity_l1) ||
        !std::isfinite(reference.conservative_face_flux_l1))
        throw std::invalid_argument("invalid Picard physical continuity evidence");
    out.continuity = sample.continuity_l1 /
        reference.conservative_face_flux_l1;
    out.combined = std::max(out.momentum/momentum_gate,
                            out.continuity/continuity_gate);
    if (!std::isfinite(out.combined))
        throw std::invalid_argument("nonfinite Picard globalization merit");
    return out;
}

inline bool picard_physical_accept(
    const PicardPhysicalMerit& before,
    const PicardPhysicalMerit& after,
    double relaxation_scale,
    const PicardPhysicalGlobalizationControls& controls)
{
    validate_picard_physical_globalization(controls);
    if (!(relaxation_scale > 0.0 && relaxation_scale <= 1.0) ||
        !std::isfinite(relaxation_scale) ||
        !std::isfinite(before.combined) || !std::isfinite(after.combined) ||
        before.combined < 0.0 || after.combined < 0.0)
        throw std::invalid_argument("invalid Picard globalization trial");
    if (before.combined <= controls.noise_floor_gate_multiple &&
        after.combined <= controls.noise_floor_gate_multiple)
        return true;
    const double limit = std::max(controls.noise_floor_gate_multiple,
        before.combined * controls.maximum_window_growth *
            (1.0 - controls.sufficient_decrease * relaxation_scale));
    return after.combined <= limit &&
        after.momentum <= std::max(controls.noise_floor_gate_multiple,
                                  controls.channel_growth_cap * before.momentum) &&
        after.continuity <= std::max(controls.noise_floor_gate_multiple,
                                    controls.channel_growth_cap * before.continuity);
}

} // namespace cfdx::physics

