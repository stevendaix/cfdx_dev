#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

// A single maximum RHS across Ux/Uy/Uz can hide an unconverged transverse
// momentum equation. Check components whose RHS has an independently resolved
// physical scale, while retaining the global absolute/RHS relative gate for
// nearly dormant components (e.g. Uz in an extruded 2D mesh).
struct MomentumResidualCertificate {
    double global_relative = 0.0;
    double max_active_component_relative = 0.0;
    std::array<double, 3> component_relative{};
    std::array<bool, 3> active{};
    bool converged = false;
};

inline MomentumResidualCertificate certify_momentum_residual_components(
    const std::array<double, 3>& raw_residuals,
    const std::array<double, 3>& rhs_scales,
    double tolerance,
    double active_rhs_fraction = 1.e-8)
{
    if (!(tolerance > 0.) || !std::isfinite(tolerance) ||
        !(active_rhs_fraction > 0.) || !(active_rhs_fraction < 1.) ||
        !std::isfinite(active_rhs_fraction))
        throw std::invalid_argument("momentum residual gate: invalid tolerance or component scale fraction");
    double global_rhs = 0.0;
    double global_residual = 0.0;
    for (std::size_t d = 0; d < 3; ++d) {
        if (!std::isfinite(raw_residuals[d]) || raw_residuals[d] < 0. ||
            !std::isfinite(rhs_scales[d]) || rhs_scales[d] < 0.)
            throw std::invalid_argument("momentum residual gate: invalid component data");
        global_rhs = std::max(global_rhs, rhs_scales[d]);
        global_residual = std::max(global_residual, raw_residuals[d]);
    }
    if (!(global_rhs > 0.))
        throw std::invalid_argument("momentum residual gate: all RHS scales are zero");
    MomentumResidualCertificate report;
    report.global_relative = global_residual / global_rhs;
    for (std::size_t d = 0; d < 3; ++d) {
        report.active[d] = rhs_scales[d] >= active_rhs_fraction * global_rhs;
        report.component_relative[d] = raw_residuals[d] /
            std::max(rhs_scales[d], active_rhs_fraction * global_rhs);
        if (report.active[d])
            report.max_active_component_relative = std::max(
                report.max_active_component_relative, report.component_relative[d]);
    }
    report.converged = std::isfinite(report.global_relative) &&
        std::isfinite(report.max_active_component_relative) &&
        report.global_relative <= tolerance &&
        report.max_active_component_relative <= tolerance;
    return report;
}

} // namespace cfdx::physics
