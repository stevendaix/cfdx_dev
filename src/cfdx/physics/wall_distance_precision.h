#pragma once

// Backwards-compatible exhaustive precision interface. In V5 a request for
// exact correction replaces the ENTIRE returned field by an explicitly named
// BVH field; partial per-point blending is prohibited. See wall_distance_service.h
// for sampled monitoring, method traits, and a production fail-closed API.
#include "cfdx/physics/wall_distance_service.h"

#include <cstddef>
#include <string>

namespace cfdx::physics {

struct WallDistanceAccuracyPolicy {
    double absolute_tolerance{0.0};
    double relative_tolerance{1e-8};
    // Historic spelling retained for callers; now means FULL-field fallback.
    bool correct_with_exact_bvh{false};
};

struct WallDistanceAccuracyReport {
    WallDistanceResult original;
    WallDistanceResult accepted;
    std::size_t fluid_samples{0};
    std::size_t violations{0};
    // Number of nodes replaced when full BVH fallback is selected, not only
    // those that failed a tolerance: avoids the old misleading interpretation.
    std::size_t corrected_samples{0};
    double original_l2_relative{0.0};
    double original_linf_absolute{0.0};
    double accepted_linf_absolute{0.0};
    // The qualified field may be EXACT BVH even if the original method failed.
    bool accuracy_qualified{false};
    bool original_algorithm_qualified{false};
    bool entire_field_replaced_by_bvh{false};
    bool uses_exact_bvh_check{true};
};

inline WallDistanceAccuracyReport compute_wall_distance_accuracy_checked(
    WallDistanceMethod method, const WallSurface& surface,
    const WallDistanceGrid& grid, std::size_t max_iterations,
    const WallDistanceAccuracyPolicy& policy = {})
{
    WallDistanceServiceOptions options;
    options.method=method;
    options.max_iterations=max_iterations;
    options.check_scope=WallDistanceCheckScope::ALL_FLUID_NODES;
    options.absolute_tolerance=policy.absolute_tolerance;
    options.relative_tolerance=policy.relative_tolerance;
    options.on_failure=policy.correct_with_exact_bvh
        ? WallDistanceFailureAction::EXPLICIT_FULL_BVH
        : WallDistanceFailureAction::REJECT;
    const auto service=audit_wall_distance_service(surface,grid,options);
    WallDistanceAccuracyReport report;
    report.original=service.raw;
    report.accepted=service.selected;
    report.fluid_samples=service.quality.fluid_nodes;
    report.violations=service.quality.geometry_tolerance_violations;
    report.original_l2_relative=service.quality.relative_l2_error_checked;
    report.original_linf_absolute=service.quality.absolute_linf_error_checked;
    report.entire_field_replaced_by_bvh=service.quality.entirely_replaced_by_bvh;
    report.corrected_samples=report.entire_field_replaced_by_bvh
        ? report.fluid_samples : 0;
    report.accepted_linf_absolute=report.entire_field_replaced_by_bvh
        ? 0.0 : report.original_linf_absolute;
    report.original_algorithm_qualified=service.quality.algorithm_qualified &&
        service.quality.geometry_within_tolerance &&
        service.quality.full_geometry_check;
    report.accuracy_qualified=service.quality.approved_for_euclidean_wall_distance;
    if (!report.accuracy_qualified) {
        report.accepted.converged=false;
        report.accepted.stopping_reason="geometry_accuracy_unqualified";
        report.accepted.convergence_status=WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH;
    }
    return report;
}

} // namespace cfdx::physics
