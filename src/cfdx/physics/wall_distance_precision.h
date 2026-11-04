#pragma once

// Opt-in wall-distance accuracy controls; deliberately does not alter the
// default solver selection or claim that a PDE residual certifies geometry.
#include "cfdx/physics/wall_distance.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// The selectable second-order fast-sweeping implementation now lives in
// wall_distance.h so compute_wall_distance(EIKONAL_SECOND_ORDER,...) and
// standalone users execute exactly the same implementation.


struct WallDistanceAccuracyPolicy {
    // Distance error is bounded at EACH fluid grid point by
    // absolute_tolerance + relative_tolerance * exact_distance.
    // Distances are measured relative to the input TRIANGULATED surface,
    // not an unprovided CAD or smooth-geometry oracle.
    double absolute_tolerance{0.0};
    double relative_tolerance{1e-8};
    bool correct_with_exact_bvh{false};
};

struct WallDistanceAccuracyReport {
    WallDistanceResult original;
    WallDistanceResult accepted;
    std::size_t fluid_samples{0};
    std::size_t violations{0};
    std::size_t corrected_samples{0};
    double original_l2_relative{0.0};
    double original_linf_absolute{0.0};
    double accepted_linf_absolute{0.0};
    bool accuracy_qualified{false};
    bool uses_exact_bvh_check{true};
};

// Exhaustive (not sampled) accuracy gate. If a grid point fails tolerance,
// correct it using the exact BVH or fail closed. Because EVERY fluid point is
// queried, this is a verified tolerance bound on the output *samples*.
// The cost necessarily includes an O(N_grid * BVH query) verification pass.
inline WallDistanceAccuracyReport compute_wall_distance_accuracy_checked(
    WallDistanceMethod method, const WallSurface& surface,
    const WallDistanceGrid& grid, std::size_t max_iterations,
    const WallDistanceAccuracyPolicy& policy = {})
{
    if (!std::isfinite(policy.absolute_tolerance) ||
        policy.absolute_tolerance < 0.0 ||
        !std::isfinite(policy.relative_tolerance) ||
        policy.relative_tolerance < 0.0)
        throw std::invalid_argument("wall distance accuracy policy: invalid tolerance");
    if (grid.points.empty() || grid.solid.size() != grid.points.size())
        throw std::invalid_argument("wall distance accuracy policy: invalid grid");
    WallDistanceAccuracyReport report;
    report.original = compute_wall_distance(method,surface,grid,max_iterations);
    report.accepted = report.original;
    if (report.original.distance.size() != grid.points.size())
        throw std::runtime_error("wall distance accuracy policy: wrong field size");
    const WallDistanceBvh bvh(surface);
    double err2=0.0, ref2=0.0;
    for (std::size_t i=0; i<grid.points.size(); ++i) {
        if (grid.solid[i]) continue;
        ++report.fluid_samples;
        const double exact = bvh.nearest_distance(grid.points[i]);
        if (!std::isfinite(exact) || exact < 0.0)
            throw std::runtime_error("wall distance accuracy policy: invalid exact reference");
        const double value = report.original.distance[i];
        const double error = std::isfinite(value) && value >= 0.0 ? std::abs(value-exact) :
            std::numeric_limits<double>::infinity();
        const double tolerance = policy.absolute_tolerance + policy.relative_tolerance * exact;
        ref2 += exact*exact;
        err2 += error*error;
        report.original_linf_absolute = std::max(report.original_linf_absolute,error);
        if (!(error <= tolerance)) {
            ++report.violations;
            if (policy.correct_with_exact_bvh) {
                report.accepted.distance[i] = exact;
                report.accepted.valid[i] = 1;
                ++report.corrected_samples;
            }
        }
        const double accepted = report.accepted.distance[i];
        report.accepted_linf_absolute = std::max(report.accepted_linf_absolute,
            std::isfinite(accepted) ? std::abs(accepted-exact) :
            std::numeric_limits<double>::infinity());
    }
    report.original_l2_relative = report.fluid_samples != 0 ?
        std::sqrt(err2/std::max(1e-30,ref2)) :
        std::numeric_limits<double>::infinity();
    report.accuracy_qualified = report.fluid_samples != 0 &&
        (report.violations == 0 ||
         (policy.correct_with_exact_bvh && report.corrected_samples == report.violations));
    if (report.corrected_samples > 0) {
        report.accepted.method = report.original.method + "+exact_bvh_corrected";
        report.accepted.stopping_reason = "geometry_accuracy_guarded";
        report.accepted.convergence_status = WallDistanceConvergenceStatus::DIRECT;
        report.accepted.converged = report.accuracy_qualified;
    }
    if (!report.accuracy_qualified) {
        report.accepted.converged = false;
        report.accepted.stopping_reason = "geometry_accuracy_unqualified";
        report.accepted.convergence_status = WallDistanceConvergenceStatus::RESIDUAL_TOO_HIGH;
    }
    return report;
}

} // namespace cfdx::physics
