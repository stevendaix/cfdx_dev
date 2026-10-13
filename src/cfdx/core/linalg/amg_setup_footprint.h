#pragma once

#include <cstddef>

namespace cfdx::core {

// Reporting-only counters for AMG setup decisions that can silently degrade
// coarse-grid quality. Nothing in the solver reads these values: they let
// tests and qualification campaigns observe how often a setup resorted to
// singleton aggregates, fallback interpolation rows, or dropped
// numerically-zero Galerkin entries.
struct AMGSetupFootprint {
    // Smoothed-aggregation aggregates containing exactly one fine point.
    std::size_t singleton_aggregates = 0;
    // Interpolation rows rebuilt as {own aggregate, 1.0} because the
    // smoothed weights summed to (numerically) zero.
    std::size_t fallback_interpolation_rows = 0;
    // Galerkin coarse entries dropped because the accumulated value was
    // non-finite or at/below the 1e-30 numeric-zero cut.
    std::size_t galerkin_dropped_entries = 0;
};

} // namespace cfdx::core
