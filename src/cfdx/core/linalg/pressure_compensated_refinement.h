#pragma once

// High/low pressure correction refinement on the ORIGINAL, unscaled CSR.
// This module does not change CFD acceptance. A successful algebraic split
// does NOT certify pressure as a single FP64 field, full Rhie-Chow face flux,
// momentum balance, or turbulence. Each residual row is independent and can
// later be executed as a deterministic parallel CPU/GPU reduction.

#include "compensated_pressure_balance.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::core {

struct CompensatedRefinementReport {
    bool valid = false;
    bool algebraic_converged = false;
    bool physical_face_flux_certified = false;
    std::size_t passes_attempted = 0;
    std::size_t passes_accepted = 0;
    double initial_relative_residual = std::numeric_limits<double>::infinity();
    double final_relative_residual = std::numeric_limits<double>::infinity();
    double final_double_only_residual = std::numeric_limits<double>::infinity();
};

// Stable pressure Laplacian multiplication: never first round (high+low) or
// compute the delicate large diagonal and off-diagonal sums separately.
// Keeping an explicit row sum also makes spurious FP64 pressure leaks visible.
inline std::pair<Vector, double> compensated_pressure_residual(
    const SparseMatrix& a, const Vector& b,
    const Vector& high, const Vector& low) {
    const std::size_t n = a.n_rows();
    if (n == 0 || n != a.n_cols() || n != b.size() ||
        n != high.size() || n != low.size() ||
        std::numeric_limits<long double>::digits <=
            std::numeric_limits<double>::digits)
        throw std::invalid_argument("pressure refinement: incompatible or unsupported dimensions");
    const auto* offsets = a.row_offsets_data();
    const auto* cols = a.columns_data();
    const auto* values = a.values_data();
    Vector r(n, 0.0);
    long double norm2 = 0.0L, bnorm2 = 0.0L;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(high(i)) || !std::isfinite(low(i)) ||
            !std::isfinite(b(i)))
            throw std::domain_error("pressure refinement: nonfinite input");
        bnorm2 += static_cast<long double>(b(i)) * b(i);
        long double row_sum = 0.0L, difference = 0.0L;
        std::size_t n_diagonal = 0;
        for (std::uint32_t k = offsets[i]; k < offsets[i + 1]; ++k) {
            const std::size_t j = cols[k];
            const double coeff = values[k];
            if (j >= n || !std::isfinite(coeff))
                throw std::domain_error("pressure refinement: nonfinite CSR");
            row_sum += static_cast<long double>(coeff);
            if (i == j) {
                ++n_diagonal;
                continue;
            }
            difference += static_cast<long double>(coeff) *
                ((static_cast<long double>(high(j)) - high(i)) +
                 (static_cast<long double>(low(j)) - low(i)));
        }
        if (n_diagonal != 1)
            throw std::domain_error("pressure refinement: invalid diagonal");
        const long double ri = static_cast<long double>(b(i)) - difference -
            row_sum * (static_cast<long double>(high(i)) + low(i));
        if (!std::isfinite(ri) ||
            std::abs(ri) > std::numeric_limits<double>::max())
            throw std::domain_error("pressure refinement: residual overflow");
        r(i) = static_cast<double>(ri);
        norm2 += ri * ri;
    }
    if (bnorm2 <= 0.0L)
        throw std::invalid_argument("pressure refinement: nonzero RHS required");
    const double relative = static_cast<double>(std::sqrt(norm2 / bnorm2));
    if (!std::isfinite(relative))
        throw std::domain_error("pressure refinement: invalid residual norm");
    return {std::move(r), relative};
}

// Callback signature: bool(const SparseMatrix&, const Vector& residual,
//                           Vector& zero_initialised_delta).
// A callback may use a loose *inner* tolerance; only the independently
// recomputed original-system residual is eligible for algebraic convergence.
// It must never change A, b, tolerances of the physical CFD gate, or FV fluxes.
template <class IncrementSolver>
CompensatedRefinementReport refine_pressure_high_low(
    const SparseMatrix& a, const Vector& b, Vector& high, Vector& low,
    double original_tolerance, std::size_t max_passes,
    IncrementSolver&& solve_increment) {
    if (!(original_tolerance > 0.0) || !std::isfinite(original_tolerance) ||
        max_passes > 8)
        throw std::invalid_argument("pressure refinement: invalid controls");
    auto [residual, relative] = compensated_pressure_residual(a, b, high, low);
    CompensatedRefinementReport result;
    result.valid = true;
    result.initial_relative_residual = relative;
    result.final_relative_residual = relative;
    for (std::size_t pass = 0; pass < max_passes &&
         relative > original_tolerance; ++pass) {
        ++result.passes_attempted;
        Vector increment(a.n_rows(), 0.0);
        if (!solve_increment(a, residual, increment)) break;
        Vector proposed_high = high;
        Vector proposed_low = low;
        bool finite = true;
        for (std::size_t i = 0; i < a.n_rows(); ++i)
            finite = compensated_add(proposed_high(i), proposed_low(i),
                                     increment(i)) && finite;
        if (!finite) break;
        const auto [proposed_residual, proposed_relative] =
            compensated_pressure_residual(a, b, proposed_high, proposed_low);
        // Stagnation is a measured failure, not an accepted iteration.
        if (!(proposed_relative < relative)) break;
        high = std::move(proposed_high);
        low = std::move(proposed_low);
        residual = proposed_residual;
        relative = proposed_relative;
        ++result.passes_accepted;
        result.final_relative_residual = relative;
    }
    result.algebraic_converged = relative <= original_tolerance;
    // The returned high/low pair is never silently rounded for certification.
    const Vector no_low(a.n_rows(), 0.0);
    result.final_double_only_residual =
        compensated_pressure_residual(a, b, high, no_low).second;
    return result;
}

} // namespace cfdx::core
