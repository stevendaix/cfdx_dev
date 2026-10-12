#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::core::numerics {

struct EquivalenceMetrics {
    std::size_t size = 0;
    double l1_abs = 0.0;
    double l2_abs = 0.0;
    double linf_abs = 0.0;
    double l2_relative = 0.0;
    double max_relative = 0.0;
};

inline EquivalenceMetrics compare_vectors(
    const std::vector<double>& reference,
    const std::vector<double>& candidate,
    double scale_floor = 1e-30)
{
    if (reference.size() != candidate.size())
        throw std::invalid_argument("backend equivalence requires equal vector sizes");

    EquivalenceMetrics m;
    m.size = reference.size();

    double sum_abs = 0.0;
    double sum_sq = 0.0;
    double ref_sq = 0.0;
    double max_ref = 0.0;

    for (std::size_t i = 0; i < reference.size(); ++i) {
        const double delta = std::abs(candidate[i] - reference[i]);
        if (!std::isfinite(delta))
            throw std::runtime_error("non-finite backend equivalence difference");
        sum_abs += delta;
        sum_sq += delta * delta;
        ref_sq += reference[i] * reference[i];
        max_ref = std::max(max_ref, std::abs(reference[i]));
    }

    m.l1_abs = sum_abs;
    m.l2_abs = std::sqrt(sum_sq);
    m.linf_abs = reference.empty() ? 0.0 : m.linf_abs;
    for (std::size_t i = 0; i < reference.size(); ++i)
        m.linf_abs = std::max(m.linf_abs, std::abs(candidate[i] - reference[i]));

    const double l2_scale = std::max(std::sqrt(ref_sq), scale_floor);
    m.l2_relative = m.l2_abs / l2_scale;

    const double max_scale = std::max(max_ref, scale_floor);
    m.max_relative = m.linf_abs / max_scale;
    return m;
}

inline bool numerically_equivalent(
    const EquivalenceMetrics& m,
    double l2_relative_tolerance,
    double linf_relative_tolerance)
{
    if (!(l2_relative_tolerance >= 0.0) ||
        !(linf_relative_tolerance >= 0.0))
        throw std::invalid_argument("equivalence tolerances must be non-negative");
    return m.l2_relative <= l2_relative_tolerance &&
           m.max_relative <= linf_relative_tolerance;
}

inline void require_numerically_equivalent(
    const std::string& label,
    const std::vector<double>& reference,
    const std::vector<double>& candidate,
    double l2_relative_tolerance,
    double linf_relative_tolerance)
{
    const auto m = compare_vectors(reference, candidate);
    if (!numerically_equivalent(m, l2_relative_tolerance, linf_relative_tolerance)) {
        throw std::runtime_error(
            label + ": numerical equivalence failed: L2rel=" +
            std::to_string(m.l2_relative) + ", Linfrel=" +
            std::to_string(m.max_relative));
    }
}

} // namespace cfdx::core::numerics
