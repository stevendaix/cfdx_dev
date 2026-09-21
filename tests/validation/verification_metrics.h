#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::verification {

struct ErrorMetrics {
    double l2 = 0.0;
    double linf = 0.0;
    double l2_relative = 0.0;
    double linf_relative = 0.0;
};

inline ErrorMetrics error_norms(const std::vector<double>& numerical,
                               const std::vector<double>& exact,
                               const std::vector<double>& weights = {})
{
    if (numerical.size() != exact.size() || numerical.empty())
        throw std::invalid_argument("error_norms: incompatible vectors");
    if (!weights.empty() && weights.size() != numerical.size())
        throw std::invalid_argument("error_norms: incompatible weights");

    double sum_w = 0.0;
    double sum_e2 = 0.0;
    double sum_exact2 = 0.0;
    double max_e = 0.0;
    double max_exact = 0.0;
    for (std::size_t i = 0; i < numerical.size(); ++i) {
        const double w = weights.empty() ? 1.0 : weights[i];
        if (!(w > 0.0) || !std::isfinite(w))
            throw std::invalid_argument("error_norms: invalid weight");
        const double e = std::abs(numerical[i] - exact[i]);
        sum_w += w;
        sum_e2 += w * e * e;
        sum_exact2 += w * exact[i] * exact[i];
        max_e = std::max(max_e, e);
        max_exact = std::max(max_exact, std::abs(exact[i]));
    }
    const double l2=std::sqrt(sum_e2/sum_w);
    const double exact_l2=std::sqrt(sum_exact2/sum_w);
    const double scale=std::max(exact_l2,std::numeric_limits<double>::min());
    const double inf_scale=std::max(max_exact,std::numeric_limits<double>::min());
    return {l2,max_e,l2/scale,max_e/inf_scale};
}

inline double observed_order(double coarse_error, double fine_error,
                             double refinement_ratio = 2.0)
{
    if (!(coarse_error > 0.0) || !(fine_error > 0.0) ||
        !(refinement_ratio > 1.0))
        return std::numeric_limits<double>::quiet_NaN();
    return std::log(coarse_error / fine_error) /
           std::log(refinement_ratio);
}

inline void require_order(const std::vector<double>& errors,
                          double expected_order,
                          double minimum_order,
                          const std::string& case_name)
{
    if (errors.size() < 2)
        throw std::invalid_argument("require_order: at least two errors required");

    for (std::size_t i = 1; i < errors.size(); ++i) {
        const double p = observed_order(errors[i - 1], errors[i]);
        if (!std::isfinite(p) || p < minimum_order) {
            throw std::runtime_error(
                case_name + ": measured order " + std::to_string(p) +
                " below required " + std::to_string(minimum_order) +
                " (expected " + std::to_string(expected_order) + ")");
        }
    }
}

} // namespace cfdx::verification
