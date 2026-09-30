#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::core {

struct StencilQuality {
    std::size_t samples = 0;
    int dimension = 0;
    int rank = 0;
    double min_pivot = 0.0;
    double max_pivot = 0.0;
    double condition_estimate = std::numeric_limits<double>::infinity();
    double scale = 0.0;
    bool full_rank = false;
    bool finite = true;

    bool valid() const {
        return finite && rank > 0 && std::isfinite(condition_estimate);
    }
    bool well_conditioned(double limit) const {
        return valid() && condition_estimate <= limit;
    }
};

inline void validate_stencil_quality(const StencilQuality& q,
                                     bool require_full_rank = true)
{
    if (!q.finite)
        throw std::runtime_error("gradient stencil contains non-finite geometry");
    if (q.rank <= 0)
        throw std::runtime_error("gradient stencil is rank deficient");
    if (require_full_rank && !q.full_rank)
        throw std::runtime_error("gradient stencil does not span the requested dimension");
}

} // namespace cfdx::core
