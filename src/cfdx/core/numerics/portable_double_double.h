#pragma once

// Portable two-component accumulation of binary64 numbers.  No assumption
// about long double width, platform endianness, or hardware FMA is needed.
// This is an accumulator, not a drop-in arbitrary-precision field type.
// Each operation uses error-free TwoSum before a final renormalization.
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cfdx::core::numerics {

struct PortableDoubleDouble {
    double high = 0.0;
    double low = 0.0;

    static std::pair<double, double> two_sum(double a, double b) {
        const double sum = a + b;
        const double b_virtual = sum - a;
        const double error = (a - (sum - b_virtual)) + (b - b_virtual);
        return {sum, error};
    }

    void add(double x) {
        if (!std::isfinite(x))
            throw std::invalid_argument("double-double: non-finite contribution");
        const auto [sum, roundoff] = two_sum(high, x);
        if (!std::isfinite(sum))
            throw std::overflow_error("double-double: accumulation overflow");
        const auto [new_high, new_low] = two_sum(sum, low + roundoff);
        if (!std::isfinite(new_high) || !std::isfinite(new_low))
            throw std::overflow_error("double-double: renormalization overflow");
        high = new_high;
        low = new_low;
    }

    // Adding two renormalized pairs retains the contribution lost in the
    // high-component sum and in the low-component sum.
    void add(const PortableDoubleDouble& x) {
        if (!std::isfinite(x.high) || !std::isfinite(x.low))
            throw std::invalid_argument("double-double: non-finite pair");
        const auto [sum, high_error] = two_sum(high, x.high);
        if (!std::isfinite(sum))
            throw std::overflow_error("double-double: high-component overflow");
        const auto [low_sum, low_error] = two_sum(low, x.low);
        const auto [merged, merged_error] = two_sum(sum, high_error + low_sum);
        const auto [new_high, new_low] = two_sum(merged, merged_error + low_error);
        if (!std::isfinite(new_high) || !std::isfinite(new_low))
            throw std::overflow_error("double-double: pair overflow");
        high = new_high;
        low = new_low;
    }

    PortableDoubleDouble operator-() const { return {-high, -low}; }

    double narrow_to_double() const {
        const double value = high + low;
        if (!std::isfinite(value))
            throw std::overflow_error("double-double: narrowing overflow");
        return value;
    }
};

} // namespace cfdx::core::numerics
