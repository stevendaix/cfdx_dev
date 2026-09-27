#pragma once
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>
namespace cfdx::core {
using Dimension = std::array<int, 7>;
inline constexpr Dimension DIMENSIONLESS{0,0,0,0,0,0,0};
struct ValueContext { double time = 0.0; double x = 0.0, y = 0.0, z = 0.0; std::size_t face_id = 0; };
class ValueProvider {
public:
    virtual ~ValueProvider() = default;
    virtual double evaluate(const ValueContext& context) const = 0;
    virtual Dimension dimension() const noexcept = 0;
};
class ConstantValueProvider final : public ValueProvider {
public:
    explicit ConstantValueProvider(double value, Dimension dimension = DIMENSIONLESS)
        : value_(value), dimension_(dimension) {}
    double evaluate(const ValueContext&) const override { return value_; }
    Dimension dimension() const noexcept override { return dimension_; }
private:
    double value_;
    Dimension dimension_;
};
}