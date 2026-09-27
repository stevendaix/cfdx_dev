#pragma once
#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace cfdx::core {

using Dimension = std::array<int, 7>;
inline constexpr Dimension DIMENSIONLESS{0,0,0,0,0,0,0};

struct ValueContext {
    double time = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    std::size_t face_id = 0;
};

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

class LinearTimeValueProvider final : public ValueProvider {
public:
    LinearTimeValueProvider(double value0, double slope, Dimension dimension = DIMENSIONLESS)
        : value0_(value0), slope_(slope), dimension_(dimension) {}
    double evaluate(const ValueContext& context) const override {
        return value0_ + slope_ * context.time;
    }
    Dimension dimension() const noexcept override { return dimension_; }
private:
    double value0_;
    double slope_;
    Dimension dimension_;
};

class SpatialLinearValueProvider final : public ValueProvider {
public:
    SpatialLinearValueProvider(double origin, double gx, double gy, double gz,
                               Dimension dimension = DIMENSIONLESS)
        : origin_(origin), gx_(gx), gy_(gy), gz_(gz), dimension_(dimension) {}
    double evaluate(const ValueContext& context) const override {
        return origin_ + gx_ * context.x + gy_ * context.y + gz_ * context.z;
    }
    Dimension dimension() const noexcept override { return dimension_; }
private:
    double origin_;
    double gx_;
    double gy_;
    double gz_;
    Dimension dimension_;
};

class TableValueProvider final : public ValueProvider {
public:
    TableValueProvider(std::vector<double> abscissa, std::vector<double> values,
                       Dimension dimension = DIMENSIONLESS)
        : x_(std::move(abscissa)), y_(std::move(values)), dimension_(dimension) {
        if (x_.empty() || x_.size() != y_.size())
            throw std::invalid_argument("TableValueProvider requires equally-sized non-empty data");
        for (std::size_t i = 1; i < x_.size(); ++i)
            if (!(x_[i] > x_[i - 1]))
                throw std::invalid_argument("TableValueProvider abscissa must be strictly increasing");
    }

    double evaluate(const ValueContext& context) const override {
        if (context.time <= x_.front()) return y_.front();
        if (context.time >= x_.back()) return y_.back();
        for (std::size_t i = 1; i < x_.size(); ++i) {
            if (context.time <= x_[i]) {
                const double t = (context.time - x_[i - 1]) / (x_[i] - x_[i - 1]);
                return y_[i - 1] + t * (y_[i] - y_[i - 1]);
            }
        }
        return y_.back();
    }

    Dimension dimension() const noexcept override { return dimension_; }
private:
    std::vector<double> x_;
    std::vector<double> y_;
    Dimension dimension_;
};

// Configuration-side provider description. It is intentionally value-based so it
// can be lowered to a backend-specific representation before numerical kernels.
struct ConstantProviderSpec {
    double value = 0.0;
    Dimension dimension = DIMENSIONLESS;
};
struct LinearTimeProviderSpec {
    double value0 = 0.0;
    double slope = 0.0;
    Dimension dimension = DIMENSIONLESS;
};
struct SpatialLinearProviderSpec {
    double origin = 0.0;
    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;
    Dimension dimension = DIMENSIONLESS;
};

using ValueProviderSpec = std::variant<ConstantProviderSpec, LinearTimeProviderSpec,
                                       SpatialLinearProviderSpec>;

struct DeviceValueProvider {
    enum class Kind : unsigned char { CONSTANT, LINEAR_TIME, SPATIAL_LINEAR };
    Kind kind = Kind::CONSTANT;
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;
    double d = 0.0;
    Dimension dimension = DIMENSIONLESS;

    double evaluate(const ValueContext& context) const noexcept {
        switch (kind) {
            case Kind::CONSTANT: return a;
            case Kind::LINEAR_TIME: return a + b * context.time;
            case Kind::SPATIAL_LINEAR: return a + b * context.x + c * context.y + d * context.z;
        }
        return 0.0;
    }
};

static_assert(std::is_trivially_copyable_v<DeviceValueProvider>,
              "DeviceValueProvider must remain a compact device-transfer representation");

inline DeviceValueProvider lower_value_provider(const ValueProviderSpec& spec) {
    return std::visit([](const auto& provider) {
        DeviceValueProvider device;
        using T = std::decay_t<decltype(provider)>;
        if constexpr (std::is_same_v<T, ConstantProviderSpec>) {
            device.kind = DeviceValueProvider::Kind::CONSTANT;
            device.a = provider.value;
            device.dimension = provider.dimension;
        } else if constexpr (std::is_same_v<T, LinearTimeProviderSpec>) {
            device.kind = DeviceValueProvider::Kind::LINEAR_TIME;
            device.a = provider.value0;
            device.b = provider.slope;
            device.dimension = provider.dimension;
        } else {
            device.kind = DeviceValueProvider::Kind::SPATIAL_LINEAR;
            device.a = provider.origin;
            device.b = provider.gx;
            device.c = provider.gy;
            device.d = provider.gz;
            device.dimension = provider.dimension;
        }
        return device;
    }, spec);
}

} // namespace cfdx::core
