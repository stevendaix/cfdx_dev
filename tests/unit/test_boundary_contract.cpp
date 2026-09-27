#include "cfdx/core/boundary/value_provider.h"
#include "common/test_harness.h"
#include <cmath>
using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("constant_provider", [] {
        ConstantValueProvider p(3.5);
        EXPECT_TRUE(p.evaluate(ValueContext{}) == 3.5);
    });
    run_case("linear_time_provider", [] {
        LinearTimeValueProvider p(2.0, 0.5);
        ValueContext c; c.time = 4.0;
        EXPECT_TRUE(std::abs(p.evaluate(c) - 4.0) < 1e-14);
    });
    run_case("spatial_linear_provider", [] {
        SpatialLinearValueProvider p(1.0, 2.0, 3.0, 4.0);
        ValueContext c; c.x = 1.0; c.y = 2.0; c.z = 3.0;
        EXPECT_TRUE(std::abs(p.evaluate(c) - 21.0) < 1e-14);
    });
    run_case("table_provider_interpolates", [] {
        TableValueProvider p({0.0, 2.0, 4.0}, {0.0, 10.0, 20.0});
        ValueContext c; c.time = 1.0;
        EXPECT_TRUE(std::abs(p.evaluate(c) - 5.0) < 1e-14);
    });
    run_case("table_provider_rejects_bad_grid", [] {
        EXPECT_THROW(TableValueProvider({0.0, 0.0}, {1.0, 2.0}), std::invalid_argument);
    });
    run_case("device_lowering_is_value_based", [] {
        ValueProviderSpec spec = SpatialLinearProviderSpec{1.0, 2.0, 3.0, 4.0, DIMENSIONLESS};
        const DeviceValueProvider d = lower_value_provider(spec);
        ValueContext c; c.x = 1.0; c.y = 2.0; c.z = 3.0;
        EXPECT_TRUE(std::abs(d.evaluate(c) - 21.0) < 1e-14);
    });
    return run_all();
}
