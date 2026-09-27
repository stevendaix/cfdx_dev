#include "cfdx/core/boundary/value_provider.h"
#include "cfdx/core/boundary/legacy_patch_field_adapter.h"
#include "cfdx/core/boundary/flow_boundary_conditions.h"
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
    run_case("legacy_fixed_value_adapter", [] {
        PatchField patch("inlet", 2, "fixedValue");
        patch.fill(5.0);
        const BoundaryConstraint c = LegacyPatchFieldAdapter::to_constraint(patch, "T");
        EXPECT_TRUE(std::holds_alternative<Dirichlet>(c.condition));
        EXPECT_TRUE(std::get<Dirichlet>(c.condition).value->evaluate(ValueContext{}) == 5.0);
    });
    run_case("legacy_zero_gradient_adapter", [] {
        PatchField patch("outlet", 2, "zeroGradient");
        const BoundaryConstraint c = LegacyPatchFieldAdapter::to_constraint(patch, "p");
        EXPECT_TRUE(std::holds_alternative<Neumann>(c.condition));
        EXPECT_TRUE(std::get<Neumann>(c.condition).gradient->evaluate(ValueContext{}) == 0.0);
    });
    run_case("legacy_varying_fixed_value_is_rejected", [] {
        PatchField patch("inlet", 2, "fixedValue");
        patch(0) = 1.0; patch(1) = 2.0;
        EXPECT_THROW(LegacyPatchFieldAdapter::to_constraint(patch, "T"), std::invalid_argument);
    });
    run_case("velocity_inlet_maps_to_field_constraints", [] {
        auto ux = std::make_shared<ConstantValueProvider>(1.0);
        auto uy = std::make_shared<ConstantValueProvider>(2.0);
        auto uz = std::make_shared<ConstantValueProvider>(0.0);
        VelocityInlet inlet(ux, uy, uz);
        Boundary b("inlet", BoundaryRole::INLET, {0});
        const auto constraints = inlet.constraints(b);
        EXPECT_TRUE(constraints.size() == 4);
        EXPECT_TRUE(constraints[0].field == "U.x");
        EXPECT_TRUE(std::holds_alternative<Dirichlet>(constraints[0].condition));
        EXPECT_TRUE(constraints[3].field == "p");
        EXPECT_TRUE(std::holds_alternative<Neumann>(constraints[3].condition));
    });
    run_case("pressure_outlet_maps_pressure_and_velocity", [] {
        auto p = std::make_shared<ConstantValueProvider>(0.0);
        PressureOutlet outlet(p);
        Boundary b("outlet", BoundaryRole::OUTLET, {1});
        const auto constraints = outlet.constraints(b);
        EXPECT_TRUE(constraints.size() == 4);
        EXPECT_TRUE(constraints[0].field == "p");
        EXPECT_TRUE(std::holds_alternative<Dirichlet>(constraints[0].condition));
        EXPECT_TRUE(constraints[1].field == "U.x");
        EXPECT_TRUE(std::holds_alternative<Neumann>(constraints[1].condition));
    });
    run_case("no_slip_is_compositional", [] {
        NoSlip wall;
        Boundary b("wall", BoundaryRole::WALL, {2});
        const auto constraints = wall.constraints(b);
        EXPECT_TRUE(constraints.size() == 3);
    });
    return run_all();
}
