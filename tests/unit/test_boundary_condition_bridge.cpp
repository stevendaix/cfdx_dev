#include "cfdx/core/boundary/flow_boundary_conditions.h"
#include "cfdx/physics/boundary_condition_bridge.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

int main()
{
    using namespace cfdx::core;
    using namespace cfdx::physics;

    Boundary boundary{"inlet", BoundaryRole::INLET, {0, 1, 2}};

    VelocityInlet inlet(
        std::make_shared<ConstantValueProvider>(2.0),
        std::make_shared<ConstantValueProvider>(-1.0),
        std::make_shared<ConstantValueProvider>(0.0));

    const auto lowered = lower_boundary_constraints(boundary, inlet);
    assert(lowered.velocity.at("inlet").type == VelocityBoundaryCondition::Type::FIXED_VALUE);
    assert(lowered.velocity.at("inlet").value.x == 2.0);
    assert(lowered.velocity.at("inlet").value.y == -1.0);
    assert(lowered.velocity.at("inlet").value.z == 0.0);
    assert(lowered.pressure.at("inlet").type == ScalarBoundaryType::ZERO_GRADIENT);

    NoSlip wall;
    const auto wall_lowered = lower_boundary_constraints(boundary, wall);
    assert(wall_lowered.velocity.at("inlet").value.x == 0.0);
    assert(wall_lowered.velocity.at("inlet").value.y == 0.0);
    assert(wall_lowered.velocity.at("inlet").value.z == 0.0);

    bool rejected = false;
    try {
        VelocityInlet varying(
            std::make_shared<SpatialLinearValueProvider>(1.0, 2.0, 0.0, 0.0),
            std::make_shared<ConstantValueProvider>(0.0),
            std::make_shared<ConstantValueProvider>(0.0));
        (void)lower_boundary_constraints(boundary, varying);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    return 0;
}
