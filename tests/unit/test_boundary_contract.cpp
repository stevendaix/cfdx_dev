#include "cfdx/core/boundary/boundary.h"
#include "cfdx/core/boundary/boundary_condition.h"
#include "cfdx/core/boundary/boundary_constraint.h"
#include "cfdx/core/boundary/boundary_validator.h"
#include "cfdx/core/boundary/mathematical_condition.h"
#include "cfdx/core/boundary/value_provider.h"
#include "common/test_harness.h"
#include <memory>
using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("boundary_is_geometry_and_role", [] {
        Boundary b("inlet", BoundaryRole::INLET, {2, 4, 8});
        EXPECT_TRUE(b.name() == "inlet");
        EXPECT_TRUE(b.role() == BoundaryRole::INLET);
        EXPECT_TRUE(b.face_ids().size() == 3);
    });
    run_case("constant_provider", [] {
        ConstantValueProvider p(3.5, Dimension{0,1,-1,0,0,0,0});
        EXPECT_TRUE(p.evaluate(ValueContext{}) == 3.5);
        EXPECT_TRUE(p.dimension()[1] == 1 && p.dimension()[2] == -1);
    });
    run_case("dirichlet_constraint", [] {
        auto p = std::make_shared<ConstantValueProvider>(2.0);
        BoundaryConstraint c("U", Dirichlet{p});
        EXPECT_TRUE(c.field == "U");
    });
    run_case("robin_validation", [] {
        auto p = std::make_shared<ConstantValueProvider>(1.0);
        EXPECT_THROW(BoundaryConstraint("T", Robin{0.0, 0.0, p}), std::invalid_argument);
    });
    run_case("periodic_validation", [] {
        EXPECT_THROW(BoundaryConstraint("U", Periodic{""}), std::invalid_argument);
    });
    run_case("well_posedness_requirements", [] {
        auto p = std::make_shared<ConstantValueProvider>(0.0);
        BoundaryValidator validator;
        validator.require_field("U");
        validator.require_field("p");
        std::vector<BoundaryConstraint> constraints;
        constraints.emplace_back("U", Dirichlet{p});
        const auto errors = validator.validate(constraints);
        EXPECT_TRUE(errors.size() == 1);
        EXPECT_TRUE(errors.front().find("p") != std::string::npos);
    });
    return run_all();
}
