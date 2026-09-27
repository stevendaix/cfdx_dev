#include "cfdx/core/boundary/boundary_validator.h"
#include "common/test_harness.h"
#include <memory>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("validator_accepts_required_valid_constraint", [] {
        BoundaryValidator validator;
        validator.require_field("p");
        std::vector<BoundaryConstraint> constraints{
            {"p", Neumann{std::make_shared<ConstantValueProvider>(0.0)}}};
        EXPECT_TRUE(validator.valid(constraints));
    });

    run_case("validator_rejects_duplicate_field", [] {
        BoundaryValidator validator;
        std::vector<BoundaryConstraint> constraints{
            {"p", Neumann{std::make_shared<ConstantValueProvider>(0.0)}},
            {"p", Neumann{std::make_shared<ConstantValueProvider>(0.0)}}};
        const auto errors = validator.validate(constraints);
        EXPECT_TRUE(!errors.empty());
    });

    run_case("validator_rejects_missing_required_field", [] {
        BoundaryValidator validator;
        validator.require_field("p");
        std::vector<BoundaryConstraint> constraints;
        const auto errors = validator.validate(constraints);
        EXPECT_TRUE(!errors.empty());
    });

    return run_all();
}
