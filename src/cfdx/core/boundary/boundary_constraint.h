#pragma once
#include "cfdx/core/boundary/mathematical_condition.h"
#include <string>
#include <utility>
namespace cfdx::core {
struct BoundaryConstraint {
    std::string field;
    MathematicalCondition condition;
    BoundaryConstraint(std::string field_name, MathematicalCondition c)
        : field(std::move(field_name)), condition(std::move(c)) {
        if (field.empty()) throw std::invalid_argument("boundary constraint requires a field");
        validate_condition(condition);
    }
};
}