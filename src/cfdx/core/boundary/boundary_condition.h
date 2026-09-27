#pragma once
#include "cfdx/core/boundary/boundary.h"
#include "cfdx/core/boundary/boundary_constraint.h"
#include <vector>
namespace cfdx::core {
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;
    virtual std::vector<BoundaryConstraint> constraints(const Boundary& boundary) const = 0;
};
}