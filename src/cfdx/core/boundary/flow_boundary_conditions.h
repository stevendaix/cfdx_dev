#pragma once
#include "cfdx/core/boundary/boundary_condition.h"
#include <memory>
#include <vector>
namespace cfdx::core {

class VelocityInlet final : public BoundaryCondition {
public:
    VelocityInlet(std::shared_ptr<const ValueProvider> ux,
                  std::shared_ptr<const ValueProvider> uy,
                  std::shared_ptr<const ValueProvider> uz)
        : ux_(std::move(ux)), uy_(std::move(uy)), uz_(std::move(uz)) {
        if (!ux_ || !uy_ || !uz_) throw std::invalid_argument("VelocityInlet requires three velocity providers");
    }

    std::vector<BoundaryConstraint> constraints(const Boundary&) const override {
        auto zero = std::make_shared<ConstantValueProvider>(0.0);
        return {
            {"U.x", Dirichlet{ux_}},
            {"U.y", Dirichlet{uy_}},
            {"U.z", Dirichlet{uz_}},
            {"p", Neumann{zero}}
        };
    }
private:
    std::shared_ptr<const ValueProvider> ux_, uy_, uz_;
};

class PressureOutlet final : public BoundaryCondition {
public:
    explicit PressureOutlet(std::shared_ptr<const ValueProvider> pressure)
        : pressure_(std::move(pressure)) {
        if (!pressure_) throw std::invalid_argument("PressureOutlet requires a pressure provider");
    }

    std::vector<BoundaryConstraint> constraints(const Boundary&) const override {
        auto zero = std::make_shared<ConstantValueProvider>(0.0);
        return {
            {"p", Dirichlet{pressure_}},
            {"U.x", FluxDependent{zero, zero}},
            {"U.y", FluxDependent{zero, zero}},
            {"U.z", FluxDependent{zero, zero}}
        };
    }
private:
    std::shared_ptr<const ValueProvider> pressure_;
};

class NoSlip final : public BoundaryCondition {
public:
    std::vector<BoundaryConstraint> constraints(const Boundary&) const override {
        auto zero = std::make_shared<ConstantValueProvider>(0.0);
        return {
            {"U.x", Dirichlet{zero}},
            {"U.y", Dirichlet{zero}},
            {"U.z", Dirichlet{zero}}
        };
    }
};

} // namespace cfdx::core
