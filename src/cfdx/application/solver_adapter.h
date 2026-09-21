#pragma once

#include "cfdx/application/simulation_controller.h"

#include <functional>
#include <stdexcept>
#include <string>

namespace cfdx::application {

class SolverAdapter {
public:
    virtual ~SolverAdapter() = default;
    virtual void validate(const CaseModel& model) = 0;
    virtual void begin(const CaseModel& model, const Checkpoint& checkpoint) = 0;
    virtual bool iterate(std::size_t iteration, double time) = 0;
    virtual void end() = 0;
};

class CallbackSolverAdapter final : public SolverAdapter {
public:
    using ValidateFn = std::function<void(const CaseModel&)>;
    using BeginFn = std::function<void(const CaseModel&, const Checkpoint&)>;
    using IterateFn = std::function<bool(std::size_t, double)>;
    using EndFn = std::function<void()>;

    ValidateFn validate_fn;
    BeginFn begin_fn;
    IterateFn iterate_fn;
    EndFn end_fn;

    void validate(const CaseModel& model) override {
        if (validate_fn) validate_fn(model);
    }

    void begin(const CaseModel& model, const Checkpoint& checkpoint) override {
        if (begin_fn) begin_fn(model, checkpoint);
    }

    bool iterate(std::size_t iteration, double time) override {
        return iterate_fn ? iterate_fn(iteration, time) : true;
    }

    void end() override {
        if (end_fn) end_fn();
    }
};

} // namespace cfdx::application
