#pragma once

#include "cfdx/application/case_model.h"

#include <cstddef>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <string>

namespace cfdx::application {

class SolverAdapter;

enum class SimulationState {
    Created, Validating, Ready, Running, Paused, Stopping, Stopped, Converged, Failed
};

struct RunTarget {
    std::size_t iterations{0};
    double end_time{-1.0};
    bool until_converged{false};
};

class SimulationController {
public:
    explicit SimulationController(CaseModel model = {});

    SimulationState state() const;
    const CaseModel& model() const noexcept { return model_; }
    const Checkpoint& checkpoint() const noexcept { return checkpoint_; }

    void validate();
    void run(RunTarget target);
    void pause();
    void stop();

    // Apply a case edit while paused/stopped. Hot edits can continue directly;
    // restart/rebuild edits are recorded and exposed through requires_restart().
    ChangeImpact edit(const std::string& key, Parameter parameter);
    bool requires_restart() const noexcept { return requires_restart_; }
    void acknowledge_restart();

    void create_checkpoint();
    void restore(const CaseSnapshot& snapshot);
    const Checkpoint& latest_checkpoint() const noexcept { return checkpoint_; }

    using IterationCallback = std::function<bool(std::size_t, double)>;
    void set_iteration_callback(IterationCallback callback) { callback_ = std::move(callback); }
    void set_solver_adapter(SolverAdapter* adapter) noexcept { solver_adapter_ = adapter; }

private:
    mutable std::mutex mutex_;
    CaseModel model_;
    Checkpoint checkpoint_;
    SimulationState state_{SimulationState::Created};
    bool stop_requested_{false};
    bool pause_requested_{false};
    bool requires_restart_{false};
    IterationCallback callback_;
    SolverAdapter* solver_adapter_{nullptr};

    void set_state(SimulationState state);
};

} // namespace cfdx::application
