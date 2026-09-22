#pragma once

#include "cfdx/application/case_model.h"

#include <cstddef>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>

namespace cfdx::application {

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
    // Return copies so callers cannot retain references into state mutated by run/edit.
    CaseModel model() const;
    Checkpoint checkpoint() const;
    CaseSnapshot snapshot() const;

    void validate();
    void run(RunTarget target);
    void pause();
    void stop();

    // Apply a case edit while paused/stopped. Hot edits can continue directly;
    // restart/rebuild edits are recorded and exposed through requires_restart().
    ChangeImpact edit(const std::string& key, Parameter parameter);
    bool requires_restart() const;
    void acknowledge_restart();

    void create_checkpoint();
    Checkpoint latest_checkpoint() const;

    using IterationCallback = std::function<bool(std::size_t, double)>;
    void set_iteration_callback(IterationCallback callback);

private:
    mutable std::mutex mutex_;
    CaseModel model_;
    Checkpoint checkpoint_;
    SimulationState state_{SimulationState::Created};
    bool stop_requested_{false};
    bool pause_requested_{false};
    bool requires_restart_{false};
    IterationCallback callback_;

    void set_state(SimulationState state);
};

} // namespace cfdx::application
