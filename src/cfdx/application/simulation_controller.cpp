#include "cfdx/application/simulation_controller.h"
#include "cfdx/application/solver_adapter.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace cfdx::application {

SimulationController::SimulationController(CaseModel model)
    : model_(std::move(model)) {}

SimulationState SimulationController::state() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return state_;
}

void SimulationController::set_state(SimulationState state) {
    std::lock_guard<std::mutex> lock(mutex_);
    state_ = state;
}

void SimulationController::validate() {
    set_state(SimulationState::Validating);
    if (solver_adapter_) solver_adapter_->validate(model_);
    if (model_.is_transient() && model_.time.step <= 0.0) {
        set_state(SimulationState::Failed);
        throw std::invalid_argument("Transient time step must be positive");
    }
    if (model_.is_transient() && model_.time.end < model_.time.start) {
        set_state(SimulationState::Failed);
        throw std::invalid_argument("Transient end time must not precede start time");
    }
    set_state(SimulationState::Ready);
}

void SimulationController::run(RunTarget target) {
    if (state() == SimulationState::Created) validate();
    if (state() != SimulationState::Ready &&
        state() != SimulationState::Paused &&
        state() != SimulationState::Stopped) {
        throw std::logic_error("Simulation can only run from Ready, Paused or Stopped");
    }

    {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_requested_ = false;
        pause_requested_ = false;
        state_ = SimulationState::Running;
    }

    const std::size_t max_iterations = target.iterations != 0
        ? target.iterations
        : (model_.is_transient() ? model_.time.max_iterations_per_step : 1);

    const double dt = model_.is_transient() ? model_.time.step : 0.0;
    if (solver_adapter_) solver_adapter_->begin(model_, checkpoint_);
    for (std::size_t i = 0; i < max_iterations; ++i) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (stop_requested_) {
                state_ = SimulationState::Stopped;
                if (solver_adapter_) solver_adapter_->end();
                return;
            }
            if (pause_requested_) {
                state_ = SimulationState::Paused;
                if (solver_adapter_) solver_adapter_->end();
                return;
            }
        }

        ++checkpoint_.iteration;
        if (model_.is_transient()) checkpoint_.time += dt;
        if (solver_adapter_ && !solver_adapter_->iterate(checkpoint_.iteration, checkpoint_.time)) {
            std::lock_guard<std::mutex> lock(mutex_);
            state_ = SimulationState::Failed;
            solver_adapter_->end();
            return;
        }
        if (callback_ && !callback_(checkpoint_.iteration, checkpoint_.time)) {
            std::lock_guard<std::mutex> lock(mutex_);
            state_ = SimulationState::Stopped;
            if (solver_adapter_) solver_adapter_->end();
            return;
        }

        if (target.end_time >= 0.0 && checkpoint_.time >= target.end_time) {
            set_state(SimulationState::Converged);
            create_checkpoint();
            if (solver_adapter_) solver_adapter_->end();
            return;
        }
    }

    if (target.until_converged) set_state(SimulationState::Converged);
    else set_state(SimulationState::Paused);
    create_checkpoint();
    if (solver_adapter_) solver_adapter_->end();
}

void SimulationController::pause() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (state_ != SimulationState::Running) return;
    pause_requested_ = true;
}

void SimulationController::stop() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (state_ != SimulationState::Running && state_ != SimulationState::Paused) return;
    stop_requested_ = true;
    pause_requested_ = false;
    if (state_ == SimulationState::Paused) state_ = SimulationState::Stopped;
}

ChangeImpact SimulationController::edit(const std::string& key, Parameter parameter) {
    const SimulationState current = state();
    if (current == SimulationState::Running || current == SimulationState::Validating) {
        throw std::logic_error("Edit the case only while paused or stopped");
    }
    const ChangeImpact impact = model_.set_parameter(key, std::move(parameter));
    if (impact != ChangeImpact::Hot) requires_restart_ = true;
    return impact;
}

void SimulationController::acknowledge_restart() {
    requires_restart_ = false;
    ++model_.numerics_revision;
}

void SimulationController::create_checkpoint() {
    checkpoint_.case_revision = model_.revision;
    checkpoint_.mesh_revision = model_.mesh_revision;
    checkpoint_.physics_revision = model_.physics_revision;
    checkpoint_.numerics_revision = model_.numerics_revision;
}

} // namespace cfdx::application
