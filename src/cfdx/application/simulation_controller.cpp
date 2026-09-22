#include "cfdx/application/simulation_controller.h"

#include <utility>

namespace cfdx::application {

SimulationController::SimulationController(CaseModel model)
    : model_(std::move(model)) {
    checkpoint_.time = model_.is_transient() ? model_.time.start : 0.0;
    create_checkpoint();
}

SimulationState SimulationController::state() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return state_;
}

CaseModel SimulationController::model() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return model_;
}

Checkpoint SimulationController::checkpoint() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return checkpoint_;
}

CaseSnapshot SimulationController::snapshot() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return {model_, checkpoint_};
}

void SimulationController::set_state(SimulationState state) {
    std::lock_guard<std::mutex> lock(mutex_);
    state_ = state;
}

void SimulationController::validate() {
    std::lock_guard<std::mutex> lock(mutex_);
    state_ = SimulationState::Validating;
    if (model_.is_transient() && model_.time.step <= 0.0) {
        state_ = SimulationState::Failed;
        throw std::invalid_argument("Transient time step must be positive");
    }
    if (model_.is_transient() && model_.time.end < model_.time.start) {
        state_ = SimulationState::Failed;
        throw std::invalid_argument("Transient end time must not precede start time");
    }
    state_ = SimulationState::Ready;
}

void SimulationController::run(RunTarget target) {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (state_ == SimulationState::Created) {
            state_ = SimulationState::Validating;
            if (model_.is_transient() && model_.time.step <= 0.0) {
                state_ = SimulationState::Failed;
                throw std::invalid_argument("Transient time step must be positive");
            }
            if (model_.is_transient() && model_.time.end < model_.time.start) {
                state_ = SimulationState::Failed;
                throw std::invalid_argument("Transient end time must not precede start time");
            }
            state_ = SimulationState::Ready;
        }
        if (state_ != SimulationState::Ready &&
            state_ != SimulationState::Paused &&
            state_ != SimulationState::Stopped) {
            throw std::logic_error("Simulation can only run from Ready, Paused or Stopped");
        }

        stop_requested_ = false;
        pause_requested_ = false;
        state_ = SimulationState::Running;
    }

    std::size_t max_iterations = 0;
    double dt = 0.0;
    bool transient = false;
    IterationCallback callback;
    {
        std::lock_guard<std::mutex> lock(mutex_);
        transient = model_.is_transient();
        max_iterations = target.iterations != 0
            ? target.iterations : (transient ? model_.time.max_iterations_per_step : 1);
        dt = transient ? model_.time.step : 0.0;
        callback = callback_;
    }

    for (std::size_t i = 0; i < max_iterations; ++i) {
        std::size_t iteration = 0;
        double time = 0.0;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (stop_requested_) {
                state_ = SimulationState::Stopped;
                return;
            }
            if (pause_requested_) {
                state_ = SimulationState::Paused;
                return;
            }

            ++checkpoint_.iteration;
            if (transient) checkpoint_.time += dt;
            iteration = checkpoint_.iteration;
            time = checkpoint_.time;
        }

        if (callback && !callback(iteration, time)) {
            std::lock_guard<std::mutex> lock(mutex_);
            state_ = SimulationState::Stopped;
            return;
        }

        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (target.end_time >= 0.0 && checkpoint_.time >= target.end_time) {
                state_ = SimulationState::Converged;
                checkpoint_.case_revision = model_.revision;
                checkpoint_.mesh_revision = model_.mesh_revision;
                checkpoint_.physics_revision = model_.physics_revision;
                checkpoint_.numerics_revision = model_.numerics_revision;
                return;
            }
        }
    }

    std::lock_guard<std::mutex> lock(mutex_);
    state_ = target.until_converged ? SimulationState::Converged
                                    : SimulationState::Paused;
    checkpoint_.case_revision = model_.revision;
    checkpoint_.mesh_revision = model_.mesh_revision;
    checkpoint_.physics_revision = model_.physics_revision;
    checkpoint_.numerics_revision = model_.numerics_revision;
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
    std::lock_guard<std::mutex> lock(mutex_);
    if (state_ == SimulationState::Running || state_ == SimulationState::Validating) {
        throw std::logic_error("Edit the case only while paused or stopped");
    }
    const ChangeImpact impact = parameter.impact;
    model_.set_parameter(key, std::move(parameter));
    // A case edit changes the model represented by the checkpoint. Keep the
    // snapshot revision metadata coherent until the next solver checkpoint.
    checkpoint_.case_revision = model_.revision;
    if (impact != ChangeImpact::Hot) requires_restart_ = true;
    return impact;
}

bool SimulationController::requires_restart() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return requires_restart_;
}

void SimulationController::acknowledge_restart() {
    std::lock_guard<std::mutex> lock(mutex_);
    requires_restart_ = false;
    ++model_.numerics_revision;
    checkpoint_.numerics_revision = model_.numerics_revision;
}

void SimulationController::create_checkpoint() {
    std::lock_guard<std::mutex> lock(mutex_);
    checkpoint_.case_revision = model_.revision;
    checkpoint_.mesh_revision = model_.mesh_revision;
    checkpoint_.physics_revision = model_.physics_revision;
    checkpoint_.numerics_revision = model_.numerics_revision;
}

Checkpoint SimulationController::latest_checkpoint() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return checkpoint_;
}

void SimulationController::set_iteration_callback(IterationCallback callback) {
    std::lock_guard<std::mutex> lock(mutex_);
    callback_ = std::move(callback);
}

} // namespace cfdx::application
