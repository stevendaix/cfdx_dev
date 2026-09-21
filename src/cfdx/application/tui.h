#pragma once

#include "cfdx/application/simulation_controller.h"

#include <sstream>
#include <string>

namespace cfdx::application {

class TuiRenderer {
public:
    static std::string render(const SimulationController& controller) {
        std::ostringstream out;
        out << "CFDX | Case: " << controller.model().name << "\n";
        out << "State: " << state_name(controller.state()) << "\n";
        out << "Iteration: " << controller.checkpoint().iteration
            << " | Time: " << controller.checkpoint().time << "\n";
        out << "[RUN] [PAUSE] [STOP] [EDIT] [CHECKPOINT]\n";
        if (controller.requires_restart()) {
            out << "WARNING: restart/rebuild required before continuing.\n";
        }
        return out.str();
    }

private:
    static const char* state_name(SimulationState state) {
        switch (state) {
        case SimulationState::Created: return "CREATED";
        case SimulationState::Validating: return "VALIDATING";
        case SimulationState::Ready: return "READY";
        case SimulationState::Running: return "RUNNING";
        case SimulationState::Paused: return "PAUSED";
        case SimulationState::Stopping: return "STOPPING";
        case SimulationState::Stopped: return "STOPPED";
        case SimulationState::Converged: return "CONVERGED";
        case SimulationState::Failed: return "FAILED";
        }
        return "UNKNOWN";
    }
};

} // namespace cfdx::application
