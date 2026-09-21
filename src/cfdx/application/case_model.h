#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace cfdx::application {

enum class SimulationMode { Steady, Transient };
enum class ChangeImpact { Hot, RequiresRestart, RequiresRebuild };

struct Parameter {
    using Value = std::variant<bool, std::int64_t, double, std::string>;
    Value value;
    std::string unit;
    ChangeImpact impact{ChangeImpact::Hot};
};

struct TimeControl {
    double start{0.0};
    double end{1.0};
    double step{1.0e-3};
    std::size_t max_iterations_per_step{20};
    bool adaptive{false};
};

struct CaseModel {
    std::string name{"untitled"};
    SimulationMode mode{SimulationMode::Steady};
    std::map<std::string, Parameter> parameters;
    TimeControl time;
    std::uint64_t revision{0};
    std::uint64_t mesh_revision{0};
    std::uint64_t physics_revision{0};
    std::uint64_t numerics_revision{0};

    ChangeImpact set_parameter(std::string key, Parameter parameter) {
        const auto old = parameters.find(key);
        const ChangeImpact impact = parameter.impact;
        parameters[std::move(key)] = std::move(parameter);
        ++revision;
        if (old != parameters.end()) {
            // Revisions are deliberately coarse-grained: the controller decides
            // whether a restart/rebuild is required from the parameter metadata.
        }
        return impact;
    }

    std::optional<Parameter> get_parameter(const std::string& key) const {
        const auto it = parameters.find(key);
        if (it == parameters.end()) return std::nullopt;
        return it->second;
    }

    bool is_transient() const noexcept { return mode == SimulationMode::Transient; }
};

struct Checkpoint {
    std::uint64_t case_revision{0};
    std::uint64_t mesh_revision{0};
    std::uint64_t physics_revision{0};
    std::uint64_t numerics_revision{0};
    std::size_t iteration{0};
    double time{0.0};
};

struct CaseSnapshot {
    CaseModel model;
    Checkpoint checkpoint;
};

} // namespace cfdx::application
