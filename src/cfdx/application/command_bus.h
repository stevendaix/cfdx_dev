#pragma once

#include "cfdx/application/checkpoint_store.h"
#include "cfdx/application/workbench.h"

#include <cctype>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::application {

struct CommandResult {
    bool ok{false};
    std::string message;
};

class CommandBus {
public:
    explicit CommandBus(Workbench& workbench) : workbench_(workbench) {}

    CommandResult execute(const std::string& line) {
        std::istringstream in(line);
        std::string command;
        in >> command;
        if (command.empty()) return {false, "Empty command"};

        if (command == "validate") {
            workbench_.controller().validate();
            return {true, "Case validated"};
        }
        if (command == "run") {
            std::size_t iterations = 0;
            in >> iterations;
            workbench_.controller().run({iterations, -1.0, false});
            return {true, "Run completed"};
        }
        if (command == "run-until") {
            double end_time = 0.0;
            if (!(in >> end_time)) return {false, "Usage: run-until <time>"};
            workbench_.controller().run({0, end_time, false});
            return {true, "Run-until completed"};
        }
        if (command == "stop") {
            workbench_.controller().stop();
            return {true, "Stop requested"};
        }
        if (command == "pause") {
            workbench_.controller().pause();
            return {true, "Pause requested"};
        }
        if (command == "checkpoint") {
            store_.save({workbench_.controller().model(), workbench_.controller().latest_checkpoint()});
            return {true, "Checkpoint saved"};
        }
        if (command == "restore") {
            if (!store_.contains()) return {false, "No checkpoint available"};
            restored_ = store_.load();
            workbench_.controller().restore(*restored_);
            workbench_.tree().reset(workbench_.controller().model());
            return {true, "Checkpoint restored"};
        }
        if (command == "set") {
            std::string key;
            std::string value;
            std::string impact;
            if (!(in >> key >> value)) return {false, "Usage: set <key> <value> [hot|restart|rebuild]"};
            in >> impact;
            const auto change = parse_impact(impact);
            Parameter parameter;
            parameter.impact = change;
            if (value == "true" || value == "false") parameter.value = (value == "true");
            else {
                try {
                    std::size_t used = 0;
                    const auto number = std::stod(value, &used);
                    if (used == value.size()) parameter.value = number;
                    else parameter.value = value;
                } catch (...) {
                    parameter.value = value;
                }
            }
            const auto actual = workbench_.properties().set(key, std::move(parameter));
            return {true, std::string("Parameter updated; impact=") + impact_name(actual)};
        }
        return {false, "Unknown command: " + command};
    }

    void set_checkpoint_store(CheckpointStore& store) { store_ = &store; }

    std::optional<CaseSnapshot> restored_snapshot() const { return restored_; }

private:
    static ChangeImpact parse_impact(const std::string& impact) {
        if (impact == "restart") return ChangeImpact::RequiresRestart;
        if (impact == "rebuild") return ChangeImpact::RequiresRebuild;
        return ChangeImpact::Hot;
    }

    static const char* impact_name(ChangeImpact impact) {
        switch (impact) {
        case ChangeImpact::Hot: return "hot";
        case ChangeImpact::RequiresRestart: return "restart";
        case ChangeImpact::RequiresRebuild: return "rebuild";
        }
        return "unknown";
    }

    Workbench& workbench_;
    MemoryCheckpointStore default_store_;
    CheckpointStore* store_{&default_store_};
    std::optional<CaseSnapshot> restored_;
};

} // namespace cfdx::application
