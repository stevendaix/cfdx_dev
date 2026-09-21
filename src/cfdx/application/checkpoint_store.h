#pragma once

#include "cfdx/application/case_model.h"

#include <cstdint>
#include <type_traits>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cfdx::application {

class CheckpointStore {
public:
    virtual ~CheckpointStore() = default;
    virtual void save(const CaseSnapshot& snapshot) = 0;
    virtual CaseSnapshot load() const = 0;
    virtual bool contains() const noexcept = 0;
};

class MemoryCheckpointStore final : public CheckpointStore {
public:
    void save(const CaseSnapshot& snapshot) override {
        snapshot_ = snapshot;
    }

    CaseSnapshot load() const override {
        if (!snapshot_) throw std::logic_error("No checkpoint is available");
        return *snapshot_;
    }

    bool contains() const noexcept override { return snapshot_.has_value(); }

private:
    std::optional<CaseSnapshot> snapshot_;
};

class TextCheckpointStore final : public CheckpointStore {
public:
    explicit TextCheckpointStore(std::string path) : path_(std::move(path)) {}

    void save(const CaseSnapshot& snapshot) override {
        std::ofstream out(path_);
        if (!out) throw std::runtime_error("Cannot create checkpoint: " + path_);
        out << std::setprecision(17);
        out << "CFDX_CHECKPOINT 2\n";
        out << std::quoted(snapshot.model.name) << "\n";
        out << static_cast<int>(snapshot.model.mode) << ' '
            << snapshot.model.time.start << ' ' << snapshot.model.time.end << ' '
            << snapshot.model.time.step << ' ' << snapshot.model.time.max_iterations_per_step << ' '
            << snapshot.model.time.adaptive << "\n";
        out << snapshot.model.revision << ' ' << snapshot.model.mesh_revision << ' '
            << snapshot.model.physics_revision << ' ' << snapshot.model.numerics_revision << "\n";
        out << snapshot.model.parameters.size() << "\n";
        for (const auto& [key, parameter] : snapshot.model.parameters) {
            out << std::quoted(key) << ' ' << std::quoted(parameter.unit) << ' '
                << static_cast<int>(parameter.impact) << ' ';
            std::visit([&out](const auto& value) {
                using T = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<T, bool>) out << "b " << value;
                else if constexpr (std::is_same_v<T, std::int64_t>) out << "i " << value;
                else if constexpr (std::is_same_v<T, double>) out << "d " << value;
                else out << "s " << std::quoted(value);
            }, parameter.value);
            out << "\n";
        }
        out << snapshot.checkpoint.case_revision << ' ' << snapshot.checkpoint.mesh_revision << ' '
            << snapshot.checkpoint.physics_revision << ' ' << snapshot.checkpoint.numerics_revision << ' '
            << snapshot.checkpoint.iteration << ' ' << snapshot.checkpoint.time << "\n";
    }

    CaseSnapshot load() const override {
        std::ifstream in(path_);
        if (!in) throw std::runtime_error("Cannot open checkpoint: " + path_);
        std::string magic;
        int version = 0;
        CaseSnapshot snapshot;
        if (!(in >> magic >> version) || magic != "CFDX_CHECKPOINT" || version != 2) {
            throw std::runtime_error("Unsupported CFDX checkpoint format");
        }
        if (!(in >> std::quoted(snapshot.model.name))) throw std::runtime_error("Invalid checkpoint name");
        int mode = 0;
        if (!(in >> mode >> snapshot.model.time.start >> snapshot.model.time.end
                 >> snapshot.model.time.step >> snapshot.model.time.max_iterations_per_step
                 >> snapshot.model.time.adaptive)) {
            throw std::runtime_error("Invalid checkpoint time control");
        }
        snapshot.model.mode = mode == 1 ? SimulationMode::Transient : SimulationMode::Steady;
        if (!(in >> snapshot.model.revision >> snapshot.model.mesh_revision
                 >> snapshot.model.physics_revision >> snapshot.model.numerics_revision)) {
            throw std::runtime_error("Invalid checkpoint revisions");
        }
        std::size_t parameter_count = 0;
        if (!(in >> parameter_count)) throw std::runtime_error("Invalid checkpoint parameter count");
        for (std::size_t n = 0; n < parameter_count; ++n) {
            std::string key;
            Parameter parameter;
            int impact = 0;
            char type = '\0';
            if (!(in >> std::quoted(key) >> std::quoted(parameter.unit) >> impact >> type)) {
                throw std::runtime_error("Invalid checkpoint parameter");
            }
            parameter.impact = static_cast<ChangeImpact>(impact);
            switch (type) {
            case 'b': {
                bool value = false;
                if (!(in >> value)) throw std::runtime_error("Invalid bool parameter");
                parameter.value = value;
                break;
            }
            case 'i': {
                std::int64_t value = 0;
                if (!(in >> value)) throw std::runtime_error("Invalid integer parameter");
                parameter.value = value;
                break;
            }
            case 'd': {
                double value = 0.0;
                if (!(in >> value)) throw std::runtime_error("Invalid double parameter");
                parameter.value = value;
                break;
            }
            case 's': {
                std::string value;
                if (!(in >> std::quoted(value))) throw std::runtime_error("Invalid string parameter");
                parameter.value = std::move(value);
                break;
            }
            default:
                throw std::runtime_error("Unknown checkpoint parameter type");
            }
            snapshot.model.parameters.emplace(std::move(key), std::move(parameter));
        }
        if (!(in >> snapshot.checkpoint.case_revision >> snapshot.checkpoint.mesh_revision
                 >> snapshot.checkpoint.physics_revision >> snapshot.checkpoint.numerics_revision
                 >> snapshot.checkpoint.iteration >> snapshot.checkpoint.time)) {
            throw std::runtime_error("Invalid checkpoint state");
        }
        return snapshot;
    }

    bool contains() const noexcept {
        std::ifstream in(path_);
        return static_cast<bool>(in);
    }

private:
    std::string path_;
};

} // namespace cfdx::application
