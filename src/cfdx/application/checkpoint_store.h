#pragma once

#include "cfdx/application/case_model.h"

#include <cstdint>
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
        out << "CFDX_CHECKPOINT 1\n";
        out << snapshot.model.name << "\n";
        out << static_cast<int>(snapshot.model.mode) << "\n";
        out << snapshot.model.revision << ' ' << snapshot.model.mesh_revision << ' '
            << snapshot.model.physics_revision << ' ' << snapshot.model.numerics_revision << "\n";
        out << snapshot.checkpoint.case_revision << ' ' << snapshot.checkpoint.mesh_revision << ' '
            << snapshot.checkpoint.physics_revision << ' ' << snapshot.checkpoint.numerics_revision << ' '
            << snapshot.checkpoint.iteration << ' ' << snapshot.checkpoint.time << "\n";
    }

    CaseSnapshot load() const override {
        std::ifstream in(path_);
        if (!in) throw std::runtime_error("Cannot open checkpoint: " + path_);
        std::string magic;
        int mode = 0;
        CaseSnapshot snapshot;
        if (!(in >> magic) || magic != "CFDX_CHECKPOINT") throw std::runtime_error("Invalid CFDX checkpoint");
        int version = 0;
        if (!(in >> version) || version != 1) throw std::runtime_error("Unsupported CFDX checkpoint version");
        in.ignore();
        std::getline(in, snapshot.model.name);
        if (!(in >> mode)) throw std::runtime_error("Invalid checkpoint mode");
        snapshot.model.mode = mode == 1 ? SimulationMode::Transient : SimulationMode::Steady;
        if (!(in >> snapshot.model.revision >> snapshot.model.mesh_revision
                 >> snapshot.model.physics_revision >> snapshot.model.numerics_revision)) {
            throw std::runtime_error("Invalid checkpoint revisions");
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
