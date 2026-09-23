#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::runtime {

enum class MemoryKind {
    Mesh,
    Geometry,
    Fields,
    LinearAlgebra,
    Preconditioner,
    Temporary,
    MPI,
    CUDAWorkspace,
    Halo,
    Output
};

struct MemoryEstimate {
    MemoryKind kind;
    std::size_t bytes = 0;
    std::string name;
};

struct MemoryPlan {
    std::vector<MemoryEstimate> components;
    std::size_t estimated_bytes = 0;
    std::size_t runtime_reserve_bytes = 0;
    std::size_t safety_margin_bytes = 0;
    std::size_t budget_bytes = 0;
    std::size_t peak_bytes = 0;
    bool overflow = false;

    bool fits() const {
        return !overflow && peak_bytes <= budget_bytes;
    }
};

struct MemoryPlannerConfig {
    std::size_t runtime_reserve_bytes = 0;
    std::size_t safety_margin_bytes = 0;
    double safety_margin_fraction = 0.0;
};

class MemoryPlanner {
public:
    explicit MemoryPlanner(MemoryPlannerConfig config = {}) : config_(config) {
        if (!(config_.safety_margin_fraction >= 0.0 &&
              config_.safety_margin_fraction <= 1.0))
            throw std::invalid_argument("safety margin fraction must be in [0,1]");
    }

    void add(MemoryKind kind, std::size_t bytes, std::string name = {}) {
        components_.push_back({kind, bytes, std::move(name)});
    }

    void clear() { components_.clear(); }

    MemoryPlan plan(std::size_t available_memory_bytes) const {
        MemoryPlan out;
        out.components = components_;
        out.budget_bytes = available_memory_bytes;
        out.runtime_reserve_bytes = config_.runtime_reserve_bytes;

        std::size_t total = 0;
        for (const auto& c : components_) {
            if (c.bytes > std::numeric_limits<std::size_t>::max() - total) {
                out.overflow = true;
                out.estimated_bytes = std::numeric_limits<std::size_t>::max();
                break;
            }
            total += c.bytes;
        }
        out.estimated_bytes = total;

        if (out.overflow) {
            out.peak_bytes = std::numeric_limits<std::size_t>::max();
            return out;
        }

        const std::size_t proportional =
            static_cast<std::size_t>(static_cast<long double>(total) *
                                     config_.safety_margin_fraction);
        out.safety_margin_bytes = proportional > config_.safety_margin_bytes
            ? proportional : config_.safety_margin_bytes;

        if (out.estimated_bytes >
            std::numeric_limits<std::size_t>::max() - out.runtime_reserve_bytes ||
            out.estimated_bytes + out.runtime_reserve_bytes >
            std::numeric_limits<std::size_t>::max() - out.safety_margin_bytes) {
            out.overflow = true;
            out.peak_bytes = std::numeric_limits<std::size_t>::max();
            return out;
        }

        out.peak_bytes = out.estimated_bytes +
                         out.runtime_reserve_bytes +
                         out.safety_margin_bytes;
        return out;
    }

private:
    MemoryPlannerConfig config_;
    std::vector<MemoryEstimate> components_;
};

inline RuntimeWorkload make_runtime_workload(const MemoryPlan& plan) {
    RuntimeWorkload workload;
    workload.estimated_bytes = plan.estimated_bytes;
    workload.runtime_reserve_bytes = plan.runtime_reserve_bytes;
    workload.safety_margin_bytes = plan.safety_margin_bytes;
    return workload;
}

} // namespace cfdx::runtime
