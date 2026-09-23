#pragma once

#include "cfdx/core/memory/memory_planner.h"
#include "cfdx/runtime/execution/execution_policy.h"

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

        cfdx::core::memory::MemoryPlanner core_planner;
        std::vector<cfdx::core::memory::BufferDescriptor> buffers;
        buffers.reserve(components_.size());

        for (std::size_t i = 0; i < components_.size(); ++i) {
            const auto& component = components_[i];
            cfdx::core::memory::BufferType type = to_buffer_type(component.kind);
            buffers.push_back({
                {static_cast<std::uint64_t>(i + 1), component.name},
                type,
                component.bytes,
                cfdx::core::memory::MemoryLocation::DEVICE,
                0,
                0,
                {}
            });
        }

        const auto core_plan = core_planner.plan(
            buffers, 0, 1, 1, 1, 1);

        out.estimated_bytes = core_plan.budget.peak_vram;
        if (!core_plan.feasible && !components_.empty()) {
            out.overflow = true;
            out.peak_bytes = std::numeric_limits<std::size_t>::max();
            return out;
        }

        const std::size_t proportional =
            static_cast<std::size_t>(
                static_cast<long double>(out.estimated_bytes) *
                config_.safety_margin_fraction);
        out.safety_margin_bytes =
            proportional > config_.safety_margin_bytes
                ? proportional : config_.safety_margin_bytes;

        if (out.estimated_bytes >
                std::numeric_limits<std::size_t>::max() -
                    out.runtime_reserve_bytes ||
            out.estimated_bytes + out.runtime_reserve_bytes >
                std::numeric_limits<std::size_t>::max() -
                    out.safety_margin_bytes) {
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
    static cfdx::core::memory::BufferType to_buffer_type(MemoryKind kind) {
        using BT = cfdx::core::memory::BufferType;
        switch (kind) {
            case MemoryKind::Mesh:
            case MemoryKind::Geometry:
            case MemoryKind::Fields:
            case MemoryKind::LinearAlgebra:
            case MemoryKind::Preconditioner:
            case MemoryKind::Temporary:
            case MemoryKind::MPI:
            case MemoryKind::CUDAWorkspace:
            case MemoryKind::Halo:
            case MemoryKind::Output:
                return BT::TEMPORARY;
        }
        return BT::TEMPORARY;
    }

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
