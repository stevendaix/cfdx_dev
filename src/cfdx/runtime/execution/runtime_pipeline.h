#pragma once

#include "cfdx/runtime/execution/execution_policy.h"
#include "cfdx/runtime/execution/memory_planner.h"
#include "cfdx/runtime/ooc/domain_decomposition.h"
#include "cfdx/runtime/ooc/ooc_executor.h"

#ifdef CFDX_HAS_MPI
#if CFDX_HAS_MPI
#include "cfdx/core/parallel/distributed_execution.h"
#endif
#endif

#include <cstddef>
#include <limits>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::runtime {

struct RuntimePipelineConfig {
    ExecutionPolicy policy = ExecutionPolicy::AUTO;
    CpuExecutionConfig cpu;
    GpuExecutionConfig gpu;
    MemoryPlannerConfig memory;
    ooc::OOCConfig ooc;
    std::size_t field_components = 1;
};

struct RuntimePipelineState {
    RuntimeDecision decision;
    MemoryPlan memory_plan;
    bool memory_plan_valid = false;
    bool ooc_prepared = false;
    std::size_t mpi_partitions = 1;
    std::size_t local_partition_cells = 0;
};

class RuntimePipeline {
public:
    explicit RuntimePipeline(RuntimePipelineConfig config = {})
        : config_(std::move(config)), ooc_(config_.ooc) {}

    RuntimePipelineState prepare(
        std::size_t n_cells,
        const std::vector<std::pair<std::uint32_t, std::uint32_t>>& adjacency,
        const RuntimeCapabilities& capabilities)
    {
        if (n_cells == 0)
            throw std::invalid_argument("RuntimePipeline: mesh has no cells");
        if (config_.field_components == 0)
            throw std::invalid_argument("RuntimePipeline: field components must be > 0");

        MemoryPlanner planner(config_.memory);
        planner.add(MemoryKind::Mesh, checked_mul(n_cells, sizeof(std::uint64_t)), "mesh");
        planner.add(
            MemoryKind::Fields,
            checked_mul(n_cells, checked_mul(config_.field_components, sizeof(double))),
            "fields");

        const std::size_t budget =
            capabilities.gpu_memory_bytes == 0
                ? std::numeric_limits<std::size_t>::max()
                : capabilities.gpu_memory_bytes;
        const auto memory_plan = planner.plan(budget);

        RuntimePipelineState state;
        state.memory_plan = memory_plan;
        state.memory_plan_valid = !memory_plan.overflow;
        state.decision = choose_execution_policy(
            config_.policy, capabilities, make_runtime_workload(memory_plan),
            config_.cpu, config_.gpu);
        require_selected_backend(
            state.decision,
            capabilities.cuda && capabilities.cuda_device_count > 0,
            capabilities.cpu);

        if (state.decision.uses_out_of_core) {
            ooc_.build_tiles(n_cells, adjacency);
            ooc_.validate_memory_budget(config_.field_components);
            state.ooc_prepared = true;
        }

        return state;
    }

    const ooc::OOCExecutor& ooc_executor() const noexcept { return ooc_; }

    static ooc::DomainPartition prepare_partition(
        std::size_t n_cells,
        std::size_t partitions,
        const std::vector<std::pair<std::uint32_t, std::uint32_t>>& adjacency)
    {
        auto partition = ooc::make_contiguous_partition(
            n_cells, partitions, adjacency);
        ooc::validate_partition(partition, n_cells, partitions);
        return partition;
    }

private:
    static std::size_t checked_mul(std::size_t a, std::size_t b)
    {
        if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
            throw std::overflow_error("RuntimePipeline size overflow");
        return a * b;
    }

    RuntimePipelineConfig config_;
    ooc::OOCExecutor ooc_;
};

} // namespace cfdx::runtime
