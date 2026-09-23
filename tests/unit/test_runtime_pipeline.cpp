#include "cfdx/runtime/execution/runtime_pipeline.h"
#include "common/test_harness.h"

#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

using namespace cfdx::testing;

int main()
{
    using namespace cfdx::runtime;

    run_case("runtime_pipeline_composes_phase5_partition_and_phase6_cpu_policy", [] {
        RuntimePipelineConfig config;
        config.policy = ExecutionPolicy::CPU;
        config.cpu.mode = CpuExecutionMode::SERIAL;
        config.cpu.deterministic = true;

        RuntimePipeline pipeline(config);
        const std::vector<std::pair<std::uint32_t, std::uint32_t>> adjacency{
            {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}};

        const auto partition = RuntimePipeline::prepare_partition(6, 2, adjacency);
        EXPECT_TRUE(partition.cells_by_partition.size() == 2);
        EXPECT_TRUE(partition.cells_by_partition[0].size() == 3);
        EXPECT_TRUE(partition.cells_by_partition[1].size() == 3);
        EXPECT_TRUE(partition.cut_edges.size() == 1);

        RuntimeCapabilities caps;
        caps.cpu = true;
        caps.cpu_threads = 8;

        const auto state = pipeline.prepare(6, adjacency, caps);
        EXPECT_TRUE(state.decision.selected == ExecutionPolicy::CPU);
        EXPECT_TRUE(state.decision.deterministic);
        EXPECT_TRUE(!state.ooc_prepared);
        EXPECT_TRUE(state.memory_plan_valid);
    });

    run_case("runtime_pipeline_gpu_ooc_prepares_phase7_tiles", [] {
        RuntimePipelineConfig config;
        config.policy = ExecutionPolicy::GPU_OUT_OF_CORE;
        config.field_components = 2;
        config.ooc.tile_cells = 2;
        config.ooc.device_memory_bytes = 64;
        config.ooc.safety_margin_bytes = 0;

        RuntimePipeline pipeline(config);
        const std::vector<std::pair<std::uint32_t, std::uint32_t>> adjacency{
            {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}};

        RuntimeCapabilities caps;
        caps.cpu = true;
        caps.cuda = true;
        caps.cuda_device_count = 1;
        caps.gpu_memory_bytes = 64;

        const auto state = pipeline.prepare(6, adjacency, caps);
        EXPECT_TRUE(state.decision.selected == ExecutionPolicy::GPU_OUT_OF_CORE);
        EXPECT_TRUE(state.decision.requires_gpu);
        EXPECT_TRUE(state.ooc_prepared);
        EXPECT_TRUE(pipeline.ooc_executor().tiles().tiles().size() == 3);
    });

    run_case("runtime_pipeline_never_hides_gpu_failure", [] {
        RuntimePipelineConfig config;
        config.policy = ExecutionPolicy::GPU;

        RuntimePipeline pipeline(config);
        RuntimeCapabilities caps;
        caps.cpu = true;
        caps.cuda = false;
        caps.cuda_device_count = 0;

        EXPECT_THROW(
            pipeline.prepare(4, {{0, 1}, {1, 2}, {2, 3}}, caps),
            std::runtime_error);
    });

    return run_all();
}
