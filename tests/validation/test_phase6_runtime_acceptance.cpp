#include "cfdx/runtime/execution/execution_policy.h"
#include "cfdx/runtime/execution/memory_planner.h"
#include "cfdx/runtime/gpu/gpu_execution.h"
#include "cfdx/runtime/gpu/gpu_kernels.h"
#include "cfdx/runtime/ooc/pinned_buffer_pool.h"
#include "common/test_harness.h"

#ifdef CFDX_ENABLE_GPU
#include "cfdx/runtime/gpu/cuda_double_buffer.h"
#include "cfdx/runtime/gpu/gpu_profiler.h"
#include "cfdx/runtime/ooc/cuda_pinned_buffer_pool.h"
#include <cuda_runtime.h>
#endif

#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::testing;

int main()
{
    using namespace cfdx::runtime;

    run_case("phase6_cpu_policy_contract", [] {
        RuntimeCapabilities caps;
        caps.cpu = true;
        caps.cuda = false;
        caps.cpu_threads = 8;

        CpuExecutionConfig cpu;
        cpu.threads = 4;
        cpu.deterministic = true;

        const auto decision =
            choose_execution_policy(ExecutionPolicy::CPU, caps, {}, cpu);

        EXPECT_TRUE(decision.selected == ExecutionPolicy::CPU);
        EXPECT_TRUE(decision.cpu_threads == 4);
        EXPECT_TRUE(decision.deterministic);
        EXPECT_TRUE(!decision.requires_gpu);
    });

    run_case("phase6_gpu_policy_is_hard_requirement", [] {
        RuntimeDecision decision;
        decision.selected = ExecutionPolicy::GPU;
        decision.requires_gpu = true;

        EXPECT_THROW(
            require_selected_backend(decision, false, true),
            std::runtime_error);

        require_selected_backend(decision, true, true);
    });

    run_case("phase6_memory_budget_is_fail_fast", [] {
        RuntimeCapabilities caps;
        caps.cpu = true;
        caps.cuda = true;
        caps.cuda_device_count = 1;
        caps.gpu_memory_bytes = 1024;

        RuntimeWorkload work;
        work.estimated_bytes = 900;
        work.safety_margin_bytes = 100;
        work.runtime_reserve_bytes = 100;

        EXPECT_THROW(
            choose_execution_policy(
                ExecutionPolicy::GPU, caps, work,
                CpuExecutionConfig{}, GpuExecutionConfig{}),
            std::runtime_error);
    });

    run_case("phase6_memory_planner_reports_peak", [] {
        MemoryPlanner planner;
        planner.add(MemoryKind::Mesh, 1024, "mesh");
        planner.add(MemoryKind::Geometry, 2048, "geometry");
        planner.add(MemoryKind::Fields, 4096, "fields");

        const auto plan = planner.plan(std::numeric_limits<std::size_t>::max());

        EXPECT_TRUE(!plan.overflow);
        EXPECT_TRUE(plan.estimated_bytes >= 7168);
        EXPECT_TRUE(plan.peak_bytes >= plan.estimated_bytes);
        EXPECT_TRUE(plan.fits());
    });

    run_case("phase6_cpu_only_gpu_path_fails_explicitly", [] {
#ifndef CFDX_ENABLE_GPU
        std::vector<double> gx, gy, gz;
        EXPECT_THROW(
            gpu::execute_gradient_cuda(
                {0.0, 1.0}, {1.0}, {0.0}, {0.0},
                {0U}, {1}, {1.0, 1.0}, gx, gy, gz),
            std::runtime_error);
#else
        EXPECT_TRUE(true);
#endif
    });

#ifdef CFDX_ENABLE_GPU
    run_case("phase6_real_cuda_pipeline_and_cpu_equivalence", [] {
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;

        const std::vector<double> phi{0.0, 1.0};
        const std::vector<double> sx{1.0};
        const std::vector<double> sy{0.0};
        const std::vector<double> sz{0.0};
        const std::vector<std::uint32_t> owner{0U};
        const std::vector<std::int64_t> neighbour{1};
        const std::vector<double> volume{1.0, 1.0};

        std::vector<double> cpu_x(2), cpu_y(2), cpu_z(2);
        const std::vector<std::size_t> owner_cpu{0U};
        gpu::gradient_gauss_reference(
            phi.data(), sx.data(), sy.data(), sz.data(),
            owner_cpu.data(), neighbour.data(), volume.data(),
            1, 2, cpu_x.data(), cpu_y.data(), cpu_z.data());

        std::vector<double> gpu_x, gpu_y, gpu_z;
        gpu::GpuExecutionMetrics metrics;
        gpu::execute_gradient_cuda_with_metrics(
            phi, sx, sy, sz, owner, neighbour, volume,
            gpu_x, gpu_y, gpu_z, 0, &metrics);

        EXPECT_TRUE(gpu_x.size() == cpu_x.size());
        for (std::size_t i = 0; i < gpu_x.size(); ++i) {
            EXPECT_NEAR(gpu_x[i], cpu_x[i], 1e-12);
            EXPECT_NEAR(gpu_y[i], cpu_y[i], 1e-12);
            EXPECT_NEAR(gpu_z[i], cpu_z[i], 1e-12);
        }
        EXPECT_TRUE(metrics.h2d_bytes > 0);
        EXPECT_TRUE(metrics.d2h_bytes > 0);
        EXPECT_TRUE(std::isfinite(metrics.h2d_ms));
        EXPECT_TRUE(std::isfinite(metrics.kernel_ms));
        EXPECT_TRUE(std::isfinite(metrics.d2h_ms));
    });

    run_case("phase6_invalid_gpu_device_is_rejected", [] {
        int devices = 0;
        if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) return;

        std::vector<double> gx, gy, gz;
        EXPECT_THROW(
            gpu::execute_gradient_cuda_with_metrics(
                {0.0, 1.0}, {1.0}, {0.0}, {0.0},
                {0U}, {1}, {1.0, 1.0},
                gx, gy, gz, devices, nullptr),
            std::invalid_argument);
    });

    run_case("phase6_pinned_pool_and_double_buffer", [] {
        int devices = 0;
        if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) return;

        ooc::CudaPinnedBufferPool pool(2048, 1024);
        EXPECT_TRUE(pool.is_pinned());
        EXPECT_TRUE(pool.size() == 2);

        auto* first = pool.acquire();
        auto* second = pool.acquire();
        EXPECT_TRUE(first != nullptr);
        EXPECT_TRUE(second != nullptr);
        EXPECT_TRUE(pool.acquire() == nullptr);

        pool.release(first);
        pool.release(second);

        gpu::CudaDoubleBuffer buffers(1024);
        EXPECT_TRUE(buffers.slot_count() == 2);
        EXPECT_TRUE(buffers.bytes() == 1024);
        EXPECT_TRUE(&buffers.buffer(0) != &buffers.buffer(1));
        EXPECT_TRUE(&buffers.stream(0) != &buffers.stream(1));
        buffers.synchronize(0);
        buffers.synchronize(1);
    });

    run_case("phase6_gpu_profiler_has_explicit_cuda_timing", [] {
        int devices = 0;
        if (cudaGetDeviceCount(&devices) != cudaSuccess || devices == 0) return;

        gpu::CudaStream stream;
        gpu::GpuProfiler profiler;
        profiler.start(stream.get());
        stream.synchronize();
        profiler.stop(stream.get());

        EXPECT_TRUE(std::isfinite(profiler.elapsed_ms()));
        EXPECT_TRUE(profiler.elapsed_ms() >= 0.0);
    });
#endif

    return 0;
}
