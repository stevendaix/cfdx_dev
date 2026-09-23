#include "cfdx/runtime/execution/execution_policy.h"
#include "cfdx/runtime/gpu/device_buffer.h"
#include "cfdx/runtime/gpu/gpu_kernels.h"
#include "cfdx/runtime/gpu/gpu_execution.h"
#include "cfdx/runtime/ooc/tile_manager.h"
#include "cfdx/runtime/ooc/working_set.h"
#include "cfdx/runtime/ooc/pinned_buffer_pool.h"
#include "cfdx/runtime/ooc/async_transfer.h"
#include "common/test_harness.h"
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::testing;

int main() {
    using namespace cfdx::runtime;

    run_case("auto_policy_selects_cpu_without_cuda", [] {
        RuntimeCapabilities caps;
        caps.cuda = false;
        RuntimeWorkload work;
        work.estimated_bytes = 100;
        auto d = choose_execution_policy(ExecutionPolicy::AUTO, caps, work);
        EXPECT_TRUE(d.selected == ExecutionPolicy::CPU);
        EXPECT_TRUE(d.cpu_threads >= 1);
    });

    run_case("explicit_cpu_policy_is_self_contained", [] {
        RuntimeCapabilities caps;
        caps.cuda = false;
        caps.cpu_threads = 8;
        CpuExecutionConfig config;
        config.threads = 4;
        config.deterministic = true;
        auto d = choose_execution_policy(ExecutionPolicy::CPU, caps, {}, config);
        EXPECT_TRUE(d.selected == ExecutionPolicy::CPU);
        EXPECT_TRUE(d.cpu_threads == 4);
        EXPECT_TRUE(d.deterministic);
        EXPECT_TRUE(!d.requires_gpu);
    });

    run_case("serial_cpu_policy_forces_one_thread", [] {
        RuntimeCapabilities caps;
        caps.cpu_threads = 32;
        CpuExecutionConfig config;
        config.mode = CpuExecutionMode::SERIAL;
        config.deterministic = true;
        auto d = choose_execution_policy(ExecutionPolicy::CPU, caps, {}, config);
        EXPECT_TRUE(d.cpu_threads == 1);
        EXPECT_TRUE(d.deterministic);
    });

    run_case("cpu_policy_rejects_impossible_thread_request", [] {
        RuntimeCapabilities caps;
        caps.cpu_threads = 4;
        CpuExecutionConfig config;
        config.threads = 8;
        EXPECT_THROW(
            choose_execution_policy(ExecutionPolicy::CPU, caps, {}, config),
            std::invalid_argument);
    });

    run_case("cpu_policy_rejects_missing_cpu_backend", [] {
        RuntimeCapabilities caps;
        caps.cpu = false;
        EXPECT_THROW(
            choose_execution_policy(ExecutionPolicy::CPU, caps, {}),
            std::runtime_error);
    });

    run_case("auto_policy_rejects_when_no_backend_exists", [] {
        RuntimeCapabilities caps;
        caps.cpu = false;
        caps.cuda = false;
        EXPECT_THROW(
            choose_execution_policy(ExecutionPolicy::AUTO, caps, {}),
            std::runtime_error);
    });

    run_case("gpu_policy_requires_cuda", [] {
        RuntimeCapabilities caps;
        EXPECT_THROW(
            choose_execution_policy(ExecutionPolicy::GPU, caps, {}),
            std::runtime_error);
    });

    run_case("gpu_policy_requires_real_device", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 0;
        EXPECT_THROW(
            choose_execution_policy(ExecutionPolicy::GPU, caps, {}),
            std::runtime_error);
    });

    run_case("gpu_policy_selects_requested_device", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 2;
        caps.gpu_memory_bytes = 1024;
        GpuExecutionConfig gpu;
        gpu.device_id = 1;
        RuntimeWorkload work;
        work.estimated_bytes = 128;
        auto d = choose_execution_policy(ExecutionPolicy::GPU, caps, work, {}, gpu);
        EXPECT_TRUE(d.selected == ExecutionPolicy::GPU);
        EXPECT_TRUE(d.requires_gpu);
        EXPECT_TRUE(d.gpu_device_id == 1);
        EXPECT_TRUE(!d.uses_out_of_core);
    });

    run_case("gpu_policy_rejects_invalid_device", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 2;
        EXPECT_THROW(
            choose_execution_policy(
                ExecutionPolicy::GPU, caps, {}, {}, GpuExecutionConfig{2, true}),
            std::invalid_argument);
    });

    run_case("explicit_gpu_never_falls_back_when_oversubscribed", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 1;
        caps.gpu_memory_bytes = 1024;
        RuntimeWorkload work;
        work.estimated_bytes = 2048;
        EXPECT_THROW(
            choose_execution_policy(ExecutionPolicy::GPU, caps, work),
            std::runtime_error);
    });

    run_case("explicit_gpu_ooc_is_distinct_policy", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 1;
        caps.gpu_memory_bytes = 1024;
        RuntimeWorkload work;
        work.estimated_bytes = 2048;
        auto d = choose_execution_policy(ExecutionPolicy::GPU_OUT_OF_CORE, caps, work);
        EXPECT_TRUE(d.selected == ExecutionPolicy::GPU_OUT_OF_CORE);
        EXPECT_TRUE(d.requires_gpu && d.uses_out_of_core);
        EXPECT_TRUE(d.gpu_device_id == 0);
    });

    run_case("auto_selects_ooc_when_full_gpu_does_not_fit", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 1;
        caps.gpu_memory_bytes = 1024;
        RuntimeWorkload work;
        work.estimated_bytes = 2048;
        auto d = choose_execution_policy(ExecutionPolicy::AUTO, caps, work);
        EXPECT_TRUE(d.selected == ExecutionPolicy::GPU_OUT_OF_CORE);
        EXPECT_TRUE(d.requires_gpu && d.uses_out_of_core);
    });

    run_case("gpu_ooc_can_be_disabled_explicitly", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 1;
        caps.gpu_memory_bytes = 1024;
        RuntimeWorkload work;
        work.estimated_bytes = 2048;
        EXPECT_THROW(
            choose_execution_policy(
                ExecutionPolicy::GPU_OUT_OF_CORE, caps, work, {}, GpuExecutionConfig{0, false}),
            std::runtime_error);
    });

    run_case("gpu_budget_reserve_overflow_is_safe", [] {
        RuntimeCapabilities caps;
        caps.cuda = true;
        caps.cuda_device_count = 1;
        caps.gpu_memory_bytes = std::numeric_limits<std::size_t>::max();
        RuntimeWorkload work;
        work.estimated_bytes = 1;
        work.safety_margin_bytes = std::numeric_limits<std::size_t>::max();
        work.runtime_reserve_bytes = 1;
        EXPECT_THROW(
            choose_execution_policy(ExecutionPolicy::GPU, caps, work),
            std::runtime_error);
    });

    run_case("cpu_gpu_gradient_numerical_equivalence", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;

        const std::vector<double> phi{0.0, 1.0};
        const std::vector<double> sx{1.0};
        const std::vector<double> sy{0.0};
        const std::vector<double> sz{0.0};
        const std::vector<std::size_t> owner_cpu{0};
        const std::vector<std::uint32_t> owner_gpu{0};
        const std::vector<std::int64_t> neighbour{1};
        const std::vector<double> volume{1.0, 1.0};

        std::vector<double> cpu_x(2), cpu_y(2), cpu_z(2);
        gpu::gradient_gauss_reference(
            phi.data(), sx.data(), sy.data(), sz.data(),
            owner_cpu.data(), neighbour.data(), volume.data(),
            1, 2, cpu_x.data(), cpu_y.data(), cpu_z.data());

        std::vector<double> gpu_x, gpu_y, gpu_z;
        gpu::execute_gradient_cuda(
            phi, sx, sy, sz, owner_gpu, neighbour, volume,
            gpu_x, gpu_y, gpu_z);

        EXPECT_TRUE(gpu_x.size() == cpu_x.size());
        EXPECT_TRUE(gpu_y.size() == cpu_y.size());
        EXPECT_TRUE(gpu_z.size() == cpu_z.size());
        for (std::size_t i = 0; i < cpu_x.size(); ++i) {
            EXPECT_NEAR(gpu_x[i], cpu_x[i], 1e-12);
            EXPECT_NEAR(gpu_y[i], cpu_y[i], 1e-12);
            EXPECT_NEAR(gpu_z[i], cpu_z[i], 1e-12);
        }
#else
        EXPECT_TRUE(true);
#endif
    });

    run_case("cpu_gpu_gradient_equivalence_is_repeatable", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;

        const std::vector<double> phi{0.25, 1.5};
        const std::vector<double> sx{0.75};
        const std::vector<double> sy{0.125};
        const std::vector<double> sz{-0.25};
        const std::vector<std::size_t> owner_cpu{0};
        const std::vector<std::uint32_t> owner_gpu{0};
        const std::vector<std::int64_t> neighbour{1};
        const std::vector<double> volume{0.5, 1.25};

        std::vector<double> cpu_x(2), cpu_y(2), cpu_z(2);
        gpu::gradient_gauss_reference(
            phi.data(), sx.data(), sy.data(), sz.data(),
            owner_cpu.data(), neighbour.data(), volume.data(),
            1, 2, cpu_x.data(), cpu_y.data(), cpu_z.data());

        std::vector<double> gx1, gy1, gz1, gx2, gy2, gz2;
        gpu::execute_gradient_cuda(
            phi, sx, sy, sz, owner_gpu, neighbour, volume,
            gx1, gy1, gz1);
        gpu::execute_gradient_cuda(
            phi, sx, sy, sz, owner_gpu, neighbour, volume,
            gx2, gy2, gz2);

        for (std::size_t i = 0; i < cpu_x.size(); ++i) {
            EXPECT_NEAR(gx1[i], cpu_x[i], 1e-12);
            EXPECT_NEAR(gy1[i], cpu_y[i], 1e-12);
            EXPECT_NEAR(gz1[i], cpu_z[i], 1e-12);
            EXPECT_NEAR(gx2[i], gx1[i], 1e-15);
            EXPECT_NEAR(gy2[i], gy1[i], 1e-15);
            EXPECT_NEAR(gz2[i], gz1[i], 1e-15);
        }
#else
        EXPECT_TRUE(true);
#endif
    });

    run_case("host_emulated_device_roundtrip", [] {
        gpu::HostEmulatedDeviceBuffer b(3*sizeof(double));
        const double src[3] = {1,2,3};
        double dst[3] = {};
        b.copy_from_host(src,sizeof(src));
        b.copy_to_host(dst,sizeof(dst));
        EXPECT_NEAR(dst[0],1.0,1e-14);
        EXPECT_NEAR(dst[2],3.0,1e-14);
    });

    run_case("tile_manager_builds_owned_and_halo_sets", [] {
        ooc::TileManager tm;
        tm.build(8,4,{{0,4},{1,4},{3,7},{4,5}});
        EXPECT_TRUE(tm.tiles().size()==2);
        EXPECT_TRUE(tm.tile(0).owned_cells.size()==4);
        EXPECT_TRUE(tm.tile(0).halo_cells.size()==2);
        EXPECT_TRUE(tm.tile(1).halo_cells.size()==3);
    });

    run_case("working_set_and_bounded_staging_pool", [] {
        ooc::TileManager tm;
        tm.build(4,2,{{1,2}});
        auto ws = ooc::make_working_set(tm.tile(0),2);
        EXPECT_TRUE(ws.values.size()==4);
        ooc::PinnedBufferPool pool(1024,512);
#ifdef CFDX_ENABLE_GPU
        EXPECT_TRUE(ooc::PinnedBufferPool::is_pinned());
#else
        EXPECT_TRUE(!ooc::PinnedBufferPool::is_pinned());
#endif
        auto* a=pool.acquire(); auto* b=pool.acquire();
        EXPECT_TRUE(a!=nullptr && b!=nullptr);
        EXPECT_TRUE(pool.acquire()==nullptr);
        pool.release(a); pool.release(b);
    });

    run_case("async_transfer_and_double_buffer", [] {
        ooc::AsyncTransfer t;
        std::vector<int> h{1,2,3}, d;
        auto f=t.h2d(h,d); f.get();
        EXPECT_TRUE(d==h);
        std::vector<int> back;
        t.d2h(d,back).get();
        EXPECT_TRUE(back==h);
        ooc::DoubleBuffer<int> db(3);
        db.acquire(0)[1]=42;
        EXPECT_TRUE(db.acquire(0)[1]==42);
        EXPECT_TRUE(db.acquire(1)[1]==0);
    });

    return run_all();
}
