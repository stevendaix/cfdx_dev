#include "cfdx/runtime/execution/execution_policy.h"
#include "cfdx/runtime/gpu/device_buffer.h"
#include "cfdx/runtime/gpu/gpu_kernels.h"
#include "cfdx/runtime/gpu/gpu_execution.h"
#include "cfdx/runtime/ooc/tile_manager.h"
#include "cfdx/runtime/ooc/working_set.h"
#include "cfdx/runtime/ooc/pinned_buffer_pool.h"
#include "cfdx/runtime/ooc/async_transfer.h"
#include "cfdx/runtime/gpu/gpu_profiler.h"
#ifdef CFDX_ENABLE_GPU
#include "cfdx/runtime/ooc/cuda_pinned_buffer_pool.h"
#include "cfdx/runtime/gpu/cuda_double_buffer.h"
#endif
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

    run_case("selected_gpu_never_falls_back_to_cpu", [] {
        RuntimeDecision d;
        d.selected = ExecutionPolicy::GPU;
        d.requires_gpu = true;
        EXPECT_THROW(require_selected_backend(d, false, true), std::runtime_error);
        bool completed = false;
        try {
            require_selected_backend(d, true, true);
            completed = true;
        } catch (...) {
        }
        EXPECT_TRUE(completed);
    });

    run_case("gpu_ooc_never_falls_back_to_cpu", [] {
        RuntimeDecision d;
        d.selected = ExecutionPolicy::GPU_OUT_OF_CORE;
        d.requires_gpu = true;
        d.uses_out_of_core = true;
        EXPECT_THROW(require_selected_backend(d, false, true), std::runtime_error);
        bool completed = false;
        try {
            require_selected_backend(d, true, true);
            completed = true;
        } catch (...) {
        }
        EXPECT_TRUE(completed);
    });

    run_case("runtime_decision_cannot_lose_gpu_requirement", [] {
        RuntimeDecision d;
        d.selected = ExecutionPolicy::CPU;
        d.requires_gpu = true;
        EXPECT_THROW(require_selected_backend(d, true, true), std::logic_error);
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

    run_case("gpu_gradient_executes_on_cuda_when_available", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;
        const std::vector<double> phi{0.0, 1.0};
        const std::vector<double> sx{1.0};
        const std::vector<double> sy{0.0};
        const std::vector<double> sz{0.0};
        const std::vector<std::uint32_t> owner{0};
        const std::vector<std::int64_t> neighbour{1};
        const std::vector<double> volume{1.0, 1.0};
        std::vector<double> gx, gy, gz;
        gpu::execute_gradient_cuda(phi, sx, sy, sz, owner, neighbour, volume, gx, gy, gz);
        EXPECT_NEAR(gx[0], 0.5, 1e-12);
        EXPECT_NEAR(gx[1], -0.5, 1e-12);
        EXPECT_NEAR(gy[0], 0.0, 1e-14);
        EXPECT_NEAR(gz[1], 0.0, 1e-14);
        const std::vector<double> face_flux{2.0};
        std::vector<double> div;
        gpu::execute_divergence_cuda(face_flux, owner, neighbour, 2, div);
        EXPECT_NEAR(div[0], 2.0, 1e-12);
        EXPECT_NEAR(div[1], -2.0, 1e-12);
#else
        std::vector<double> gx, gy, gz;
        EXPECT_THROW(
            gpu::execute_gradient_cuda({}, {}, {}, {}, {}, {}, {}, gx, gy, gz),
            std::runtime_error);
#endif
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

    run_case("cuda_pinned_pool_reuses_and_enforces_capacity", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;
        ooc::CudaPinnedBufferPool pool(1024, 512);
        EXPECT_TRUE(pool.size() == 2);
        EXPECT_TRUE(pool.buffer_bytes() == 512);
        auto* a = pool.acquire();
        auto* b = pool.acquire();
        EXPECT_TRUE(a != nullptr && b != nullptr);
        EXPECT_TRUE(pool.acquire() == nullptr);
        pool.release(a);
        EXPECT_TRUE(pool.acquire() == a);
        pool.release(b);
        pool.release(a);
#else
        EXPECT_TRUE(true);
#endif
    });

    run_case("cuda_pinned_pool_rejects_foreign_buffer", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;
        ooc::CudaPinnedBufferPool pool(512, 512);
        ooc::CudaPinnedBufferPool other(512, 512);
        auto* b = other.acquire();
        EXPECT_THROW(pool.release(b), std::invalid_argument);
        other.release(b);
#else
        EXPECT_TRUE(true);
#endif
    });

    run_case("cuda_async_h2d_d2h_roundtrip", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;

        const double source[4] = {1.0, -2.0, 3.5, 8.0};
        double destination[4] = {};
        gpu::CudaDeviceBuffer device(sizeof(source));
        gpu::CudaStream stream;
        gpu::async_copy_h2d(device, source, sizeof(source), stream.get());
        stream.synchronize();
        gpu::async_copy_d2h(device, destination, sizeof(destination), stream.get());
        stream.synchronize();

        for (std::size_t i = 0; i < 4; ++i)
            EXPECT_NEAR(destination[i], source[i], 0.0);
#else
        EXPECT_TRUE(true);
#endif
    });

    run_case("cuda_double_buffer_has_two_independent_streams_and_buffers", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;
        gpu::CudaDoubleBuffer db(256);
        EXPECT_TRUE(gpu::CudaDoubleBuffer::slot_count() == 2);
        EXPECT_TRUE(db.bytes() == 256);
        EXPECT_TRUE(db.buffer(0).data() != db.buffer(1).data());
        EXPECT_TRUE(db.stream(0).get() != db.stream(1).get());
        const int src0[4] = {1,2,3,4};
        const int src1[4] = {5,6,7,8};
        gpu::async_copy_h2d(db.buffer(0), src0, sizeof(src0), db.stream(0).get());
        gpu::async_copy_h2d(db.buffer(1), src1, sizeof(src1), db.stream(1).get());
        db.synchronize(0);
        db.synchronize(1);
        int out0[4] = {}, out1[4] = {};
        gpu::async_copy_d2h(db.buffer(0), out0, sizeof(out0), db.stream(0).get());
        gpu::async_copy_d2h(db.buffer(1), out1, sizeof(out1), db.stream(1).get());
        db.synchronize(0);
        db.synchronize(1);
        for (int i = 0; i < 4; ++i) {
            EXPECT_TRUE(out0[i] == src0[i]);
            EXPECT_TRUE(out1[i] == src1[i]);
        }
#else
        EXPECT_TRUE(true);
#endif
    });

    run_case("cuda_profiler_measures_stream_work", [] {
#ifdef CFDX_ENABLE_GPU
        int devices = 0;
        const auto status = cudaGetDeviceCount(&devices);
        if (status != cudaSuccess || devices == 0) return;
        gpu::CudaStream stream;
        gpu::GpuProfiler profiler;
        profiler.start(stream.get());
        const auto launch = cudaDeviceSynchronize();
        EXPECT_TRUE(launch == cudaSuccess);
        profiler.stop(stream.get());
        EXPECT_TRUE(profiler.elapsed_ms() >= 0.0);
#else
        gpu::GpuProfiler profiler;
        EXPECT_THROW(profiler.start(), std::runtime_error);
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
