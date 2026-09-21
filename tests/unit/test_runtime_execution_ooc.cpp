#include "cfdx/runtime/execution/execution_policy.h"
#include "cfdx/runtime/gpu/device_buffer.h"
#include "cfdx/runtime/gpu/gpu_kernels.h"
#include "cfdx/runtime/ooc/tile_manager.h"
#include "cfdx/runtime/ooc/working_set.h"
#include "cfdx/runtime/ooc/pinned_buffer_pool.h"
#include "cfdx/runtime/ooc/async_transfer.h"
#include "common/test_harness.h"
#include <cstdint>
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
    });

    run_case("gpu_policy_requires_cuda", [] {
        RuntimeCapabilities caps;
        EXPECT_THROW(choose_execution_policy(ExecutionPolicy::GPU, caps, {}), std::runtime_error);
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
        EXPECT_TRUE(tm.tile(0).halo_cells.size()==1);
        EXPECT_TRUE(tm.tile(1).halo_cells.size()==1);
    });

    run_case("working_set_and_bounded_staging_pool", [] {
        ooc::TileManager tm;
        tm.build(4,2,{{1,2}});
        auto ws = ooc::make_working_set(tm.tile(0),2);
        EXPECT_TRUE(ws.values.size()==4);
        ooc::PinnedBufferPool pool(1024,512);
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
