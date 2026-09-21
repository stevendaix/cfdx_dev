#include "cfdx/runtime/ooc/ooc_executor.h"
#include "common/test_harness.h"
#include <cstddef>

using namespace cfdx::testing;

int main()
{
    run_case("ooc_pipeline_executes_every_tile", [] {
        using namespace cfdx::runtime::ooc;
        OOCConfig cfg;
        cfg.tile_cells = 2;
        cfg.staging_buffers = 2;
        cfg.staging_buffer_bytes = 256;
        OOCExecutor executor(cfg);
        executor.build_tiles(6, {{0,2},{1,2},{2,3},{3,4},{4,5}});
        std::size_t loaded = 0;
        std::size_t computed = 0;
        executor.for_each_tile_pipelined(
            [&loaded](const Tile&, WorkingSet& ws, auto&) {
                for (double& value : ws.values) value = 1.0;
                ++loaded;
            },
            [&computed](const WorkingSet& ws) {
                EXPECT_TRUE(!ws.values.empty());
                ++computed;
            });
        EXPECT_TRUE(loaded == executor.tiles().tiles().size());
        EXPECT_TRUE(computed == loaded);
    });
    return run_all();
}
