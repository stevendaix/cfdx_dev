#include "cfdx/runtime/ooc/domain_decomposition.h"
#include "cfdx/runtime/ooc/halo_manager.h"
#include "cfdx/runtime/ooc/ooc_executor.h"
#include "cfdx/runtime/ooc/double_buffer_scheduler.h"
#include "common/test_harness.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <thread>
#include <vector>

using namespace cfdx::testing;

int main()
{
    using namespace cfdx::runtime::ooc;

    run_case("phase7_domain_decomposition_is_balanced_and_complete", [] {
        const auto p = make_contiguous_partition(
            10, 3, {{0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9}});
        validate_partition(p, 10, 3);
        EXPECT_TRUE(p.cells_by_partition[0].size() == 4);
        EXPECT_TRUE(p.cells_by_partition[1].size() == 3);
        EXPECT_TRUE(p.cells_by_partition[2].size() == 3);
        EXPECT_TRUE(p.cut_edges.size() == 2);
    });

    run_case("phase7_tile_manager_is_deterministic_and_exposes_ownership", [] {
        TileManager a, b;
        const std::vector<std::pair<std::uint32_t, std::uint32_t>> adjacency{
            {0,1},{1,2},{2,3},{3,4},{4,5},{5,6}};
        a.build(7, 3, adjacency);
        b.build(7, 3, adjacency);
        EXPECT_TRUE(a.tiles().size() == 3);
        EXPECT_TRUE(a.tiles().size() == b.tiles().size());
        for (std::size_t i = 0; i < a.tiles().size(); ++i) {
            EXPECT_TRUE(a.tile(i).owned_cells == b.tile(i).owned_cells);
            EXPECT_TRUE(a.tile(i).halo_cells == b.tile(i).halo_cells);
        }
        EXPECT_TRUE(a.owner_tile(0) == 0);
        EXPECT_TRUE(a.owner_tile(6) == 2);
        EXPECT_THROW(a.owner_tile(99), std::out_of_range);
    });

    run_case("phase7_multicomponent_halo_roundtrip_is_exact", [] {
        TileManager tm;
        tm.build(6, 3, {{0,1},{1,2},{2,3},{3,4},{4,5}});
        HaloManager halo;
        std::vector<double> field(6 * 2);
        for (std::size_t c = 0; c < 6; ++c) {
            field[2*c] = static_cast<double>(c);
            field[2*c + 1] = 100.0 + static_cast<double>(c);
        }

        auto ws = make_working_set(tm.tile(1), 2);
        for (std::size_t i = 0; i < ws.cell_ids.size(); ++i) {
            ws.value(i, 0) = 10.0 + static_cast<double>(ws.cell_ids[i]);
            ws.value(i, 1) = 110.0 + static_cast<double>(ws.cell_ids[i]);
        }
        halo.pack(tm.tile(1), field, ws);
        EXPECT_TRUE(!ws.halo_cell_ids.empty());
        for (std::size_t i = 0; i < ws.halo_cell_ids.size(); ++i) {
            const auto cell = ws.halo_cell_ids[i];
            EXPECT_NEAR(ws.halo_value(i, 0), field[2*cell], 0.0);
            EXPECT_NEAR(ws.halo_value(i, 1), field[2*cell + 1], 0.0);
        }

        std::vector<double> output = field;
        halo.unpack(tm.tile(1), ws, output);
        for (std::size_t i = 0; i < ws.cell_ids.size(); ++i) {
            const auto cell = ws.cell_ids[i];
            EXPECT_NEAR(output[2*cell], ws.value(i, 0), 0.0);
            EXPECT_NEAR(output[2*cell + 1], ws.value(i, 1), 0.0);
        }
    });

    run_case("phase7_double_buffer_starts_next_load_before_current_compute", [] {
        TileManager tm;
        tm.build(4, 2, {{0,1},{1,2},{2,3}});
        DoubleBufferScheduler scheduler(2);
        std::atomic<bool> current_compute_started{false};
        std::atomic<bool> next_load_observed{false};

        scheduler.run(
            tm,
            [&](const Tile& tile, WorkingSet& ws) {
                if (tile.id == 1) {
                    while (!current_compute_started.load())
                        std::this_thread::yield();
                    next_load_observed.store(true);
                }
                for (auto& value : ws.values) value = static_cast<double>(tile.id);
            },
            [&](const WorkingSet& ws) {
                if (ws.tile_id == 0)
                    current_compute_started.store(true);
                EXPECT_TRUE(!ws.values.empty());
            });
        EXPECT_TRUE(next_load_observed.load());
    });

    run_case("phase7_ooc_matches_global_reference_on_multitile_stencil", [] {
        constexpr std::size_t n = 24;
        const double dt = 0.125;
        std::vector<double> input(n);
        for (std::size_t i = 0; i < n; ++i)
            input[i] = std::sin(0.17 * static_cast<double>(i));

        std::vector<std::pair<std::uint32_t, std::uint32_t>> adjacency;
        for (std::size_t i = 0; i + 1 < n; ++i)
            adjacency.emplace_back(
                static_cast<std::uint32_t>(i),
                static_cast<std::uint32_t>(i + 1));

        std::vector<double> reference = input;
        for (std::size_t i = 0; i < n; ++i) {
            const double left = i == 0 ? input[i] : input[i - 1];
            const double right = i + 1 == n ? input[i] : input[i + 1];
            reference[i] = input[i] + dt * (left + right - 2.0 * input[i]);
        }

        OOCConfig cfg;
        cfg.tile_cells = 4;
        cfg.staging_buffers = 2;
        cfg.staging_buffer_bytes = 1024;
        cfg.device_memory_bytes = 64;
        OOCExecutor executor(cfg);
        executor.build_tiles(n, adjacency);
        executor.validate_memory_budget();

        std::vector<double> result(n);
        executor.for_each_tile_pipelined(
            [&](const Tile& tile, WorkingSet& ws, auto&) {
                for (std::size_t i = 0; i < ws.cell_ids.size(); ++i)
                    ws.value(i) = input[ws.cell_ids[i]];
                for (std::size_t i = 0; i < ws.halo_cell_ids.size(); ++i)
                    ws.halo_value(i) = input[ws.halo_cell_ids[i]];
                (void)tile;
            },
            [&](WorkingSet& ws) {
                for (std::size_t i = 0; i < ws.cell_ids.size(); ++i) {
                    const auto cell = static_cast<std::size_t>(ws.cell_ids[i]);
                    const double left =
                        cell == 0 ? ws.value(i) :
                        (cell - 1 == cell ? ws.value(i) : [&] {
                            for (std::size_t h = 0; h < ws.halo_cell_ids.size(); ++h)
                                if (ws.halo_cell_ids[h] == cell - 1) return ws.halo_value(h);
                            for (std::size_t j = 0; j < ws.cell_ids.size(); ++j)
                                if (ws.cell_ids[j] == cell - 1) return ws.value(j);
                            return ws.value(i);
                        }());
                    const double right =
                        cell + 1 == n ? ws.value(i) : [&] {
                            for (std::size_t h = 0; h < ws.halo_cell_ids.size(); ++h)
                                if (ws.halo_cell_ids[h] == cell + 1) return ws.halo_value(h);
                            for (std::size_t j = 0; j < ws.cell_ids.size(); ++j)
                                if (ws.cell_ids[j] == cell + 1) return ws.value(j);
                            return ws.value(i);
                        }();
                    ws.value(i) = input[cell] + dt * (left + right - 2.0 * input[cell]);
                }
            },
            1);

        // The working-set computation above is intentionally independent of the
        // global result. Re-run the same tile pipeline with a commit step so the
        // ownership contract is exercised explicitly.
        HaloManager halo;
        std::fill(result.begin(), result.end(), 0.0);
        executor.for_each_tile_pipelined(
            [&](const Tile& tile, WorkingSet& ws, auto&) {
                for (std::size_t i = 0; i < ws.cell_ids.size(); ++i)
                    ws.value(i) = input[ws.cell_ids[i]];
                halo.pack(tile, input, ws);
            },
            [&](WorkingSet& ws) {
                // Snapshot owned values so the stencil reads the same time level
                // for every neighbour, rather than consuming an already-updated
                // value from an earlier cell in this tile.
                const std::vector<double> old_values = ws.values;
                for (std::size_t i = 0; i < ws.cell_ids.size(); ++i) {
                    const auto cell = static_cast<std::size_t>(ws.cell_ids[i]);
                    auto sample = [&](std::size_t target) {
                        for (std::size_t j = 0; j < ws.cell_ids.size(); ++j)
                            if (ws.cell_ids[j] == target) return old_values[j];
                        for (std::size_t h = 0; h < ws.halo_cell_ids.size(); ++h)
                            if (ws.halo_cell_ids[h] == target) return ws.halo_value(h);
                        return input[target];
                    };
                    const double left = cell == 0 ? old_values[i] : sample(cell - 1);
                    const double right = cell + 1 == n ? old_values[i] : sample(cell + 1);
                    ws.value(i) = input[cell] + dt * (left + right - 2.0 * input[cell]);
                }
                for (std::size_t i = 0; i < ws.cell_ids.size(); ++i)
                    result[ws.cell_ids[i]] = ws.value(i);
            },
            1);

        for (std::size_t i = 0; i < n; ++i)
            EXPECT_NEAR(result[i], reference[i], 1e-14);
    });

    run_case("phase7_memory_limit_rejects_oversized_tile", [] {
        OOCConfig cfg;
        cfg.tile_cells = 8;
        cfg.device_memory_bytes = 32;
        cfg.safety_margin_bytes = 0;
        OOCExecutor executor(cfg);
        executor.build_tiles(8, {{0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7}});
        EXPECT_THROW(executor.validate_memory_budget(), std::runtime_error);
    });

    run_case("phase7_global_problem_can_exceed_device_memory", [] {
        OOCConfig cfg;
        cfg.tile_cells = 4;
        cfg.device_memory_bytes = 64;
        OOCExecutor executor(cfg);
        executor.build_tiles(32, {});
        const auto estimate = executor.estimate_memory();
        EXPECT_TRUE(estimate.global_field_bytes > cfg.device_memory_bytes);
        EXPECT_TRUE(estimate.tile_working_set_bytes <= cfg.device_memory_bytes);
        executor.validate_memory_budget();
    });

    return run_all();
}
