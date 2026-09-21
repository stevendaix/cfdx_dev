#pragma once

#include "async_transfer.h"
#include "pinned_buffer_pool.h"
#include "tile_manager.h"
#include "working_set.h"
#include <cstddef>
#include <future>
#include <vector>

namespace cfdx::runtime::ooc {

struct OOCConfig {
    std::size_t tile_cells = 10000;
    std::size_t staging_buffers = 2;
    std::size_t staging_buffer_bytes = 4 * 1024 * 1024;
};

class OOCExecutor {
public:
    explicit OOCExecutor(const OOCConfig& config)
        : config_(config),
          pool_(config.staging_buffers * config.staging_buffer_bytes,
                config.staging_buffer_bytes) {}

    void build_tiles(std::size_t n_cells,
                     const std::vector<std::pair<std::uint32_t,std::uint32_t>>& adjacency) {
        tiles_.build(n_cells, config_.tile_cells, adjacency);
    }

    const TileManager& tiles() const noexcept { return tiles_; }
    PinnedBufferPool& staging_pool() noexcept { return pool_; }

    template<class Loader, class Compute>
    void for_each_tile(Loader load, Compute compute) {
        for (const auto& tile : tiles_.tiles()) {
            auto* buffer = pool_.acquire();
            if (!buffer) throw std::runtime_error("OOCExecutor: staging pool exhausted");
            WorkingSet ws = make_working_set(tile);
            load(tile, ws, *buffer);
            compute(ws);
            pool_.release(buffer);
        }
    }

    // Pipeline tile loading with computation: while tile N is computed, tile
    // N+1 is loaded on a worker thread. This gives the CPU OOC backend real
    // transfer/compute overlap without requiring CUDA.
    template<class Loader, class Compute>
    void for_each_tile_pipelined(Loader load, Compute compute) {
        const auto& tiles = tiles_.tiles();
        if (tiles.empty()) return;

        struct Pending {
            const Tile* tile = nullptr;
            typename PinnedBufferPool::Buffer* buffer = nullptr;
            WorkingSet ws;
            std::future<void> load_future;
        };

        Pending pending;
        auto start_load = [&](const Tile& tile) {
            auto* buffer = pool_.acquire();
            if (!buffer) throw std::runtime_error("OOCExecutor: staging pool exhausted");
            Pending next;
            next.tile = &tile;
            next.buffer = buffer;
            next.ws = make_working_set(tile);
            next.load_future = std::async(
                std::launch::async,
                [&load, &tile, &next]() { load(tile, next.ws, *next.buffer); });
            return next;
        };

        pending = start_load(tiles.front());
        for (std::size_t i = 0; i < tiles.size(); ++i) {
            pending.load_future.get();

            Pending next;
            if (i + 1 < tiles.size()) {
                next = start_load(tiles[i + 1]);
            }

            compute(pending.ws);
            pool_.release(pending.buffer);
            pending = std::move(next);
        }
    }

private:
    OOCConfig config_;
    TileManager tiles_;
    PinnedBufferPool pool_;
};

} // namespace cfdx::runtime::ooc
