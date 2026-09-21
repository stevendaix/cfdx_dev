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
        std::size_t slot = 0;
        for (const auto& tile : tiles_.tiles()) {
            auto* buffer = pool_.acquire();
            if (!buffer) throw std::runtime_error("OOCExecutor: staging pool exhausted");
            WorkingSet ws = make_working_set(tile);
            load(tile, ws, *buffer);
            compute(ws);
            pool_.release(buffer);
            ++slot;
        }
    }

private:
    OOCConfig config_;
    TileManager tiles_;
    PinnedBufferPool pool_;
};

} // namespace cfdx::runtime::ooc
