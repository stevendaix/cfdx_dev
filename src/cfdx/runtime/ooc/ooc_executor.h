#pragma once

#include "async_transfer.h"
#include "double_buffer_scheduler.h"
#include "halo_manager.h"
#include "pinned_buffer_pool.h"
#include "tile_manager.h"
#include "working_set.h"
#include <algorithm>
#include <cstddef>
#include <future>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::runtime::ooc {

struct OOCConfig {
    std::size_t tile_cells = 10000;
    std::size_t staging_buffers = 2;
    std::size_t staging_buffer_bytes = 4 * 1024 * 1024;
    std::size_t device_memory_bytes = 0;
    std::size_t safety_margin_bytes = 0;
};

struct OOCMemoryEstimate {
    std::size_t global_field_bytes = 0;
    std::size_t tile_owned_bytes = 0;
    std::size_t tile_halo_bytes = 0;
    std::size_t tile_working_set_bytes = 0;
};

class OOCExecutor {
public:
    explicit OOCExecutor(const OOCConfig& config)
        : config_(config),
          pool_(checked_capacity(config.staging_buffers, config.staging_buffer_bytes),
                config.staging_buffer_bytes) {}

    void build_tiles(std::size_t n_cells,
                     const std::vector<std::pair<std::uint32_t,std::uint32_t>>& adjacency)
    {
        tiles_.build(n_cells, config_.tile_cells, adjacency);
        validate_tiles();
    }

    const TileManager& tiles() const noexcept { return tiles_; }
    PinnedBufferPool& staging_pool() noexcept { return pool_; }

    OOCMemoryEstimate estimate_memory(std::size_t field_components = 1) const
    {
        if (field_components == 0)
            throw std::invalid_argument("OOCExecutor: components must be > 0");
        if (tiles_.tiles().empty())
            return {};

        std::size_t max_owned = 0;
        std::size_t max_halo = 0;
        std::size_t total_owned = 0;
        for (const auto& tile : tiles_.tiles()) {
            max_owned = std::max(max_owned, tile.owned_cells.size());
            max_halo = std::max(max_halo, tile.halo_cells.size());
            total_owned = checked_add(total_owned, tile.owned_cells.size());
        }
        const auto scalar = sizeof(double);
        OOCMemoryEstimate e;
        e.global_field_bytes = checked_mul(total_owned, checked_mul(field_components, scalar));
        e.tile_owned_bytes = checked_mul(max_owned, checked_mul(field_components, scalar));
        e.tile_halo_bytes = checked_mul(max_halo, checked_mul(field_components, scalar));
        e.tile_working_set_bytes = checked_add(e.tile_owned_bytes, e.tile_halo_bytes);
        return e;
    }

    void validate_memory_budget(std::size_t field_components = 1) const
    {
        if (config_.device_memory_bytes == 0) return;
        const auto e = estimate_memory(field_components);
        if (config_.safety_margin_bytes > config_.device_memory_bytes)
            throw std::runtime_error("OOCExecutor: safety margin exceeds device memory");
        const auto usable = config_.device_memory_bytes - config_.safety_margin_bytes;
        if (e.tile_working_set_bytes > usable)
            throw std::runtime_error("OOCExecutor: one tile working set exceeds the device memory budget");
        if (e.global_field_bytes <= usable)
            throw std::runtime_error(
                "OOCExecutor: configured workload fits in device memory; GPU-OOC is not required");
    }

    template<class Loader, class Compute>
    void for_each_tile(Loader load, Compute compute, std::size_t field_components = 1)
    {
        validate_memory_budget_if_configured(field_components);
        for (const auto& tile : tiles_.tiles()) {
            auto* buffer = pool_.acquire();
            if (!buffer) throw std::runtime_error("OOCExecutor: staging pool exhausted");
            WorkingSet ws = make_working_set(tile, field_components);
            try {
                load(tile, ws, *buffer);
                compute(ws);
                pool_.release(buffer);
            } catch (...) {
                pool_.release(buffer);
                throw;
            }
        }
    }

    template<class Loader, class Compute>
    void for_each_tile_pipelined(Loader load, Compute compute,
                                 std::size_t field_components = 1)
    {
        validate_memory_budget_if_configured(field_components);
        const auto& tiles = tiles_.tiles();
        if (tiles.empty()) return;

        struct Pending {
            const Tile* tile = nullptr;
            typename PinnedBufferPool::Buffer* buffer = nullptr;
            std::shared_ptr<WorkingSet> ws;
            std::future<void> load_future;

            void release(PinnedBufferPool& pool) noexcept
            {
                if (buffer) {
                    try { pool.release(buffer); } catch (...) {}
                    buffer = nullptr;
                }
                tile = nullptr;
                ws.reset();
                load_future = std::future<void>{};
            }
        };

        Pending pending;
        auto start_load = [&](const Tile& tile) {
            auto* buffer = pool_.acquire();
            if (!buffer) throw std::runtime_error("OOCExecutor: staging pool exhausted");
            Pending next;
            next.tile = &tile;
            next.buffer = buffer;
            next.ws = std::make_shared<WorkingSet>(make_working_set(tile, field_components));
            auto* next_buffer = next.buffer;
            auto next_ws = next.ws;
            auto* tile_ptr = next.tile;
            try {
                next.load_future = std::async(
                    std::launch::async,
                    [&load, tile_ptr, next_ws, next_buffer]() {
                        load(*tile_ptr, *next_ws, *next_buffer);
                    });
            } catch (...) {
                pool_.release(next.buffer);
                throw;
            }
            return next;
        };

        pending = start_load(tiles.front());
        for (std::size_t i = 0; i < tiles.size(); ++i) {
            Pending next;
            try {
                pending.load_future.get();
                if (i + 1 < tiles.size())
                    next = start_load(tiles[i + 1]);

                compute(*pending.ws);
                pending.release(pool_);
                pending = std::move(next);
            } catch (...) {
                pending.release(pool_);
                next.release(pool_);
                throw;
            }
        }
    }

private:
    static std::size_t checked_add(std::size_t a, std::size_t b)
    {
        if (b > std::numeric_limits<std::size_t>::max() - a)
            throw std::overflow_error("OOCExecutor size overflow");
        return a + b;
    }

    static std::size_t checked_mul(std::size_t a, std::size_t b)
    {
        if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
            throw std::overflow_error("OOCExecutor size overflow");
        return a * b;
    }

    static std::size_t checked_capacity(std::size_t buffers, std::size_t bytes)
    {
        if (buffers == 0 || bytes == 0)
            throw std::invalid_argument("OOCExecutor: staging capacities must be > 0");
        return checked_mul(buffers, bytes);
    }

    void validate_tiles() const
    {
        if (tiles_.tiles().empty())
            throw std::invalid_argument("OOCExecutor: tile decomposition is empty");
        for (const auto& tile : tiles_.tiles()) {
            for (const auto cell : tile.owned_cells)
                if (tiles_.owner_tile(cell) != tile.id)
                    throw std::logic_error("OOCExecutor: inconsistent tile ownership");
        }
    }

    void validate_memory_budget_if_configured(std::size_t components) const
    {
        if (config_.device_memory_bytes != 0)
            validate_memory_budget(components);
    }

    OOCConfig config_;
    TileManager tiles_;
    PinnedBufferPool pool_;
};

} // namespace cfdx::runtime::ooc
