#pragma once

#include "tile_manager.h"
#include "working_set.h"
#include <array>
#include <cstddef>
#include <future>
#include <memory>
#include <utility>

namespace cfdx::runtime::ooc {

class DoubleBufferScheduler {
public:
    explicit DoubleBufferScheduler(std::size_t capacity)
        : capacity_(capacity) {}

    template<class Loader, class Compute>
    void run(const TileManager& tiles, Loader load, Compute compute,
             std::size_t field_components = 1)
    {
        if (capacity_ == 0 || tiles.tiles().empty()) return;

        struct Pending {
            const Tile* tile = nullptr;
            std::shared_ptr<WorkingSet> ws;
            std::future<void> load_future;
        };

        std::array<Pending, 2> pending;
        std::size_t current = 0;

        auto start_load = [&](const Tile& tile, std::size_t slot) {
            Pending p;
            p.tile = &tile;
            p.ws = std::make_shared<WorkingSet>(make_working_set(tile, field_components));
            auto ws = p.ws;
            p.load_future = std::async(
                std::launch::async,
                [&load, &tile, ws]() { load(tile, *ws); });
            pending[slot] = std::move(p);
        };

        start_load(tiles.tile(0), current);

        for (std::size_t i = 0; i < tiles.tiles().size(); ++i) {
            pending[current].load_future.get();

            if (i + 1 < tiles.tiles().size()) {
                const std::size_t next = current ^ 1U;
                start_load(tiles.tile(i + 1), next);
            }

            compute(*pending[current].ws);
            pending[current] = Pending{};
            current ^= 1U;
        }
    }

private:
    std::size_t capacity_ = 0;
};

} // namespace cfdx::runtime::ooc
