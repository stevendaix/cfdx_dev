#pragma once

#include "async_transfer.h"
#include "tile_manager.h"
#include "working_set.h"
#include <cstddef>
#include <vector>

namespace cfdx::runtime::ooc {

class DoubleBufferScheduler {
public:
    explicit DoubleBufferScheduler(std::size_t capacity)
        : buffers_{WorkingSet{}, WorkingSet{}}, capacity_(capacity) {}

    template<class Loader, class Compute>
    void run(const TileManager& tiles, Loader load, Compute compute) {
        if (capacity_ == 0) return;
        std::size_t slot = 0;
        for (const auto& tile : tiles.tiles()) {
            auto& ws = buffers_[slot & 1U];
            ws = make_working_set(tile);
            load(tile, ws);
            compute(ws);
            slot++;
        }
    }

private:
    WorkingSet buffers_[2];
    std::size_t capacity_ = 0;
};

} // namespace cfdx::runtime::ooc
