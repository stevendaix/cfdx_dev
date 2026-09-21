#pragma once

#include "tile_manager.h"
#include "working_set.h"
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cfdx::runtime::ooc {

class HaloManager {
public:
    // Build a deterministic global-cell -> value lookup and pack each tile's
    // halo values. The same representation works for CPU staging and can be
    // replaced by device buffers without changing the halo contract.
    void pack(const Tile& tile,
              const std::vector<double>& global_values,
              WorkingSet& ws) const
    {
        if (ws.halo_values.size() != tile.halo_cells.size())
            ws.halo_values.resize(tile.halo_cells.size());
        for (std::size_t i = 0; i < tile.halo_cells.size(); ++i) {
            const auto cell = static_cast<std::size_t>(tile.halo_cells[i]);
            if (cell >= global_values.size())
                throw std::out_of_range("HaloManager: halo cell out of range");
            ws.halo_values[i] = global_values[cell];
        }
    }

    // Scatter owned values back into the global field. Halo values are never
    // written to owned storage; ownership is therefore deterministic.
    void unpack(const Tile& tile,
                const WorkingSet& ws,
                std::vector<double>& global_values) const
    {
        if (ws.values.size() < tile.owned_cells.size())
            throw std::invalid_argument("HaloManager: working set has too few values");
        if (global_values.empty() && !tile.owned_cells.empty())
            throw std::invalid_argument("HaloManager: empty global field");
        for (std::size_t i = 0; i < tile.owned_cells.size(); ++i) {
            const auto cell = static_cast<std::size_t>(tile.owned_cells[i]);
            if (cell >= global_values.size())
                throw std::out_of_range("HaloManager: owned cell out of range");
            global_values[cell] = ws.values[i];
        }
    }
};

} // namespace cfdx::runtime::ooc
