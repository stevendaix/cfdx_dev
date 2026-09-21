#pragma once

#include "tile_manager.h"
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::ooc {

struct WorkingSet {
    std::uint64_t tile_id = 0;
    std::vector<std::uint32_t> cell_ids;
    std::vector<double> values;
    std::vector<double> halo_values;
};

inline WorkingSet make_working_set(const Tile& tile, std::size_t field_components = 1) {
    if (field_components == 0) throw std::invalid_argument("WorkingSet: components must be > 0");
    WorkingSet ws;
    ws.tile_id = tile.id;
    ws.cell_ids = tile.owned_cells;
    ws.values.resize(tile.owned_cells.size() * field_components);
    ws.halo_values.resize(tile.halo_cells.size() * field_components);
    return ws;
}

} // namespace cfdx::runtime::ooc
