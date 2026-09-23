#pragma once

#include "tile_manager.h"
#include "working_set.h"
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::ooc {

class HaloManager {
public:
    void pack(const Tile& tile,
              const std::vector<double>& global_values,
              WorkingSet& ws) const
    {
        if (ws.field_components == 0)
            throw std::invalid_argument("HaloManager: working set has zero components");
        if (global_values.size() % ws.field_components != 0)
            throw std::invalid_argument("HaloManager: global field size is not component-aligned");
        if (ws.halo_cell_ids != tile.halo_cells)
            ws.halo_cell_ids = tile.halo_cells;
        ws.halo_values.resize(
            checked_size_mul(tile.halo_cells.size(), ws.field_components));

        const std::size_t global_cells = global_values.size() / ws.field_components;
        for (std::size_t i = 0; i < tile.halo_cells.size(); ++i) {
            const auto cell = static_cast<std::size_t>(tile.halo_cells[i]);
            if (cell >= global_cells)
                throw std::out_of_range("HaloManager: halo cell out of range");
            for (std::size_t c = 0; c < ws.field_components; ++c)
                ws.halo_value(i, c) = global_values[cell * ws.field_components + c];
        }
    }

    void unpack(const Tile& tile,
                const WorkingSet& ws,
                std::vector<double>& global_values) const
    {
        if (ws.field_components == 0)
            throw std::invalid_argument("HaloManager: working set has zero components");
        if (ws.cell_ids != tile.owned_cells)
            throw std::invalid_argument("HaloManager: working-set ownership mismatch");
        if (ws.values.size() != checked_size_mul(tile.owned_cells.size(), ws.field_components))
            throw std::invalid_argument("HaloManager: working-set value size mismatch");
        if (global_values.size() % ws.field_components != 0)
            throw std::invalid_argument("HaloManager: global field size is not component-aligned");

        const std::size_t global_cells = global_values.size() / ws.field_components;
        for (std::size_t i = 0; i < tile.owned_cells.size(); ++i) {
            const auto cell = static_cast<std::size_t>(tile.owned_cells[i]);
            if (cell >= global_cells)
                throw std::out_of_range("HaloManager: owned cell out of range");
            for (std::size_t c = 0; c < ws.field_components; ++c)
                global_values[cell * ws.field_components + c] = ws.value(i, c);
        }
    }
};

} // namespace cfdx::runtime::ooc
