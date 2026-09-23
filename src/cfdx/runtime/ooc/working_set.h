#pragma once

#include "tile_manager.h"
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime::ooc {

inline std::size_t checked_size_mul(std::size_t a, std::size_t b)
{
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b)
        throw std::overflow_error("WorkingSet size overflow");
    return a * b;
}

struct WorkingSet {
    std::uint64_t tile_id = 0;
    std::size_t field_components = 1;
    std::vector<std::uint32_t> cell_ids;
    std::vector<std::uint32_t> halo_cell_ids;
    std::vector<double> values;
    std::vector<double> halo_values;

    double& value(std::size_t local_cell, std::size_t component = 0)
    {
        return values.at(checked_size_mul(local_cell, field_components) + component);
    }

    const double& value(std::size_t local_cell, std::size_t component = 0) const
    {
        return values.at(checked_size_mul(local_cell, field_components) + component);
    }

    double& halo_value(std::size_t local_halo, std::size_t component = 0)
    {
        return halo_values.at(checked_size_mul(local_halo, field_components) + component);
    }

    const double& halo_value(std::size_t local_halo, std::size_t component = 0) const
    {
        return halo_values.at(checked_size_mul(local_halo, field_components) + component);
    }
};

inline WorkingSet make_working_set(const Tile& tile, std::size_t field_components = 1)
{
    if (field_components == 0)
        throw std::invalid_argument("WorkingSet: components must be > 0");

    WorkingSet ws;
    ws.tile_id = tile.id;
    ws.field_components = field_components;
    ws.cell_ids = tile.owned_cells;
    ws.halo_cell_ids = tile.halo_cells;
    ws.values.resize(checked_size_mul(tile.owned_cells.size(), field_components));
    ws.halo_values.resize(checked_size_mul(tile.halo_cells.size(), field_components));
    return ws;
}

} // namespace cfdx::runtime::ooc
