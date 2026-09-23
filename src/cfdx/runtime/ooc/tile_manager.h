#pragma once

#include "domain_decomposition.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <limits>
#include <unordered_map>
#include <vector>

namespace cfdx::runtime::ooc {

struct Tile {
    std::uint64_t id = 0;
    std::vector<std::uint32_t> owned_cells;
    std::vector<std::uint32_t> halo_cells;
    std::vector<std::uint64_t> neighbor_tiles;
};

class TileManager {
public:
    void clear() { tiles_.clear(); cell_to_tile_.clear(); }

    void build(std::size_t n_cells, std::size_t target_cells,
               const std::vector<std::pair<std::uint32_t, std::uint32_t>>& adjacency)
    {
        if (target_cells == 0)
            throw std::invalid_argument("TileManager: target_cells must be > 0");
        clear();

        const std::size_t max = std::numeric_limits<std::size_t>::max();
        if (n_cells > 0 && target_cells > max - (n_cells - 1))
            throw std::overflow_error("TileManager: tile count calculation overflow");
        const std::size_t ntiles = (n_cells + target_cells - 1) / target_cells;
        const auto partition = make_contiguous_partition(n_cells, ntiles, adjacency);
        validate_partition(partition, n_cells, ntiles);

        tiles_.resize(ntiles);
        for (std::size_t t = 0; t < ntiles; ++t) {
            tiles_[t].id = static_cast<std::uint64_t>(t);
            tiles_[t].owned_cells = partition.cells_by_partition[t];
            for (const auto cell : tiles_[t].owned_cells)
                cell_to_tile_[cell] = t;
        }

        for (const auto& [a, b] : adjacency) {
            const auto ta = cell_to_tile_.at(a);
            const auto tb = cell_to_tile_.at(b);
            if (ta == tb) continue;
            tiles_[ta].halo_cells.push_back(b);
            tiles_[tb].halo_cells.push_back(a);
            tiles_[ta].neighbor_tiles.push_back(static_cast<std::uint64_t>(tb));
            tiles_[tb].neighbor_tiles.push_back(static_cast<std::uint64_t>(ta));
        }

        for (auto& tile : tiles_) {
            std::sort(tile.halo_cells.begin(), tile.halo_cells.end());
            tile.halo_cells.erase(
                std::unique(tile.halo_cells.begin(), tile.halo_cells.end()),
                tile.halo_cells.end());
            std::sort(tile.neighbor_tiles.begin(), tile.neighbor_tiles.end());
            tile.neighbor_tiles.erase(
                std::unique(tile.neighbor_tiles.begin(), tile.neighbor_tiles.end()),
                tile.neighbor_tiles.end());
        }
    }

    const std::vector<Tile>& tiles() const noexcept { return tiles_; }
    const Tile& tile(std::size_t i) const { return tiles_.at(i); }

    std::size_t owner_tile(std::uint32_t cell) const
    {
        const auto it = cell_to_tile_.find(cell);
        if (it == cell_to_tile_.end())
            throw std::out_of_range("TileManager: cell is not assigned to a tile");
        return it->second;
    }

private:
    std::vector<Tile> tiles_;
    std::unordered_map<std::uint32_t, std::size_t> cell_to_tile_;
};

} // namespace cfdx::runtime::ooc
