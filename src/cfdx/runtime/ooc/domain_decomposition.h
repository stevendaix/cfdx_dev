#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::runtime::ooc {

struct DomainPartition {
    std::vector<std::uint32_t> cell_owner;
    std::vector<std::vector<std::uint32_t>> cells_by_partition;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> cut_edges;
};

inline DomainPartition make_contiguous_partition(
    std::size_t n_cells,
    std::size_t partitions,
    const std::vector<std::pair<std::uint32_t, std::uint32_t>>& adjacency)
{
    if (n_cells == 0) throw std::invalid_argument("domain partition requires at least one cell");
    if (partitions == 0) throw std::invalid_argument("domain partition requires at least one partition");
    partitions = std::min(partitions, n_cells);

    DomainPartition result;
    result.cell_owner.resize(n_cells);
    result.cells_by_partition.resize(partitions);

    const std::size_t base = n_cells / partitions;
    const std::size_t remainder = n_cells % partitions;
    std::size_t cursor = 0;
    for (std::size_t p = 0; p < partitions; ++p) {
        const std::size_t count = base + (p < remainder ? 1U : 0U);
        auto& cells = result.cells_by_partition[p];
        cells.reserve(count);
        for (std::size_t i = 0; i < count; ++i) {
            const auto cell = static_cast<std::uint32_t>(cursor++);
            result.cell_owner[cell] = static_cast<std::uint32_t>(p);
            cells.push_back(cell);
        }
    }

    for (const auto& [a, b] : adjacency) {
        if (a >= n_cells || b >= n_cells || a == b)
            throw std::invalid_argument("domain partition adjacency contains an invalid edge");
        if (result.cell_owner[a] != result.cell_owner[b])
            result.cut_edges.emplace_back(a, b);
    }
    return result;
}

inline void validate_partition(
    const DomainPartition& partition,
    std::size_t n_cells,
    std::size_t expected_partitions)
{
    if (partition.cell_owner.size() != n_cells)
        throw std::invalid_argument("domain partition owner map has the wrong size");
    if (partition.cells_by_partition.size() != expected_partitions)
        throw std::invalid_argument("domain partition count mismatch");

    std::vector<std::uint8_t> seen(n_cells, 0);
    for (std::size_t p = 0; p < partition.cells_by_partition.size(); ++p) {
        for (const auto cell : partition.cells_by_partition[p]) {
            if (cell >= n_cells || seen[cell] != 0)
                throw std::invalid_argument("domain partition does not assign each cell exactly once");
            if (partition.cell_owner[cell] != p)
                throw std::invalid_argument("domain partition owner map disagrees with cell lists");
            seen[cell] = 1;
        }
    }
    if (std::find(seen.begin(), seen.end(), 0) != seen.end())
        throw std::invalid_argument("domain partition leaves cells unassigned");
}

} // namespace cfdx::runtime::ooc
