#include "cfdx/core/memory/memory_planner.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace cfdx::core::memory {
namespace {

struct ActiveAllocation {
    BufferReuseOptimizer::Allocation allocation;
    int death_op;
};
// The planner provides deterministic topology/field/solver memory estimates;
// runtime allocation tracking is handled by MemoryLedger.
// The planner provides deterministic topology/field/solver memory estimates;
// runtime allocation tracking is handled by MemoryLedger.

size_t saturating_add(size_t a, size_t b, bool& overflow) {
    if (b > std::numeric_limits<size_t>::max() - a) {
        overflow = true;
        return std::numeric_limits<size_t>::max();
    }
    return a + b;
}

size_t saturating_mul(size_t a, size_t b, bool& overflow) {
    if (a != 0 && b > std::numeric_limits<size_t>::max() / a) {
        overflow = true;
        return std::numeric_limits<size_t>::max();
    }
    return a * b;
}

const char* location_name(MemoryLocation location) {
    switch (location) {
        case MemoryLocation::HOST: return "HOST";
        case MemoryLocation::DEVICE: return "DEVICE";
        case MemoryLocation::MANAGED: return "MANAGED";
    }
    return "UNKNOWN";
}

} // namespace

std::vector<BufferReuseOptimizer::MemoryPool>
BufferReuseOptimizer::optimize(
    const LifetimeAnalyzer::LifetimeResult& lifetimes) const {
    std::vector<MemoryPool> pools;

    std::vector<BufferDescriptor> buffers = lifetimes.buffers;
    std::sort(buffers.begin(), buffers.end(),
              [](const BufferDescriptor& a, const BufferDescriptor& b) {
                  if (a.birth_op != b.birth_op) return a.birth_op < b.birth_op;
                  if (a.size_bytes != b.size_bytes) return a.size_bytes > b.size_bytes;
                  return a.id.id < b.id.id;
              });

    const MemoryLocation locations[] = {
        MemoryLocation::HOST, MemoryLocation::DEVICE, MemoryLocation::MANAGED
    };

    for (const MemoryLocation location : locations) {
        MemoryPool pool{location, 0, {}};
        std::vector<ActiveAllocation> active;
        struct FreeBlock { size_t offset; size_t size; };
        std::vector<FreeBlock> free_blocks;

        for (const auto& buffer : buffers) {
            if (buffer.location != location || buffer.size_bytes == 0) continue;

            active.erase(
                std::remove_if(
                    active.begin(), active.end(),
                    [&](const ActiveAllocation& a) {
                        if (a.death_op < buffer.birth_op) {
                            free_blocks.push_back(
                                {a.allocation.offset, a.allocation.size_bytes});
                            return true;
                        }
                        return false;
                    }),
                active.end());

            // Best-fit reuse minimizes fragmentation while retaining a
            // deterministic allocation order.
            size_t best = free_blocks.size();
            for (size_t i = 0; i < free_blocks.size(); ++i) {
                if (free_blocks[i].size >= buffer.size_bytes &&
                    (best == free_blocks.size() ||
                     free_blocks[i].size < free_blocks[best].size)) {
                    best = i;
                }
            }

            size_t offset = 0;
            if (best != free_blocks.size()) {
                offset = free_blocks[best].offset;
                const size_t remaining =
                    free_blocks[best].size - buffer.size_bytes;
                if (remaining == 0) {
                    free_blocks.erase(free_blocks.begin() +
                                      static_cast<std::ptrdiff_t>(best));
                } else {
                    free_blocks[best].offset += buffer.size_bytes;
                    free_blocks[best].size = remaining;
                }
            } else {
                offset = pool.total_size;
                pool.total_size += buffer.size_bytes;
            }

            BufferReuseOptimizer::Allocation allocation{
                buffer.id, offset, buffer.size_bytes, location};
            pool.allocations.push_back(allocation);
            active.push_back({allocation, buffer.death_op});
        }

        if (!pool.allocations.empty()) {
            pools.push_back(std::move(pool));
        }
    }

    return pools;
}

// MemoryPlanner
// MemoryPlanner
// --------------------------------------------
MemoryPlanner::Plan MemoryPlanner::plan(
    const std::vector<BufferDescriptor>& buffers,
    int total_operations,
    size_t num_cells,
    size_t num_faces,
    int num_fields,
    int num_solver_vectors) const {
    (void)buffers;
    (void)total_operations;
    (void)buffers;
    (void)total_operations;
    Plan result;

    if (total_operations < 0 || num_fields < 0 || num_solver_vectors < 0) {
        result.error_message = "Invalid negative planner parameter";
        return result;
    }

    bool overflow = false;

    if (!buffers.empty()) {
        LifetimeAnalyzer::LifetimeResult lifetimes;
        lifetimes.buffers = buffers;
        lifetimes.total_operations = total_operations;

        for (const auto& buffer : buffers) {
            if (buffer.birth_op < 0 || buffer.death_op < buffer.birth_op ||
                buffer.birth_op > total_operations ||
                buffer.death_op > total_operations) {
                result.error_message = "Invalid buffer lifetime for " +
                                       buffer.id.name;
                return result;
            }
            if (buffer.size_bytes == 0) continue;

            switch (buffer.type) {
                case BufferType::TOPOLOGY:
                    result.budget.topology_bytes =
                        saturating_add(result.budget.topology_bytes,
                                       buffer.size_bytes, overflow);
                    break;
                case BufferType::GEOMETRY:
                    result.budget.geometry_bytes =
                        saturating_add(result.budget.geometry_bytes,
                                       buffer.size_bytes, overflow);
                    break;
                case BufferType::FIELD:
                    result.budget.fields_bytes =
                        saturating_add(result.budget.fields_bytes,
                                       buffer.size_bytes, overflow);
                    break;
                case BufferType::SOLVER_VECTOR:
                    result.budget.solver_vectors_bytes =
                        saturating_add(result.budget.solver_vectors_bytes,
                                       buffer.size_bytes, overflow);
                    break;
                case BufferType::TEMPORARY:
                    result.budget.temporaries_bytes =
                        saturating_add(result.budget.temporaries_bytes,
                                       buffer.size_bytes, overflow);
                    break;
                case BufferType::AMG_LEVEL:
                    result.budget.amg_bytes =
                        saturating_add(result.budget.amg_bytes,
                                       buffer.size_bytes, overflow);
                    break;
                case BufferType::HALO:
                    result.budget.halo_bytes =
                        saturating_add(result.budget.halo_bytes,
                                       buffer.size_bytes, overflow);
                    break;
            }
        }

        result.pools = BufferReuseOptimizer{}.optimize(lifetimes);

        for (const auto& pool : result.pools) {
            if (pool.location == MemoryLocation::HOST ||
                pool.location == MemoryLocation::MANAGED) {
                result.budget.peak_ram =
                    saturating_add(result.budget.peak_ram,
                                   pool.total_size, overflow);
            }
            if (pool.location == MemoryLocation::DEVICE ||
                pool.location == MemoryLocation::MANAGED) {
                result.budget.peak_vram =
                    saturating_add(result.budget.peak_vram,
                                   pool.total_size, overflow);
            }
        }
    } else {
        // Backward-compatible estimate for callers that have not yet built an
        // execution graph. These estimates are deliberately conservative and
        // are replaced by descriptor-based lifetimes as soon as buffers exist.
        result.budget.topology_bytes =
            saturating_add(
                saturating_mul(num_faces, sizeof(uint32_t), overflow),
                saturating_add(
                    saturating_mul(num_faces, sizeof(uint32_t), overflow),
                    saturating_mul(num_cells, sizeof(uint32_t), overflow),
                    overflow),
                overflow);

        result.budget.geometry_bytes =
            saturating_add(
                saturating_mul(num_cells, sizeof(double), overflow),
                saturating_mul(num_faces, sizeof(double) * 3, overflow),
                overflow);

        result.budget.fields_bytes =
            saturating_mul(
                saturating_mul(num_cells, static_cast<size_t>(num_fields), overflow),
                sizeof(double), overflow);

        result.budget.solver_vectors_bytes =
            saturating_mul(
                saturating_mul(num_cells,
                               static_cast<size_t>(num_solver_vectors), overflow),
                sizeof(double), overflow);

        result.budget.amg_bytes =
            saturating_mul(num_cells, sizeof(double) * 3, overflow);

        result.budget.peak_ram =
            saturating_add(result.budget.topology_bytes,
                           result.budget.geometry_bytes, overflow);
        result.budget.peak_ram =
            saturating_add(result.budget.peak_ram,
                           result.budget.fields_bytes, overflow);

        result.budget.peak_vram =
            saturating_add(result.budget.solver_vectors_bytes,
                           result.budget.amg_bytes, overflow);
    }

    const size_t traffic_per_face = 6 * sizeof(double);
    result.budget.bytes_per_cell_per_iteration =
        num_cells > 0
            ? static_cast<double>(
                  saturating_mul(num_faces, traffic_per_face, overflow)) /
                  static_cast<double>(num_cells)
            : 0.0;

    const size_t total_peak =
        saturating_add(result.budget.peak_ram, result.budget.peak_vram, overflow);
    result.budget.bytes_per_cell =
        num_cells > 0
            ? static_cast<double>(total_peak) / static_cast<double>(num_cells)
            : 0.0;

    result.feasible = !overflow && num_cells > 0 &&
                      result.budget.bytes_per_cell > 0.0;
    result.error_message = result.feasible
        ? "Memory plan valid"
        : (overflow ? "Memory budget overflow" : "Invalid/empty mesh budget");

    return result;
}

void MemoryLedger::allocate(const BufferID& id, void* ptr, size_t size,
                            MemoryLocation loc) {
    deallocate(id);
    allocations_[id] = {id, ptr, size, loc, std::chrono::steady_clock::now()};
    current_usage_[loc] += size;
    peak_usage_[loc] = std::max(peak_usage_[loc], current_usage_[loc]);
}

void MemoryLedger::deallocate(const BufferID& id) {
    const auto it = allocations_.find(id);
    if (it == allocations_.end()) return;
    auto usage = current_usage_.find(it->second.location);
    if (usage != current_usage_.end()) {
        usage->second = usage->second >= it->second.size_bytes
            ? usage->second - it->second.size_bytes
            : 0;
    }
    allocations_.erase(it);
}

size_t MemoryLedger::currentUsage(MemoryLocation loc) const {
    const auto it = current_usage_.find(loc);
    return it == current_usage_.end() ? 0 : it->second;
}

size_t MemoryLedger::peakUsage(MemoryLocation loc) const {
    const auto it = peak_usage_.find(loc);
    return it == peak_usage_.end() ? 0 : it->second;
}

void MemoryLedger::printReport() const {
    std::cout << "========== Memory Ledger Report ==========\n";
    std::cout << "Active allocations: " << allocations_.size() << "\n";
    for (const auto& [loc, usage] : current_usage_) {
        std::cout << "Location " << location_name(loc)
                  << " | Current: " << std::setw(10) << usage
                  << " bytes | Peak: " << std::setw(10) << peakUsage(loc)
                  << " bytes\n";
    }
    std::cout << "==========================================\n";
}

RCMReorderer::ReorderingResult RCMReorderer::reorder(
    const std::vector<uint32_t>& owner,
    const std::vector<uint32_t>& neighbour,
    size_t num_cells) {
    ReorderingResult result{};

    if (num_cells == 0) return result;

    std::vector<int> degree(num_cells, 0);
    const size_t nfaces = std::min(owner.size(), neighbour.size());
    for (size_t face = 0; face < nfaces; ++face) {
        const uint32_t o = owner[face];
        const uint32_t n = neighbour[face];
        if (o < num_cells) ++degree[o];
        if (n != 0xFFFFFFFFU && n < num_cells) ++degree[n];
    }

    uint32_t start = 0;
    int min_degree = degree[0];
    for (size_t i = 1; i < num_cells; ++i) {
        if (degree[i] < min_degree) {
            min_degree = degree[i];
            start = static_cast<uint32_t>(i);
        }
    }

    std::vector<bool> visited(num_cells, false);
    std::vector<uint32_t> ordering;
    ordering.reserve(num_cells);
    std::queue<uint32_t> q;
    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        const uint32_t current = q.front();
        q.pop();
        ordering.push_back(current);

        std::vector<uint32_t> neighbours;
        for (size_t f = 0; f < nfaces; ++f) {
            const uint32_t o = owner[f];
            const uint32_t n = neighbour[f];
            if (o == current && n != 0xFFFFFFFFU && n < num_cells &&
                !visited[n]) {
                neighbours.push_back(n);
            } else if (n == current && o < num_cells && !visited[o]) {
                neighbours.push_back(o);
            }
        }

        std::sort(neighbours.begin(), neighbours.end(),
                  [&degree](uint32_t a, uint32_t b) {
                      if (degree[a] != degree[b]) return degree[a] < degree[b];
                      return a < b;
                  });

        for (const uint32_t nb : neighbours) {
            if (!visited[nb]) {
                visited[nb] = true;
                q.push(nb);
            }
        }
    }

    // A disconnected graph is valid; append unvisited components in
    // deterministic minimum-degree order.
    while (ordering.size() < num_cells) {
        uint32_t component_start = 0;
        int component_degree = std::numeric_limits<int>::max();
        for (size_t i = 0; i < num_cells; ++i) {
            if (!visited[i] &&
                (degree[i] < component_degree ||
                 (degree[i] == component_degree &&
                  i < component_start))) {
                component_start = static_cast<uint32_t>(i);
                component_degree = degree[i];
            }
        }
        q.push(component_start);
        visited[component_start] = true;
        while (!q.empty()) {
            const uint32_t current = q.front();
            q.pop();
            ordering.push_back(current);
            std::vector<uint32_t> neighbours;
            for (size_t f = 0; f < nfaces; ++f) {
                const uint32_t o = owner[f];
                const uint32_t n = neighbour[f];
                if (o == current && n != 0xFFFFFFFFU && n < num_cells &&
                    !visited[n]) {
                    neighbours.push_back(n);
                } else if (n == current && o < num_cells && !visited[o]) {
                    neighbours.push_back(o);
                }
            }
            std::sort(neighbours.begin(), neighbours.end(),
                      [&degree](uint32_t a, uint32_t b) {
                          if (degree[a] != degree[b]) return degree[a] < degree[b];
                          return a < b;
                      });
            for (const uint32_t nb : neighbours) {
                if (!visited[nb]) {
                    visited[nb] = true;
                    q.push(nb);
                }
            }
        }
    }

    std::reverse(ordering.begin(), ordering.end());
    result.old_to_new.resize(num_cells);
    result.new_to_old.resize(num_cells);
    for (size_t new_idx = 0; new_idx < ordering.size(); ++new_idx) {
        const uint32_t old_idx = ordering[new_idx];
        result.old_to_new[old_idx] = static_cast<uint32_t>(new_idx);
        result.new_to_old[new_idx] = old_idx;
    }

    result.bandwidth_before = computeBandwidth(owner, neighbour, {});
    result.bandwidth_after = computeBandwidth(
        owner, neighbour, result.old_to_new);
    result.profile_before = computeProfile(owner, neighbour, {});
    result.profile_after = computeProfile(
        owner, neighbour, result.old_to_new);
    return result;
}

double RCMReorderer::computeBandwidth(
    const std::vector<uint32_t>& owner,
    const std::vector<uint32_t>& neighbour,
    const std::vector<uint32_t>& mapping) {
    size_t max_bw = 0;
    const size_t nfaces = std::min(owner.size(), neighbour.size());
    for (size_t f = 0; f < nfaces; ++f) {
        const uint32_t o = owner[f];
        const uint32_t n = neighbour[f];
        if (n == 0xFFFFFFFFU) continue;
        const uint32_t o_new = mapping.empty() ? o : mapping[o];
        const uint32_t n_new = mapping.empty() ? n : mapping[n];
        const size_t bw = o_new > n_new
            ? static_cast<size_t>(o_new - n_new)
            : static_cast<size_t>(n_new - o_new);
        max_bw = std::max(max_bw, bw);
    }
    return static_cast<double>(max_bw);
}

double RCMReorderer::computeProfile(
    const std::vector<uint32_t>& owner,
    const std::vector<uint32_t>& neighbour,
    const std::vector<uint32_t>& mapping) {
    std::vector<std::vector<uint32_t>> adjacency;
    uint32_t max_cell = 0;
    bool has_cell = false;
    const size_t nfaces = std::min(owner.size(), neighbour.size());
    for (size_t f = 0; f < nfaces; ++f) {
        if (owner[f] != 0xFFFFFFFFU) {
            max_cell = std::max(max_cell, owner[f]);
            has_cell = true;
        }
        if (neighbour[f] != 0xFFFFFFFFU) {
            max_cell = std::max(max_cell, neighbour[f]);
            has_cell = true;
        }
    }
    if (!has_cell) return 0.0;
    adjacency.resize(static_cast<size_t>(max_cell) + 1);

    for (size_t f = 0; f < nfaces; ++f) {
        const uint32_t o = owner[f];
        const uint32_t n = neighbour[f];
        if (n == 0xFFFFFFFFU || o == 0xFFFFFFFFU) continue;
        adjacency[o].push_back(n);
        adjacency[n].push_back(o);
    }

    size_t profile = 0;
    for (size_t cell = 0; cell < adjacency.size(); ++cell) {
        const uint32_t cell_new = mapping.empty()
            ? static_cast<uint32_t>(cell)
            : mapping[cell];
        for (const uint32_t nb : adjacency[cell]) {
            const uint32_t nb_new = mapping.empty() ? nb : mapping[nb];
            if (nb_new < cell_new) profile += cell_new - nb_new;
        }
    }
    return static_cast<double>(profile);
}

// Mesh Reordering
// --------------------------------------------
MeshReorderer::ReorderingResult MeshReorderer::reorder(
    const std::vector<uint32_t>& owner,
    const std::vector<uint32_t>& neighbour,
    int num_partitions,
    MeshOrdering policy) {
    (void)num_partitions;
    ReorderingResult result{};
    result.selected_strategy = policy;

    size_t n = 0;
    const size_t nfaces = std::min(owner.size(), neighbour.size());
    for (size_t f = 0; f < nfaces; ++f) {
        if (owner[f] != 0xFFFFFFFFU)
            n = std::max(n, static_cast<size_t>(owner[f]) + 1);
        if (neighbour[f] != 0xFFFFFFFFU)
            n = std::max(n, static_cast<size_t>(neighbour[f]) + 1);
    }

    result.old_to_new.resize(n);
    result.new_to_old.resize(n);
    for (size_t i = 0; i < n; ++i) {
        result.old_to_new[i] = static_cast<uint32_t>(i);
        result.new_to_old[i] = static_cast<uint32_t>(i);
    }

    if (policy == MeshOrdering::RCM && n > 0) {
        RCMReorderer rcm;
        const auto r = rcm.reorder(owner, neighbour, n);
        result.old_to_new = r.old_to_new;
        result.new_to_old = r.new_to_old;
    }

    result.estimated_speedup = 0.0;
    return result;
}

void MeshReorderer::applyReordering(
    std::vector<uint32_t>& owner,
    std::vector<uint32_t>& neighbour,
    const std::vector<uint32_t>& old_to_new) const {
    for (size_t i = 0; i < owner.size(); ++i) {
        if (owner[i] != 0xFFFFFFFFU && owner[i] < old_to_new.size())
            owner[i] = old_to_new[owner[i]];
        if (neighbour[i] != 0xFFFFFFFFU && neighbour[i] < old_to_new.size())
            neighbour[i] = old_to_new[neighbour[i]];
    }
}

} // namespace cfdx::core::memory
