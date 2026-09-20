#include "cfdx/core/memory/memory_planner.h"
#include <queue>
#include <algorithm>

using namespace cfdx::core::memory;

// ============================================
// Memory Planner Implementation — M0.10-A1-T01/T02/T03/T04/T05
// ============================================
// This file provides skeleton/stub implementations for the Memory Planner
// and Mesh Reordering components defined in memory_planner.md (v4 spec).
// Full production implementations require ExecutionGraph integration,
// mesh connectivity analysis, and hardware profiling.

// --------------------------------------------
// MemoryPlanner (stub)
// --------------------------------------------
MemoryPlanner::Plan MemoryPlanner::plan(
    const std::vector<BufferDescriptor>& buffers,
    int total_operations,
    size_t num_cells,
    size_t num_faces,
    int num_fields,
    int num_solver_vectors) {
    Plan result;

    // Budget prediction (v4 — 6 KPIs réels)
    result.budget.topology_bytes =
        num_faces * sizeof(uint32_t) * 2 +
        num_cells * sizeof(uint32_t) +
        num_faces * sizeof(uint32_t);

    result.budget.geometry_bytes =
        num_cells * sizeof(double) +
        num_faces * sizeof(double) * 3;  // face area + centre + normal approximatif

    result.budget.fields_bytes =
        num_cells * num_fields * sizeof(float);

    result.budget.solver_vectors_bytes =
        num_cells * num_solver_vectors * sizeof(float);

    result.budget.amg_bytes =
        (num_cells > 0) ? static_cast<size_t>(num_cells * sizeof(double) * 2.5) : 0;

    result.budget.temporaries_bytes = 0;  // À calculer via buffer reuse
    result.budget.halo_bytes = 0;         // Dépend du partitionnement MPI

    result.budget.peak_ram =
        result.budget.topology_bytes + result.budget.geometry_bytes +
        result.budget.fields_bytes + result.budget.halo_bytes;

    result.budget.peak_vram =
        result.budget.solver_vectors_bytes + result.budget.amg_bytes +
        result.budget.temporaries_bytes;

    // 6 KPIs (v4 spec)
    result.budget.bytes_per_cell =
        (num_cells > 0) ? static_cast<double>(result.budget.peak_ram + result.budget.peak_vram)
                        / static_cast<double>(num_cells) : 0.0;

    // Traffic : 4 reads + 2 writes per face
    size_t traffic_face = 4 * sizeof(float) + 2 * sizeof(float);
    result.budget.bytes_per_cell_per_iteration =
        (num_cells > 0) ? static_cast<double>(num_faces * traffic_face)
                        / static_cast<double>(num_cells) : 0.0;

    // Vérification faisabilité : budget doit être positif et cohérent
    result.feasible = (num_cells > 0) && (result.budget.bytes_per_cell > 0);
    result.error_message = result.feasible ? "Budget valide" : "Erreur budget";
    return result;
}

// --------------------------------------------
// Memory Ledger (stub)
// --------------------------------------------
void MemoryLedger::allocate(const BufferID& id, void* ptr, size_t size, MemoryLocation loc) {
    (void)id; (void)ptr; (void)size; (void)loc;
}
void MemoryLedger::deallocate(const BufferID& id) { (void)id; }
size_t MemoryLedger::currentUsage(MemoryLocation loc) const { return 0; }
size_t MemoryLedger::peakUsage(MemoryLocation loc) const { return 0; }
void MemoryLedger::printReport() const {}

// ============================================
// RCM Reordering — REAL ALGORITHM (v4 A3-T01)
// ============================================

RCMReorderer::ReorderingResult RCMReorderer::reorder(
    const std::vector<uint32_t>& owner,
    const std::vector<uint32_t>& neighbour,
    size_t num_cells) {
    ReorderingResult result;

    // 1. Compute degrees from connectivity (owner + neighbour)
    std::vector<int> degree(num_cells, 0);
    for (size_t face = 0; face < owner.size(); ++face) {
        uint32_t o = owner[face];
        uint32_t n = neighbour[face];
        if (o < num_cells) degree[o]++;
        if (n != 0xFFFFFFFF && n < num_cells) degree[o]++;
    }

    // 2. Find start node (minimum degree)
    uint32_t start = 0;
    int min_degree = degree[0];
    for (size_t i = 1; i < num_cells; ++i) {
        if (degree[i] < min_degree) {
            min_degree = degree[i];
            start = static_cast<uint32_t>(i);
        }
    }

    // 3. BFS (Cuthill-McKee) with neighbours sorted by degree
    std::vector<bool> visited(num_cells, false);
    std::vector<uint32_t> ordering;
    std::queue<uint32_t> q;
    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
        uint32_t current = q.front(); q.pop();
        ordering.push_back(current);

        // Collect neighbours of this cell from faces
        std::vector<uint32_t> neighbours;
        for (size_t f = 0; f < owner.size(); ++f) {
            uint32_t o = owner[f];
            uint32_t n = neighbour[f];
            if (o == current && n != 0xFFFFFFFF && n < num_cells) {
                if (!visited[n]) neighbours.push_back(n);
            }
            if (n == current && o < num_cells) {
                if (!visited[o]) neighbours.push_back(o);
            }
        }
        std::sort(neighbours.begin(), neighbours.end(),
            [&degree](uint32_t a, uint32_t b) { return degree[a] < degree[b]; });

        for (uint32_t nb : neighbours) {
            if (!visited[nb]) {
                visited[nb] = true;
                q.push(nb);
            }
        }
    }

    // 4. Reverse (RCM = Reverse Cuthill-McKee)
    std::reverse(ordering.begin(), ordering.end());

    // 5. Build mappings
    result.old_to_new.resize(num_cells);
    result.new_to_old.resize(num_cells);
    for (uint32_t new_idx = 0; new_idx < ordering.size(); ++new_idx) {
        uint32_t old_idx = ordering[new_idx];
        result.old_to_new[old_idx] = new_idx;
        result.new_to_old[new_idx] = old_idx;
    }

    // 6. Compute bandwidth before/after
    result.bandwidth_before = computeBandwidth(owner, neighbour, std::vector<uint32_t>{});
    result.bandwidth_after = computeBandwidth(owner, neighbour, result.old_to_new);
    result.profile_before = computeProfile(owner, neighbour, std::vector<uint32_t>{});
    result.profile_after = computeProfile(owner, neighbour, result.old_to_new);

    return result;
}

// Helper: bandwidth = max |owner_new - neighbour_new|
double RCMReorderer::computeBandwidth(const std::vector<uint32_t>& owner,
                                        const std::vector<uint32_t>& neighbour,
                                        const std::vector<uint32_t>& mapping) {
    size_t max_bw = 0;
    for (size_t f = 0; f < owner.size(); ++f) {
        uint32_t o = owner[f];
        uint32_t n = neighbour[f];
        if (n == 0xFFFFFFFF) continue;  // Boundary face
        uint32_t o_new = mapping.empty() ? o : mapping[o];
        uint32_t n_new = mapping.empty() ? n : mapping[n];
        size_t bw = std::abs(static_cast<int>(o_new) - static_cast<int>(n_new));
        max_bw = std::max(max_bw, bw);
    }
    return static_cast<double>(max_bw);
}

// Helper: profile = sum of distances for each cell's neighbours
// (only count neighbours with index < current)
double RCMReorderer::computeProfile(const std::vector<uint32_t>& owner,
                                     const std::vector<uint32_t>& neighbour,
                                     const std::vector<uint32_t>& mapping) {
    size_t profile = 0;
    std::map<uint32_t, std::vector<uint32_t>> cell_neighbours;
    for (size_t f = 0; f < owner.size(); ++f) {
        uint32_t o = owner[f];
        uint32_t n = neighbour[f];
        if (n == 0xFFFFFFFF) continue;
        if (o < cell_neighbours.size()) cell_neighbours[o].push_back(n);
        else {
            while (cell_neighbours.size() <= o) cell_neighbours[o];
            cell_neighbours[o].push_back(n);
        }
        if (n < cell_neighbours.size()) cell_neighbours[n].push_back(o);
        else {
            while (cell_neighbours.size() <= n) cell_neighbours[n];
            cell_neighbours[n].push_back(o);
        }
    }
    for (auto& kv : cell_neighbours) {
        uint32_t cell_new = mapping.empty() ? kv.first : mapping[kv.first];
        for (uint32_t nb_old : kv.second) {
            uint32_t nb_new = mapping.empty() ? nb_old : mapping[nb_old];
            if (nb_new < cell_new) {
                profile += cell_new - nb_new;
            }
        }
    }
    return static_cast<double>(profile);
}

// --------------------------------------------
// Mesh Reordering (stub skeleton)
// --------------------------------------------
MeshReorderer::ReorderingResult MeshReorderer::reorder(
    const std::vector<uint32_t>& owner,
    const std::vector<uint32_t>& neighbour,
    int num_partitions,
    MeshOrdering policy) {
    ReorderingResult result;
    result.selected_strategy = policy;
    result.estimated_speedup = 1.0;
    size_t n = owner.size();
    result.old_to_new.resize(n);
    result.new_to_old.resize(n);
    for (size_t i = 0; i < n; ++i) {
        result.old_to_new[i] = static_cast<uint32_t>(i);
        result.new_to_old[i] = static_cast<uint32_t>(i);
    }
    return result;
}

void MeshReorderer::applyReordering(
    std::vector<uint32_t>& owner,
    std::vector<uint32_t>& neighbour,
    const std::vector<uint32_t>& old_to_new) const {
    for (size_t i = 0; i < owner.size(); ++i) {
        if (owner[i] < old_to_new.size()) {
            owner[i] = old_to_new[owner[i]];
        }
        if (neighbour[i] < old_to_new.size()) {
            neighbour[i] = old_to_new[neighbour[i]];
        }
    }
}
