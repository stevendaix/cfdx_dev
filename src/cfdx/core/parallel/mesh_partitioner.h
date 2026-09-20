// M0.11-T02 — Mesh Partitioning
//
// Simple mesh partitioning for parallel execution.
// Supports METIS (if available) and a simple geometric fallback.

#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/parallel/mpi_utils.h"
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <numeric>
#include <random>

namespace cfdx {
namespace core {
namespace parallel {

// Partition result: cell -> rank mapping
struct Partition {
    std::vector<int> cell_rank;           // size = n_cells, value = rank
    std::vector<int> face_owner_rank;     // size = n_faces, value = owner rank
    std::vector<int> face_ghost_rank;     // size = n_faces, value = ghost rank (or -1)
    int n_parts = 0;
};

// Simple geometric partitioning using coordinate-based space-filling curve (Hilbert/Z-order)
inline Partition partition_geometric(const Mesh& m, int n_parts, MPI_Comm comm = MPI_COMM_WORLD) {
    int rank = mpi_rank(comm);
    int size = mpi_size(comm);
    
    if (n_parts != size) {
        if (rank == 0) {
            mpi_print_rank0("[partition_geometric] n_parts must equal mpi_size, adjusting\n");
        }
        n_parts = size;
    }
    
    const std::size_t n_cells = m.n_cells();
    const std::size_t n_faces = m.n_faces();
    
    Partition part;
    part.n_parts = n_parts;
    part.cell_rank.resize(n_cells);
    part.face_owner_rank.resize(n_faces, -1);
    part.face_ghost_rank.resize(n_faces, -1);
    
    if (n_cells == 0) return part;
    
    // Get cell centres for geometric partitioning
    std::vector<Vec3> cell_centres(n_cells);
    const PointCloud& pts = m.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = m.faces().vertices_data();
    const auto* face_offsets = m.faces().offsets_data();
    const auto* cell_faces = m.cells().faces_data();
    const auto* cell_offsets = m.cells().offsets_data();
    
    // Compute face centres first
    std::vector<Vec3> face_centres(n_faces);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const VertexIndex off = face_offsets[f];
        const VertexIndex n = face_offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        face_centres[f] = fg.centre;
    }
    
    // Compute cell centres as average of face centres
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        Vec3 centre;
        for (Offset k = 0; k < n; ++k) {
            const FaceIndex f = cell_faces[off + k];
            centre = centre + face_centres[f];
        }
        centre = centre * (1.0 / n);
        cell_centres[c] = centre;
    }
    
    // Find bounding box
    Vec3 min_corner = cell_centres[0];
    Vec3 max_corner = cell_centres[0];
    for (const auto& c : cell_centres) {
        min_corner.x = std::min(min_corner.x, c.x);
        min_corner.y = std::min(min_corner.y, c.y);
        min_corner.z = std::min(min_corner.z, c.z);
        max_corner.x = std::max(max_corner.x, c.x);
        max_corner.y = std::max(max_corner.y, c.y);
        max_corner.z = std::max(max_corner.z, c.z);
    }
    
    // Normalize coordinates to [0, 1]^3 and compute Morton codes (Z-order)
    struct CellKey {
        std::uint64_t morton;
        std::size_t index;
    };
    
    std::vector<CellKey> keys(n_cells);
    Vec3 range = {max_corner.x - min_corner.x, max_corner.y - min_corner.y, max_corner.z - min_corner.z};
    
    for (std::size_t c = 0; c < n_cells; ++c) {
        double nx = (range.x > 0) ? (cell_centres[c].x - min_corner.x) / range.x : 0.5;
        double ny = (range.y > 0) ? (cell_centres[c].y - min_corner.y) / range.y : 0.5;
        double nz = (range.z > 0) ? (cell_centres[c].z - min_corner.z) / range.z : 0.5;
        
        std::uint32_t ix = static_cast<std::uint32_t>(std::clamp(nx * 1023.0, 0.0, 1023.0));
        std::uint32_t iy = static_cast<std::uint32_t>(std::clamp(ny * 1023.0, 0.0, 1023.0));
        std::uint32_t iz = static_cast<std::uint32_t>(std::clamp(nz * 1023.0, 0.0, 1023.0));
        
        // Interleave bits (Morton code)
        std::uint64_t morton = 0;
        for (int i = 0; i < 10; ++i) {
            morton |= (std::uint64_t)((ix >> i) & 1) << (3 * i);
            morton |= (std::uint64_t)((iy >> i) & 1) << (3 * i + 1);
            morton |= (std::uint64_t)((iz >> i) & 1) << (3 * i + 2);
        }
        
        keys[c] = {morton, c};
    }
    
    // Sort by Morton code
    std::sort(keys.begin(), keys.end(), [](const CellKey& a, const CellKey& b) {
        return a.morton < b.morton;
    });
    
    // Assign cells to ranks evenly
    const std::size_t cells_per_part = (n_cells + n_parts - 1) / n_parts;
    for (std::size_t i = 0; i < n_cells; ++i) {
        int target_rank = static_cast<int>(i / cells_per_part);
        if (target_rank >= n_parts) target_rank = n_parts - 1;
        part.cell_rank[keys[i].index] = target_rank;
    }
    
    // Set face owner ranks
    const FaceOwnership& own = m.ownership();
    for (std::size_t f = 0; f < n_faces; ++f) {
        part.face_owner_rank[f] = part.cell_rank[own.owner(f)];
    }
    
    // Determine ghost faces (faces shared with other ranks)
    for (std::size_t f = 0; f < n_faces; ++f) {
        const int owner_rank = part.face_owner_rank[f];
        const int neighbour = own.neighbour(f);
        
        if (neighbour >= 0) {
            const int neighbour_rank = part.cell_rank[neighbour];
            if (neighbour_rank != owner_rank) {
                part.face_ghost_rank[f] = neighbour_rank;
            }
        }
    }
    
    return part;
}

// Build halo exchange plan from partition
struct HaloPlan {
    std::vector<std::vector<int>> send_faces;  // send_faces[rank] = list of face indices to send
    std::vector<std::vector<int>> recv_faces;  // recv_faces[rank] = list of face indices to receive
};

inline HaloPlan build_halo_plan(const Mesh& m, const Partition& part, MPI_Comm comm = MPI_COMM_WORLD) {
    int rank = mpi_rank(comm);
    int size = mpi_size(comm);
    
    const std::size_t n_faces = m.n_faces();
    const FaceOwnership& own = m.ownership();
    
    HaloPlan plan;
    plan.send_faces.resize(size);
    plan.recv_faces.resize(size);
    
    for (std::size_t f = 0; f < n_faces; ++f) {
        const int owner_rank = part.face_owner_rank[f];
        const int ghost_rank = part.face_ghost_rank[f];
        
        if (ghost_rank >= 0) {
            // Face is shared between owner_rank and ghost_rank
            if (owner_rank == rank) {
                plan.send_faces[ghost_rank].push_back(static_cast<int>(f));
            }
            if (ghost_rank == rank) {
                plan.recv_faces[owner_rank].push_back(static_cast<int>(f));
            }
        }
    }
    
    return plan;
}

// Exchange halo data for a face field
inline void exchange_halo_faces(
    Field<double, Location::FACE>& field,
    const HaloPlan& plan,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    int rank = mpi_rank(comm);
    int size = mpi_size(comm);
    const std::size_t dim = field.dimension();
    
    // Prepare send buffers
    std::vector<std::vector<double>> send_buffers(size);
    for (int r = 0; r < size; ++r) {
        const auto& faces = plan.send_faces[r];
        send_buffers[r].resize(faces.size() * dim);
        for (std::size_t i = 0; i < faces.size(); ++i) {
            int f = faces[i];
            double* src = field.component_data(0) + f * dim;
            std::copy(src, src + dim, send_buffers[r].data() + i * dim);
        }
    }
    
    // Exchange using MPI_Sendrecv
    for (int r = 0; r < size; ++r) {
        if (r == rank) continue;
        
        int send_count = static_cast<int>(plan.send_faces[r].size() * dim);
        int recv_count = static_cast<int>(plan.recv_faces[r].size() * dim);
        
        std::vector<double> recv_buffer(recv_count);
        
        MPI_Status status;
        MPI_Sendrecv(
            send_buffers[r].data(), send_count, MPI_DOUBLE, r, 0,
            recv_buffer.data(), recv_count, MPI_DOUBLE, r, 0,
            comm, &status
        );
        
        // Copy received data back to field
        const auto& recv_faces = plan.recv_faces[r];
        for (std::size_t i = 0; i < recv_faces.size(); ++i) {
            int f = recv_faces[i];
            double* dst = field.component_data(0) + f * dim;
            std::copy(recv_buffer.data() + i * dim, recv_buffer.data() + (i + 1) * dim, dst);
        }
    }
}

// Exchange halo data for a cell field
inline void exchange_halo_cells(
    Field<double, Location::CELL>& field,
    const HaloPlan& plan,
    MPI_Comm comm = MPI_COMM_WORLD)
{
    // For cell fields, we need to exchange cells that are adjacent to ghost faces
    // This is more complex - for now, implement basic face-based exchange
    // and derive cell values from face interpolation if needed
    (void)field; (void)plan; (void)comm;
    // TODO: Implement proper cell halo exchange
}

}  // namespace parallel
}  // namespace core
}  // namespace cfdx
