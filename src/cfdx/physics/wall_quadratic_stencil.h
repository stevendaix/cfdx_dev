#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/wall_face_reconstruction_quadratic.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Build a one-sided wall stencil from the actual cell-face adjacency graph.
// The wall sample carries the exact prescribed velocity; wall_pressure must
// be supplied by the caller using the pressure boundary-condition contract.
inline std::vector<WallQuadraticSample> build_wall_quadratic_stencil(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    std::size_t wall_face,
    const std::vector<double>& cell_pressure,
    const std::vector<cfdx::core::Vec3>& cell_velocity,
    const cfdx::core::Vec3& wall_velocity,
    double wall_pressure,
    std::size_t max_rings = 3)
{
    using cfdx::core::FaceOwnership;
    using cfdx::core::PatchType;
    using cfdx::core::Vec3;

    if (wall_face >= mesh.n_faces() || geometry.face_centres.size() != mesh.n_faces() ||
        geometry.cell_centres.size() != mesh.n_cells() ||
        cell_pressure.size() != mesh.n_cells() ||
        cell_velocity.size() != mesh.n_cells() || max_rings == 0 ||
        !std::isfinite(wall_pressure) ||
        !std::isfinite(wall_velocity.x) || !std::isfinite(wall_velocity.y) ||
        !std::isfinite(wall_velocity.z))
        throw std::invalid_argument("wall quadratic stencil: invalid input");

    const auto neighbour = mesh.ownership().neighbour(wall_face);
    const auto owner = mesh.ownership().owner(wall_face);
    if (neighbour != FaceOwnership::BOUNDARY || owner >= mesh.n_cells())
        throw std::invalid_argument("wall quadratic stencil: face is not a boundary face");

    if (geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument("wall quadratic stencil: missing face-patch map");
    const std::size_t patch_id = geometry.face_patch[wall_face];
    if (patch_id >= mesh.boundary().n_patches() ||
        mesh.boundary().patch(patch_id).type != PatchType::WALL)
        throw std::invalid_argument("wall quadratic stencil: face does not belong to a WALL patch");

    std::vector<WallQuadraticSample> samples;
    samples.reserve(1 + mesh.n_cells());
    samples.push_back({geometry.face_centres[wall_face], wall_pressure, wall_velocity});

    // Breadth-first traversal gives a deterministic, topological ring stencil.
    // Do not infer neighbours from Euclidean proximity: that can cross walls.
    std::vector<std::size_t> depth(
        mesh.n_cells(), std::numeric_limits<std::size_t>::max());
    std::queue<std::size_t> frontier;
    depth[owner] = 0;
    frontier.push(owner);
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();

    while (!frontier.empty()) {
        const std::size_t cell = frontier.front();
        frontier.pop();
        if (!std::isfinite(cell_pressure[cell]) ||
            !std::isfinite(cell_velocity[cell].x) ||
            !std::isfinite(cell_velocity[cell].y) ||
            !std::isfinite(cell_velocity[cell].z) ||
            !std::isfinite(geometry.cell_centres[cell].x) ||
            !std::isfinite(geometry.cell_centres[cell].y) ||
            !std::isfinite(geometry.cell_centres[cell].z))
            throw std::invalid_argument("wall quadratic stencil: non-finite cell sample");

        samples.push_back({geometry.cell_centres[cell], cell_pressure[cell],
                           cell_velocity[cell]});
        if (depth[cell] >= max_rings) continue;

        for (std::size_t k = cell_offsets[cell]; k < cell_offsets[cell + 1]; ++k) {
            const std::size_t face = cell_faces[k];
            const std::size_t face_owner = mesh.ownership().owner(face);
            const std::int64_t face_neighbour = mesh.ownership().neighbour(face);
            if (face_neighbour < 0) continue;
            const std::size_t other = face_owner == cell
                ? static_cast<std::size_t>(face_neighbour) : face_owner;
            if (other >= mesh.n_cells())
                throw std::invalid_argument("wall quadratic stencil: invalid adjacent cell");
            if (depth[other] <= max_rings) continue;
            depth[other] = depth[cell] + 1;
            frontier.push(other);
        }
    }
    return samples;
}

} // namespace cfdx::physics
