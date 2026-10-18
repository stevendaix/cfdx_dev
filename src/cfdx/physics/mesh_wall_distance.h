#pragma once

#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/wall_distance.h"

#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// Extract the triangulated surface of WALL boundary patches from the solver mesh.
// If patch_names is empty, all and only PatchType::WALL patches are selected.
// Explicitly requested patches must exist and be typed as walls.
inline WallSurface wall_surface_from_mesh(
    const cfdx::core::Mesh& mesh,
    const std::vector<std::string>& patch_names = {})
{
    if (mesh.n_points() == 0 || mesh.n_faces() == 0 || mesh.n_cells() == 0)
        throw std::invalid_argument("wall_surface_from_mesh: mesh is empty");
    if (!mesh.faces().is_consistent() ||
        !mesh.faces().indices_valid(mesh.n_points()))
        throw std::invalid_argument("wall_surface_from_mesh: invalid face connectivity");

    std::set<std::string> requested;
    for (const auto& name : patch_names) {
        if (name.empty() || !requested.insert(name).second)
            throw std::invalid_argument(
                "wall_surface_from_mesh: patch names must be non-empty and unique");
    }

    WallSurface surface;
    surface.points.reserve(mesh.n_points());
    const auto* px = mesh.points().x_data();
    const auto* py = mesh.points().y_data();
    const auto* pz = mesh.points().z_data();
    for (std::size_t i = 0; i < mesh.n_points(); ++i) {
        const WallDistanceVec3 point{px[i], py[i], pz[i]};
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z))
            throw std::invalid_argument(
                "wall_surface_from_mesh: mesh contains non-finite point coordinates");
        surface.points.push_back(point);
    }

    std::set<std::string> found;
    for (const auto& patch : mesh.boundary().patches()) {
        const bool selected = patch_names.empty()
            ? patch.type == cfdx::core::PatchType::WALL
            : requested.count(patch.name) != 0;
        if (!selected)
            continue;
        if (patch.type != cfdx::core::PatchType::WALL)
            throw std::invalid_argument(
                "wall_surface_from_mesh: selected patch '" + patch.name +
                "' is not typed as a wall");
        found.insert(patch.name);
        for (const auto face_id : patch.face_ids) {
            const std::size_t face = static_cast<std::size_t>(face_id);
            if (face >= mesh.n_faces())
                throw std::invalid_argument(
                    "wall_surface_from_mesh: patch face index is out of range");
            if (mesh.ownership().neighbour(face) >= 0)
                throw std::invalid_argument(
                    "wall_surface_from_mesh: wall patch contains an internal face");
            const std::size_t count = mesh.faces().face_size(face);
            if (count < 3)
                throw std::invalid_argument(
                    "wall_surface_from_mesh: wall face has fewer than three vertices");
            const std::size_t offset =
                static_cast<std::size_t>(mesh.faces().face_offset(face));
            const auto* vertices = mesh.faces().vertices_data();
            const auto v0 = static_cast<std::size_t>(vertices[offset]);
            for (std::size_t j = 1; j + 1 < count; ++j) {
                const auto v1 = static_cast<std::size_t>(vertices[offset + j]);
                const auto v2 = static_cast<std::size_t>(vertices[offset + j + 1]);
                const auto ab = surface.points[v1] - surface.points[v0];
                const auto ac = surface.points[v2] - surface.points[v0];
                const WallDistanceVec3 cross{
                    ab.y * ac.z - ab.z * ac.y,
                    ab.z * ac.x - ab.x * ac.z,
                    ab.x * ac.y - ab.y * ac.x};
                const double area2 = wd_norm2(cross);
                if (!std::isfinite(area2) || area2 <= 0.0)
                    throw std::invalid_argument(
                        "wall_surface_from_mesh: degenerate triangle in wall face");
                surface.triangles.push_back({{v0, v1, v2}});
            }
        }
    }

    if (!patch_names.empty() && found != requested) {
        for (const auto& name : requested) {
            if (found.count(name) == 0)
                throw std::invalid_argument(
                    "wall_surface_from_mesh: requested wall patch '" + name +
                    "' does not exist");
        }
    }
    if (surface.triangles.empty())
        throw std::invalid_argument(
            "wall_surface_from_mesh: no usable wall faces; SST wall distance is undefined");
    return surface;
}

// Compute the exact nearest-wall distance at each finite-volume cell centre.
// The BVH avoids the O(N_cells * N_wall_triangles) brute-force scan. A missing
// or zero distance is rejected instead of being silently floored for SST.
inline WallDistanceResult compute_mesh_cell_wall_distance(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const std::vector<std::string>& wall_patch_names = {})
{
    if (geometry.cell_centres.size() != mesh.n_cells())
        throw std::invalid_argument(
            "compute_mesh_cell_wall_distance: cell-centre count does not match mesh");

    const WallSurface surface = wall_surface_from_mesh(mesh, wall_patch_names);
    const WallDistanceBvh bvh(surface);
    WallDistanceResult result;
    result.method = "mesh_exact_bvh";
    result.distance.resize(mesh.n_cells());
    result.valid.assign(mesh.n_cells(), 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const auto& centre = geometry.cell_centres[cell];
        const WallDistanceVec3 point{centre.x, centre.y, centre.z};
        if (!std::isfinite(point.x) || !std::isfinite(point.y) ||
            !std::isfinite(point.z))
            throw std::invalid_argument(
                "compute_mesh_cell_wall_distance: non-finite cell centre");
        const double distance = bvh.nearest_distance(point);
        if (!std::isfinite(distance) || distance <= 0.0)
            throw std::domain_error(
                "compute_mesh_cell_wall_distance: wall distance must be finite and positive");
        result.distance[cell] = distance;
    }
    result.converged = true;
    result.stopping_reason = "exact_cell_centre_distance";
    result.convergence_status = WallDistanceConvergenceStatus::DIRECT;
    return result;
}

}  // namespace cfdx::physics
