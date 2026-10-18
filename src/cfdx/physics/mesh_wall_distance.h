#pragma once

// Method and production-coupling requirements: docs/audits/09-turbulence-rans.md.
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
        !mesh.faces().indices_valid(mesh.n_points()) ||
        mesh.ownership().size() != mesh.n_faces() ||
        !mesh.boundary().is_consistent(mesh.n_faces()))
        throw std::invalid_argument(
            "wall_surface_from_mesh: invalid face, ownership, or patch connectivity");

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
            std::vector<std::size_t> face_vertices;
            face_vertices.reserve(count);
            for (std::size_t j = 0; j < count; ++j)
                face_vertices.push_back(
                    static_cast<std::size_t>(vertices[offset + j]));

            const auto cross = [](const WallDistanceVec3& a,
                                  const WallDistanceVec3& b) {
                return WallDistanceVec3{
                    a.y * b.z - a.z * b.y,
                    a.z * b.x - a.x * b.z,
                    a.x * b.y - a.y * b.x};
            };
            const auto& origin = surface.points[face_vertices[0]];
            WallDistanceVec3 normal{};
            double normal2 = 0.0;
            double scale = 0.0;
            for (std::size_t j = 1; j + 1 < count; ++j) {
                const auto a = surface.points[face_vertices[j]] - origin;
                const auto b = surface.points[face_vertices[j + 1]] - origin;
                normal = cross(a, b);
                normal2 = wd_norm2(normal);
                if (normal2 > 0.0 && std::isfinite(normal2))
                    break;
            }
            for (const auto vertex : face_vertices)
                scale = std::max(
                    scale, wd_norm(surface.points[vertex] - origin));
            if (!(scale > 0.0) || !std::isfinite(scale) ||
                !(normal2 > 0.0) || !std::isfinite(normal2))
                throw std::invalid_argument(
                    "wall_surface_from_mesh: degenerate wall polygon");
            const double normal_length = std::sqrt(normal2);
            const WallDistanceVec3 unit_normal{
                normal.x / normal_length,
                normal.y / normal_length,
                normal.z / normal_length};
            const double planar_tolerance = 1e-10 * scale;
            const double turn_tolerance = 1e-12 * scale * scale;
            for (const auto vertex : face_vertices) {
                const auto delta = surface.points[vertex] - origin;
                if (std::abs(wd_dot(delta, unit_normal)) > planar_tolerance)
                    throw std::invalid_argument(
                        "wall_surface_from_mesh: non-planar wall polygon");
            }

            int turn_sign = 0;
            for (std::size_t j = 0; j < count; ++j) {
                const auto& a = surface.points[face_vertices[j]];
                const auto& b = surface.points[face_vertices[(j + 1) % count]];
                const auto& c = surface.points[face_vertices[(j + 2) % count]];
                const double turn = wd_dot(
                    cross(b - a, c - b), unit_normal);
                if (std::abs(turn) <= turn_tolerance)
                    continue;
                const int sign = turn > 0.0 ? 1 : -1;
                if (turn_sign != 0 && sign != turn_sign)
                    throw std::invalid_argument(
                        "wall_surface_from_mesh: non-convex or self-intersecting wall polygon");
                turn_sign = sign;
            }

            const auto v0 = face_vertices[0];
            for (std::size_t j = 1; j + 1 < count; ++j) {
                const auto v1 = face_vertices[j];
                const auto v2 = face_vertices[j + 1];
                const double area2 = wd_norm2(cross(
                    surface.points[v1] - surface.points[v0],
                    surface.points[v2] - surface.points[v0]));
                if (!std::isfinite(area2) ||
                    area2 <= turn_tolerance * turn_tolerance)
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
