#pragma once

#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::runtime {

struct MeshMotionResult {
    double min_volume = 0.0;
    double max_volume = 0.0;
    double min_volume_ratio = 0.0;
    bool valid = false;
};

inline std::vector<double> cell_volumes(const cfdx::core::Mesh& mesh)
{
    using namespace cfdx::core;
    std::vector<Vec3> fc(mesh.n_faces());
    std::vector<Vec3> sf(mesh.n_faces());
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto off = mesh.faces().offsets_data()[f];
        const auto count = mesh.faces().offsets_data()[f + 1] - off;
        const auto g = compute_face_geometry(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(),
            mesh.faces().vertices_data(), off, count);
        fc[f] = g.centre;
        sf[f] = g.Sf;
    }
    std::vector<double> volumes(mesh.n_cells());
    std::vector<Vec3> cell_centres(mesh.n_cells());
    compute_area_weighted_cell_centres(mesh, fc.data(), sf.data(), cell_centres.data());
    orient_mesh_face_vectors(mesh, fc, cell_centres, sf);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto off = mesh.cells().offsets_data()[c];
        const auto count = mesh.cells().offsets_data()[c + 1] - off;
        const auto cg = compute_cell_geometry(mesh, fc.data(), sf.data(), mesh.cells().faces_data() + off, c, count);
        if (!(cg.signed_volume > 0.0))
            throw std::runtime_error("cell_volumes: inverted cell orientation");
        volumes[c] = cg.volume;
    }
    return volumes;
}

inline MeshMotionResult apply_point_displacement(
    cfdx::core::Mesh& mesh,
    const std::vector<cfdx::core::Vec3>& displacement,
    double minimum_volume_ratio = 1e-8)
{
    if (displacement.size() != mesh.n_points())
        throw std::invalid_argument("mesh motion displacement size mismatch");
    if (!(minimum_volume_ratio > 0.0 && minimum_volume_ratio < 1.0))
        throw std::invalid_argument("invalid minimum volume ratio");

    const auto old_volume = cell_volumes(mesh);
    for (std::size_t i = 0; i < mesh.n_points(); ++i) {
        if (!std::isfinite(displacement[i].x) ||
            !std::isfinite(displacement[i].y) ||
            !std::isfinite(displacement[i].z))
            throw std::invalid_argument("mesh motion contains non-finite displacement");
    }

    for (std::size_t i = 0; i < mesh.n_points(); ++i) {
        const auto p = mesh.points();
        mesh.points().set(
            i, p.x(i) + displacement[i].x,
            p.y(i) + displacement[i].y,
            p.z(i) + displacement[i].z);
    }

    const auto new_volume = cell_volumes(mesh);
    MeshMotionResult result;
    result.min_volume = new_volume.empty()
        ? 0.0 : *std::min_element(new_volume.begin(), new_volume.end());
    result.max_volume = new_volume.empty()
        ? 0.0 : *std::max_element(new_volume.begin(), new_volume.end());
    result.min_volume_ratio = 1.0;
    result.valid = true;

    for (std::size_t c = 0; c < new_volume.size(); ++c) {
        if (!(new_volume[c] > 0.0) || !std::isfinite(new_volume[c]) ||
            !(old_volume[c] > 0.0)) {
            result.valid = false;
            result.min_volume_ratio = 0.0;
            break;
        }
        result.min_volume_ratio =
            std::min(result.min_volume_ratio, new_volume[c] / old_volume[c]);
        if (new_volume[c] / old_volume[c] < minimum_volume_ratio)
            result.valid = false;
    }
    return result;
}

struct TopologyChange {
    std::vector<std::size_t> old_to_new_cell;
    std::vector<std::size_t> new_to_old_cell;

    void validate(std::size_t old_cells, std::size_t new_cells) const
    {
        if (old_to_new_cell.size() != old_cells ||
            new_to_old_cell.size() != new_cells)
            throw std::invalid_argument("invalid dynamic-mesh topology map");
        for (const auto i : old_to_new_cell)
            if (i != static_cast<std::size_t>(-1) && i >= new_cells)
                throw std::invalid_argument("old-to-new topology index out of range");
        for (const auto i : new_to_old_cell)
            if (i != static_cast<std::size_t>(-1) && i >= old_cells)
                throw std::invalid_argument("new-to-old topology index out of range");
    }
};

} // namespace cfdx::runtime
