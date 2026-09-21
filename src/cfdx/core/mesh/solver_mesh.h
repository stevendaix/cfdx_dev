#pragma once
#include "cfdx/core/mesh/mesh.h"
#include <cstdint>
#include <vector>

namespace cfdx::core {

// Compact runtime representation for solver kernels. Construction-only
// connectivity (point coordinates and face-vertex lists) can be released
// after this object has been built for a fixed mesh.
struct SolverMesh {
    std::vector<std::uint32_t> owner;
    std::vector<std::int32_t> neighbour;
    std::vector<Vec3> face_area;
    std::vector<Vec3> cell_centres;
    std::vector<double> cell_volume;

    void clear_construction_data() {
        owner.shrink_to_fit();
        neighbour.shrink_to_fit();
    }

    std::size_t n_cells() const noexcept { return cell_centres.size(); }
    std::size_t n_faces() const noexcept { return owner.size(); }
};

inline SolverMesh build_solver_mesh(
    const Mesh& mesh,
    const std::vector<Vec3>& face_area,
    const std::vector<Vec3>& cell_centres,
    const std::vector<double>& cell_volume) {
    SolverMesh out;
    out.owner.resize(mesh.n_faces());
    out.neighbour.resize(mesh.n_faces());
    out.face_area = face_area;
    out.cell_centres = cell_centres;
    out.cell_volume = cell_volume;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        out.owner[f] = static_cast<std::uint32_t>(mesh.ownership().owner(f));
        out.neighbour[f] = static_cast<std::int32_t>(mesh.ownership().neighbour(f));
    }
    return out;
}

} // namespace cfdx::core
