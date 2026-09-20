#pragma once
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <vector>
#include <cmath>
#include <algorithm>
namespace cfdx { namespace core {
inline void assemble_laplacian_csr(
    const Mesh& mesh, SparseMatrix& A, Vector& b,
    const std::vector<bool>& is_dirichlet_cell) {
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();
    A = SparseMatrix(n_cells, n_cells);
    b = Vector(n_cells, 0.0);
    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = mesh.faces().vertices_data();
    const auto* f_offsets = mesh.faces().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const CellConnectivity& cells = mesh.cells();
    const auto* c_faces = cells.faces_data();
    const auto* c_offsets = cells.offsets_data();
    std::vector<Vec3> face_centres(n_faces);
    std::vector<Vec3> face_Sf(n_faces);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const auto off = f_offsets[f];
        const auto n = f_offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }
    std::vector<Vec3> cell_centres(n_cells);
    for (std::size_t c = 0; c < n_cells; ++c) {
        const auto off = c_offsets[c];
        const auto n = c_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry(
            face_centres.data(), face_Sf.data(), c_faces + off, n);
        cell_centres[c] = cg.centre;
    }
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        const int neighbour = own.neighbour(f);
        if (is_dirichlet_cell[owner]) continue;
        const Vec3 Sf = face_Sf[f];
        const double mag_Sf2 = Sf.x * Sf.x + Sf.y * Sf.y + Sf.z * Sf.z;
        if (neighbour >= 0) {
            const std::size_t nb = static_cast<std::size_t>(neighbour);
            const Vec3 d = cell_centres[nb] - cell_centres[owner];
            const double d_dot_Sf = d.x * Sf.x + d.y * Sf.y + d.z * Sf.z;
            const double D_f = mag_Sf2 / std::max(std::abs(d_dot_Sf), 1e-12);
            if (is_dirichlet_cell[nb]) {
                A.push_back(owner, owner, D_f);
            } else {
                A.push_back(owner, owner, D_f);
                A.push_back(owner, nb, -D_f);
            }
        } else {
        }
    }
    for (std::size_t c = 0; c < n_cells; ++c) {
        if (is_dirichlet_cell[c]) {
            A.push_back(c, c, 1.0);
            b(c) = 0.0;
        }
    }
    A.finalize();
}
} } // namespace cfdx::core
