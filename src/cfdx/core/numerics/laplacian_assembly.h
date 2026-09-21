#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx { namespace core {

// Assemble a cell-centred Laplacian.
//
// is_dirichlet_cell identifies cells whose equation is replaced by a strong
// Dirichlet constraint. When dirichlet_values is supplied, its value for each
// constrained cell is imposed in the right-hand side. An empty vector keeps
// the historical homogeneous Dirichlet behaviour.
inline void assemble_laplacian_csr(
    const Mesh& mesh,
    SparseMatrix& A,
    Vector& b,
    const std::vector<bool>& is_dirichlet_cell,
    const std::vector<double>& dirichlet_values = {}) {

    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();

    if (n_cells == 0 || n_faces == 0) {
        throw std::runtime_error(
            "assemble_laplacian_csr: mesh is empty (n_cells=" +
            std::to_string(n_cells) + ", n_faces=" +
            std::to_string(n_faces) + ")");
    }
    if (mesh.cells().offsets_data() == nullptr ||
        mesh.cells().n_cells() != n_cells) {
        throw std::runtime_error(
            "assemble_laplacian_csr: cell connectivity is invalid or incomplete. "
            "Expected n_cells=" + std::to_string(n_cells));
    }
    if (mesh.ownership().size() != n_faces) {
        throw std::runtime_error(
            "assemble_laplacian_csr: ownership size mismatch. Expected " +
            std::to_string(n_faces) + ", got " +
            std::to_string(mesh.ownership().size()));
    }
    if (is_dirichlet_cell.size() != n_cells) {
        throw std::runtime_error(
            "assemble_laplacian_csr: is_dirichlet_cell size mismatch. Expected " +
            std::to_string(n_cells) + ", got " +
            std::to_string(is_dirichlet_cell.size()));
    }
    if (!dirichlet_values.empty() && dirichlet_values.size() != n_cells) {
        throw std::runtime_error(
            "assemble_laplacian_csr: dirichlet_values size mismatch. Expected " +
            std::to_string(n_cells) + ", got " +
            std::to_string(dirichlet_values.size()));
    }
    for (std::size_t c = 0; c < n_cells; ++c) {
        if (!dirichlet_values.empty() &&
            is_dirichlet_cell[c] &&
            !std::isfinite(dirichlet_values[c])) {
            throw std::runtime_error(
                "assemble_laplacian_csr: Dirichlet value must be finite");
        }
    }

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

        if (n == 2) {
            const std::size_t n0 = verts[off];
            const std::size_t n1 = verts[off + 1];
            face_centres[f] = Vec3{
                (px[n0] + px[n1]) * 0.5,
                (py[n0] + py[n1]) * 0.5,
                (pz[n0] + pz[n1]) * 0.5
            };
            const double dx = px[n1] - px[n0];
            const double dy = py[n1] - py[n0];
            face_Sf[f] = Vec3{-dy, dx, 0.0};
        } else if (n >= 3) {
            const FaceGeometry fg =
                compute_face_geometry(px, py, pz, verts, off, n);
            face_centres[f] = fg.centre;
            face_Sf[f] = fg.Sf;
        } else {
            throw std::runtime_error(
                "assemble_laplacian_csr: face " + std::to_string(f) +
                " has invalid vertex count: " + std::to_string(n));
        }
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

        // The Dirichlet row is replaced by an identity equation below.
        if (is_dirichlet_cell[owner]) continue;

        const Vec3 Sf = face_Sf[f];
        const double mag_Sf2 =
            Sf.x * Sf.x + Sf.y * Sf.y + Sf.z * Sf.z;

        if (neighbour >= 0) {
            const std::size_t nb = static_cast<std::size_t>(neighbour);
            const Vec3 d = cell_centres[nb] - cell_centres[owner];
            const double d_dot_Sf =
                d.x * Sf.x + d.y * Sf.y + d.z * Sf.z;
            const double D_f =
                mag_Sf2 / std::max(std::abs(d_dot_Sf), 1e-12);

            A.push_back(owner, owner, D_f);
            if (!is_dirichlet_cell[nb]) {
                A.push_back(owner, nb, -D_f);
            }
        } else {
            // Boundary faces currently represent homogeneous Dirichlet
            // diffusion. Patch-specific non-zero values belong in the
            // boundary-condition layer, not in this mesh-only assembler.
            const Vec3 d = face_centres[f] - cell_centres[owner];
            const double d_dot_Sf =
                d.x * Sf.x + d.y * Sf.y + d.z * Sf.z;
            const double D_f =
                mag_Sf2 / std::max(std::abs(d_dot_Sf), 1e-12);
            A.push_back(owner, owner, D_f);
        }
    }

    for (std::size_t c = 0; c < n_cells; ++c) {
        if (!is_dirichlet_cell[c]) continue;

        // Strong Dirichlet enforcement avoids the 1e15 penalty used before,
        // which unnecessarily damages the conditioning of the linear system.
        A.push_back(c, c, 1.0);
        b(c) = dirichlet_values.empty() ? 0.0 : dirichlet_values[c];
    }

    A.finalize();
}

}} // namespace cfdx::core
