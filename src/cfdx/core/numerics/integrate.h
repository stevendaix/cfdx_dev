// M0.7-T05 — Surface / Volume integration
//
// Spécification CFDX v0.7 §27 :
//   surface_integrate  : Σ_f φ_f * Sf_f  (par cellule, ou globale)
//   volume_integrate   : Σ_c φ_c * V_c   (globale)
//
//   - surface_integrate calcule, pour chaque cellule, la somme des contributions
//     de ses faces. Le résultat est un champ scalaire par cellule.
//   - volume_integrate calcule l'intégrale de volume d'un champ scalaire cellulaire.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/mesh/index_types.h"
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

// Intègre un champ de faces par cellule : Σ_{f ∈ faces(c)} φ_f * Sf_f.
//
// Args:
//   face_field : champ scalaire (dim=1) sur les faces.
//   mesh       : maillage.
//
// Retourne un Field<double, CELL> de dimension 3 (le résultat est un vecteur
// 3D par cellule, puisque Sf_f est un vecteur).
inline Field<double, Location::CELL> surface_integrate(
    const Field<double, Location::FACE>& face_field,
    const Mesh& mesh)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();

    if (face_field.size() != n_faces) {
        throw std::runtime_error("surface_integrate: field size != n_faces");
    }
    if (face_field.dimension() != 1) {
        throw std::runtime_error("surface_integrate: field must be scalar (dim=1)");
    }

    Field<double, Location::CELL> result(n_cells, face_field.name() + "_surf_int", "m^3/s", 3);

    // Géométrie des faces.
    std::vector<Vec3> face_Sf(n_faces);

    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = mesh.faces().vertices_data();
    const auto* offsets = mesh.faces().offsets_data();

    for (std::size_t f = 0; f < n_faces; ++f) {
        const VertexIndex off = offsets[f];
        const VertexIndex n = offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        face_Sf[f] = fg.Sf;
    }

    const FaceOwnership& own = mesh.ownership();
    std::vector<Vec3> provisional_centres(n_cells);
    compute_area_weighted_cell_centres(mesh, face_centres.data(), face_Sf.data(), provisional_centres.data());
    orient_mesh_face_vectors(mesh, face_centres, provisional_centres, face_Sf);

    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    const double* phi = face_field.component_data(0);
    double* r_x = result.component_data(0);
    double* r_y = result.component_data(1);
    double* r_z = result.component_data(2);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;

        Vec3 sum;
        for (Offset k = 0; k < n; ++k) {
            const std::size_t f = cell_faces[off + k];
            const bool is_owner = (own.owner(f) == c);
            const Vec3 Sf_cell = is_owner ? face_Sf[f] : face_Sf[f] * (-1.0);
            sum = sum + Sf_cell * phi[f];
        }
        r_x[c] = sum.x;
        r_y[c] = sum.y;
        r_z[c] = sum.z;
    }

    return result;
}

// Intègre un champ scalaire cellulaire sur le volume total.
//
// Args:
//   cell_field : champ scalaire (dim=1) sur les cellules.
//   mesh       : maillage.
//
// Retourne un double (somme pondérée par le volume).
inline double volume_integrate(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh)
{
    const std::size_t n_cells = mesh.n_cells();

    if (cell_field.size() != n_cells) {
        throw std::runtime_error("volume_integrate: field size != n_cells");
    }
    if (cell_field.dimension() != 1) {
        throw std::runtime_error("volume_integrate: field must be scalar (dim=1)");
    }

    // Géométrie des cellules.
    std::vector<double> cell_volume(n_cells, 0.0);

    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();

    // Géométrie des faces (centres + Sf).
    std::vector<Vec3> face_centres(mesh.n_faces());
    std::vector<Vec3> face_Sf(mesh.n_faces());

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const VertexIndex off = face_offsets[f];
        const VertexIndex n = face_offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry(mesh, face_centres.data(), face_Sf.data(), cell_faces + off, c, n);
        cell_volume[c] = cg.volume;
    }

    const double* phi = cell_field.component_data(0);
    double sum = 0.0;
    for (std::size_t c = 0; c < n_cells; ++c) {
        sum += phi[c] * cell_volume[c];
    }

    return sum;
}

}  // namespace core
}  // namespace cfdx
