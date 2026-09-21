// M0.7-T01 — Gauss gradient
//
// Spécification CFDX v0.7 §28 :
//   Première implémentation : Gauss linear.
//
// Pour un champ scalaire φ (Field<double, CELL>) :
//   ∇φ_c = (1/V_c) Σ_f φ_f Sf_f
//
//   - φ_f est interpolé géométriquement à partir des cellules owner/neighbour.
//     φ_f = (1-w) * φ_owner + w * φ_neighbour
//     w = |Cf - C_owner| / |C_neighbour - C_owner|  (pondération géométrique)
//   - Pour les faces de frontière, on utilise φ_f = φ_owner (zero-gradient).
//
// Le résultat est un Field<double, CELL> de dimension 3 (le gradient est un vecteur
// 3D par cellule). Le stockage est SoA (x, y, z contigus).

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/mesh/index_types.h"
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx {
namespace core {

// Calcule le poids d'interpolation géométrique pour une face.
// w = |Cf - C_owner| / |C_neighbour - C_owner|
// Pour face frontière : w = 0 (zero gradient)
inline double geometric_interpolation_weight(
    const Vec3& face_centre,
    const Vec3& owner_centre,
    const Vec3* neighbour_centre)
{
    if (!neighbour_centre) return 0.0;  // Boundary face

    const Vec3 d_owner = face_centre - owner_centre;
    const Vec3 d_neigh = *neighbour_centre - owner_centre;
    const double dist_owner = d_owner.mag();
    const double dist_total = d_neigh.mag();

    if (dist_total < 1e-15) return 0.5;  // Degenerate: fallback to midpoint
    return std::min(1.0, std::max(0.0, dist_owner / dist_total));
}

// Calcule le gradient d'un champ scalaire cellulaire par la méthode Gauss.
//
inline Field<double, Location::CELL> compute_gradient_gauss(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();

    if (cell_field.size() != n_cells) {
        throw std::runtime_error("compute_gradient_gauss: field size != n_cells");
    }
    if (cell_field.dimension() != 1) {
        throw std::runtime_error("compute_gradient_gauss: field must be scalar (dim=1)");
    }

    Field<double, Location::CELL> grad(n_cells, cell_field.name() + "_grad", "1/s", 3);

    // 1. Géométrie des faces : centre + Sf.
    std::vector<Vec3> face_centres(n_faces);
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
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    // 2. Establish global owner -> neighbour orientation, then compute
    // cell-local geometry from the oriented face vectors.
    std::vector<Vec3> provisional_centres(n_cells);
    compute_area_weighted_cell_centres(mesh, face_centres.data(), face_Sf.data(), provisional_centres.data());
    orient_mesh_face_vectors(mesh, face_centres, provisional_centres, face_Sf);

    // 2. Géométrie des cellules : centre + volume.
    std::vector<double> cell_volume(n_cells, 0.0);
    std::vector<Vec3> cell_centre(n_cells);

    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry(mesh, face_centres.data(), face_Sf.data(), cell_faces + off, c, n);
        cell_centre[c] = cg.centre;
        cell_volume[c] = cg.volume;
    }

    // 3. Interpolation cellule → face avec pondération géométrique.
    std::vector<double> face_field_values(n_faces, 0.0);
    const FaceOwnership& own = mesh.ownership();
    const double* cell_values = cell_field.component_data(0);

    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        const bool is_internal = (own.neighbour(f) >= 0);
        const std::size_t neighbour = is_internal
            ? static_cast<std::size_t>(own.neighbour(f))
            : owner;

        const Vec3* neigh_ptr = is_internal ? &cell_centre[neighbour] : nullptr;
        const double w = geometric_interpolation_weight(face_centres[f], cell_centre[owner], neigh_ptr);

        face_field_values[f] = (1.0 - w) * cell_values[owner] + w * cell_values[neighbour];
    }

    // 4. Boucle de calcul du gradient :
    //    ∇φ_c = (1/V_c) Σ_f φ_f Sf_f
    //
    //    Pour chaque face f de la cellule c, la contribution est φ_f * Sf_f
    //    où Sf_f est orienté vers l'extérieur de la cellule.
    //
    //    Pour une face interne :
    //      - si c est le owner, Sf est déjà vers l'extérieur.
    //      - si c est le voisin, Sf est vers l'intérieur → on utilise −Sf.
    //
    //    Pour une face de frontière, c est le owner → Sf vers l'extérieur.

    double* gx = grad.component_data(0);
    double* gy = grad.component_data(1);
    double* gz = grad.component_data(2);
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;

        Vec3 sum;
        for (Offset k = 0; k < n; ++k) {
            const std::size_t f = cell_faces[off + k];
            const double phi_f = face_field_values[f];

            // Déterminer si Sf est orienté vers l'extérieur de cette cellule.
            const bool is_owner = (own.owner(f) == c);
            Vec3 Sf_cell = is_owner ? face_Sf[f] : face_Sf[f] * (-1.0);

            sum = sum + Sf_cell * phi_f;
        }
        const double inv_vol = (cell_volume[c] > 0.0) ? 1.0 / cell_volume[c] : 0.0;
        gx[c] = sum.x * inv_vol;
        gy[c] = sum.y * inv_vol;
        gz[c] = sum.z * inv_vol;
    }

    return grad;
}

}  // namespace core
}  // namespace cfdx