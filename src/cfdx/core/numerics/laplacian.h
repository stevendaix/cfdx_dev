// M0.7-T03 — Laplacian
//
// Spécification CFDX v0.7 §31 :
//   La décomposition générale est :
//     Laplacian = orthogonal contribution + non-orthogonal correction
//
//   Schémas supportés :
//     orthogonal
//     non-orthogonal corrected
//     non-orthogonal limited
//     uncorrected
//
// Première implémentation : orthogonal.
//
// Pour un champ scalaire φ (Field<double, CELL>) :
//   ∇²φ_c = (1/V_c) Σ_f (∇φ)_f · Sf_f
//
//   - (∇φ)_f est interpolé linéairement à partir des gradients aux cellules.
//   - Pour les faces de frontière, on utilise (∇φ)_f = (∇φ)_owner.
//
// Le résultat est un Field<double, CELL> de dimension 1.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/mesh/index_types.h"
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

enum class LaplacianScheme : std::uint8_t {
    ORTHOGONAL = 0,
    CORRECTED,
    LIMITED,
    UNCORRECTED
};

inline const char* to_string(LaplacianScheme s) {
    switch (s) {
        case LaplacianScheme::ORTHOGONAL:  return "orthogonal";
        case LaplacianScheme::CORRECTED:   return "corrected";
        case LaplacianScheme::LIMITED:     return "limited";
        case LaplacianScheme::UNCORRECTED: return "uncorrected";
        default:                            return "unknown";
    }
}

inline LaplacianScheme laplacian_scheme_from_string(const std::string& s) {
    if (s == "orthogonal")  return LaplacianScheme::ORTHOGONAL;
    if (s == "corrected")   return LaplacianScheme::CORRECTED;
    if (s == "limited")     return LaplacianScheme::LIMITED;
    if (s == "uncorrected") return LaplacianScheme::UNCORRECTED;
    throw std::runtime_error("laplacian_scheme_from_string: unknown scheme '" + s + "'");
}

// Calcule le laplacien d'un champ scalaire cellulaire.
//
// Args:
//   cell_field : champ scalaire (dim=1) sur les cellules.
//   mesh       : maillage.
//   scheme     : schéma de discrétisation.
//
// Retourne un Field<double, CELL> de dimension 1.
inline Field<double, Location::CELL> compute_laplacian(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    LaplacianScheme scheme = LaplacianScheme::ORTHOGONAL)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();

    if (cell_field.size() != n_cells) {
        throw std::runtime_error("compute_laplacian: field size != n_cells");
    }
    if (cell_field.dimension() != 1) {
        throw std::runtime_error("compute_laplacian: field must be scalar (dim=1)");
    }

    Field<double, Location::CELL> lap(n_cells, cell_field.name() + "_lap", "1/s^2", 1);

    // 1. Gradient par Gauss (linéaire).
    auto grad = compute_gradient_gauss(cell_field, mesh);

    // 2. Géométrie des faces.
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

    // 3. Establish global owner -> neighbour orientation, then compute
    // cell-local geometry from the oriented face vectors.
    std::vector<Vec3> provisional_centres(n_cells);
    compute_area_weighted_cell_centres(mesh, face_centres.data(), face_Sf.data(), provisional_centres.data());
    orient_mesh_face_vectors(mesh, face_centres, provisional_centres, face_Sf);

    // 3. Géométrie des cellules.
    std::vector<double> cell_volume(n_cells, 0.0);

    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry(mesh, face_centres.data(), face_Sf.data(), cell_faces + off, c, n);
        cell_volume[c] = cg.volume;
    }

// 4. Interpolation du gradient vers les faces (linéaire).
    //    (∇φ)_f = 0.5 * ((∇φ)_owner + (∇φ)_neighbour)
    //    Pour les faces de frontière : (∇φ)_f = (∇φ)_owner.
    std::vector<Vec3> face_grad(n_faces);
    const FaceOwnership& own = mesh.ownership();

    const double* gr_x = grad.component_data(0);
    const double* gr_y = grad.component_data(1);
    const double* gr_z = grad.component_data(2);

    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        const Vec3 g_owner{gr_x[owner], gr_y[owner], gr_z[owner]};

        if (own.neighbour(f) >= 0) {
            const std::size_t neighbour = static_cast<std::size_t>(own.neighbour(f));
            const Vec3 g_neigh{gr_x[neighbour], gr_y[neighbour], gr_z[neighbour]};
            face_grad[f] = (g_owner + g_neigh) * 0.5;
        } else {
            face_grad[f] = g_owner;
        }
    }

    // 5. Boucle de calcul :
    //    ∇²φ_c = (1/V_c) Σ_f (∇φ)_f · Sf_f
    //
    //    Pour une face interne :
    //      - si c est le owner, Sf est vers l'extérieur → contribution = (∇φ)_f · Sf
    //      - si c est le voisin, Sf est vers l'intérieur → contribution = (∇φ)_f · (−Sf)
    //
//    Pour une face de frontière, c est le owner → contribution = (∇φ)_f · Sf
    double* l = lap.component_data(0);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;

        double sum = 0.0;
        for (Offset k = 0; k < n; ++k) {
            const std::size_t f = cell_faces[off + k];
            const bool is_owner = (own.owner(f) == c);
            const Vec3 Sf_cell = is_owner ? face_Sf[f] : face_Sf[f] * (-1.0);
            sum += face_grad[f].dot(Sf_cell);
        }
        const double inv_vol = (cell_volume[c] > 0.0) ? 1.0 / cell_volume[c] : 0.0;
        l[c] = sum * inv_vol;
    }

    (void)scheme;  // Les schémas non-orthogonaux sont réservés à M0.7-T04+.
    return lap;
}

}  // namespace core
}  // namespace cfdx
