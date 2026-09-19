// M0.7-T01 — Gauss gradient
//
// Spécification CFDX v0.7 §28 :
//   Première implémentation : Gauss linear.
//
// Pour un champ scalaire φ (Field<double, CELL>) :
//   ∇φ_c = (1/V_c) Σ_f φ_f Sf_f
//
//   - φ_f est interpolé linéairement à partir des cellules owner/neighbour.
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
#include "cfdx/core/numerics/interpolation.h"
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

// Calcule le gradient d'un champ scalaire cellulaire par la méthode Gauss.
//
// Args:
//   cell_field : champ scalaire (dim=1) sur les cellules.
//   mesh       : maillage (doit contenir la géométrie des faces — calculée à partir des points).
//
// Retourne un Field<double, CELL> de dimension 3 : ∇φ = (∂φ/∂x, ∂φ/∂y, ∂φ/∂z).
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

    // 1. Interpolation cellule → face (linéaire).
    auto face_field = interpolate_cell_to_face(cell_field, mesh, InterpScheme::LINEAR);

    // 2. Géométrie des faces : centre + Sf.
    std::vector<Vec3> face_centres(n_faces);
    std::vector<Vec3> face_Sf(n_faces);

    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = mesh.faces().vertices_data();
    const auto* offsets = mesh.faces().offsets_data();

    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::uint32_t off = offsets[f];
        const std::uint32_t n = offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    // 3. Géométrie des cellules : volume.
    std::vector<double> cell_volume(n_cells, 0.0);
    std::vector<Vec3> cell_centre(n_cells);

    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::uint32_t off = cell_offsets[c];
        const std::uint32_t n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry(
            face_centres.data(), face_Sf.data(), cell_faces + off, n);
        cell_centre[c] = cg.centre;
        cell_volume[c] = cg.volume;
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
    const FaceOwnership& own = mesh.ownership();

    double* g = grad.data();
    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::uint32_t off = cell_offsets[c];
        const std::uint32_t n = cell_offsets[c + 1] - off;

        Vec3 sum;
        for (std::uint32_t k = 0; k < n; ++k) {
            const std::size_t f = cell_faces[off + k];
            const double phi_f = face_field(f);

            // Déterminer si Sf est orienté vers l'extérieur de cette cellule.
            const bool is_owner = (own.owner(f) == c);
            Vec3 Sf_cell = is_owner ? face_Sf[f] : face_Sf[f] * (-1.0);

            sum = sum + Sf_cell * phi_f;
        }
        const double inv_vol = (cell_volume[c] > 0.0) ? 1.0 / cell_volume[c] : 0.0;
        const std::size_t base = c * 3;
        g[base + 0] = sum.x * inv_vol;
        g[base + 1] = sum.y * inv_vol;
        g[base + 2] = sum.z * inv_vol;
    }

    return grad;
}

}  // namespace core
}  // namespace cfdx