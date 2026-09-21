// M0.7-T02 — Divergence
//
// Spécification CFDX v0.7 §30 :
//   div(phi) = (1/V) * Σ_f phi_f
//
//   - phi_f est un champ scalaire sur les faces (flux massique).
//   - Le résultat est un champ scalaire par cellule.
//   - div_c = (1/V_c) * Σ_{f ∈ faces(c)} phi_f
//
//   La contribution d'une face interne est comptée une seule fois (elle est
//   partagée entre owner et neighbour). Le signe est déjà inclus dans phi_f :
//     - si c est l'owner, contribution = +phi_f
//     - si c est le voisin, contribution = −phi_f
//
//   Pour une face de frontière, c est le owner → contribution = +phi_f.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx {
namespace core {

// Calcule la divergence d'un champ de flux par face.
//
// Args:
//   face_field : champ scalaire (dim=1) sur les faces (flux).
//   mesh       : maillage.
//
// Retourne un Field<double, CELL> de dimension 1 (divergence scalaire).
//   div_c = (1/V_c) * Σ_{f ∈ faces(c)} phi_f
inline Field<double, Location::CELL> compute_divergence(
    const Field<double, Location::FACE>& face_field,
    const Mesh& mesh)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();

    if (face_field.size() != n_faces) {
        throw std::runtime_error("compute_divergence: field size != n_faces");
    }
    if (face_field.dimension() != 1) {
        throw std::runtime_error("compute_divergence: field must be scalar (dim=1)");
    }

    Field<double, Location::CELL> div(n_cells, face_field.name() + "_div", "1/s", 1);

    // Compute face geometry (centres, Sf)
    std::vector<Vec3> face_centres(n_faces);
    std::vector<Vec3> face_Sf(n_faces);

    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = mesh.faces().vertices_data();
    const auto* offsets = mesh.faces().offsets_data();

    for (std::size_t f = 0; f < n_faces; ++f) {
        const Offset off = offsets[f];
        const Offset n = offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    // Compute cell volumes
    std::vector<double> cell_volume(n_cells, 0.0);

    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry_oriented(
            face_centres.data(), face_Sf.data(), cell_faces + off, n,
            static_cast<CellIndex>(c), mesh.ownership());
        cell_volume[c] = cg.volume;
    }

    const FaceOwnership& own = mesh.ownership();
    const double* phi = face_field.component_data(0);
    double* d = div.component_data(0);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;

        double sum = 0.0;
        for (Offset k = 0; k < n; ++k) {
            const std::size_t f = cell_faces[off + k];
            const double phi_f = phi[f];

            // Si c est le owner, contribution = +phi_f.
            // Si c est le voisin, contribution = −phi_f.
            const bool is_owner = (own.owner(f) == c);
            sum += is_owner ? phi_f : -phi_f;
        }
        if (!(cell_volume[c] > 0.0) || !std::isfinite(cell_volume[c]))
            throw std::runtime_error("compute_divergence: non-positive cell volume");
        const double inv_vol = 1.0 / cell_volume[c];
        d[c] = sum * inv_vol;
    }

    return div;
}

}  // namespace core
}  // namespace cfdx