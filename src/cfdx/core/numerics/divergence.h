// M0.7-T02 — Divergence
//
// Spécification CFDX v0.7 §30 :
//   div(phi) = Σ_f phi_f
//
//   - phi_f est un champ scalaire sur les faces (flux massique).
//   - Le résultat est un champ scalaire par cellule.
//
//   Pour chaque cellule c :
//     div_c = Σ_{f ∈ faces(c)} phi_f
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
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

// Calcule la divergence d'un champ de flux par face.
//
// Args:
//   face_field : champ scalaire (dim=1) sur les faces (flux).
//   mesh       : maillage.
//
// Retourne un Field<double, CELL> de dimension 1 (divergence scalaire).
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

    const FaceOwnership& own = mesh.ownership();
    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    const double* phi = face_field.data();
    double* d = div.data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const std::uint32_t off = cell_offsets[c];
        const std::uint32_t n = cell_offsets[c + 1] - off;

        double sum = 0.0;
        for (std::uint32_t k = 0; k < n; ++k) {
            const std::size_t f = cell_faces[off + k];
            const double phi_f = phi[f];

            // Si c est le owner, contribution = +phi_f.
            // Si c est le voisin, contribution = −phi_f.
            const bool is_owner = (own.owner(f) == c);
            sum += is_owner ? phi_f : -phi_f;
        }
        d[c] = sum;
    }

    return div;
}

}  // namespace core
}  // namespace cfdx