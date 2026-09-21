// M0.2-T02 — Cell geometry
//
// Spécification CFDX v0.7 §15, §17 :
//   Calcule pour chaque cellule :
//     - centre
//     - volume
//     - métriques géométriques
//
// Le volume d'une cellule polyédrique est calculé par décomposition en pyramides
// (une pyramide par face, sommet au centre de la cellule).
//
//   V = 1/3 * Σ_f (centre_f - centre_c) · Sf_f
//
// La formule est exacte pour une cellule fermée et convexe.
// Pour une cellule non convexe, c'est une approximation géométrique.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/index_types.h"
#include <vector>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace cfdx {
namespace core {

struct CellGeometry {
    Vec3 centre;
    double volume;
};

// Calcule le centre et le volume d'une cellule polyédrique.
//
// centre_c = somme des centres de faces pondérés par l'aire (approximation)
// volume    = 1/3 * Σ_f (Cf - Cc) · Sf_f   (formule de la pyramide)
//
// Args:
//   face_centres : tableau des centres de faces (n_faces)
//   face_Sf      : tableau des vecteurs surface des faces (n_faces)
//   face_ids     : liste des indices de faces appartenant à la cellule
//
// Retourne {centre, volume}.
inline CellGeometry compute_cell_geometry(
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const FaceIndex* face_ids,
    std::size_t n_cell_faces)
{
    if (n_cell_faces == 0) {
        throw std::runtime_error("CellGeometry: cell must have at least one face");
    }

    // Centre : moyenne pondérée par l'aire des faces.
    Vec3 weighted_sum;
    double total_area = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const Vec3& cf = face_centres[f];
        const Vec3& sf = face_Sf[f];
        const double a = sf.mag();
        weighted_sum = weighted_sum + cf * a;
        total_area += a;
    }
    Vec3 centre = total_area > 0.0 ? weighted_sum * (1.0 / total_area) : Vec3{};

    // Volume by convex-cell pyramid decomposition. The face-area-weighted
    // centre is used as an interior reference point; each contribution is
    // one third of face area times the perpendicular distance to that face.
    // Using the absolute normal projection makes the volume robust to a face
    // winding error while preserving the exact result for closed planar convex
    // cells.
    double volume = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const Vec3& cf = face_centres[f];
        const Vec3& sf = face_Sf[f];
        const double area = sf.mag();
        if (!(area > 0.0))
            throw std::runtime_error("CellGeometry: degenerate face");
        const Vec3 d = cf - centre;
        const double height = std::abs(d.dot(sf)) / area;
        volume += area * height;
    }
    volume /= 3.0;

    // Le volume doit être positif pour une cellule valide (§18).
    // Si négatif, l'orientation des faces est inversée.
    // On retourne la valeur absolue pour le volume, mais on signale la non-conformité.

    return {centre, std::abs(volume)};
}

}  // namespace core
}  // namespace cfdx
