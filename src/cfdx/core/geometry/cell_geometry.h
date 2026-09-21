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
    // Signed volume before taking the absolute value. A negative value is an
    // inverted cell and must be rejected by the mesh validator.
    double signed_volume;
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
inline CellGeometry compute_cell_geometry_raw(
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const FaceIndex* face_ids,
    std::size_t n_cell_faces)
{
    if (n_cell_faces == 0) {
        throw std::runtime_error("CellGeometry: cell must have at least one face");
    }

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

    double signed_volume = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const Vec3& cf = face_centres[f];
        const Vec3& sf = face_Sf[f];
        const double area = sf.mag();
        if (!(area > 0.0))
            throw std::runtime_error("CellGeometry: degenerate face");
        signed_volume += (cf - centre).dot(sf);
    }
    signed_volume /= 3.0;
    return {centre, std::abs(signed_volume), signed_volume};
}

// Compute geometric cell measures from face data whose normals use the global
// mesh face orientation. A shared face has one owner-oriented normal, so it
// must be flipped for the neighbouring cell. The local orientation is inferred
// from the face-centre to cell-centre direction and is therefore independent of
// the global owner/neighbour convention.
inline CellGeometry compute_cell_geometry(
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const FaceIndex* face_ids,
    std::size_t n_cell_faces)
{
    if (n_cell_faces == 0) {
        throw std::runtime_error("CellGeometry: cell must have at least one face");
    }

    Vec3 weighted_sum;
    double total_area = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const Vec3& cf = face_centres[f];
        const Vec3& sf = face_Sf[f];
        const double a = sf.mag();
        if (!(a > 0.0))
            throw std::runtime_error("CellGeometry: degenerate face");
        weighted_sum = weighted_sum + cf * a;
        total_area += a;
    }
    if (!(total_area > 0.0))
        throw std::runtime_error("CellGeometry: cell has zero total face area");

    const Vec3 centre = weighted_sum * (1.0 / total_area);

    double volume = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const Vec3& cf = face_centres[f];
        const Vec3& sf_global = face_Sf[f];
        const double area = sf_global.mag();
        Vec3 sf = sf_global;
        if ((cf - centre).dot(sf) < 0.0)
            sf = sf * -1.0;
        volume += (cf - centre).dot(sf);
    }
    volume /= 3.0;

    if (!(volume > 0.0) || !std::isfinite(volume))
        throw std::runtime_error("CellGeometry: unable to construct positive cell volume");

    // This value is the locally outward-oriented geometric volume. Validation
    // of mesh winding must use compute_cell_geometry_raw().
    return {centre, volume, volume};
}

}  // namespace core
}  // namespace cfdx
