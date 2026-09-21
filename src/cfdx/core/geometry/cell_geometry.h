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
#include "cfdx/core/mesh/ownership.h"
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

// Compute cell geometry when the face list uses global owner-oriented vectors.
// Internal-face vectors are reversed for neighbour cells so every contribution
// is expressed with the local outward normal. Boundary faces remain owner-outward.
inline CellGeometry compute_cell_geometry_oriented(
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const std::size_t* face_owners,
    const int* face_neighbours,
    const FaceIndex* face_ids,
    std::size_t cell_id,
    std::size_t n_cell_faces)
{
    if (n_cell_faces == 0)
        throw std::runtime_error("CellGeometry: cell must have at least one face");

    Vec3 weighted_sum;
    double total_area = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        if (face_owners[f] != cell_id &&
            (face_neighbours[f] < 0 || static_cast<std::size_t>(face_neighbours[f]) != cell_id))
            throw std::runtime_error("CellGeometry: cell-face ownership mismatch");
        const Vec3 sf = face_owners[f] == cell_id ? face_Sf[f] : face_Sf[f] * -1.0;
        const double area = sf.mag();
        if (!(area > 0.0) || !std::isfinite(area))
            throw std::runtime_error("CellGeometry: degenerate face");
        weighted_sum = weighted_sum + face_centres[f] * area;
        total_area += area;
    }
    if (!(total_area > 0.0) || !std::isfinite(total_area))
        throw std::runtime_error("CellGeometry: invalid total face area");
    const Vec3 centre = weighted_sum * (1.0 / total_area);

    double signed_volume = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const Vec3 sf = face_owners[f] == cell_id ? face_Sf[f] : face_Sf[f] * -1.0;
        signed_volume += (face_centres[f] - centre).dot(sf);
    }
    signed_volume /= 3.0;
    if (!std::isfinite(signed_volume) || std::abs(signed_volume) <= 0.0)
        throw std::runtime_error("CellGeometry: invalid cell volume");
    return {centre, std::abs(signed_volume), signed_volume};
}

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

    // Signed pyramid decomposition. The sign is an invariant of the face
    // orientation and must not be erased here: validation uses it to reject
    // globally inverted cells. The absolute value is exposed separately as
    // the geometric measure for kernels that require a positive volume.
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


inline CellGeometry compute_cell_geometry(
    const FaceOwnership& ownership,
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const FaceIndex* face_ids,
    std::size_t cell_id,
    std::size_t n_cell_faces)
{
    if (ownership.size() == 0)
        throw std::runtime_error("CellGeometry: empty face ownership");

    Vec3 weighted_sum;
    double total_area = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const std::size_t owner = ownership.owner(f);
        const std::int64_t neighbour = ownership.neighbour(f);
        if (owner != cell_id &&
            (neighbour < 0 || static_cast<std::size_t>(neighbour) != cell_id))
            throw std::runtime_error("CellGeometry: cell-face ownership mismatch");
        const Vec3 sf = owner == cell_id ? face_Sf[f] : face_Sf[f] * -1.0;
        const double area = sf.mag();
        if (!(area > 0.0) || !std::isfinite(area))
            throw std::runtime_error("CellGeometry: degenerate face");
        weighted_sum = weighted_sum + face_centres[f] * area;
        total_area += area;
    }
    if (!(total_area > 0.0) || !std::isfinite(total_area))
        throw std::runtime_error("CellGeometry: invalid total face area");
    const Vec3 centre = weighted_sum * (1.0 / total_area);

    double signed_volume = 0.0;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        const std::size_t owner = ownership.owner(f);
        const Vec3 sf = owner == cell_id ? face_Sf[f] : face_Sf[f] * -1.0;
        signed_volume += (face_centres[f] - centre).dot(sf);
    }
    signed_volume /= 3.0;
    if (!std::isfinite(signed_volume) || signed_volume == 0.0)
        throw std::runtime_error("CellGeometry: invalid cell volume");
    return {centre, std::abs(signed_volume), signed_volume};
}

inline CellGeometry compute_cell_geometry(
    const Mesh& mesh,
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const FaceIndex* face_ids,
    std::size_t cell_id,
    std::size_t n_cell_faces)
{
    return compute_cell_geometry(mesh.ownership(), face_centres, face_Sf,
                                 face_ids, cell_id, n_cell_faces);
}

}  // namespace core
}  // namespace cfdx
