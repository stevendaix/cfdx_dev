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
inline CellGeometry compute_cell_geometry(
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const FaceIndex* face_ids,
    std::size_t n_cell_faces)
{
    if (n_cell_faces == 0) {
        throw std::runtime_error("CellGeometry: cell must have at least one face");
    }

    // Detect 2D mesh: all face centres have z ≈ 0 and all Sf vectors
    // have only z-component (x≈0, y≈0). This occurs when 2D elements
    // are represented as cells with a single face (the cell itself).
    bool is_2d = true;
    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        if (std::abs(face_centres[f].z) > 1e-12) { is_2d = false; break; }
        const Vec3& sf = face_Sf[f];
        if (std::abs(sf.x) > 1e-12 || std::abs(sf.y) > 1e-12) { is_2d = false; break; }
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

    double signed_volume = 0.0;

    if (is_2d) {
        // 2D mesh: the pyramid formula gives zero because (Cf - Cc)·Sf = 0
        // (all z=0, Sf only has z-component). Use face areas as cell area.
        // For the current topology, each cell has one face (the cell itself),
        // so the sum of face areas equals the cell area.
        for (std::size_t k = 0; k < n_cell_faces; ++k) {
            const FaceIndex f = face_ids[k];
            const double area = face_Sf[f].mag();
            if (!(area > 0.0))
                throw std::runtime_error("CellGeometry: degenerate face");
            signed_volume += area;
        }
    } else {
        // 3D mesh: signed pyramid decomposition.
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
    }

    return {centre, std::abs(signed_volume), signed_volume};
}


inline CellGeometry compute_cell_geometry_oriented(
    const Vec3* face_centres,
    const Vec3* face_Sf,
    const FaceIndex* face_ids,
    std::size_t n_cell_faces,
    CellIndex cell,
    const FaceOwnership& ownership)
{
    if (n_cell_faces == 0)
        throw std::runtime_error("CellGeometry: cell must have at least one face");

    std::vector<Vec3> local_centres(n_cell_faces);
    std::vector<Vec3> local_Sf(n_cell_faces);
    std::vector<FaceIndex> local_ids(n_cell_faces);

    for (std::size_t k = 0; k < n_cell_faces; ++k) {
        const FaceIndex f = face_ids[k];
        if (f >= ownership.size())
            throw std::runtime_error("CellGeometry: face ownership index out of range");

        const CellIndex owner = ownership.owner(f);
        const std::int64_t neighbour = ownership.neighbour(f);
        if (owner != cell && neighbour != static_cast<std::int64_t>(cell))
            throw std::runtime_error("CellGeometry: face is not attached to requested cell");

        local_centres[k] = face_centres[f];
        local_Sf[k] = (owner == cell) ? face_Sf[f] : face_Sf[f] * (-1.0);
        local_ids[k] = static_cast<FaceIndex>(k);
    }

    return compute_cell_geometry(
        local_centres.data(), local_Sf.data(), local_ids.data(), n_cell_faces);
}

}  // namespace core
}  // namespace cfdx
