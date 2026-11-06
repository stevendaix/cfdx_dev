// M0.2-T01 — Face geometry
//
// Spécification CFDX v0.7 §15, §17 :
//   La géométrie est dérivée de la topologie.
//   /mesh/geometry/
//     face_centres, face_area_vectors, face_areas,
//     delta_coeffs, non_orthogonal_correction,
//     skewness, non_orthogonality
//
//   La topologie reste la source de vérité (§15).
//
// Pour une face, le moteur géométrique calcule (§17) :
//   - centre
//   - aire
//   - vecteur surface Sf
//   - normale
//   - triangulation d'intégration si nécessaire

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/index_types.h"
#include <vector>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cfdx {
namespace core {

// Calcule le centre et le vecteur surface d'une face à partir de ses sommets.
// Supporte les faces non planaires par triangulation en éventail (fan triangulation)
// autour du premier sommet (§16).
//
// Pour une face planaire, le résultat est exact.
// Pour une face non planaire, la triangulation est une vue d'intégration géométrique
// et non une modification de la topologie (§16).
struct FaceGeometry {
    Vec3 centre;
    Vec3 Sf;       // vecteur surface (normale * aire)
    double area;   // aire scalaire
    Vec3 normal;   // normale unitaire (Sf / area)
    double triangulated_area = 0.0; // sum of physical triangle areas
    double warpage_ratio = 0.0; // 1 - |integrated Sf| / triangulated_area
};

// Calcule la géométrie d'une face à partir de ses sommets (positions).
// Les sommets doivent être dans l'ordre géométrique (consécutifs le long de la face).
inline FaceGeometry compute_face_geometry(
    const double* px,
    const double* py,
    const double* pz,
    const VertexIndex* vertices,
    VertexIndex offset,
    VertexIndex n_verts)
{
    if (n_verts < 3) {
        throw std::runtime_error("FaceGeometry: face must have at least 3 vertices");
    }

    // Integrate the face polygon through a triangle fan.  The integrated
    // vector Sf is the conservative FVM surface vector; its magnitude may be
    // smaller than the true triangulated area on a warped face.
    const VertexIndex v0 = vertices[offset];
    const Vec3 p0{px[v0], py[v0], pz[v0]};
    Vec3 Sf{}, area_weighted_centroid{};
    double triangulated_area = 0.0;
    for (VertexIndex k = 1; k + 1 < n_verts; ++k) {
        const VertexIndex va = vertices[offset + k];
        const VertexIndex vb = vertices[offset + k + 1];
        const Vec3 pa{px[va], py[va], pz[va]};
        const Vec3 pb{px[vb], py[vb], pz[vb]};
        const Vec3 triangle_Sf = (pa - p0).cross(pb - p0) * 0.5;
        const double triangle_area = triangle_Sf.mag();
        if (!std::isfinite(triangle_area))
            throw std::runtime_error("FaceGeometry: non-finite triangle area");
        Sf += triangle_Sf;
        area_weighted_centroid += (p0 + pa + pb) * (triangle_area / 3.0);
        triangulated_area += triangle_area;
    }
    const double area = Sf.mag();
    if (!(triangulated_area > 0.0) || !(area > 0.0) || !std::isfinite(area))
        throw std::runtime_error("FaceGeometry: degenerate or cancelling face fan");
    const Vec3 centre = area_weighted_centroid * (1.0 / triangulated_area);
    const Vec3 normal = Sf * (1.0 / area);
    const double warpage_ratio = std::max(0.0,
        1.0 - std::min(1.0, area / triangulated_area));
    return {centre, Sf, area, normal, triangulated_area, warpage_ratio};
}

// Surcharge de commodité : accepte des std::vector (API utilisée par les tests).
inline FaceGeometry compute_face_geometry(
    const std::vector<double>& px,
    const std::vector<double>& py,
    const std::vector<double>& pz,
    const VertexIndex* vertices,
    VertexIndex offset,
    VertexIndex n_verts)
{
    return compute_face_geometry(px.data(), py.data(), pz.data(), vertices, offset, n_verts);
}

// Calcule la distance (signée) entre le centre d'une face et le centre de la cellule owner.
// Utilisée pour le calcul du delta coefficient (§15).
inline double face_cell_distance(const Vec3& face_centre, const Vec3& cell_centre) {
    return (face_centre - cell_centre).mag();
}

// --- Face orientation utilities (§17) ---

// Vérifie et corrige l'orientation du vecteur surface Sf pour une face.
// Invariant : Sf · (C_neighbour - C_owner) > 0 pour face interne
//             Sf · (Cf - C_owner) > 0 pour face frontière
inline void ensure_face_orientation(Vec3& Sf, const Vec3& face_centre,
                                     const Vec3& owner_centre,
                                     const Vec3* neighbour_centre = nullptr) {
    Vec3 d;
    if (neighbour_centre) {
        d = *neighbour_centre - owner_centre;
    } else {
        d = face_centre - owner_centre;
    }
    if (Sf.dot(d) < 0.0) {
        Sf = Sf * (-1.0);
    }
}

// Calcule la géométrie d'une face et assure l'orientation correcte du vecteur surface.
// Nécessite les centres des cellules owner/neighbour pour déterminer l'orientation.
inline FaceGeometry compute_face_geometry_oriented(
    const double* px,
    const double* py,
    const double* pz,
    const VertexIndex* vertices,
    VertexIndex offset,
    VertexIndex n_verts,
    const Vec3& owner_centre,
    const Vec3* neighbour_centre = nullptr)
{
    FaceGeometry fg = compute_face_geometry(px, py, pz, vertices, offset, n_verts);
    ensure_face_orientation(fg.Sf, fg.centre, owner_centre, neighbour_centre);
    fg.normal = fg.area > 0.0 ? fg.Sf.normalized() : Vec3{0, 0, 0};
    return fg;
}

}  // namespace core
}  // namespace cfdx
