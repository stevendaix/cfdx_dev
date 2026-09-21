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
#include "cfdx/core/mesh/mesh.h"
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

    // Centre : moyenne arithmétique des sommets (valable pour une face planaire ;
    // pour une face non planaire c'est une approximation acceptable du centre
    // de gravité géométrique).
    Vec3 centre;
    for (VertexIndex k = 0; k < n_verts; ++k) {
        const VertexIndex v = vertices[offset + k];
        centre.x += px[v];
        centre.y += py[v];
        centre.z += pz[v];
    }
    centre.x /= static_cast<double>(n_verts);
    centre.y /= static_cast<double>(n_verts);
    centre.z /= static_cast<double>(n_verts);

    // Triangulation en éventail (fan) autour du premier sommet.
    Vec3 Sf;
    const VertexIndex v0 = vertices[offset];
    const Vec3 p0{px[v0], py[v0], pz[v0]};
    for (VertexIndex k = 1; k + 1 < n_verts; ++k) {
        const VertexIndex va = vertices[offset + k];
        const VertexIndex vb = vertices[offset + k + 1];
        const Vec3 pa{px[va], py[va], pz[va]};
        const Vec3 pb{px[vb], py[vb], pz[vb]};
        const Vec3 ea = pa - p0;
        const Vec3 eb = pb - p0;
        Sf = Sf + ea.cross(eb);
    }
    Sf = Sf * 0.5;

    const double area = Sf.mag();
    Vec3 normal = area > 0.0 ? Sf.normalized() : Vec3{0, 0, 0};

    return {centre, Sf, area, normal};
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


inline void orient_mesh_face_vectors(
    const Mesh& mesh,
    const std::vector<Vec3>& face_centres,
    const std::vector<Vec3>& cell_centres,
    std::vector<Vec3>& face_Sf)
{
    if (face_centres.size() != mesh.n_faces() ||
        face_Sf.size() != mesh.n_faces() ||
        cell_centres.size() != mesh.n_cells())
        throw std::invalid_argument("orient_mesh_face_vectors: geometry size mismatch");

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t owner = mesh.ownership().owner(f);
        if (owner >= mesh.n_cells())
            throw std::runtime_error("orient_mesh_face_vectors: invalid owner");
        const int neighbour = mesh.ownership().neighbour(f);
        if (neighbour >= 0 && static_cast<std::size_t>(neighbour) >= mesh.n_cells())
            throw std::runtime_error("orient_mesh_face_vectors: invalid neighbour");
        const Vec3 d = neighbour >= 0
            ? cell_centres[static_cast<std::size_t>(neighbour)] - cell_centres[owner]
            : face_centres[f] - cell_centres[owner];
        const double projection = face_Sf[f].dot(d);
        if (projection < 0.0) {
            face_Sf[f] = face_Sf[f] * -1.0;
        } else if (!(projection > 0.0) || !std::isfinite(projection)) {
            throw std::runtime_error("orient_mesh_face_vectors: face orientation is undefined");
        }
    }
}

}  // namespace core
}  // namespace cfdx
