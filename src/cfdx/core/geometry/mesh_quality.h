// M0.2-T03/T04 — Mesh quality metrics (skewness, non-orthogonality)
//
// Spécification CFDX v0.7 §15, §17, §19 :
//   - skewness
//   - non-orthogonality
//
// La topologie reste la source de vérité (§15).
// Ces métriques sont dérivées de la géométrie (§17).

#pragma once

#include "face_geometry.h"
#include <vector>
#include <cstddef>
#include <cmath>
#include <stdexcept>

namespace cfdx {
namespace core {

struct FaceQuality {
    double skewness = 0.0;
    double non_orthogonality = 0.0;  // en radians
    double non_orthogonality_deg = 0.0;
};

// Calcule le skewness et la non-orthogonalité d'une face.
//
// Non-orthogonalité (angle entre Sf et le vecteur C_cell → C_face) :
//   cos(θ) = (Cf - Cc) · Sf / (|Cf - Cc| * |Sf|)
//   non_orthogonality = acos(cos(θ))   [radians]
//
// Skewness (approximation basée sur la projection) :
//   Le skewness mesure à quel point la normale de la face s'écarte de la
//   direction idéale (Cf - Cc). Pour une face parfaitement orthogonale,
//   skewness = 0.
//   skewness = | (Cf - Cc) - ((Cf - Cc) · n) * n | / |Cf - Cc|
//   où n = Sf / |Sf| est la normale unitaire de la face.
inline FaceQuality compute_face_quality(
    const Vec3& face_centre,
    const Vec3& cell_centre,
    const Vec3& Sf)
{
    FaceQuality q;

    const Vec3 d = face_centre - cell_centre;
    const double d_mag = d.mag();
    const double Sf_mag = Sf.mag();

    if (d_mag <= 0.0 || Sf_mag <= 0.0) {
        return q;  // valeurs par défaut (0)
    }

    // Non-orthogonalité.
    const double cos_theta = d.dot(Sf) / (d_mag * Sf_mag);
    // Éviter les erreurs numériques hors de [-1, 1].
    double clamped = cos_theta;
    if (clamped > 1.0) clamped = 1.0;
    if (clamped < -1.0) clamped = -1.0;
    q.non_orthogonality = std::acos(clamped);
    q.non_orthogonality_deg = q.non_orthogonality * 180.0 / std::acos(-1.0);

    // Skewness : composante de d perpendiculaire à la normale, normalisée.
    const Vec3 n = Sf.normalized();
    const double d_parallel = d.dot(n);
    const Vec3 d_perp = d - n * d_parallel;
    const double d_perp_mag = d_perp.mag();
    q.skewness = d_perp_mag / d_mag;

    return q;
}

}  // namespace core
}  // namespace cfdx