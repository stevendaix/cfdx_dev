// M0.2-T03/T04 — Mesh quality metrics (skewness, non-orthogonality)
//
// Spécification CFDX v0.7 §15, §17, §19 :
//   - skewness
//   - non-orthogonality
//
// La topologie reste la source de vérité (§15).
// Ces métriques sont dérivées de la géométrie (§17).

#pragma once

#include "cfdx/core/field/field.h"
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

// Calcule la qualité d'une face frontière à partir du segment cellule -> face.
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
        return q;
    }

    const double cos_theta = d.dot(Sf) / (d_mag * Sf_mag);
    const double clamped = std::max(-1.0, std::min(1.0, cos_theta));
    q.non_orthogonality = std::acos(clamped);
    q.non_orthogonality_deg = q.non_orthogonality * 180.0 / std::acos(-1.0);

    const Vec3 n = Sf.normalized();
    const double d_parallel = d.dot(n);
    const Vec3 d_perp = d - n * d_parallel;
    q.skewness = d_perp.mag() / d_mag;
    return q;
}

// Calcule la qualité d'une face interne à partir de la ligne des centres
// owner -> neighbour. C'est cette direction qui définit l'orthogonalité FV
// interne; utiliser uniquement owner -> face surestime ou masque la
// non-orthogonalité sur un maillage skewed.
inline FaceQuality compute_face_quality(
    const Vec3& face_centre,
    const Vec3& owner_centre,
    const Vec3& neighbour_centre,
    const Vec3& Sf)
{
    FaceQuality q;

    const Vec3 d = neighbour_centre - owner_centre;
    const double d_mag = d.mag();
    const double Sf_mag = Sf.mag();
    if (d_mag <= 0.0 || Sf_mag <= 0.0) {
        return q;
    }

    const double cos_theta = d.dot(Sf) / (d_mag * Sf_mag);
    const double clamped = std::max(-1.0, std::min(1.0, cos_theta));
    q.non_orthogonality = std::acos(clamped);
    q.non_orthogonality_deg = q.non_orthogonality * 180.0 / std::acos(-1.0);

    // Skewness: distance between the face centre and the projection of the
    // face centre onto the owner-neighbour line, normalized by |d|.
    const Vec3 owner_to_face = face_centre - owner_centre;
    const double projection = owner_to_face.dot(d) / (d_mag * d_mag);
    const Vec3 projected = owner_centre + d * projection;
    q.skewness = (face_centre - projected).mag() / d_mag;
    return q;
}

}  // namespace core
}  // namespace cfdx
