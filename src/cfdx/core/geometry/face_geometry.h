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

#include <vector>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace cfdx {
namespace core {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double mag() const { return std::sqrt(x * x + y * y + z * z); }
    double mag2() const { return x * x + y * y + z * z; }

    Vec3 normalized() const {
        double m = mag();
        if (m == 0.0) return {0.0, 0.0, 0.0};
        return {x / m, y / m, z / m};
    }
};

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
    const std::vector<double>& px,
    const std::vector<double>& py,
    const std::vector<double>& pz,
    const std::uint32_t* vertices,
    std::uint32_t offset,
    std::uint32_t n_verts)
{
    if (n_verts < 3) {
        throw std::runtime_error("FaceGeometry: face must have at least 3 vertices");
    }

    // Centre : moyenne arithmétique des sommets (valable pour une face planaire ;
    // pour une face non planaire c'est une approximation acceptable du centre
    // de gravité géométrique).
    Vec3 centre;
    for (std::uint32_t k = 0; k < n_verts; ++k) {
        const std::uint32_t v = vertices[offset + k];
        centre.x += px[v];
        centre.y += py[v];
        centre.z += pz[v];
    }
    centre.x /= n_verts;
    centre.y /= n_verts;
    centre.z /= n_verts;

    // Triangulation en éventail (fan) autour du premier sommet.
    // Aire vectorielle = 1/2 * Σ (P_i × P_{i+1})  (formule de shoelace 3D).
    // On utilise le premier sommet comme origine.
    Vec3 Sf;
    const std::uint32_t v0 = vertices[offset];
    const Vec3 p0{px[v0], py[v0], pz[v0]};
    for (std::uint32_t k = 1; k + 1 < n_verts; ++k) {
        const std::uint32_t va = vertices[offset + k];
        const std::uint32_t vb = vertices[offset + k + 1];
        const Vec3 pa{px[va], py[va], pz[va]};
        const Vec3 pb{px[vb], py[vb], pz[vb]};
        // Triangle (p0, pa, pb) — aire vectorielle = 1/2 * ( (pa-p0) × (pb-p0) )
        const Vec3 ea = pa - p0;
        const Vec3 eb = pb - p0;
        Sf = Sf + ea.cross(eb);
    }
    Sf = Sf * 0.5;

    const double area = Sf.mag();
    Vec3 normal = area > 0.0 ? Sf.normalized() : Vec3{0, 0, 0};

    return {centre, Sf, area, normal};
}

// Calcule la distance (signée) entre le centre d'une face et le centre de la cellule owner.
// Utilisée pour le calcul du delta coefficient (§15).
inline double face_cell_distance(const Vec3& face_centre, const Vec3& cell_centre) {
    return (face_centre - cell_centre).mag();
}

}  // namespace core
}  // namespace cfdx