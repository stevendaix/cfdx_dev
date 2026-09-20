// M0.7-T04 — Flux
//
// Spécification CFDX v0.7 §27, §30 :
//   flux = U_f · Sf_f
//
//   - U_f est un champ de vitesse aux faces (Field<double, FACE>, dim=3).
//   - Sf_f est le vecteur surface de la face f.
//   - Le résultat est un champ scalaire par face (flux massique).
//
//   Le flux est utilisé par la divergence (§30).

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/mesh/index_types.h"
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

// Calcule le flux massique phi_f = U_f · Sf_f pour chaque face.
//
// Args:
//   face_velocity : champ de vitesse aux faces (Field<double, FACE>, dim=3).
//   mesh         : maillage.
//
// Retourne un Field<double, FACE> de dimension 1 (flux scalaire).
inline Field<double, Location::FACE> compute_flux(
    const Field<double, Location::FACE>& face_velocity,
    const Mesh& mesh)
{
    const std::size_t n_faces = mesh.n_faces();

    if (face_velocity.size() != n_faces) {
        throw std::runtime_error("compute_flux: field size != n_faces");
    }
    if (face_velocity.dimension() != 3) {
        throw std::runtime_error("compute_flux: field must have dimension 3");
    }

    Field<double, Location::FACE> phi(n_faces, face_velocity.name() + "_flux", "m^3/s", 1);

    // Géométrie des faces.
    std::vector<Vec3> face_Sf(n_faces);

    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* verts = mesh.faces().vertices_data();
    const auto* offsets = mesh.faces().offsets_data();

    for (std::size_t f = 0; f < n_faces; ++f) {
        const VertexIndex off = offsets[f];
        const VertexIndex n = offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, verts, off, n);
        face_Sf[f] = fg.Sf;
    }

    const double* ux = face_velocity.component_data(0);
    const double* uy = face_velocity.component_data(1);
    const double* uz = face_velocity.component_data(2);
    double* p = phi.component_data(0);

    for (std::size_t f = 0; f < n_faces; ++f) {
        const Vec3& Sf = face_Sf[f];
        p[f] = ux[f] * Sf.x + uy[f] * Sf.y + uz[f] * Sf.z;
    }

    return phi;
}

}  // namespace core
}  // namespace cfdx
