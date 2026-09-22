// M0.7-T01 — Gauss gradient
//
// Spécification CFDX v0.7 §28 :
//   Première implémentation : Gauss linear.
//
// Pour un champ scalaire φ (Field<double, CELL>) :
//   ∇φ_c = (1/V_c) Σ_f φ_f Sf_f
//
//   - φ_f est interpolé géométriquement à partir des cellules owner/neighbour.
//     φ_f = (1-w) * φ_owner + w * φ_neighbour
//     w = |Cf - C_owner| / |C_neighbour - C_owner|  (pondération géométrique)
//   - Pour les faces de frontière, on utilise φ_f = φ_owner (zero-gradient).
//
// Le résultat est un Field<double, CELL> de dimension 3 (le gradient est un vecteur
// 3D par cellule). Le stockage est SoA (x, y, z contigus).

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/mesh/index_types.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx {
namespace core {

// Calcule le poids d'interpolation géométrique pour une face.
// w = |Cf - C_owner| / |C_neighbour - C_owner|
// Pour face frontière : w = 0 (zero gradient)
inline double geometric_interpolation_weight(
    const Vec3& face_centre,
    const Vec3& owner_centre,
    const Vec3* neighbour_centre)
{
    if (!neighbour_centre) return 0.0;  // Boundary face

    const Vec3 d_owner = face_centre - owner_centre;
    const Vec3 d_neigh = *neighbour_centre - owner_centre;
    const double dist_owner = d_owner.mag();
    const double dist_total = d_neigh.mag();

    if (!std::isfinite(dist_total) || dist_total < 1e-15)
        throw std::runtime_error("geometric_interpolation_weight: degenerate cell-centre distance");
    return std::min(1.0, std::max(0.0, dist_owner / dist_total));
}

// Calcule le gradient d'un champ scalaire cellulaire par la méthode Gauss.
//
inline Field<double, Location::CELL> compute_gradient_gauss(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    const GeometryCache& geometry)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();

    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("compute_gradient_gauss: invalid geometry cache");
    if (cell_field.size() != n_cells) {
        throw std::runtime_error("compute_gradient_gauss: field size != n_cells");
    }
    if (cell_field.dimension() != 1) {
        throw std::runtime_error("compute_gradient_gauss: field must be scalar (dim=1)");
    }

    Field<double, Location::CELL> grad(n_cells, cell_field.name() + "_grad", "1/s", 3);

    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* cell_values = cell_field.component_data(0);

    std::vector<double> face_field_values(n_faces, 0.0);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        if (owner >= n_cells)
            throw std::runtime_error("compute_gradient_gauss: owner index out of range");
        const std::int64_t neighbour = own.neighbour(f);
        if (neighbour >= 0) {
            const std::size_t nb = static_cast<std::size_t>(neighbour);
            if (nb >= n_cells)
                throw std::runtime_error("compute_gradient_gauss: neighbour index out of range");
            const double w = geometric_interpolation_weight(
                geometry.face_centres[f], geometry.cell_centres[owner],
                &geometry.cell_centres[nb]);
            face_field_values[f] =
                (1.0 - w) * cell_values[owner] + w * cell_values[nb];
        } else {
            face_field_values[f] = cell_values[owner];
        }
    }

    double* gx = grad.component_data(0);
    double* gy = grad.component_data(1);
    double* gz = grad.component_data(2);
    for (std::size_t c = 0; c < n_cells; ++c) {
        Vec3 sum;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const Vec3 Sf_cell = (own.owner(f) == c)
                ? geometry.face_Sf[f]
                : geometry.face_Sf[f] * (-1.0);
            sum = sum + Sf_cell * face_field_values[f];
        }
        const double volume = geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error("compute_gradient_gauss: non-positive cell volume");
        const double inv_volume = 1.0 / volume;
        gx[c] = sum.x * inv_volume;
        gy[c] = sum.y * inv_volume;
        gz[c] = sum.z * inv_volume;
    }
    return grad;
}

inline Field<double, Location::CELL> compute_gradient_gauss(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return compute_gradient_gauss(cell_field, mesh, geometry);
}

}  // namespace core
}  // namespace cfdx
