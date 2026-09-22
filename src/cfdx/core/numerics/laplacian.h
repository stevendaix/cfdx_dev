// M0.7-T03 — Laplacian
//
// Spécification CFDX v0.7 §31 :
//   La décomposition générale est :
//     Laplacian = orthogonal contribution + non-orthogonal correction
//
//   Schémas supportés :
//     orthogonal
//     non-orthogonal corrected
//     non-orthogonal limited
//     uncorrected
//
// Première implémentation : orthogonal.
//
// Pour un champ scalaire φ (Field<double, CELL>) :
//   ∇²φ_c = (1/V_c) Σ_f (∇φ)_f · Sf_f
//
//   - (∇φ)_f est interpolé linéairement à partir des gradients aux cellules.
//   - Pour les faces de frontière, on utilise (∇φ)_f = (∇φ)_owner.
//
// Le résultat est un Field<double, CELL> de dimension 1.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/mesh/index_types.h"
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include "cfdx/core/geometry/geometry_cache.h"

namespace cfdx {
namespace core {

enum class LaplacianScheme : std::uint8_t {
    ORTHOGONAL = 0,
    CORRECTED,
    LIMITED,
    UNCORRECTED
};

inline const char* to_string(LaplacianScheme s) {
    switch (s) {
        case LaplacianScheme::ORTHOGONAL:  return "orthogonal";
        case LaplacianScheme::CORRECTED:   return "corrected";
        case LaplacianScheme::LIMITED:     return "limited";
        case LaplacianScheme::UNCORRECTED: return "uncorrected";
        default:                            return "unknown";
    }
}

inline LaplacianScheme laplacian_scheme_from_string(const std::string& s) {
    if (s == "orthogonal")  return LaplacianScheme::ORTHOGONAL;
    if (s == "corrected")   return LaplacianScheme::CORRECTED;
    if (s == "limited")     return LaplacianScheme::LIMITED;
    if (s == "uncorrected") return LaplacianScheme::UNCORRECTED;
    throw std::runtime_error("laplacian_scheme_from_string: unknown scheme '" + s + "'");
}

// Calcule le laplacien d'un champ scalaire cellulaire.
//
// Args:
//   cell_field : champ scalaire (dim=1) sur les cellules.
//   mesh       : maillage.
//   scheme     : schéma de discrétisation.
//
// Retourne un Field<double, CELL> de dimension 1.
inline Field<double, Location::CELL> compute_laplacian(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    const GeometryCache& geometry,
    LaplacianScheme scheme = LaplacianScheme::ORTHOGONAL)
{
    const std::size_t n_cells = mesh.n_cells();
    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("compute_laplacian: invalid geometry cache");
    if (cell_field.size() != n_cells)
        throw std::runtime_error("compute_laplacian: field size != n_cells");
    if (cell_field.dimension() != 1)
        throw std::runtime_error("compute_laplacian: field must be scalar (dim=1)");
    if (scheme == LaplacianScheme::CORRECTED || scheme == LaplacianScheme::LIMITED)
        throw std::runtime_error("compute_laplacian: requested non-orthogonal scheme is not implemented");

    Field<double, Location::CELL> lap(
        n_cells, cell_field.name() + "_lap", cell_field.metadata().unit + "/m^2", 1);
    const auto grad = compute_gradient_gauss(cell_field, mesh, geometry);
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* gx = grad.component_data(0);
    const double* gy = grad.component_data(1);
    const double* gz = grad.component_data(2);

    std::vector<Vec3> face_grad(geometry.face_Sf.size());
    for (std::size_t f = 0; f < geometry.face_Sf.size(); ++f) {
        const std::size_t owner = own.owner(f);
        if (owner >= n_cells)
            throw std::runtime_error("compute_laplacian: owner index out of range");
        const Vec3 g_owner{gx[owner], gy[owner], gz[owner]};
        const std::int64_t neighbour = own.neighbour(f);
        if (neighbour >= 0) {
            const std::size_t nb = static_cast<std::size_t>(neighbour);
            if (nb >= n_cells)
                throw std::runtime_error("compute_laplacian: neighbour index out of range");
            const Vec3 g_neighbour{gx[nb], gy[nb], gz[nb]};
            face_grad[f] = (g_owner + g_neighbour) * 0.5;
        } else {
            face_grad[f] = g_owner;
        }
    }

    double* out = lap.component_data(0);
    for (std::size_t c = 0; c < n_cells; ++c) {
        double sum = 0.0;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            if (own.owner(f) != c && own.neighbour(f) != static_cast<std::int64_t>(c))
                throw std::runtime_error("compute_laplacian: face is not attached to cell");
            const Vec3 Sf_cell = (own.owner(f) == c)
                ? geometry.face_Sf[f]
                : geometry.face_Sf[f] * (-1.0);
            sum += face_grad[f].dot(Sf_cell);
        }
        const double volume = geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error("compute_laplacian: non-positive cell volume");
        out[c] = sum / volume;
    }
    return lap;
}

inline Field<double, Location::CELL> compute_laplacian(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    LaplacianScheme scheme = LaplacianScheme::ORTHOGONAL)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return compute_laplacian(cell_field, mesh, geometry, scheme);
}

}  // namespace core
}  // namespace cfdx
