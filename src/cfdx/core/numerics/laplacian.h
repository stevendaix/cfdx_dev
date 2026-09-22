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
// Pour un champ scalaire φ (Field<double, CELL>) :
//   For ORTHOGONAL/UNCORRECTED:
//     ∇²φ_c = (1/V_c) Σ_internal_f (|Sf|/|d|) (φ_N - φ_P)
//   where d is the owner-neighbour centre distance.
//   Boundary faces have zero contribution because no BoundaryField value is supplied.
//
//   For CORRECTED:
//     Sf = Sf_orth + Sf_corr
//     Sf_orth = d (Sf·d)/|d|²
//     F_f = (Sf·d)/|d|² (φ_N-φ_P) + Sf_corr·grad(phi)_f
//   where grad(phi)_f is the arithmetic mean of the cached Gauss gradients
//   in the two adjacent cells.
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
    if (scheme == LaplacianScheme::LIMITED)
        throw std::runtime_error("compute_laplacian: limited non-orthogonal scheme is not implemented");

    Field<double, Location::CELL> lap(
        n_cells, cell_field.name() + "_lap", cell_field.metadata().unit + "/m^2", 1);
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* phi = cell_field.component_data(0);
    double* out = lap.component_data(0);
    Field<double, Location::CELL> gradients;
    if (scheme == LaplacianScheme::CORRECTED)
        gradients = compute_gradient_gauss(cell_field, mesh, geometry);

    for (std::size_t c = 0; c < n_cells; ++c) {
        double sum = 0.0;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const std::size_t owner = own.owner(f);
            if (owner >= n_cells)
                throw std::runtime_error("compute_laplacian: owner index out of range");
            const std::int64_t neighbour = own.neighbour(f);
            if (owner != c && neighbour != static_cast<std::int64_t>(c))
                throw std::runtime_error("compute_laplacian: face is not attached to cell");

            // Boundary faces have zero normal gradient because this Module-0
            // operator has no BoundaryField value.
            if (neighbour < 0)
                continue;

            const std::size_t nb = static_cast<std::size_t>(neighbour);
            if (nb >= n_cells)
                throw std::runtime_error("compute_laplacian: neighbour index out of range");

            const Vec3 dvec = geometry.cell_centres[nb] - geometry.cell_centres[owner];
            const Vec3 Sf = geometry.face_Sf[f];
            const double d2 = dvec.dot(dvec);
            const double d = std::sqrt(d2);
            const double area = Sf.mag();
            if (!(d > 1e-14) || !(area > 0.0) ||
                !std::isfinite(d) || !std::isfinite(area))
                throw std::runtime_error("compute_laplacian: invalid internal-face geometry");

            const double orth_dot = Sf.dot(dvec);
            const Vec3 Sf_orth = dvec * (orth_dot / d2);
            const double orth_conductance = orth_dot / d2;

            double contribution = orth_conductance * (phi[nb] - phi[owner]);

            if (scheme == LaplacianScheme::CORRECTED) {
                const Vec3 Sf_corr = Sf - Sf_orth;
                const Vec3 grad_face{
                    0.5 * (gradients(owner, 0) + gradients(nb, 0)),
                    0.5 * (gradients(owner, 1) + gradients(nb, 1)),
                    0.5 * (gradients(owner, 2) + gradients(nb, 2))
                };
                contribution += Sf_corr.dot(grad_face);
            }

            sum += (owner == c) ? contribution : -contribution;
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
