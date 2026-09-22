// M0.7-T03 — Laplacian
//
// Current production implementation: orthogonal / uncorrected two-point
// Gauss evaluation. Corrected and limited non-orthogonal schemes are explicit
// API values but are rejected until their quantitative MMS gates are implemented.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/mesh/index_types.h"
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>

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
        default:                           return "unknown";
    }
}

inline LaplacianScheme laplacian_scheme_from_string(const std::string& s) {
    if (s == "orthogonal")  return LaplacianScheme::ORTHOGONAL;
    if (s == "corrected")   return LaplacianScheme::CORRECTED;
    if (s == "limited")     return LaplacianScheme::LIMITED;
    if (s == "uncorrected") return LaplacianScheme::UNCORRECTED;
    throw std::runtime_error("laplacian_scheme_from_string: unknown scheme '" + s + "'");
}

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
    if (scheme == LaplacianScheme::CORRECTED || scheme == LaplacianScheme::LIMITED) {
        throw std::runtime_error(
            "compute_laplacian: scheme '" + std::string(to_string(scheme)) +
            "' is not implemented; use ORTHOGONAL/UNCORRECTED or implement the non-orthogonal MMS gate first");
    }

    Field<double, Location::CELL> lap(
        n_cells, cell_field.name() + "_lap", cell_field.metadata().unit + "/m^2", 1);
    const auto grad = compute_gradient_gauss(cell_field, mesh, geometry);
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* gx = grad.component_data(0);
    const double* gy = grad.component_data(1);
    const double* gz = grad.component_data(2);
    double* out = lap.component_data(0);

    for (std::size_t c = 0; c < n_cells; ++c) {
        double sum = 0.0;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const Vec3 gf{gx[c], gy[c], gz[c]};
            const Vec3 Sf_cell = (own.owner(f) == c)
                ? geometry.face_Sf[f]
                : geometry.face_Sf[f] * (-1.0);
            sum += gf.dot(Sf_cell);
        }
        out[c] = sum / geometry.cell_volumes[c];
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
