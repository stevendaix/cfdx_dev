// M0.7-T02 — Divergence
//
// div(phi)_c = (1/V_c) * sum_f s_cf * phi_f
// where s_cf = +1 for owner and -1 for neighbour.
// Geometry is supplied through GeometryCache so repeated evaluations do not
// rebuild face/cell geometry.

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx {
namespace core {

inline Field<double, Location::CELL> compute_divergence(
    const Field<double, Location::FACE>& face_field,
    const Mesh& mesh,
    const GeometryCache& geometry)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();
    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("compute_divergence: invalid geometry cache");
    if (face_field.size() != n_faces)
        throw std::runtime_error("compute_divergence: field size != n_faces");
    if (face_field.dimension() != 1)
        throw std::runtime_error("compute_divergence: field must be scalar (dim=1)");

    Field<double, Location::CELL> div(
        n_cells, face_field.name() + "_div", "1/s", 1);
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* phi = face_field.component_data(0);
    double* d = div.component_data(0);

    for (std::size_t c = 0; c < n_cells; ++c) {
        double sum = 0.0;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const std::size_t owner = own.owner(f);
            if (owner != c && own.neighbour(f) != static_cast<std::int64_t>(c))
                throw std::runtime_error("compute_divergence: face is not attached to cell");
            sum += (owner == c) ? phi[f] : -phi[f];
        }
        const double volume = geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error("compute_divergence: non-positive cell volume");
        d[c] = sum / volume;
    }
    return div;
}

inline Field<double, Location::CELL> compute_divergence(
    const Field<double, Location::FACE>& face_field,
    const Mesh& mesh)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return compute_divergence(face_field, mesh, geometry);
}

}  // namespace core
}  // namespace cfdx
