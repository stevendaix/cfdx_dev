#pragma once

#include "cfdx/core/numerics/gradient.h"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::core {

// Gauss gradient with explicit Dirichlet values on boundary faces.
//
// boundary_face_values must have one entry per mesh face. Values on internal
// faces are ignored; values on boundary faces must be finite and represent the
// prescribed/reconstructed face value. Unlike compute_gradient_gauss(), this
// routine never silently substitutes the owner-cell value at a boundary.
//
// This API is deliberately Dirichlet-only. A prescribed Neumann condition
// cannot in general be represented by a face value without a reconstruction;
// callers must construct a consistent face value or use a dedicated Neumann
// gradient formulation.
inline Field<double, Location::CELL> compute_gradient_gauss_boundary_values(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    const GeometryCache& geometry,
    const std::vector<double>& boundary_face_values)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();
    if (!is_valid(geometry, mesh))
        throw std::invalid_argument(
            "compute_gradient_gauss_boundary_values: invalid geometry cache");
    if (cell_field.size() != n_cells || cell_field.dimension() != 1)
        throw std::invalid_argument(
            "compute_gradient_gauss_boundary_values: field must be scalar and match mesh cells");
    if (boundary_face_values.size() != n_faces)
        throw std::invalid_argument(
            "compute_gradient_gauss_boundary_values: boundary value array must match mesh faces");

    const auto& ownership = mesh.ownership();
    const double* values = cell_field.component_data(0);
    for (std::size_t c = 0; c < n_cells; ++c) {
        if (!std::isfinite(values[c]))
            throw std::invalid_argument(
                "compute_gradient_gauss_boundary_values: non-finite cell value");
    }

    std::vector<double> face_values(n_faces, 0.0);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = ownership.owner(f);
        if (owner >= n_cells)
            throw std::runtime_error(
                "compute_gradient_gauss_boundary_values: owner index out of range");
        const std::int64_t neighbour = ownership.neighbour(f);
        if (neighbour >= 0) {
            const auto nb = static_cast<std::size_t>(neighbour);
            if (nb >= n_cells)
                throw std::runtime_error(
                    "compute_gradient_gauss_boundary_values: neighbour index out of range");
            const double weight = geometric_interpolation_weight(
                geometry.face_centres[f], geometry.cell_centres[owner],
                &geometry.cell_centres[nb]);
            face_values[f] = (1.0 - weight) * values[owner] + weight * values[nb];
        } else {
            const double boundary_value = boundary_face_values[f];
            if (!std::isfinite(boundary_value))
                throw std::invalid_argument(
                    "compute_gradient_gauss_boundary_values: missing/non-finite Dirichlet value on boundary face "
                    + std::to_string(f));
            face_values[f] = boundary_value;
        }
    }

    Field<double, Location::CELL> gradient(
        n_cells, cell_field.name() + "_grad_boundary",
        cell_field.metadata().unit + "/m", 3);
    double* gx = gradient.component_data(0);
    double* gy = gradient.component_data(1);
    double* gz = gradient.component_data(2);
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    for (std::size_t c = 0; c < n_cells; ++c) {
        Vec3 sum;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const Vec3 sf = ownership.owner(f) == c
                ? geometry.face_Sf[f] : geometry.face_Sf[f] * (-1.0);
            sum = sum + sf * face_values[f];
        }
        const double volume = geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error(
                "compute_gradient_gauss_boundary_values: non-positive cell volume");
        gx[c] = sum.x / volume;
        gy[c] = sum.y / volume;
        gz[c] = sum.z / volume;
        if (!std::isfinite(gx[c]) || !std::isfinite(gy[c]) || !std::isfinite(gz[c]))
            throw std::runtime_error(
                "compute_gradient_gauss_boundary_values: non-finite gradient");
    }
    return gradient;
}

} // namespace cfdx::core
