#pragma once

#include "cfdx/core/numerics/gradient.h"
#include <array>
#include <cmath>
#include <map>
#include <map>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Return the SST strain-rate magnitude |S| = sqrt(2 S_ij S_ij) from the
// nine components of grad(U), ordered as du_i/dx_j with j varying fastest.
inline double sst_strain_rate_magnitude(const std::array<double, 9>& grad_u)
{
    for (const double value : grad_u) {
        if (!std::isfinite(value))
            throw std::invalid_argument("sst_strain_rate_magnitude: non-finite velocity gradient");
    }

    const double sxx = grad_u[0];
    const double syy = grad_u[4];
    const double szz = grad_u[8];
    const double sxy = 0.5 * (grad_u[1] + grad_u[3]);
    const double sxz = 0.5 * (grad_u[2] + grad_u[6]);
    const double syz = 0.5 * (grad_u[5] + grad_u[7]);
    const double norm_squared = 2.0 * (
        sxx * sxx + syy * syy + szz * szz +
        2.0 * (sxy * sxy + sxz * sxz + syz * syz));
    const double magnitude = std::sqrt(norm_squared);
    if (!std::isfinite(magnitude))
        throw std::invalid_argument("sst_strain_rate_magnitude: strain-rate magnitude overflow");
    return magnitude;
}

// Compute the cell-centred SST strain-rate magnitude from the production
// velocity field using the baseline Gauss gradient operator. The gradient
// operator's existing boundary treatment is preserved, not replaced here.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
compute_sst_strain_rate(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Mesh& mesh)
{
    using namespace cfdx::core;
    if (velocity.dimension() != 3 || velocity.size() != mesh.n_cells())
        throw std::invalid_argument("compute_sst_strain_rate: velocity must be a 3-component cell field");

    const GeometryCache geometry = make_geometry_cache(mesh);
    Field<double, Location::CELL> u(mesh.n_cells(), "u", "m/s", 1);
    Field<double, Location::CELL> v(mesh.n_cells(), "v", "m/s", 1);
    Field<double, Location::CELL> w(mesh.n_cells(), "w", "m/s", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        for (std::size_t component = 0; component < 3; ++component) {
            if (!std::isfinite(velocity.component_data(component)[cell]))
                throw std::invalid_argument("compute_sst_strain_rate: non-finite velocity");
        }
        u(cell) = velocity.component_data(0)[cell];
        v(cell) = velocity.component_data(1)[cell];
        w(cell) = velocity.component_data(2)[cell];
    }
    const auto grad_u = compute_gradient_gauss(u, mesh, geometry);
    const auto grad_v = compute_gradient_gauss(v, mesh, geometry);
    const auto grad_w = compute_gradient_gauss(w, mesh, geometry);

    Field<double, Location::CELL> strain(
        mesh.n_cells(), "sst_strain_rate", "1/s", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const std::array<double, 9> grad{
            grad_u.component_data(0)[cell], grad_u.component_data(1)[cell],
            grad_u.component_data(2)[cell], grad_v.component_data(0)[cell],
            grad_v.component_data(1)[cell], grad_v.component_data(2)[cell],
            grad_w.component_data(0)[cell], grad_w.component_data(1)[cell],
            grad_w.component_data(2)[cell]};
        strain(cell) = sst_strain_rate_magnitude(grad);
    }
    return strain;
}


// Compute SST strain rate with explicit fixed-value velocity boundaries. The
// map is keyed by mesh patch name and contains the prescribed (Ux, Uy, Uz).
// Patches absent from the map use zero-gradient extrapolation, matching the
// production velocity-boundary contract; EMPTY extrusion faces are excluded.
// This overload is intended for production coupling, where boundary values must
// influence the strain invariant used by the SST closure.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
compute_sst_strain_rate_with_boundary_values(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Mesh& mesh,
    const std::map<std::string, std::array<double, 3>>& fixed_boundary_values)
{
    using namespace cfdx::core;
    if (velocity.dimension() != 3 || velocity.size() != mesh.n_cells())
        throw std::invalid_argument(
            "compute_sst_strain_rate_with_boundary_values: velocity must be a 3-component cell field");

    for (const auto& [patch_name, values] : fixed_boundary_values) {
        (void)patch_name;
        for (const double value : values) {
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "compute_sst_strain_rate_with_boundary_values: non-finite boundary velocity");
        }
    }
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        for (std::size_t component = 0; component < 3; ++component) {
            if (!std::isfinite(velocity.component_data(component)[cell]))
                throw std::invalid_argument(
                    "compute_sst_strain_rate_with_boundary_values: non-finite velocity");
        }
    }

    const GeometryCache geometry = make_geometry_cache(mesh);
    const std::size_t n_faces = mesh.n_faces();
    const std::size_t n_patches = mesh.boundary().n_patches();
    std::vector<std::size_t> face_patch(n_faces, n_patches);
    std::map<std::string, PatchType> patch_types;
    for (std::size_t patch = 0; patch < n_patches; ++patch) {
        const auto& boundary_patch = mesh.boundary().patch(patch);
        if (!patch_types.emplace(boundary_patch.name, boundary_patch.type).second)
            throw std::invalid_argument(
                "compute_sst_strain_rate_with_boundary_values: duplicate boundary patch name");
        for (const auto face : boundary_patch.face_ids) {
            if (face >= n_faces)
                throw std::invalid_argument(
                    "compute_sst_strain_rate_with_boundary_values: patch face out of range");
            if (face_patch[face] != n_patches)
                throw std::invalid_argument(
                    "compute_sst_strain_rate_with_boundary_values: face belongs to multiple patches");
            face_patch[face] = patch;
        }
    }
    for (const auto& [patch_name, values] : fixed_boundary_values) {
        (void)values;
        const auto found = patch_types.find(patch_name);
        if (found == patch_types.end())
            throw std::invalid_argument(
                "compute_sst_strain_rate_with_boundary_values: fixed-value boundary references an unknown patch");
        if (found->second == PatchType::EMPTY)
            throw std::invalid_argument(
                "compute_sst_strain_rate_with_boundary_values: EMPTY patches cannot prescribe velocity values");
    }

    const auto gradient_component = [&](std::size_t component) {
        Field<double, Location::CELL> scalar(
            mesh.n_cells(), "velocity_component", "m/s", 1);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
            scalar(cell) = velocity.component_data(component)[cell];

        std::vector<double> face_values(n_faces, 0.0);
        for (std::size_t face = 0; face < n_faces; ++face) {
            const std::size_t owner = mesh.ownership().owner(face);
            if (owner >= mesh.n_cells())
                throw std::invalid_argument(
                    "compute_sst_strain_rate_with_boundary_values: face owner out of range");
            const auto neighbour_raw = mesh.ownership().neighbour(face);
            if (neighbour_raw >= 0) {
                const std::size_t neighbour = static_cast<std::size_t>(neighbour_raw);
                if (neighbour >= mesh.n_cells())
                    throw std::invalid_argument(
                        "compute_sst_strain_rate_with_boundary_values: face neighbour out of range");
                const double weight = geometric_interpolation_weight(
                    geometry.face_centres[face], geometry.cell_centres[owner],
                    &geometry.cell_centres[neighbour]);
                face_values[face] = (1.0 - weight) * scalar(owner) +
                                    weight * scalar(neighbour);
                continue;
            }

            const std::size_t patch = face_patch[face];
            if (patch >= n_patches)
                throw std::invalid_argument(
                    "compute_sst_strain_rate_with_boundary_values: boundary face has no patch");
            if (mesh.boundary().patch(patch).type == PatchType::EMPTY) {
                face_values[face] = scalar(owner);
                continue;
            }
            const auto found = fixed_boundary_values.find(
                mesh.boundary().patch(patch).name);
            face_values[face] = found == fixed_boundary_values.end()
                ? scalar(owner)
                : found->second[component];
        }

        Field<double, Location::CELL> gradient(
            mesh.n_cells(), "velocity_component_gradient", "1/s", 3);
        const auto* cell_faces = mesh.cells().faces_data();
        const auto* cell_offsets = mesh.cells().offsets_data();
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
            Vec3 sum{0.0, 0.0, 0.0};
            for (Offset offset = cell_offsets[cell];
                 offset < cell_offsets[cell + 1]; ++offset) {
                const std::size_t face = cell_faces[offset];
                const std::size_t patch = face_patch[face];
                if (mesh.ownership().neighbour(face) < 0 &&
                    patch < n_patches &&
                    mesh.boundary().patch(patch).type == PatchType::EMPTY)
                    continue;
                const Vec3 area = mesh.ownership().owner(face) == cell
                    ? geometry.face_Sf[face]
                    : geometry.face_Sf[face] * -1.0;
                sum += area * face_values[face];
            }
            const double volume = geometry.cell_volumes[cell];
            if (!(volume > 0.0) || !std::isfinite(volume))
                throw std::invalid_argument(
                    "compute_sst_strain_rate_with_boundary_values: invalid cell volume");
            gradient.component_data(0)[cell] = sum.x / volume;
            gradient.component_data(1)[cell] = sum.y / volume;
            gradient.component_data(2)[cell] = sum.z / volume;
        }
        return gradient;
    };

    const auto grad_u = gradient_component(0);
    const auto grad_v = gradient_component(1);
    const auto grad_w = gradient_component(2);
    Field<double, Location::CELL> strain(
        mesh.n_cells(), "sst_strain_rate", "1/s", 1);
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const std::array<double, 9> grad{
            grad_u.component_data(0)[cell], grad_u.component_data(1)[cell],
            grad_u.component_data(2)[cell], grad_v.component_data(0)[cell],
            grad_v.component_data(1)[cell], grad_v.component_data(2)[cell],
            grad_w.component_data(0)[cell], grad_w.component_data(1)[cell],
            grad_w.component_data(2)[cell]};
        strain(cell) = sst_strain_rate_magnitude(grad);
    }
    return strain;
}

} // namespace cfdx::physics
