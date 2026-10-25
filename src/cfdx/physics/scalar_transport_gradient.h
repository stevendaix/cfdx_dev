#pragma once

#include "cfdx/physics/finite_volume_transport.h"

namespace cfdx::physics {

// Explicit face data is an all-faces contract for each named patch. Entries
// outside that patch may be NaN; no value on an actual patch face may be absent.
// Inspect mesh face_ids, not just the supplied vectors or geometry lookup.
inline void validate_transport_boundary_face_values(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryFaceValues& values,
    bool require_positive = false)
{
    if (geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument("boundary face values: geometry shape mismatch");
    for (const auto& [name, data] : values.values) {
        const auto index = mesh.boundary().find(name);
        if (index == mesh.boundary().n_patches() || data.size() != mesh.n_faces())
            throw std::invalid_argument("boundary face values: unknown patch or invalid global-face shape: " + name);
        const auto& patch = mesh.boundary().patch(index);
        if (patch.face_ids.empty() || patch.type == cfdx::core::PatchType::EMPTY)
            throw std::invalid_argument("boundary face values: empty patch: " + name);
        for (const auto face : patch.face_ids) {
            if (face >= mesh.n_faces() || geometry.face_patch[face] != index ||
                mesh.ownership().neighbour(face) >= 0 ||
                mesh.ownership().owner(face) >= mesh.n_cells())
                throw std::invalid_argument("boundary face values: invalid mesh membership: " + name);
            if (!std::isfinite(data[face]) || (require_positive && data[face] <= 0.0))
                throw std::invalid_argument("boundary face values: missing or invalid value at face " +
                                            std::to_string(face) + " on patch " + name);
        }
    }
}

// Green--Gauss reconstruction with exactly the scalar transport boundary
// states: prescribed value, prescribed outward normal gradient, or owner value
// for zero gradient. Per-face values override patch values. Internal faces use
// the same longitudinal geometric interpolation as the core Gauss operator.
// This does not add a skewness correction or claim affine exactness on skew meshes.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
scalar_gradient_with_transport_boundaries(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& bcs = {},
    const ScalarBoundaryFaceValues* face_values = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells(), nf = mesh.n_faces();
    if (field.size() != nc || field.dimension() != 1 ||
        geometry.cell_centres.size() != nc || geometry.cell_volumes.size() != nc ||
        geometry.face_centres.size() != nf || geometry.face_area_vectors.size() != nf ||
        geometry.face_patch.size() != nf)
        throw std::invalid_argument("scalar transport gradient: field or geometry shape mismatch");
    if (face_values)
        validate_transport_boundary_face_values(mesh, geometry, *face_values);
    std::vector<Vec3> sums(nc);
    for (std::size_t cell = 0; cell < nc; ++cell) {
        if (!std::isfinite(field(cell)) || !(geometry.cell_volumes[cell] > 0.0) ||
            !std::isfinite(geometry.cell_volumes[cell]))
            throw std::invalid_argument("scalar transport gradient: invalid field or cell volume");
    }
    for (std::size_t face = 0; face < nf; ++face) {
        const std::size_t owner = mesh.ownership().owner(face);
        const auto neighbour = mesh.ownership().neighbour(face);
        if (owner >= nc || (neighbour >= 0 && static_cast<std::size_t>(neighbour) >= nc))
            throw std::invalid_argument("scalar transport gradient: invalid face ownership");
        double value = field(owner);
        if (neighbour >= 0) {
            const auto nb = static_cast<std::size_t>(neighbour);
            const double weight = geometric_interpolation_weight(
                geometry.face_centres[face], geometry.cell_centres[owner],
                &geometry.cell_centres[nb]);
            value = (1.0 - weight) * value + weight * field(nb);
        } else {
            const auto patch = geometry.face_patch[face];
            if (patch < mesh.boundary().n_patches()) {
                const auto& boundary = mesh.boundary().patch(patch);
                if (boundary.type == PatchType::EMPTY)
                    continue;
                ScalarBoundaryCondition bc;
                const auto entry = bcs.find(boundary.name);
                if (entry != bcs.end()) bc = entry->second;
                if (face_values) {
                    const auto data = face_values->values.find(boundary.name);
                    if (data != face_values->values.end())
                        bc = {ScalarBoundaryType::FIXED_VALUE, data->second[face], 0.0};
                }
                if (bc.type == ScalarBoundaryType::FIXED_VALUE)
                    value = bc.value;
                else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT)
                    value += bc.gradient * boundary_normal_distance(geometry, face, owner);
            }
        }
        const auto contribution = geometry.face_area_vectors[face] * value;
        if (!std::isfinite(contribution.x) || !std::isfinite(contribution.y) ||
            !std::isfinite(contribution.z))
            throw std::invalid_argument("scalar transport gradient: non-finite face contribution");
        sums[owner] += contribution;
        if (neighbour >= 0)
            sums[static_cast<std::size_t>(neighbour)] += contribution * -1.0;
    }
    Field<double, Location::CELL> result(nc, field.name() + "_grad_bc", "", 3);
    for (std::size_t cell = 0; cell < nc; ++cell) {
        result(cell, 0) = sums[cell].x / geometry.cell_volumes[cell];
        result(cell, 1) = sums[cell].y / geometry.cell_volumes[cell];
        result(cell, 2) = sums[cell].z / geometry.cell_volumes[cell];
    }
    return result;
}

} // namespace cfdx::physics
