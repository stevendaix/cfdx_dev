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

// Boundary-constrained weighted LS for skewed/polyhedral FV cells.
//
// Unlike cell-centre Gauss interpolation, this operator fits the actual
// face-centre Dirichlet data and normal Neumann constraints, so it reconstructs
// affine fields on skew meshes. EXPLICIT opt-in: the historical Gauss operator
// above remains the default until each SST/viscous validation is complete.
// A 2-D embedded mesh must request intrinsic_dimension=2; an invalid or
// deficient stencil is rejected rather than silently returning a gradient.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
scalar_gradient_with_transport_boundaries_wls(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& bcs = {},
    const ScalarBoundaryFaceValues* face_values = nullptr,
    int intrinsic_dimension = 3,
    double max_condition = 1e10)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells(), nf = mesh.n_faces();
    // Finite condition-number gate avoids returning an amplified gradient on
    // nearly singular stencils. Callers may choose a stricter threshold.
    if ((intrinsic_dimension != 2 && intrinsic_dimension != 3) ||
        !(max_condition > 0.0) || !std::isfinite(max_condition) ||
        field.size() != nc || field.dimension() != 1 ||
        geometry.cell_centres.size() != nc ||
        geometry.cell_volumes.size() != nc ||
        geometry.face_centres.size() != nf ||
        geometry.face_area_vectors.size() != nf ||
        geometry.face_patch.size() != nf ||
        !mesh.boundary().is_consistent(nf) ||
        !mesh.faces().is_consistent() ||
        mesh.ownership().size() != nf)
        throw std::invalid_argument("scalar WLS gradient: invalid field, dimension, or mesh geometry");
    if (face_values)
        validate_transport_boundary_face_values(mesh, geometry, *face_values);

    Field<double, Location::CELL> result(nc,
        field.name() + "_grad_bc_wls", field.metadata().unit + "/m", 3);
    const auto* cell_faces = mesh.cells().faces_data();
    const auto* offsets = mesh.cells().offsets_data();
    for (std::size_t cell = 0; cell < nc; ++cell) {
        const auto centre = geometry.cell_centres[cell];
        const double phi = field(cell);
        if (!std::isfinite(phi) || !std::isfinite(centre.x) ||
            !std::isfinite(centre.y) || !std::isfinite(centre.z) ||
            !(geometry.cell_volumes[cell] > 0.0) ||
            !std::isfinite(geometry.cell_volumes[cell]))
            throw std::invalid_argument("scalar WLS gradient: non-finite cell geometry or field");
        std::vector<Vec3> sample_centres;
        std::vector<double> sample_values;
        for (Offset o = offsets[cell]; o < offsets[cell+1]; ++o) {
            const std::size_t face = cell_faces[o];
            if (face >= nf)
                throw std::invalid_argument("scalar WLS gradient: invalid cell face index");
            const std::size_t owner = mesh.ownership().owner(face);
            const auto other = mesh.ownership().neighbour(face);
            if (owner >= nc || (other >= 0 && static_cast<std::size_t>(other) >= nc))
                throw std::invalid_argument("scalar WLS gradient: invalid ownership");
            if (other >= 0) {
                const std::size_t neighbour = owner == cell
                    ? static_cast<std::size_t>(other) : owner;
                if ((owner != cell && static_cast<std::size_t>(other) != cell) ||
                    neighbour == cell)
                    throw std::invalid_argument("scalar WLS gradient: cell-face connectivity mismatch");
                sample_centres.push_back(geometry.cell_centres[neighbour]);
                sample_values.push_back(field(neighbour));
                continue;
            }
            if (owner != cell)
                throw std::invalid_argument("scalar WLS gradient: boundary face owned by another cell");
            const std::size_t patch_id = geometry.face_patch[face];
            if (patch_id >= mesh.boundary().n_patches())
                throw std::invalid_argument("scalar WLS gradient: boundary face without a patch");
            const auto& patch = mesh.boundary().patch(patch_id);
            if (patch.type == PatchType::EMPTY)
                continue;
            ScalarBoundaryCondition bc;
            const auto it = bcs.find(patch.name);
            if (it != bcs.end()) bc = it->second;
            if (face_values) {
                const auto by_patch = face_values->values.find(patch.name);
                if (by_patch != face_values->values.end())
                    bc = {ScalarBoundaryType::FIXED_VALUE,
                          by_patch->second[face], 0.0};
            }
            if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                if (!std::isfinite(bc.value))
                    throw std::invalid_argument("scalar WLS gradient: non-finite Dirichlet face value");
                sample_centres.push_back(geometry.face_centres[face]);
                sample_values.push_back(bc.value);
                continue;
            }
            const Vec3 sf = geometry.face_area_vectors[face];
            const double area = sf.mag();
            if (!(area > 0.0) || !std::isfinite(area))
                throw std::invalid_argument("scalar WLS gradient: invalid boundary face area");
            const Vec3 n = sf*(1.0/area);
            const double dn = (geometry.face_centres[face] - centre).dot(n);
            if (!(dn > 0.0) || !std::isfinite(dn))
                throw std::invalid_argument("scalar WLS gradient: non-positive outward boundary distance");
            if (bc.type != ScalarBoundaryType::ZERO_GRADIENT &&
                bc.type != ScalarBoundaryType::FIXED_GRADIENT)
                throw std::invalid_argument("scalar WLS gradient: unknown boundary condition");
            const double derivative = bc.type == ScalarBoundaryType::FIXED_GRADIENT
                ? bc.gradient : 0.0;
            if (!std::isfinite(derivative))
                throw std::invalid_argument("scalar WLS gradient: non-finite Neumann gradient");
            // A normal ghost point encodes only the prescribed normal
            // derivative; no fictitious tangential component is imposed.
            sample_centres.push_back(centre + n * (2.0*dn));
            sample_values.push_back(phi + 2.0*dn*derivative);
        }
        StencilQuality quality;
        const Vec3 g = weighted_least_squares_gradient_robust(
            centre, phi, sample_centres, sample_values,
            GradientWeighting::INVERSE_DISTANCE_SQUARED,
            intrinsic_dimension, &quality);
        if (!quality.full_rank || !quality.well_conditioned(max_condition))
            throw std::runtime_error("scalar WLS gradient: ill-conditioned stencil at cell "
                                     + std::to_string(cell));
        result(cell,0) = g.x;
        result(cell,1) = g.y;
        result(cell,2) = g.z;
    }
    return result;
}

} // namespace cfdx::physics
