#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/flux.h"
#include "cfdx/core/numerics/interpolation.h"
#include <cstddef>
#include <stdexcept>

namespace cfdx::core {

inline Field<double, Location::CELL> compute_convection(
    const Field<double, Location::CELL>& scalar,
    const Field<double, Location::FACE>& face_flux,
    const Mesh& mesh,
    InterpScheme scheme = InterpScheme::UPWIND)
{
    if (scalar.dimension() != 1) {
        throw std::runtime_error("compute_convection: scalar field required");
    }
    if (face_flux.dimension() != 1 || face_flux.size() != mesh.n_faces()) {
        throw std::runtime_error("compute_convection: invalid face flux");
    }
    if (scalar.size() != mesh.n_cells()) {
        throw std::runtime_error("compute_convection: field size != n_cells");
    }

    const auto face_value = interpolate_cell_to_face(scalar, mesh, scheme, &face_flux);
    const auto& ownership = mesh.ownership();
    const auto& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();
    const double* phi = face_value.component_data(0);
    const double* flux = face_flux.component_data(0);

    Field<double, Location::CELL> result(
        mesh.n_cells(), scalar.name() + "_convection",
        scalar.metadata().unit + "/s", 1);
    double* out = result.component_data(0);

    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        double sum = 0.0;
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        for (Offset k = 0; k < n; ++k) {
            const std::size_t f = cell_faces[off + k];
            const double signed_flux = ownership.owner(f) == c ? flux[f] : -flux[f];
            sum += signed_flux * phi[f];
        }
        out[c] = sum;
    }
    return result;
}

inline Field<double, Location::CELL> compute_vector_convection(
    const Field<double, Location::CELL>& vector,
    const Field<double, Location::FACE>& face_flux,
    const Mesh& mesh,
    InterpScheme scheme = InterpScheme::UPWIND)
{
    if (vector.dimension() != 3) {
        throw std::runtime_error("compute_vector_convection: 3-component field required");
    }
    Field<double, Location::CELL> result(
        mesh.n_cells(), vector.name() + "_convection",
        vector.metadata().unit + "/s", 3);
    for (std::size_t d = 0; d < 3; ++d) {
        Field<double, Location::CELL> component(
            mesh.n_cells(), vector.name() + "_component", vector.metadata().unit, 1);
        double* dst = component.component_data(0);
        const double* src = vector.component_data(d);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) dst[c] = src[c];
        const auto conv = compute_convection(component, face_flux, mesh, scheme);
        double* out = result.component_data(d);
        const double* in = conv.component_data(0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) out[c] = in[c];
    }
    return result;
}

}  // namespace cfdx::core
