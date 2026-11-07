#pragma once

// 2-D planar (XY) NACA qualification on one extruded hexahedral layer.
// This is a validation guard, NOT an axisymmetric or reduced-order solver.
// It intentionally rejects multi-layer 3-D verification meshes and any
// front/back patch that is not EMPTY.  It does not replace a Jacobian MMS.

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

struct NacaPlanarGeometry {
    double z_min = 0.0;
    double z_max = 0.0;
    double span = 0.0;
    std::size_t cells = 0;
};

struct NacaPlanarState {
    double max_abs_uz = 0.0;
    double max_abs_empty_face_flux = 0.0;
    std::size_t checked_empty_faces = 0;
};

inline NacaPlanarGeometry require_naca_planar_geometry(const cfdx::core::Mesh& mesh) {
    using namespace cfdx::core;
    if (!mesh.n_cells() || !mesh.n_points() || !mesh.n_faces() ||
        !mesh.points().is_valid())
        throw std::invalid_argument("NACA2D: missing or non-finite geometry");

    NacaPlanarGeometry result;
    result.z_min = std::numeric_limits<double>::infinity();
    result.z_max = -std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < mesh.n_points(); ++i) {
        const double z = mesh.points().z(i);
        result.z_min = std::min(result.z_min, z);
        result.z_max = std::max(result.z_max, z);
    }
    result.span = result.z_max - result.z_min;
    result.cells = mesh.n_cells();
    if (!std::isfinite(result.span) || !(result.span > 0.0))
        throw std::invalid_argument("NACA2D: extrusion span must be finite and positive");

    // The tolerance scales with the span but never admits an intermediate
    // z-layer. Avoid interpreting a near-zero span as a planar extrusion.
    const double absolute_scale = std::max({1.0, std::abs(result.z_min),
                                           std::abs(result.z_max)});
    const double tolerance = std::max(1e-12 * result.span, 1e-14 * absolute_scale);
    if (!(result.span > 4.0 * tolerance))
        throw std::invalid_argument("NACA2D: extrusion span is not resolvable");
    for (std::size_t i = 0; i < mesh.n_points(); ++i) {
        const double z = mesh.points().z(i);
        if (std::min(std::abs(z - result.z_min),
                     std::abs(z - result.z_max)) > tolerance)
            throw std::invalid_argument("NACA2D: intermediate z plane: not a 2-D EMPTY mesh");
    }

    std::vector<std::uint8_t> front_count(mesh.n_cells(), 0);
    std::vector<std::uint8_t> back_count(mesh.n_cells(), 0);
    const auto inspect = [&](const char* name, double z_plane,
                             std::vector<std::uint8_t>& count) {
        const auto patch_id = mesh.boundary().find(name);
        if (patch_id >= mesh.boundary().n_patches())
            throw std::invalid_argument(std::string("NACA2D: missing ") + name + " patch");
        const auto& patch = mesh.boundary().patch(patch_id);
        if (patch.type != PatchType::EMPTY || patch.face_ids.empty())
            throw std::invalid_argument(std::string("NACA2D: ") + name + " must be EMPTY");
        for (const auto face : patch.face_ids) {
            if (face >= mesh.n_faces() || mesh.ownership().neighbour(face) >= 0)
                throw std::invalid_argument("NACA2D: invalid EMPTY boundary face");
            const auto owner = mesh.ownership().owner(face);
            if (owner >= mesh.n_cells() || ++count[owner] != 1)
                throw std::invalid_argument("NACA2D: more than one EMPTY face per cell");
            const auto start = mesh.faces().face_offset(face);
            const auto length = mesh.faces().face_size(face);
            if (length < 3)
                throw std::invalid_argument("NACA2D: degenerate EMPTY face");
            for (std::size_t j = 0; j < length; ++j) {
                const auto vertex = mesh.faces().vertices()[start + j];
                if (vertex >= mesh.n_points() ||
                    std::abs(mesh.points().z(vertex) - z_plane) > tolerance)
                    throw std::invalid_argument("NACA2D: EMPTY face not on extrusion plane");
            }
        }
    };
    inspect("front", result.z_min, front_count);
    inspect("back", result.z_max, back_count);
    for (std::size_t i = 0; i < mesh.n_cells(); ++i)
        if (front_count[i] != 1 || back_count[i] != 1)
            throw std::invalid_argument("NACA2D: expected exactly one z-layer per cell");
    return result;
}

inline NacaPlanarState require_naca_planar_state(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& authoritative_mass_flux,
    double u_ref, double flux_ref,
    double uz_relative_gate = 1e-10,
    double flux_relative_gate = 1e-12) {
    using namespace cfdx::core;
    // Check the geometric contract again: a zero W on a wrong mesh is not 2-D.
    (void)require_naca_planar_geometry(mesh);
    if (!std::isfinite(u_ref) || !(u_ref > 0.0) ||
        !std::isfinite(flux_ref) || !(flux_ref > 0.0) ||
        !std::isfinite(uz_relative_gate) || !(uz_relative_gate >= 0.0) ||
        !std::isfinite(flux_relative_gate) || !(flux_relative_gate >= 0.0) ||
        velocity.size() != mesh.n_cells() || velocity.dimension() != 3 ||
        authoritative_mass_flux.size() != mesh.n_faces())
        throw std::invalid_argument("NACA2D: invalid field shapes or reference scales");
    NacaPlanarState result;
    for (std::size_t i = 0; i < mesh.n_cells(); ++i) {
        for (std::size_t component = 0; component < 3; ++component)
            if (!std::isfinite(velocity.component_data(component)[i]))
                throw std::runtime_error("NACA2D: non-finite velocity");
        result.max_abs_uz = std::max(result.max_abs_uz,
                                    std::abs(velocity.component_data(2)[i]));
    }
    if (result.max_abs_uz > uz_relative_gate * u_ref)
        throw std::runtime_error("NACA2D: nonzero out-of-plane velocity");
    for (const char* patch_name : {"front", "back"}) {
        const auto& patch = mesh.boundary().patch(mesh.boundary().find(patch_name));
        for (const auto face : patch.face_ids) {
            const double flux = authoritative_mass_flux(face);
            if (!std::isfinite(flux))
                throw std::runtime_error("NACA2D: non-finite EMPTY face flux");
            result.max_abs_empty_face_flux =
                std::max(result.max_abs_empty_face_flux, std::abs(flux));
            ++result.checked_empty_faces;
        }
    }
    if (result.max_abs_empty_face_flux > flux_relative_gate * flux_ref)
        throw std::runtime_error("NACA2D: nonzero EMPTY face mass flux");
    return result;
}

}  // namespace cfdx::physics
