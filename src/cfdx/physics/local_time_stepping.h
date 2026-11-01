#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct LocalTimeStepControls {
    double cfl = 1.0;
    double dt_min = 1e-12;
    double dt_max = 1e12;
    double growth_limit = 2.0;
};

inline void validate_local_time_step_controls(const LocalTimeStepControls& c)
{
    if (!std::isfinite(c.cfl) || !(c.cfl > 0.0) ||
        !std::isfinite(c.dt_min) || !(c.dt_min > 0.0) ||
        !std::isfinite(c.dt_max) || c.dt_max < c.dt_min ||
        !std::isfinite(c.growth_limit) || c.growth_limit < 1.0)
        throw std::invalid_argument("invalid local time-step controls");
}

// The independent, cell-volume-aware form permits geometry to be cached by
// callers performing multiple pseudo-time iterations.  For mass flux [kg/s],
// CFL_c = dt_c * sum_f(abs(m_dot_f)) / (rho_c * V_c).
// Only the completed result is written to dt (fail-closed on invalid data).
inline void compute_local_time_step(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& rho,
    const std::vector<double>& cell_volumes,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& dt,
    const LocalTimeStepControls& controls = {})
{
    validate_local_time_step_controls(controls);
    const std::size_t n = mesh.n_cells();
    if (n == 0 || mass_flux.size() != mesh.n_faces() ||
        mass_flux.dimension() != 1 || rho.size() != n ||
        rho.dimension() != 1 || cell_volumes.size() != n ||
        dt.dimension() != 1 ||
        !mesh.cells().is_consistent() ||
        !mesh.cells().face_ids_valid(mesh.n_faces()))
        throw std::invalid_argument("invalid local time-step field or topology");

    const auto& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* offsets = cells.offsets_data();
    std::vector<double> next(n, controls.dt_max);

    for (std::size_t cell = 0; cell < n; ++cell) {
        const double volume = cell_volumes[cell];
        const double density = rho(cell);
        if (!std::isfinite(volume) || !(volume > 0.0) ||
            !std::isfinite(density) || !(density > 0.0))
            throw std::invalid_argument("invalid local CFL density or cell volume");

        double sum_flux = 0.0;
        for (auto k = offsets[cell]; k < offsets[cell + 1]; ++k) {
            const double flux = mass_flux(cell_faces[k]);
            if (!std::isfinite(flux))
                throw std::invalid_argument("non-finite local CFL face mass flux");
            sum_flux += std::abs(flux);
        }
        if (!std::isfinite(sum_flux))
            throw std::invalid_argument("non-finite local CFL flux sum");
        if (sum_flux > 0.0) {
            // Volume is essential for physical dimensional consistency: the
            // previous rho/sum_flux form was only correct for V == 1 m^3.
            const double raw_dt = (controls.cfl * density * volume) / sum_flux;
            if (!std::isfinite(raw_dt) || !(raw_dt > 0.0))
                throw std::invalid_argument("non-representable local CFL time step");
            next[cell] = std::clamp(raw_dt, controls.dt_min, controls.dt_max);
        }
    }
    if (dt.size() != n) dt.resize(n);
    for (std::size_t cell = 0; cell < n; ++cell)
        dt(cell) = next[cell];
}

// Backward-compatible convenience API.  Computing cell geometry at every
// iteration is more costly than passing cached volumes to the overload above.
inline void compute_local_time_step(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& rho,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& dt,
    const LocalTimeStepControls& controls = {})
{
    validate_local_time_step_controls(controls);
    const auto& faces = mesh.faces();
    const auto& cells = mesh.cells();
    if (mesh.n_cells() == 0 || mesh.n_faces() == 0 ||
        !faces.is_consistent() || !faces.indices_valid(mesh.n_points()) ||
        !cells.is_consistent() || !cells.face_ids_valid(mesh.n_faces()) ||
        mesh.ownership().size() != mesh.n_faces() ||
        !mesh.ownership().is_consistent(mesh.n_cells()))
        throw std::invalid_argument("invalid mesh for local CFL cell geometry");

    std::vector<cfdx::core::Vec3> centres(mesh.n_faces());
    std::vector<cfdx::core::Vec3> area_vectors(mesh.n_faces());
    const auto* face_offsets = faces.offsets_data();
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto fg = cfdx::core::compute_face_geometry(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(),
            faces.vertices_data(), face_offsets[f], face_offsets[f + 1] - face_offsets[f]);
        if (!std::isfinite(fg.area) || !(fg.area > 0.0) ||
            !std::isfinite(fg.centre.x) || !std::isfinite(fg.centre.y) ||
            !std::isfinite(fg.centre.z) || !std::isfinite(fg.Sf.x) ||
            !std::isfinite(fg.Sf.y) || !std::isfinite(fg.Sf.z))
            throw std::invalid_argument("invalid face geometry for local CFL");
        centres[f] = fg.centre;
        area_vectors[f] = fg.Sf;
    }

    const auto* offsets = cells.offsets_data();
    const auto* cell_faces = cells.faces_data();
    std::vector<double> volumes(mesh.n_cells());
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const auto cg = cfdx::core::compute_cell_geometry_oriented(
            centres.data(), area_vectors.data(), cell_faces + offsets[cell],
            offsets[cell + 1] - offsets[cell],
            static_cast<cfdx::core::CellIndex>(cell), mesh.ownership());
        if (!std::isfinite(cg.signed_volume) || !(cg.signed_volume > 0.0))
            throw std::invalid_argument("invalid signed cell volume for local CFL");
        volumes[cell] = cg.volume;
    }
    compute_local_time_step(mesh, mass_flux, rho, volumes, dt, controls);
}

inline double global_pseudo_time_step(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& local_dt)
{
    if (local_dt.empty() || local_dt.dimension() != 1)
        throw std::invalid_argument("invalid local time-step field");
    double result = local_dt(0);
    for (std::size_t i = 0; i < local_dt.size(); ++i) {
        const double v = local_dt(i);
        if (!std::isfinite(v) || !(v > 0.0))
            throw std::invalid_argument("non-finite or non-positive local time step");
        result = std::min(result, v);
    }
    return result;
}

} // namespace cfdx::physics
