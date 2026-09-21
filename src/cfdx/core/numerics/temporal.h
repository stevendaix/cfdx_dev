// M0.13 — Temporal Discretization
//
// Spécification CFDX v0.7 §33-35 :
//   Schémas d'intégration temporelle pour ∂φ/∂t + ∇·F = S
//
//   Schémas supportés :
//     - Euler explicite : φ^{n+1} = φ^n + Δt * RHS(φ^n)
//     - Euler implicite : φ^{n+1} = φ^n + Δt * RHS(φ^{n+1})
//     - Crank-Nicolson : φ^{n+1} = φ^n + 0.5*Δt * (RHS(φ^n) + RHS(φ^{n+1}))
//     - BDF2 : (3φ^{n+1} - 4φ^n + φ^{n-1}) / (2Δt) = RHS(φ^{n+1})
//
//   Local time stepping : Δt_local = CFL * h / |u| pour convergence régime permanent
//   Contrôle adaptatif : Δt^{n+1} = Δt^n * min(max_factor, max(min_factor, CFL_target / CFL_current))

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <cstddef>
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

namespace cfdx {
namespace core {

enum class TimeScheme : std::uint8_t {
    EULER_EXPLICIT = 0,
    EULER_IMPLICIT,
    CRANK_NICOLSON,
    BDF2
};

inline const char* to_string(TimeScheme s) {
    switch (s) {
        case TimeScheme::EULER_EXPLICIT:  return "euler_explicit";
        case TimeScheme::EULER_IMPLICIT:  return "euler_implicit";
        case TimeScheme::CRANK_NICOLSON:  return "crank_nicolson";
        case TimeScheme::BDF2:            return "bdf2";
        default:                           return "unknown";
    }
}

inline TimeScheme time_scheme_from_string(const std::string& s) {
    if (s == "euler_explicit")  return TimeScheme::EULER_EXPLICIT;
    if (s == "euler_implicit")  return TimeScheme::EULER_IMPLICIT;
    if (s == "crank_nicolson")  return TimeScheme::CRANK_NICOLSON;
    if (s == "bdf2")            return TimeScheme::BDF2;
    throw std::runtime_error("time_scheme_from_string: unknown scheme '" + s + "'");
}

// RHS function type: computes spatial residual RHS = -∇·F + S
// Input: current field φ, Output: residual (size = n_cells, dim = field dimension)
using RhsFunction = std::function<void(const Field<double, Location::CELL>&, Field<double, Location::CELL>&)>;

// Time integration context - holds history for multi-step schemes
struct TimeIntegrationContext {
    std::size_t n_cells = 0;
    std::size_t field_dim = 0;

    // For BDF2: stores φ^{n-1} and φ^n
    Field<double, Location::CELL> phi_prev;  // φ^{n-1}
    Field<double, Location::CELL> phi_curr;  // φ^n

    bool has_prev = false;  // true after first step

    TimeIntegrationContext() = default;

    TimeIntegrationContext(std::size_t n, std::size_t dim, const std::string& name)
        : n_cells(n), field_dim(dim),
          phi_prev(n, name + "_prev", "1", dim),
          phi_curr(n, name + "_curr", "1", dim)
    {}

    void initialize(const Field<double, Location::CELL>& phi0) {
        phi_curr = phi0;  // copy
        has_prev = false;
    }

    // Shift history: φ^{n-1} = φ^n, φ^n = φ^{n+1}
    void shift(const Field<double, Location::CELL>& phi_new) {
        phi_prev = phi_curr;
        phi_curr = phi_new;
        has_prev = true;
    }
};

// Advance field in time by one step
// Returns the new field φ^{n+1}
inline Field<double, Location::CELL> advance_time(
    const Field<double, Location::CELL>& phi,
    double dt,
    const RhsFunction& rhs_func,
    TimeScheme scheme,
    TimeIntegrationContext* ctx = nullptr)
{
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::invalid_argument("advance_time: dt must be finite and strictly positive");
    }

    const std::size_t n_cells = phi.size();
    const std::size_t dim = phi.dimension();
    Field<double, Location::CELL> phi_new(
        n_cells, "phi_new", phi.metadata().unit, dim);

    Field<double, Location::CELL> rhs(
        n_cells, "rhs", phi.metadata().unit + "/s", dim);
    Field<double, Location::CELL> rhs_iter(
        n_cells, "rhs_iter", phi.metadata().unit + "/s", dim);

    rhs_func(phi, rhs);

    constexpr int max_iterations = 100;
    constexpr double tolerance = 1e-12;

    auto solve_fixed_point = [&](auto&& predictor) {
        predictor();
        for (int iteration = 0; iteration < max_iterations; ++iteration) {
            rhs_func(phi_new, rhs_iter);
            double max_delta = 0.0;
            for (std::size_t d = 0; d < dim; ++d) {
                const double* old_data = phi_new.component_data(d);
                const double* phi_data = phi.component_data(d);
                const double* rhs_data = rhs_iter.component_data(d);
                double* new_data = phi_new.component_data(d);
                (void)old_data;
                for (std::size_t cell = 0; cell < n_cells; ++cell) {
                    double candidate = 0.0;
                    if (scheme == TimeScheme::EULER_IMPLICIT) {
                        candidate = phi_data[cell] + dt * rhs_data[cell];
                    } else if (scheme == TimeScheme::CRANK_NICOLSON) {
                        const double* rhs_old = rhs.component_data(d);
                        candidate = phi_data[cell] +
                            0.5 * dt * (rhs_old[cell] + rhs_data[cell]);
                    } else {
                        const double* prev = ctx->phi_prev.component_data(d);
                        candidate = (4.0 * phi_data[cell] - prev[cell] +
                            2.0 * dt * rhs_data[cell]) / 3.0;
                    }
                    max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                    new_data[cell] = candidate;
                }
            }
            if (max_delta <= tolerance) {
                return;
            }
        }
        throw std::runtime_error("advance_time: implicit iteration did not converge");
    };

    switch (scheme) {
        case TimeScheme::EULER_EXPLICIT:
            for (std::size_t d = 0; d < dim; ++d) {
                const double* phi_data = phi.component_data(d);
                const double* rhs_data = rhs.component_data(d);
                double* new_data = phi_new.component_data(d);
                for (std::size_t cell = 0; cell < n_cells; ++cell) {
                    new_data[cell] = phi_data[cell] + dt * rhs_data[cell];
                }
            }
            break;

        case TimeScheme::EULER_IMPLICIT:
            solve_fixed_point([&] {
                for (std::size_t d = 0; d < dim; ++d) {
                    const double* phi_data = phi.component_data(d);
                    double* new_data = phi_new.component_data(d);
                    for (std::size_t cell = 0; cell < n_cells; ++cell) {
                        new_data[cell] = phi_data[cell];
                    }
                }
            });
            break;

        case TimeScheme::CRANK_NICOLSON:
            solve_fixed_point([&] {
                for (std::size_t d = 0; d < dim; ++d) {
                    const double* phi_data = phi.component_data(d);
                    const double* rhs_data = rhs.component_data(d);
                    double* new_data = phi_new.component_data(d);
                    for (std::size_t cell = 0; cell < n_cells; ++cell) {
                        new_data[cell] = phi_data[cell] + dt * rhs_data[cell];
                    }
                }
            });
            break;

        case TimeScheme::BDF2:
            if (!ctx || !ctx->has_prev) {
                // Bootstrap BDF2 with a genuinely implicit Euler step.
                for (std::size_t d = 0; d < dim; ++d) {
                    const double* phi_data = phi.component_data(d);
                    double* new_data = phi_new.component_data(d);
                    for (std::size_t cell = 0; cell < n_cells; ++cell) {
                        new_data[cell] = phi_data[cell];
                    }
                }
                for (int iteration = 0; iteration < max_iterations; ++iteration) {
                    rhs_func(phi_new, rhs_iter);
                    double max_delta = 0.0;
                    for (std::size_t d = 0; d < dim; ++d) {
                        const double* phi_data = phi.component_data(d);
                        const double* rhs_data = rhs_iter.component_data(d);
                        double* new_data = phi_new.component_data(d);
                        for (std::size_t cell = 0; cell < n_cells; ++cell) {
                            const double candidate = phi_data[cell] + dt * rhs_data[cell];
                            max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                            new_data[cell] = candidate;
                        }
                    }
                    if (max_delta <= tolerance) break;
                    if (iteration == max_iterations - 1) {
                        throw std::runtime_error("advance_time: implicit Euler bootstrap did not converge");
                    }
                }
            } else {
                solve_fixed_point([&] {});
            }
            break;
    }

    if (ctx) {
        ctx->shift(phi_new);
    }
    return phi_new;
}


// Compute stable time step based on CFL condition
// CFL = |u| * Δt / h ≤ CFL_max
// For convection: Δt = CFL_max * h / (|u| + c) where c = speed of sound (if compressible)
inline double compute_cfl_time_step(
    const Mesh& mesh,
    const Field<double, Location::CELL>& velocity,  // dim=3
    double cfl_max = 0.5,
    double speed_of_sound = 0.0)
{
    const std::size_t n_cells = mesh.n_cells();
    const double* ux = velocity.component_data(0);
    const double* uy = velocity.component_data(1);
    const double* uz = velocity.component_data(2);

    // Compute cell characteristic length h = V^(1/3) (or similar)
    // Need cell volumes
    std::vector<Vec3> face_centres(mesh.n_faces());
    std::vector<Vec3> face_Sf(mesh.n_faces());

    const PointCloud& pts = mesh.points();
    const double* px = pts.x_data();
    const double* py = pts.y_data();
    const double* pz = pts.z_data();
    const auto* face_verts = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const Offset off = face_offsets[f];
        const Offset n = face_offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(px, py, pz, face_verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    std::vector<Vec3> provisional_centres(n_cells);
    compute_area_weighted_cell_centres(mesh, face_centres.data(), face_Sf.data(), provisional_centres.data());
    orient_mesh_face_vectors(mesh, face_centres, provisional_centres, face_Sf);

    std::vector<double> cell_volume(n_cells, 0.0);
    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        const CellGeometry cg = compute_cell_geometry(mesh, face_centres.data(), face_Sf.data(), cell_faces + off, c, n);
        cell_volume[c] = cg.volume;
    }

    double dt_min = 1e30;
    for (std::size_t c = 0; c < n_cells; ++c) {
        if (cell_volume[c] <= 0) continue;

        const double h = std::cbrt(cell_volume[c]);  // characteristic length
        const double vel_mag = std::sqrt(ux[c]*ux[c] + uy[c]*uy[c] + uz[c]*uz[c]);
        const double wave_speed = vel_mag + speed_of_sound;

        if (wave_speed > 1e-12) {
            double dt_local = cfl_max * h / wave_speed;
            dt_min = std::min(dt_min, dt_local);
        }
    }

    return dt_min;
}

// Adaptive time step controller
struct AdaptiveTimeStepper {
    double dt = 1e-3;          // current time step
    double dt_min = 1e-12;     // minimum allowed
    double dt_max = 1.0;       // maximum allowed
    double cfl_target = 0.5;   // target CFL
    double cfl_current = 0.0;  // current CFL (computed externally)
    double safety_factor = 0.9;
    double growth_factor = 1.2;
    double shrink_factor = 0.5;

    // Update time step based on current CFL
    void update(double cfl_achieved) {
        cfl_current = cfl_achieved;
        if (cfl_achieved > 0) {
            double ratio = cfl_target / cfl_achieved;
            double new_dt = dt * safety_factor * ratio;
            // Limit growth/shrink
            if (new_dt > dt * growth_factor) new_dt = dt * growth_factor;
            if (new_dt < dt * shrink_factor) new_dt = dt * shrink_factor;
            dt = std::clamp(new_dt, dt_min, dt_max);
        }
    }

    // Called when step is rejected (e.g., non-convergence)
    void reject() {
        dt = std::max(dt * shrink_factor, dt_min);
    }
};

// Local time stepping for steady-state acceleration
// Δt_local = CFL * h / |u| for each cell
inline std::vector<double> compute_local_time_steps(
    const Mesh& mesh,
    const Field<double, Location::CELL>& velocity,
    double cfl = 0.5)
{
    const std::size_t n_cells = mesh.n_cells();
    const double* ux = velocity.component_data(0);
    const double* uy = velocity.component_data(1);
    const double* uz = velocity.component_data(2);

    // Compute cell volumes
    std::vector<Vec3> face_centres(mesh.n_faces());
    std::vector<Vec3> face_Sf(mesh.n_faces());

    const auto* face_verts = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();
    (void)face_verts; // used in compute_face_geometry

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const Offset off = face_offsets[f];
        const Offset n = face_offsets[f + 1] - off;
        const FaceGeometry fg = compute_face_geometry(
            mesh.points().x_data(), mesh.points().y_data(), mesh.points().z_data(), face_verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    std::vector<double> cell_volume(n_cells, 0.0);
    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();
    (void)cell_faces; // used in compute_cell_geometry

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        (void)n; // used in compute_cell_geometry
        const CellGeometry cg = compute_cell_geometry(mesh, face_centres.data(), face_Sf.data(), cell_faces + off, c, n);
        cell_volume[c] = cg.volume;
    }

    // Large initial dt for cells with zero velocity
    const double dt_max = 1e30;
    std::vector<double> dt_local(n_cells, dt_max);
    for (std::size_t c = 0; c < n_cells; ++c) {
        if (cell_volume[c] <= 0) continue;
        const double h = std::cbrt(cell_volume[c]);
        const double vel_mag = std::sqrt(ux[c]*ux[c] + uy[c]*uy[c] + uz[c]*uz[c]);
        if (vel_mag > 1e-12) {
            dt_local[c] = cfl * h / vel_mag;
        }
    }
    return dt_local;
}

// Apply local time stepping: φ^{n+1}_c = φ^n_c + Δt_c * RHS_c
inline void apply_local_time_stepping(
    const Field<double, Location::CELL>& phi,
    const std::vector<double>& dt_local,
    const RhsFunction& rhs_func,
    Field<double, Location::CELL>& phi_new)
{
    const std::size_t n_cells = phi.size();
    const std::size_t dim = phi.dimension();

    Field<double, Location::CELL> rhs(n_cells, "rhs", phi.metadata().unit + "/s", dim);
    rhs_func(phi, rhs);

    for (std::size_t c = 0; c < n_cells; ++c) {
        for (std::size_t d = 0; d < dim; ++d) {
            const double* phi_data = phi.component_data(d);
            const double* rhs_data = rhs.component_data(d);
            double* new_data = phi_new.component_data(d);
            new_data[c] = phi_data[c] + dt_local[c] * rhs_data[c];
        }
    }
}

}  // namespace core
}  // namespace cfdx