#pragma once

#include "cfdx/physics/sst_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {


struct SSTProductionState {
    cfdx::core::Field<double, cfdx::core::Location::CELL> k;
    cfdx::core::Field<double, cfdx::core::Location::CELL> omega;
    cfdx::core::Field<double, cfdx::core::Location::CELL> nu_t;
    cfdx::core::Field<double, cfdx::core::Location::CELL> wall_distance;
    cfdx::core::Field<double, cfdx::core::Location::CELL> strain_rate;

    explicit SSTProductionState(std::size_t n_cells)
        : k(n_cells, "k", "m2/s2", 1),
          omega(n_cells, "omega", "1/s", 1),
          nu_t(n_cells, "nu_t", "m2/s", 1),
          wall_distance(n_cells, "wall_distance", "m", 1),
          strain_rate(n_cells, "sst_strain_rate", "1/s", 1)
    {
        if (n_cells == 0)
            throw std::invalid_argument("SST production state requires at least one cell");
    }
};

// Validate the production-owned fields before a state is used by transport or
// momentum assembly. No flooring or automatic initialization is performed.
inline void validate_sst_production_state(const SSTProductionState& state,
                                          std::size_t expected_cells)
{
    const auto valid_shape = [expected_cells](const auto& field) {
        return field.size() == expected_cells && field.dimension() == 1;
    };
    if (expected_cells == 0 ||
        !valid_shape(state.k) || !valid_shape(state.omega) ||
        !valid_shape(state.nu_t) || !valid_shape(state.wall_distance) ||
        !valid_shape(state.strain_rate))
        throw std::invalid_argument("SST production state has an invalid field shape");

    for (std::size_t cell = 0; cell < expected_cells; ++cell) {
        if (!std::isfinite(state.k(cell)) || state.k(cell) < 0.0 ||
            !std::isfinite(state.omega(cell)) || state.omega(cell) <= 0.0 ||
            !std::isfinite(state.nu_t(cell)) || state.nu_t(cell) < 0.0 ||
            !std::isfinite(state.wall_distance(cell)) || state.wall_distance(cell) <= 0.0 ||
            !std::isfinite(state.strain_rate(cell)) || state.strain_rate(cell) < 0.0)
            throw std::invalid_argument(
                "SST production state contains a non-finite or physically invalid cell value");
    }
}

// A value snapshot is deliberately explicit: the production outer iteration
// must include this snapshot alongside U, p and face mass flux in its own
// accept/reject decision. This helper does not claim to provide that coupling.
inline SSTProductionState snapshot_sst_production_state(
    const SSTProductionState& state)
{
    return state;
}

inline void restore_sst_production_state(
    SSTProductionState& state,
    const SSTProductionState& snapshot)
{
    state = snapshot;
}

struct SSTProductionStepResult {
    TurbulenceTransportResult transport;
    bool committed = false;
    std::string reason;
};

// Refresh kinematic turbulent viscosity from one consistent SST iterate.
// Inputs are validated rather than floored here: positivity enforcement belongs
// to the transport solve, and invalid states must not be hidden by the closure.
inline void refresh_sst_production_eddy_viscosity(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& F2,
    const TurbulenceTransportControls& controls,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& nu_t)
{
    const std::size_t n = k.size();
    if (controls.model != TurbulenceModel::SST)
        throw std::invalid_argument("SST production closure requires the SST model");
    if (n == 0 || k.dimension() != 1 || omega.size() != n ||
        omega.dimension() != 1 || strain_rate.size() != n ||
        strain_rate.dimension() != 1 || F2.size() != n ||
        F2.dimension() != 1 || nu_t.size() != n || nu_t.dimension() != 1)
        throw std::invalid_argument("SST production closure field shape mismatch");
    if (!(controls.a1 > 0.0) || !std::isfinite(controls.a1))
        throw std::invalid_argument("SST production closure requires finite positive a1");

    std::vector<double> refreshed(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double ki = k(i), wi = omega(i), si = strain_rate(i), f2 = F2(i);
        if (!std::isfinite(ki) || ki < 0.0 ||
            !std::isfinite(wi) || wi <= 0.0 ||
            !std::isfinite(si) || si < 0.0 ||
            !std::isfinite(f2) || f2 < 0.0 || f2 > 1.0)
            throw std::invalid_argument("SST production closure received an invalid cell state");
        const double denominator = std::max(controls.a1 * wi, si * f2);
        const double value = controls.a1 * ki / denominator;
        if (!std::isfinite(value) || value < 0.0)
            throw std::domain_error("SST production closure produced invalid nu_t");
        refreshed[i] = value;
    }
    for (std::size_t i = 0; i < n; ++i)
        nu_t(i) = refreshed[i];
}

// Execute a dynamic-blending SST transport update transactionally. k, omega,
// and nu_t are committed together only after transport convergence and a valid
// post-update closure. Any non-converged result restores all three fields.
inline SSTProductionStepResult advance_sst_transport_transactionally(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& nu_t,
    const ScalarBoundaryConditions& k_bcs = {},
    const ScalarBoundaryConditions& omega_bcs = {},
    std::size_t max_iterations = 100,
    double tolerance = 1e-8)
{
    const std::size_t n = mesh.n_cells();
    if (controls.model != TurbulenceModel::SST)
        throw std::invalid_argument("transactional SST step requires the SST model");
    if (n == 0 || k.size() != n || omega.size() != n ||
        strain_rate.size() != n || wall_distance.size() != n ||
        nu_t.size() != n || mass_flux.size() != mesh.n_faces())
        throw std::invalid_argument("transactional SST step field size mismatch");
    if (max_iterations == 0 || !(tolerance > 0.0) || !std::isfinite(tolerance))
        throw std::invalid_argument("transactional SST step requires valid iteration controls");
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(wall_distance(i)) || wall_distance(i) <= 0.0 ||
            !std::isfinite(strain_rate(i)) || strain_rate(i) < 0.0 ||
            !std::isfinite(k(i)) || k(i) < 0.0 ||
            !std::isfinite(omega(i)) || omega(i) <= 0.0 ||
            !std::isfinite(nu_t(i)) || nu_t(i) < 0.0)
            throw std::invalid_argument("transactional SST step received an invalid initial field");
    }

    const auto old_k = k;
    const auto old_omega = omega;
    const auto old_nu_t = nu_t;
    SSTProductionStepResult step;
    const auto rollback = [&]() {
        k = old_k;
        omega = old_omega;
        nu_t = old_nu_t;
    };

    try {
        step.transport = solve_sst_transport_dynamic_blending(
            mesh, geometry, mass_flux, k, omega, strain_rate, wall_distance,
            controls, k_bcs, omega_bcs, max_iterations, tolerance);
        if (!step.transport.converged ||
            step.transport.nonfinite_bound_events != 0 ||
            !std::isfinite(step.transport.k_residual) ||
            !std::isfinite(step.transport.second_residual)) {
            rollback();
            step.reason = "SST transport did not satisfy convergence and finite-state gates";
            return step;
        }

        const auto cross = cell_gradient_dot(k, omega, mesh);
        cfdx::core::Field<double, cfdx::core::Location::CELL> F2(
            n, "sst_F2", "1", 1);
        for (std::size_t i = 0; i < n; ++i) {
            const auto blend = compute_sst_blending(
                k(i), omega(i), wall_distance(i), controls.molecular_viscosity,
                controls.beta_star, cross[i], controls.density,
                controls.sst_sigma_w2);
            F2(i) = blend.second;
        }
        refresh_sst_production_eddy_viscosity(
            k, omega, strain_rate, F2, controls, nu_t);
        step.committed = true;
        step.reason = "SST transport and eddy-viscosity refresh committed";
        return step;
    } catch (...) {
        rollback();
        throw;
    }
}

}  // namespace cfdx::physics
