#pragma once

#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_case_contract.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
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
// momentum assembly. Model conventions are governed by docs/audits/09-turbulence-rans.md.
// No flooring or automatic initialization is performed.
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

// A single accept/reject boundary for the carrier and SST state. This mirrors
// the production outer-iteration transaction contract: U, p, the authoritative
// face mass flux, and all SST-owned fields must either be accepted together or
// restored together. The production loop must explicitly call begin(), then
// commit() or rollback(); this class does not enable SST by itself.
class SSTCoupledStateTransaction {
public:
    SSTCoupledStateTransaction(
        cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
        cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure,
        cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
        SSTProductionState& sst_state)
        : velocity_(velocity), pressure_(pressure),
          mass_flux_(mass_flux), sst_state_(sst_state)
    {
        if (velocity_.dimension() != 3 || pressure_.dimension() != 1 ||
            velocity_.size() == 0 || velocity_.size() != pressure_.size() ||
            mass_flux_.dimension() != 1)
            throw std::invalid_argument(
                "SST coupled transaction received incompatible carrier fields");
        validate_sst_production_state(sst_state_, velocity_.size());
    }

    void begin()
    {
        if (active_)
            throw std::logic_error("SST coupled transaction is already active");
        velocity_snapshot_ = velocity_;
        pressure_snapshot_ = pressure_;
        mass_flux_snapshot_ = mass_flux_;
        sst_snapshot_ = snapshot_sst_production_state(sst_state_);
        active_ = true;
    }

    void rollback()
    {
        if (!active_)
            throw std::logic_error(
                "cannot rollback SST coupled state without an active transaction");
        velocity_ = velocity_snapshot_;
        pressure_ = pressure_snapshot_;
        mass_flux_ = mass_flux_snapshot_;
        restore_sst_production_state(sst_state_, sst_snapshot_);
        active_ = false;
    }

    void commit()
    {
        if (!active_)
            throw std::logic_error(
                "cannot commit SST coupled state without an active transaction");
        active_ = false;
    }

    bool active() const noexcept { return active_; }

private:
    cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity_;
    cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_;
    cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux_;
    SSTProductionState& sst_state_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> velocity_snapshot_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> pressure_snapshot_;
    cfdx::core::Field<double, cfdx::core::Location::FACE> mass_flux_snapshot_;
    SSTProductionState sst_snapshot_{1};
    bool active_ = false;
};

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
    double tolerance = 1e-8,
    const ScalarBoundaryFaceValues* omega_face_values = nullptr)
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
            controls, k_bcs, omega_bcs, max_iterations, tolerance,
            omega_face_values);
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


struct SSTRuntimeBoundaryContract {
    ScalarBoundaryConditions k;
    ScalarBoundaryConditions omega;
    // WALL_OMEGA_SST is retained as an explicit wall-treatment obligation.
    // It must not be coerced to fixed-value or zero-gradient: the current
    // scalar assembler has no SST-specific near-wall omega operator.
    std::vector<std::string> wall_omega_sst_patches;
};

// Translate persisted case metadata only where the runtime scalar assembler
// has an equivalent operator. Wall omega remains explicit and unconsumed until
// a qualified SST wall treatment is implemented.
inline SSTRuntimeBoundaryContract make_sst_runtime_boundary_contract(
    const cfdx::io::CaseSetup& setup)
{
    cfdx::physics::validate_sst_scalar_boundary_contract(setup);
    SSTRuntimeBoundaryContract runtime;
    for (const auto& patch : setup.boundary_conditions) {
        for (const auto& condition : patch.sst_scalar_conditions) {
            auto& target = condition.field == "k" ? runtime.k : runtime.omega;
            switch (condition.type) {
                case cfdx::io::SSTScalarBoundaryType::FIXED_VALUE:
                    target.emplace(patch.patch_name, ScalarBoundaryCondition{
                        ScalarBoundaryType::FIXED_VALUE, condition.value, 0.0});
                    break;
                case cfdx::io::SSTScalarBoundaryType::ZERO_GRADIENT:
                    target.emplace(patch.patch_name, ScalarBoundaryCondition{
                        ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0});
                    break;
                case cfdx::io::SSTScalarBoundaryType::WALL_K_ZERO:
                    if (condition.field != "k")
                        throw std::invalid_argument(
                            "WALL_K_ZERO is only valid for the SST k field");
                    target.emplace(patch.patch_name, ScalarBoundaryCondition{
                        ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0});
                    break;
                case cfdx::io::SSTScalarBoundaryType::WALL_OMEGA_SST:
                    if (condition.field != "omega")
                        throw std::invalid_argument(
                            "WALL_OMEGA_SST is only valid for the SST omega field");
                    runtime.wall_omega_sst_patches.push_back(patch.patch_name);
                    break;
                case cfdx::io::SSTScalarBoundaryType::UNKNOWN:
                    throw std::invalid_argument(
                        "unknown SST scalar boundary type reached runtime adapter");
            }
        }
    }
    return runtime;
}


inline double sst_wall_omega_value(
    double molecular_kinematic_viscosity, double beta1, double wall_normal_distance)
{
    if (!std::isfinite(molecular_kinematic_viscosity) ||
        molecular_kinematic_viscosity <= 0.0 ||
        !std::isfinite(beta1) || beta1 <= 0.0 ||
        !std::isfinite(wall_normal_distance) || wall_normal_distance <= 0.0)
        throw std::invalid_argument(
            "SST wall omega requires finite positive viscosity, beta1, and wall distance");
    const double value = 6.0 * molecular_kinematic_viscosity /
        (beta1 * wall_normal_distance * wall_normal_distance);
    if (!std::isfinite(value) || value <= 0.0)
        throw std::domain_error("SST wall omega value is non-finite or non-positive");
    return value;
}

// Construct per-face omega boundary values for resolved-wall SST treatment.
// The wall distance is the geometric normal distance from the boundary face
// to its owner-cell centre. No characteristic-length fallback is allowed here:
// a bad wall geometry must fail closed instead of silently changing the model.
inline ScalarBoundaryFaceValues make_sst_wall_omega_face_values(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const std::vector<std::string>& wall_omega_sst_patches,
    const TurbulenceTransportControls& controls)
{
    if (controls.model != TurbulenceModel::SST)
        throw std::invalid_argument("SST wall omega operator requires the SST model");
    if (!std::isfinite(controls.molecular_viscosity) ||
        controls.molecular_viscosity <= 0.0 ||
        !std::isfinite(controls.beta1) || controls.beta1 <= 0.0)
        throw std::invalid_argument("SST wall omega operator requires positive viscosity and beta1");
    const std::size_t nf = mesh.n_faces();
    if (geometry.face_centres.size() != nf ||
        geometry.face_area_vectors.size() != nf ||
        geometry.cell_centres.size() != mesh.n_cells() ||
        geometry.face_patch.size() != nf)
        throw std::invalid_argument("SST wall omega geometry shape mismatch");

    ScalarBoundaryFaceValues values;
    const double unset = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, bool> requested;
    for (const auto& name : wall_omega_sst_patches) {
        if (name.empty() || !requested.emplace(name, true).second)
            throw std::invalid_argument("SST wall omega patch names must be non-empty and unique");
        bool found_patch = false;
        for (std::size_t patch = 0; patch < mesh.boundary().n_patches(); ++patch) {
            if (mesh.boundary().patch(patch).name == name) {
                found_patch = true;
                values.values.emplace(name, std::vector<double>(nf, unset));
                break;
            }
        }
        if (!found_patch)
            throw std::invalid_argument("SST wall omega patch is absent from mesh: " + name);
    }

    std::map<std::string, std::size_t> face_counts;
    for (std::size_t face = 0; face < nf; ++face) {
        const std::size_t patch_index = geometry.face_patch[face];
        if (patch_index >= mesh.boundary().n_patches())
            continue;
        const auto& patch_name = mesh.boundary().patch(patch_index).name;
        const auto it = values.values.find(patch_name);
        if (it == values.values.end())
            continue;
        if (mesh.ownership().neighbour(face) >= 0)
            throw std::invalid_argument("SST wall omega patch contains an internal face: " + patch_name);
        const std::size_t owner = mesh.ownership().owner(face);
        if (owner >= mesh.n_cells())
            throw std::invalid_argument("SST wall omega face has invalid owner cell");
        const auto& sf = geometry.face_area_vectors[face];
        const double area = sf.mag();
        if (!(area > 0.0) || !std::isfinite(area))
            throw std::invalid_argument("SST wall omega face has invalid area");
        const auto delta = geometry.face_centres[face] - geometry.cell_centres[owner];
        const double distance = std::abs(delta.dot(sf)) / area;
        it->second[face] = sst_wall_omega_value(
            controls.molecular_viscosity, controls.beta1, distance);
        ++face_counts[patch_name];
    }
    for (const auto& [name, _] : requested) {
        if (face_counts[name] == 0)
            throw std::invalid_argument("SST wall omega patch has no boundary faces in mesh: " + name);
    }
    return values;
}

// Fail closed before transport assembly until every declared wall-omega
// condition has a corresponding production-qualified wall operator.
inline void require_supported_sst_runtime_boundaries(
    const SSTRuntimeBoundaryContract& runtime,
    bool wall_omega_operator_available = false)
{
    if (!runtime.wall_omega_sst_patches.empty() && !wall_omega_operator_available)
        throw std::invalid_argument(
            "SST production transport is blocked: WALL_OMEGA_SST requires the "
            "qualified per-face wall operator");
}

}  // namespace cfdx::physics
