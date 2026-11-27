#pragma once
// CoupledJfnkDriver: experimental, unpreconditioned JFNK prototype for the
// 6N U+p+k+omega system (Option D). It is not a production solver and must
// not be used to qualify a production SIMPLE/COUPLED state. Its residual and
// boundary/operator contracts are intentionally incomplete; see audit issue
// #1207. A future production integration requires exact operator equivalence,
// pressure gauge handling, wall/SST traces, and a physical block preconditioner.
//
// The wall-operator and pseudo-transient arguments are retained for API
// compatibility but are NOT consumed by this prototype. Passing either
// feature is rejected below; silently ignoring them would solve a different
// discrete problem from the production SIMPLE path.
#include "cfdx/physics/coupled_jfnk_residual_v6n.h"
#include "cfdx/physics/coupled_newton_krylov.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/sst_wall_operator_v17.h"
#include "cfdx/physics/sst_openfoam_wall_triplet_v15.h"
#include "cfdx/physics/sst_coupled_pseudo_time_v30.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace cfdx::physics {

inline std::vector<double> coupled_jfnk_fixed_scales(
    const PhysicalResidualV51& initial_residual, std::size_t cell_count)
{
    if (cell_count == 0 ||
        cell_count > std::numeric_limits<std::size_t>::max() / 6 ||
        initial_residual.signed_raw.size() != 6 * cell_count ||
        initial_residual.normalized.size() != 6 * cell_count ||
        initial_residual.normalization_scales.size() != 6 * cell_count)
        throw std::invalid_argument(
            "JFNK residual scales: incomplete initial physical residual");
    for (double scale : initial_residual.normalization_scales)
        if (!(scale > 0.0) || !std::isfinite(scale))
            throw std::invalid_argument(
                "JFNK residual scales: nonpositive or nonfinite equation scale");
    return initial_residual.normalization_scales;
}

inline std::vector<double> coupled_jfnk_variable_scales(
    const std::vector<double>& state, std::size_t cell_count)
{
    if (cell_count == 0 ||
        cell_count > std::numeric_limits<std::size_t>::max() / 6 ||
        state.size() != 6 * cell_count)
        throw std::invalid_argument("JFNK variable scales: invalid six-field state shape");
    std::array<double, 6> block_scales{};
    for (std::size_t block = 0; block < 6; ++block) {
        for (std::size_t cell = 0; cell < cell_count; ++cell) {
            const double value = state[block * cell_count + cell];
            if (!std::isfinite(value))
                throw std::invalid_argument("JFNK variable scales: nonfinite state");
            block_scales[block] = std::max(block_scales[block], std::abs(value));
        }
        // A zero field still needs a finite perturbation scale. The unit
        // fallback is used only for the finite-difference coordinate system,
        // never as a residual/acceptance scale.
        if (block_scales[block] == 0.0)
            block_scales[block] = 1.0;
    }
    std::vector<double> scales(6 * cell_count);
    for (std::size_t block = 0; block < 6; ++block)
        std::fill(scales.begin() + static_cast<std::ptrdiff_t>(block * cell_count),
                  scales.begin() + static_cast<std::ptrdiff_t>((block + 1) * cell_count),
                  block_scales[block]);
    return scales;
}

inline bool coupled_jfnk_nonlinear_converged(
    const CoupledNewtonKrylovResult& step_result, double tolerance)
{
    return step_result.accepted && std::isfinite(tolerance) &&
        tolerance > 0.0 && std::isfinite(step_result.physical_merit_after) &&
        step_result.physical_merit_after <= tolerance;
}

struct CoupledJfnkDriverResult {
    bool accepted = false;
    bool production_committed = false;
    int newton_iterations = 0;
    int krylov_iterations = 0;
    int residual_evaluations = 0;
    double initial_merit = std::numeric_limits<double>::infinity();
    double final_merit = std::numeric_limits<double>::infinity();
    double linear_relative_residual = std::numeric_limits<double>::infinity();
    std::string reason;
};

// Optional wall operator staged on the trial state. When present, the
// monolithic residual consumes the SAME omega face values, Spalding nut and
// owner omega constraint used by the SIMPLE path, so the JFNK operator and
// the SIMPLE operator agree at the wall. nullptr leaves the legacy
// wall-free behaviour unchanged.
struct CoupledJfnkWallOperator {
    // Reserved API fields. The current residual does not consume these.
    const ScalarBoundaryFaceValues* omega_face_values = nullptr;
    const SSTWallProductionV17* wall_production = nullptr;
    const std::vector<double>* momentum_face_mu = nullptr;
};

inline bool coupled_jfnk_has_unsupported_extensions(
    const CoupledJfnkWallOperator& wall_op,
    const SSTV30ContinuationState* ptc_state,
    double ptc_max_diagonal_fraction)
{
    return wall_op.omega_face_values != nullptr ||
        wall_op.wall_production != nullptr ||
        wall_op.momentum_face_mu != nullptr ||
        ptc_state != nullptr ||
        ptc_max_diagonal_fraction != std::numeric_limits<double>::infinity();
}

inline CoupledJfnkDriverResult solve_coupled_jfnk_6n(
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& Ux,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& Uy,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& Uz,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& omega,
    SSTProductionState& sst_state,
    const TurbulenceTransportControls& sst_controls,
    const ScalarBoundaryConditions& ubc_x,
    const ScalarBoundaryConditions& ubc_y,
    const ScalarBoundaryConditions& ubc_z,
    const ScalarBoundaryConditions& pbc,
    const ScalarBoundaryConditions& k_bc,
    const ScalarBoundaryConditions& omega_bc,
    const std::map<std::string, std::array<double, 3>>& fixed_velocity_values,
    const std::vector<std::string>& wall_patch_names,
    double density,
    double kinematic_viscosity,
    const CoupledNewtonKrylovControls& jfnk_controls,
    int max_newton_steps = 10,
    const CoupledJfnkWallOperator& wall_op = {},
    SSTV30ContinuationState* ptc_state = nullptr,
    double ptc_max_diagonal_fraction = std::numeric_limits<double>::infinity(),
    double nonlinear_tolerance = 1e-8)
{
    CoupledJfnkDriverResult result;
    using namespace cfdx::core;
    if (!(nonlinear_tolerance > 0.0) || !std::isfinite(nonlinear_tolerance)) {
        result.reason = "JFNK driver: invalid nonlinear convergence tolerance";
        return result;
    }
    const std::size_t nc = mesh.n_cells();
    if (nc == 0) {
        result.reason = "JFNK driver: empty mesh";
        return result;
    }
    if (coupled_jfnk_has_unsupported_extensions(
            wall_op, ptc_state, ptc_max_diagonal_fraction)) {
        result.reason =
            "JFNK prototype does not implement wall operators or pseudo-transient continuation";
        return result;
    }
    const std::size_t n = 6 * nc;

    // Build FvGeometry once; captured by the residual callback.
    const auto geometry = build_fv_geometry(mesh);

    // Pack [Ux, Uy, Uz, p, k, omega] into the monolithic state vector.
    auto pack = [&]() -> std::vector<double> {
        std::vector<double> x(n, 0.0);
        for (std::size_t i = 0; i < nc; ++i) {
            x[0 * nc + i] = Ux(i);
            x[1 * nc + i] = Uy(i);
            x[2 * nc + i] = Uz(i);
            x[3 * nc + i] = p(i);
            x[4 * nc + i] = k(i);
            x[5 * nc + i] = omega(i);
        }
        return x;
    };

    // Unpack the accepted state back into the production fields.
    auto unpack = [&](const std::vector<double>& x) {
        for (std::size_t i = 0; i < nc; ++i) {
            Ux(i) = x[0 * nc + i];
            Uy(i) = x[1 * nc + i];
            Uz(i) = x[2 * nc + i];
            p(i) = x[3 * nc + i];
            k(i) = x[4 * nc + i];
            omega(i) = x[5 * nc + i];
        }
    };

    // Build a 3-component velocity field from the separate components.
    auto build_velocity = [&]() -> Field<double, Location::CELL> {
        Field<double, Location::CELL> U(nc, "U_jfnk", "m/s", 3);
        for (std::size_t i = 0; i < nc; ++i) {
            U.component_data(0)[i] = Ux(i);
            U.component_data(1)[i] = Uy(i);
            U.component_data(2)[i] = Uz(i);
        }
        return U;
    };

    // Freeze equation-specific residual scales from the initial accepted
    // physical state. The same divisors are used by Jv and every line-search
    // candidate; no state-dependent renormalization enters the Newton merit.
    const auto initial_physics = CoupledJfnkResidualV6N::evaluate(
        mesh, geometry, Ux, Uy, Uz, p, k, omega,
        sst_controls, ubc_x, ubc_y, ubc_z, pbc, k_bc, omega_bc,
        fixed_velocity_values, wall_patch_names, density, kinematic_viscosity);
    const auto scales = coupled_jfnk_fixed_scales(initial_physics, nc);

    // This experimental driver is deliberately unpreconditioned.
    // Wrapping identity operators in SixFieldBlockPreconditionerV8 does not
    // create a physical block preconditioner and would misrepresent the solver.
    // A real block action must be derived from the coupled momentum/pressure
    // and SST operators before this driver can be qualified for production.

    // Residual callback: unpack the trial state, assemble the full 6N
    // residual from scratch, return signed_raw.
    const auto residual_cb = [&](const std::vector<double>& x) -> std::vector<double> {
        for (std::size_t i = 0; i < n; ++i)
            if (!std::isfinite(x[i]))
                throw std::runtime_error("JFNK residual: nonfinite trial state");
        Field<double, Location::CELL> Ux_t(nc, "Ux_t", "m/s", 1);
        Field<double, Location::CELL> Uy_t(nc, "Uy_t", "m/s", 1);
        Field<double, Location::CELL> Uz_t(nc, "Uz_t", "m/s", 1);
        Field<double, Location::CELL> p_t(nc, "p_t", "Pa", 1);
        Field<double, Location::CELL> k_t(nc, "k_t", "m2/s2", 1);
        Field<double, Location::CELL> w_t(nc, "omega_t", "1/s", 1);
        for (std::size_t i = 0; i < nc; ++i) {
            Ux_t(i) = x[0 * nc + i];
            Uy_t(i) = x[1 * nc + i];
            Uz_t(i) = x[2 * nc + i];
            p_t(i) = x[3 * nc + i];
            k_t(i) = x[4 * nc + i];
            w_t(i) = x[5 * nc + i];
        }
        const auto phys = CoupledJfnkResidualV6N::evaluate(
            mesh, geometry,
            Ux_t, Uy_t, Uz_t, p_t, k_t, w_t,
            sst_controls, ubc_x, ubc_y, ubc_z, pbc, k_bc, omega_bc,
            fixed_velocity_values, wall_patch_names,
            density, kinematic_viscosity);
        return phys.signed_raw;
    };

    // Admissibility: k and omega must remain positive.
    const auto admissible = [&](const std::vector<double>& x) -> bool {
        for (std::size_t i = 0; i < nc; ++i) {
            if (x[4 * nc + i] < 0.0 || x[5 * nc + i] <= 0.0)
                return false;
        }
        return true;
    };

    // Run JFNK iterations from the current state.
    std::vector<double> state = pack();
    const auto variable_scales = coupled_jfnk_variable_scales(state, nc);
    result.initial_merit = std::numeric_limits<double>::infinity();
    bool converged = false;
    for (int step = 0; step < max_newton_steps; ++step) {
        ++result.newton_iterations;
        const auto step_result = coupled_newton_krylov_step(
            residual_cb, state, scales, jfnk_controls, admissible,
            CoupledRightPreconditioner{}, variable_scales);
        result.residual_evaluations +=
            static_cast<int>(step_result.residual_evaluations);
        result.linear_relative_residual = step_result.linear_true_relative_residual;
        if (step == 0)
            result.initial_merit = step_result.physical_merit_before;
        result.final_merit = step_result.physical_merit_after;
        if (!step_result.accepted) {
            // A rejected Newton step means JFNK failed; caller falls back to SIMPLE.
            break;
        }
        state = step_result.state;
        if (coupled_jfnk_nonlinear_converged(step_result, nonlinear_tolerance)) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        result.reason = "JFNK failed to reach nonlinear tolerance " +
            std::to_string(nonlinear_tolerance) + " within " +
            std::to_string(max_newton_steps) + " Newton steps";
        return result;
    }

    // Acceptance: unpack back to production fields and refresh SST derived
    // state (wall distance, strain rate, F2, nu_t) from the accepted carrier.
    unpack(state);
    const auto U = build_velocity();
    refresh_sst_production_derived_state_from_mesh(
        sst_state, mesh, sst_controls, U, fixed_velocity_values, wall_patch_names);
    result.accepted = true;
    result.production_committed = true;
    result.reason = "JFNK reached the nonlinear residual tolerance and committed production state";
    return result;
}

} // namespace cfdx::physics
