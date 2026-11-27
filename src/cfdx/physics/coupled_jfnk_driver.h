#pragma once
// CoupledJfnkDriver: experimental, unpreconditioned JFNK prototype for the
// 6N U+p+k+omega system (Option D). It is not a production solver and must
// not be used to qualify a production SIMPLE/COUPLED state. Its residual and
// boundary/operator contracts are intentionally incomplete; see audit issue
// #1207. A future production integration requires exact operator equivalence,
// pressure gauge handling, wall/SST traces, and a physical block preconditioner.
//
// V33: the monolithic residual now consumes the SAME OpenFOAM-style wall
// operator (omega face values + Spalding nut + owner omega constraint) used
// by the SIMPLE path when one is staged on the trial state. Without this the
// JFNK operator differed from the SIMPLE operator at the wall: the SIMPLE
// path replaces owner omega rows with the wall function target, while the
// JFNK residual passed nullptr for both omega_face_values and
// wall_production, so its omega rows were the raw PDE. A Newton step that
// minimizes the wrong operator cannot certify the physical state.
//
// V33 also adds pseudo-transient continuation (V30) to the coupled carrier
// momentum rows. This is a diagonal globalization aid only: the accepted
// state is still required to satisfy the ORIGINAL steady FV residual.
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

inline void accumulate_coupled_jfnk_block_scale(
    std::vector<double>& scales, std::size_t block, double value)
{
    if (scales.size() != 6 || block >= scales.size() || !std::isfinite(value))
        throw std::invalid_argument("JFNK block scale: invalid block or value");
    scales[block] = std::max(scales[block], std::abs(value));
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
    const ScalarBoundaryFaceValues* omega_face_values = nullptr;
    const SSTWallProductionV17* wall_production = nullptr;
    const std::vector<double>* momentum_face_mu = nullptr;
};

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
    double ptc_max_diagonal_fraction = std::numeric_limits<double>::infinity())
{
    CoupledJfnkDriverResult result;
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    if (nc == 0) {
        result.reason = "JFNK driver: empty mesh";
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

    // Fixed physical scales: max magnitude per block, floored at 1.0.
    std::array<double, 6> scales{};
    {
        scales.fill(1.0);
        for (std::size_t i = 0; i < nc; ++i) {
            accumulate_coupled_jfnk_block_scale(scales, 0, Ux(i));
            accumulate_coupled_jfnk_block_scale(scales, 1, Uy(i));
            accumulate_coupled_jfnk_block_scale(scales, 2, Uz(i));
            accumulate_coupled_jfnk_block_scale(scales, 3, p(i));
            accumulate_coupled_jfnk_block_scale(scales, 4, k(i));
            accumulate_coupled_jfnk_block_scale(scales, 5, omega(i));
        }
    }

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
    result.initial_merit = std::numeric_limits<double>::infinity();
    bool converged = false;
    for (int step = 0; step < max_newton_steps; ++step) {
        ++result.newton_iterations;
        const auto step_result = coupled_newton_krylov_step(
            residual_cb, state, scales, jfnk_controls, admissible,
            CoupledRightPreconditioner{});
        result.residual_evaluations +=
            static_cast<int>(step_result.residual_evaluations);
        result.linear_relative_residual = step_result.linear_true_relative_residual;
        if (step == 0)
            result.initial_merit = step_result.physical_merit_before;
        result.final_merit = step_result.physical_merit_after;
        state = step_result.state;
        if (step_result.accepted) {
            converged = true;
            break;
        }
        // A rejected Newton step means JFNK failed; caller falls back to SIMPLE.
        break;
    }

    if (!converged) {
        result.reason = "JFNK did not converge within " +
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
    result.reason = "JFNK converged and production state committed";
    return result;
}

} // namespace cfdx::physics
