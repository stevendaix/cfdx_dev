#pragma once

// BFS convergence accelerator: PETSc/OpenFOAM-inspired settings and helpers
// for converging the backward-facing step with SST turbulence.
//
// PETSc-inspired:
//   - NativeFieldSplit with Full Schur for the coupled U-p system
//   - Stronger coupling than segregated SIMPLE for separated flows
//
// OpenFOAM-inspired:
//   - PIMPLE-style multiple pressure correctors with under-relaxation
//   - Adaptive relaxation driven by residual history
//   - Full-order convection for the final operator

#include "cfdx/core/linalg/linear_solver_models.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/solver_control.h"
#include <cmath>

namespace cfdx::physics {

// Recommended controls for a converged BFS/SST calculation. These are
// opt-in recommendations: the caller applies them explicitly. All physical
// convergence tolerances are unchanged from the caller-provided criteria.
inline IncompressibleSolverControls bfs_converged_controls(
    double alpha_u = 0.6, double alpha_p = 0.25,
    std::size_t max_outer = 2000, double linear_tol = 1e-10)
{
    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::COUPLED;
    c.coupling.alpha_u = alpha_u;
    c.coupling.alpha_p = alpha_p;
    c.coupling.coupled_max_iterations = 2000;
    c.coupling.coupled_linear_tolerance = linear_tol;
    c.coupling.schur_model = CoupledSchurModel::PCD;
    c.coupled_linear_solver.krylov = cfdx::core::KrylovModel::FGMRES;
    c.coupled_linear_solver.preconditioner = cfdx::core::PreconditionerModel::PCD;
    c.coupled_linear_solver.gmres_restart = 512;
    c.linear_max_iterations = 500;
    c.linear_tolerance = linear_tol;
    c.convergence.max_iterations = max_outer;
    c.convergence.absolute_tolerance = 1e-10;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-10;
    c.convergence.require_turbulence = true;
    c.require_post_sst_momentum_reassembly = true;
    c.require_global_nonlinear_certificate_v8 = true;
    c.use_bounded_convection = true;
    c.convection_scheme = ConvectionScheme::UPWIND;
    c.adaptive_relaxation.enabled = true;
    c.adaptive_relaxation.min_alpha_u = 0.2;
    c.adaptive_relaxation.max_alpha_u = 0.9;
    c.adaptive_relaxation.min_alpha_p = 0.1;
    c.adaptive_relaxation.max_alpha_p = 0.5;
    c.adaptive_relaxation.improvement_threshold = 0.05;
    c.adaptive_relaxation.degradation_threshold = 0.10;
    c.adaptive_relaxation.increase_step = 0.005;
    c.adaptive_relaxation.decrease_step = 0.02;
    c.adaptive_relaxation.adaptation_window = 4;
    c.adaptive_relaxation.required_trend_windows = 2;
    c.pressure_continuity_line_search.enabled = true;
    c.pressure_continuity_line_search.affine_flux_trials = true;
    return c;
}

// PETSc PCFIELDSPLIT-style request for the coupled solve: velocity block +
// pressure block with Full Schur factorization. Stronger than BlockLocal PCD
// for separated flows where the momentum-pressure feedback loop dominates.
inline cfdx::core::LinearSolverRequest bfs_fieldsplit_request()
{
    cfdx::core::LinearSolverRequest req;
    req.krylov = cfdx::core::KrylovModel::FGMRES;
    req.preconditioner = cfdx::core::PreconditionerModel::NativeFieldSplit;
    req.gmres_restart = 512;
    req.allow_fallback = false;
    return req;
}

// OpenFOAM-style PIMPLE for steady state: multiple pressure correctors per
// outer iteration with under-relaxation, converging the pressure-velocity
// coupling without the SIMPLE one-corrector-per-iteration limit.
inline void bfs_pimple_settings(IncompressibleSolverControls& c,
                                int n_pressure_correctors = 4,
                                double alpha_u = 0.5, double alpha_p = 0.2)
{
    c.algorithm = PressureVelocityAlgorithm::PIMPLE;
    c.coupling.alpha_u = alpha_u;
    c.coupling.alpha_p = alpha_p;
    c.coupling.n_pressure_correctors = n_pressure_correctors;
    c.coupling.piso_neighbor_correction = true;
    c.adaptive_relaxation.enabled = true;
    c.pressure_continuity_line_search.enabled = true;
    c.pressure_continuity_line_search.affine_flux_trials = true;
}

} // namespace cfdx::physics
