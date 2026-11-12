#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/cell_overlapping_schwarz.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/linear_solver_context.h"
#include "cfdx/core/linalg/mgr_preconditioner.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/pressure_roundoff_leakage.h"
#include "cfdx/core/linalg/pressure_compensated_refinement.h"
#include "cfdx/core/linalg/pressure_conservative_face_operator.h"
#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/physics/initial_face_flux_contract.h"
#include "cfdx/physics/advanced_convergence.h"
#include "cfdx/physics/picard_physical_globalization.h"
#include "cfdx/physics/tvd_deferred_defect_memory.h"
#include "cfdx/physics/tvd_secant_implicit_correction.h"
#include "cfdx/physics/nonlinear_multisecant_acceleration.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/cell_convergence_tracker.h"
#include "cfdx/physics/momentum_viscous_stress.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/pressure_compensated_assembly.h"
#include "cfdx/physics/pcd_fp_equilibration.h"
#include "cfdx/physics/pcd_fp_strong_component_anchor.h"
#include "cfdx/physics/pcd_fp_capture.h"
#include "cfdx/physics/pcd_pressure_operator_capture_v38.h"
#include "cfdx/physics/pcd_schur_row_audit_v40.h"
#include "cfdx/physics/pcd_fp_bisided_basis_v13.h"
#include "cfdx/physics/pcd_fp_deflated_quality_v14.h"
#include "cfdx/physics/pcd_fp_basis_registry_v15.h"
#include "cfdx/physics/pcd_fp_refinement_v16.h"
#include "cfdx/physics/coupled_block_residual_v16.h"
#include "cfdx/physics/coupled_monolithic_ilut_v16.h"
#include "cfdx/physics/pcd_blocklocal_recovery_v42.h"
#include "cfdx/physics/pressure_velocity_system.h"
#include "cfdx/physics/pressure_continuity_line_search.h"
#include "cfdx/physics/convection_continuation.h"
#include "cfdx/physics/pressure_affine_flux_trials.h"
#include "cfdx/physics/solver_control.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/sst_outer_ptc_v31.h"
#include "cfdx/physics/sst_outer_ptc_balance_v32.h"
#include "cfdx/physics/sst_openfoam_wall_audit_v7.h"
#include "cfdx/physics/sst_openfoam_wall_coupling_v8.h"
#include "cfdx/physics/sst_openfoam_wall_triplet_v15.h"
#include "cfdx/physics/sst_wall_operator_v17.h"
#include "cfdx/physics/momentum_residual_certification.h"
#include "cfdx/physics/nonlinear_face_residual_certificate.h"
#include "cfdx/physics/global_nonlinear_residual_v8.h"
#include "cfdx/physics/global_nonlinear_acceptance_v13.h"
#include "cfdx/physics/global_nonlinear_contraction_v14.h"
#include "cfdx/physics/outer_field_period_v43.h"
#include "cfdx/physics/coupled_outer_step_v43.h"
#include "cfdx/physics/coupled_phase_escape_v45.h"
#include "cfdx/physics/post_sst_momentum_hotspot_v45.h"
#include "cfdx/physics/post_sst_face_forensics_v46.h"
#include "cfdx/physics/coupled_full_period_v44.h"
#include "cfdx/physics/coupled_linear_retry_policy_v37.h"
#include "cfdx/physics/momentum_face_budget.h"
#include "cfdx/physics/momentum_same_state_audit.h"
#include "cfdx/physics/timestep_control.h"
#include "cfdx/core/numerics/conservation.h"
#include "cfdx/core/numerics/temporal.h"
#include "cfdx/cfdx/core/utils/convergence_monitor.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <numeric>
#include <stdexcept>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <functional>
#include <fstream>
#include <cstdlib>

namespace cfdx::physics {

struct VelocityBoundaryCondition {
    enum class Type { FIXED_VALUE, ZERO_GRADIENT };
    Type type = Type::ZERO_GRADIENT;
    cfdx::core::Vec3 value{0.0, 0.0, 0.0};
};

using VelocityBoundaryConditions = std::map<std::string, VelocityBoundaryCondition>;
// Resolved from the authoritative preceding-iterate phi. The same three
// component BCs are shared by the momentum matrix, face transport oracle and
// the Rhie-Chow H/A predictor. A per-patch zero-gradient approximation is NOT
// sufficient when a nominal outlet has local inflow.
using VelocityFaceConditions = std::array<ScalarBoundaryFaceConditions, 3>;

// A Dirichlet entry removes the pressure nullspace only when it is actually
// attached to a non-impermeable boundary face. An unused/misspelled entry (or
// one on EMPTY/SYMMETRY) must not suppress gauge fixing and silently create
// a singular pressure system. Shared by SIMPLE, COUPLED, and nested PCD.
inline bool has_effective_fixed_pressure_boundary(
    const cfdx::core::Mesh& mesh,
    const ScalarBoundaryConditions& pressure_bcs)
{
    bool found = false;
    for (const auto& [name, bc] : pressure_bcs) {
        if (bc.type != ScalarBoundaryType::FIXED_VALUE) continue;
        const std::size_t id = mesh.boundary().find(name);
        if (id >= mesh.boundary().n_patches())
            throw std::invalid_argument(
                "fixed pressure boundary references unknown patch: " + name);
        const auto& patch = mesh.boundary().patch(id);
        if (patch.empty() || patch.type == cfdx::core::PatchType::EMPTY ||
            patch.type == cfdx::core::PatchType::SYMMETRY)
            throw std::invalid_argument(
                "fixed pressure boundary has no effective Dirichlet face: " + name);
        found = true;
    }
    return found;
}

enum class IncompressibleProbeField {
    PRESSURE,
    U_X,
    U_Y,
    U_Z,
    U_MAGNITUDE,
};

struct IncompressiblePointProbe {
    std::string name;
    cfdx::core::Vec3 location{0.0, 0.0, 0.0};
    IncompressibleProbeField field = IncompressibleProbeField::PRESSURE;
};

struct IncompressibleProbeSample {
    std::string name;
    std::size_t iteration = 0;
    double time = 0.0;
    double value = 0.0;
};

struct DiagnosticsControls {
    bool coupled_matrix_summary = false;
    bool freeze_state_probe = false;
    long debug_cell = -1;
    bool debug_cell_auto = false;
    // Bounded progress trace for long validation solves. Disabled by default.
    bool iteration_trace = false;
    std::size_t iteration_trace_frequency = 50;
    // Accepted-state spatial diagnostics only; empty path keeps allocation disabled.
    std::string cell_convergence_trace_path;
    std::string cell_convergence_mesh_source;
    std::size_t cell_convergence_every = 5;
    std::string cell_convergence_full_map_path;
    // V17: topology/patch mass audit on the AUTHORITATIVE Rhie-Chow flux.
    // Diagnostic only; explicitly enabled in case settings.
    bool flux_patch_audit = false;
};

struct V17PseudoMomentumControls {
    // Steady momentum globalization, NOT a physical transient time-step.
    // Adds M_P (U-U_outer_old) to the integrated momentum equation only.
    bool enabled = false;
    double initial_cfl = 0.3;
    double minimum_cfl = 0.03;
    double maximum_cfl = 3.0;
    double growth = 1.25;
    double shrink = 0.5;
    // Optional safety cap: M_p <= fraction * min(A_uu, A_vv, A_ww).
    // The original uncapped scheme is preserved when this is infinite.
    // A tiny pseudo-CFL on a stretched mesh can otherwise make M_p orders
    // of magnitude larger than the *physical* predictor diagonal, freezing
    // the nonlinear flow while the pressure equation appears solved.
    double maximum_diagonal_fraction = std::numeric_limits<double>::infinity();
};

inline void validate_v17_pseudo_momentum_controls(const V17PseudoMomentumControls& c) {
    if (!std::isfinite(c.initial_cfl) || !std::isfinite(c.minimum_cfl) ||
        !std::isfinite(c.maximum_cfl) || !std::isfinite(c.growth) ||
        !std::isfinite(c.shrink) || !(c.minimum_cfl > 0.0) ||
        !(c.minimum_cfl <= c.initial_cfl) ||
        !(c.initial_cfl <= c.maximum_cfl) || !(c.growth >= 1.0) ||
        !(c.shrink > 0.0 && c.shrink <= 1.0) ||
        !(c.maximum_diagonal_fraction > 0.0) ||
        std::isnan(c.maximum_diagonal_fraction))
        throw std::invalid_argument("V17 pseudo-momentum controls: invalid CFL/diagonal cap");
}

struct IncompressibleTransientControls {
    bool enabled = false;
    cfdx::core::TimeScheme scheme = cfdx::core::TimeScheme::BDF2;
    double dt = 0.0;
    double dt_previous = 0.0;
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* previous = nullptr;
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* previous_previous = nullptr;
    bool history_valid = false;
};

// ---------------------------------------------------------------------------
// Initialization contract (N7)
// ---------------------------------------------------------------------------
// The initial state is part of the numerical method: it fixes the denominator
// of every relative residual reported by the solve, because
// result.reference_momentum_residual is taken from the first iteration. The
// strategy is therefore declared explicitly and validated rather than inferred
// from whether a restart path happens to be non-empty.
//
// Automatic selection is deliberately NOT implemented. Choosing a strategy by
// probing the filesystem or mesh quality would change the initial state, and
// therefore the relative-residual reference, without the resulting convergence
// being attributable to a declared numerical method.
enum class InitializationMode {
    // Caller-provided fields are used unchanged. The caller is responsible for
    // having produced a meaningful state.
    Provided,
    // Fields are reset to a uniform state before the solve.
    Uniform,
    // Fields are overwritten from a CFDX-DAT checkpoint. Requires
    // controls.initialization.restart_path to be non-empty.
    Restart
};

inline const char* to_string(InitializationMode mode) {
    switch (mode) {
        case InitializationMode::Provided: return "provided";
        case InitializationMode::Uniform: return "uniform";
        case InitializationMode::Restart: return "restart";
    }
    return "unknown";
}

// Parses an initialization mode name. Returns false for an unrecognised name so
// that callers reject it instead of silently defaulting.
inline bool parse_initialization_mode(const std::string& name, InitializationMode& out) {
    if (name == "provided") { out = InitializationMode::Provided; return true; }
    if (name == "uniform") { out = InitializationMode::Uniform; return true; }
    if (name == "restart") { out = InitializationMode::Restart; return true; }
    return false;
}

struct InitializationControls {
    InitializationMode mode = InitializationMode::Provided;
    // Uniform state applied when mode == Uniform.
    cfdx::core::Vec3 uniform_velocity{0.0, 0.0, 0.0};
    double uniform_pressure = 0.0;
    // Checkpoint path applied when mode == Restart.
    std::string restart_path;
};

// Read-only view of the exact coupled saddle-point system handed to the
// Krylov solver, exposed to the optional evidence hook in
// IncompressibleSolverControls. The hook fires at the Krylov call point -
// after assembly, gauge treatment, finite checks and preconditioner setup -
// so the aliased matrix, RHS and initial guess ARE the exact solve arguments
// (same objects, no copy, no further mutation). The references are only
// valid for the duration of the callback: a consumer that needs the data
// beyond the call must copy it.
struct CoupledAssemblySnapshot {
    const cfdx::core::SparseMatrix* matrix = nullptr;
    const cfdx::core::Vector* rhs = nullptr;
    // Production initial guess ([U_old; p_old] with the reference pressure
    // component initialized to the configured value), i.e. x passed to Krylov.
    const cfdx::core::Vector* initial_guess = nullptr;
    // Momentum-weighted Rhie-Chow pressure Schur diagonal, in both forms:
    // gauge-adjusted (reference cell = 1.0 only when the pure-Neumann gauge
    // row is replaced) and raw (as assembled, before gauge treatment). With a
    // fixed-pressure boundary these arrays must be identical. Diagnostic
    // evidence channels; production preconditioning consumes the adjusted form.
    const std::vector<double>* pressure_schur_diagonal = nullptr;
    const std::vector<double>* pressure_schur_diagonal_raw = nullptr;
    // Production PCD pressure operators of this assembly, exactly as passed
    // to PcdSchurApproximation after the applicable boundary/gauge treatment.
    // Null when the resolved coupled plan is not PCD.
    const cfdx::core::SparseMatrix* pcd_pressure_mass = nullptr;
    const cfdx::core::SparseMatrix* pcd_pressure_laplacian = nullptr;
    const cfdx::core::SparseMatrix* pcd_pressure_convection_diffusion = nullptr;
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* velocity_old = nullptr;
    std::size_t reference_cell = 0;
    bool has_fixed_pressure = false;
    // Optional, read-only momentum and frozen-pressure evidence.  Valid only
    // inside coupled_assembly_callback; an external audit must not retain it.
    // These are the very scalar equations handed to the monolithic assembler.
    std::array<const ScalarEquation*, 3> momentum_equations{{nullptr, nullptr, nullptr}};
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* pressure_old = nullptr;
};

struct SSTProductionRuntime {
    // State must be explicitly initialized from case data or restored from a
    // DAT checkpoint with the helpers in sst_production_coupling.h. The
    // production solver advances this state transactionally when non-null.
    // In this incompressible SST path the carrier pressure and pressure BCs
    // represent p*=p_static+2*rho*k/3. The isotropic Reynolds stress is thus
    // absorbed in the pressure gradient. Convert explicitly before reporting
    // static pressure/Cp or prescribing static-pressure boundary data.
    SSTProductionState* state = nullptr;
    TurbulenceTransportControls transport;
    SSTV30ContinuationState coupled_ptc_v30;
    SSTV31OuterPTC outer_ptc_v31;
    SSTV32CouplingBudget outer_ptc_balance_v32;
    ScalarBoundaryConditions k_boundary_conditions;
    ScalarBoundaryConditions omega_boundary_conditions;
    // Optional per-face omega data, including the existing geometric wall-omega
    // contract. Kept owned by the runtime so the pointer passed to transport
    // remains valid for the duration of the update.
    ScalarBoundaryFaceValues omega_face_values;
    std::vector<std::string> wall_patch_names;
    // Physical WALL_OMEGA_SST obligations from the case contract. This is
    // distinct from the list of walls used to compute nearest wall distance.
    std::vector<std::string> wall_omega_sst_patches;
    std::size_t max_transport_iterations = 100;
    double transport_tolerance = 1e-8;
    // Experimental outer-coupling schedule. Finite and unclipped intermediate
    // SST predictors may advance; NOT a relaxation of the final SST gate.
    bool allow_bounded_provisional_transport = false;
    // V15 staged reference. Refuses incomplete per-face/cell values, promoted
    // only on a fully accepted carrier+SST outer transaction. NO physical
    // wall-function substitution; all three production operators must still
    // be integrated/qualified before enabling such a substitution.
    // Experimental opt-in: consume per-face nut in momentum and dynamic
    // omega/G wall fields in SST. This is not verified OpenFOAM equivalence.
    bool active_openfoam_wall_operator_v17 = false;
    // V18 additional, disabled by default: replace omega rows at owners,
    // including area-weighted shared/corner cells. The boundary face-only
    // V17 branch is otherwise unchanged. Not OpenFOAM equivalent by itself.
    bool omega_wall_owner_constraints_v18 = false;
    // V19 explicit OpenFOAM Foundation corner policy: equal contributions
    // per wall FACE. V15/V17 retain their area-weighted policy otherwise.
    bool omega_wall_corner_face_count_v19 = false;
    bool stage_openfoam_wall_triplet_v15 = false;
    std::optional<SstOpenFoamWallTripletV15> committed_wall_triplet_v15;
};

struct IncompressibleSolverControls {
    InitializationControls initialization;
    // Opt-in only. The production RANS selector remains fail-closed; supplying
    // this runtime explicitly opts the caller into the SST-coupled path.
    SSTProductionRuntime* sst_runtime = nullptr;
    // V7 strict opt-in: recompute the physical carrier momentum residual
    // AFTER the accepted candidate SST transport changed nu_t. A pressure
    // solve using the previous turbulent viscosity cannot certify this state.
    bool require_post_sst_momentum_reassembly = false;
    // V8 opt-in physical residual certificate shared by SIMPLE, SIMPLEC,
    // COUPLED BlockLocal/PCD/LSC/BFBT/SIMPLE/SIMPLEC. Do not inspect a
    // preconditioner residual to certify the actual U-p-k-omega state.
    // Reassembles momentum after SST and both SST equations on the returned
    // state; all existing acceptance gates remain mandatory.
    bool require_global_nonlinear_certificate_v8 = false;
    // V13 opt-in: reject a complete U,p,k,omega trial when a reassembled
    // physical component grows severely relative to the last accepted state.
    // Requires a steady SST solve, V8 and post-SST residual reconstruction.
    // Convergence thresholds are NOT changed.
    bool joint_globalization_v13_experimental = false;
    // V14 stronger opt-in: require a strict decrease of the independent
    // post-SST global merit between COMMITTED states at identical operator.
    // Never accept by the numerical linear solver's reported residual.
    bool monotone_global_merit_v14_experimental = false;
    // V43 experimental nonlinear correction for the monolithic COUPLED
    // carrier update. 0 disables it; 0<beta<=1 interpolates U,p between the
    // previous accepted state and the certified linear candidate. Retries use
    // the shared nonlinear retry shrink; the authoritative Rhie-Chow flux is
    // recomputed at the interpolated state (never blended as an independent
    // face field). Physical residuals are reassembled AFTER SST as usual.
    // The damped candidate is not claimed to solve the frozen block A*x=b.
    double coupled_outer_step_v43_experimental = 0.0;
    // V43 separate opt-in fixed-point composite: when all committed fields
    // return to their state n-4 for eight consecutive measurements, seed the
    // NEXT full coupled/SST nonlinear trial from the four-phase mean. The
    // average is an unaccepted predictor only; it is never exported as a
    // converged flow or used as the authoritative final Rhie-Chow flux.
    bool coupled_period_four_mean_v43_experimental = false;
    // V44 experimental six-field version of the phase predictor. Unlike V43,
    // this averages k and omega as well as U/p/phi, then RECOMPUTES nu_t from
    // the SST closure for the predicted U/k/omega. The output is an initial
    // guess, never a nonlinear certificate. Disabled by default.
    bool coupled_period_four_full_state_v44_experimental = false;
    // V45 opt-in continuation, only AFTER a four-cycle of committed full
    // fields is detected by V43 and V44's coherent phase predictor is used.
    // For the following accepted steps, damp the COUPLED carrier candidate
    // while SST is re-solved and the existing V7/V8 merit is reassembled.
    // Zero disables; default equations, physical criteria and Schur untouched.
    double coupled_period_four_damping_v45_experimental = 0.0;
    std::size_t coupled_period_four_damping_steps_v45 = 12;
    // V31: one CFL for momentum and SST, updated ONLY on committed full
    // U,p,k,omega physical residuals. Pressure remains algebraic.
    bool outer_ptc_v31_experimental = false;
    // V32: cap *both* momentum and SST pseudo masses by a relative FV
    // diagonal budget; failed retries reduce only the unpublished trial.
    bool outer_ptc_balance_v32_experimental = false;
    // Separate controlled A/B: share CFL and the accepted outer ledger,
    // but distinguish which physical operator needs the pseudo damping.
    bool outer_ptc_balance_momentum_v32 = true;
    bool outer_ptc_balance_sst_v32 = false;
    PressureVelocityAlgorithm algorithm = PressureVelocityAlgorithm::SIMPLE;
    IncompressibleTransientControls transient;
    V17PseudoMomentumControls pseudo_momentum;
    CouplingControls coupling;
    ConvergenceCriteria convergence;
    // Opt-in for verification/qualification: require an independently
    // reassembled local face momentum balance on the accepted state before
    // declaring nonlinear convergence. Shared by SIMPLE/SIMPLEC/COUPLED;
    // disabled by default to preserve existing solver behaviour.
    bool require_independent_face_residual_for_convergence = false;
    // Optional fail-closed independent SIMPLEC rAtU-to-rAU face identity.
    // Keep separate from the nonlinear physical conservation gates.
    bool require_simplec_same_state_flux_identity = false;
    // Opt-in for a long validation run: STAGNATED remains diagnostic and
    // MUST NOT count as converged. Divergence, maximum-iteration, component
    // residual, conservation and true-face momentum gates remain unchanged.
    bool keep_iterating_after_stagnation = false;
    std::size_t linear_max_iterations = 1000;
    double linear_tolerance = 1e-10;
    // Optional independent pressure linear budget; zero inherits the global one.
    // The nonlinear mass and momentum conservation gates are unaffected.
    double pressure_linear_tolerance = 0.0;
    std::size_t pressure_linear_max_iterations = 0;
    // Experimental CPU-side assembly: accumulate the same orthogonal face
    // coefficients and RHS contributions in extended precision, then round
    // each final CSR diagonal/RHS once. OFF by default. This changes neither
    // FV equations nor linear/nonlinear convergence gates, but may change
    // rounding and thus requires independent NACA/BFS qualification.
    // It does NOT carry a two-part pressure field into Rhie-Chow fluxes.
    bool pressure_compensated_assembly = false;
    // Validation-only replay on the ACTUAL pressure matrix, after an FP64
    // linear failure. Native AMG/CG plus a high/low correction can certify
    // the algebraic CSR without declaring a full FV pressure correction.
    // It does not overwrite the physical pressure, fluxes or solver status.
    bool pressure_compensated_refinement_diagnostic = false;
    // V6 experimental actual SIMPLE production path. The full nonorthogonal
    // Rhie-Chow face operator is the physical pressure residual; V18 sparse
    // CSR serves only as an inexact Jacobian/preconditioner. Carry a high/low
    // cell pressure through successive accepted iterations, use it in cell
    // velocity gradients and full-step transport flux. Never enable by default.
    bool pressure_high_low_production = false;
    // V7: exact split-pressure restart requires DAT6 U/p/p_low/phi on the
    // same ordered topology. A DAT5 or unpaired flux is refused outright.
    ConvergenceAccelerationControls acceleration;
    double density = 1.0;
    double kinematic_viscosity = 1.0e-3;
    double turbulent_viscosity = 0.0;
    // Optional cellwise kinematic eddy viscosity (m^2/s). Replaces the scalar
    // value when supplied; the production RANS guard remains active.
    std::vector<double> cell_turbulent_viscosity;
    // Explicit V23 hybrid startup: a projected FACE mass flux produced on the
    // same mesh. Not a restart and never re-created from a relaxed pressure.
    // The conservative flux is only a startup state; all final residual and
    // physical gates are unchanged. Empty means legacy interpolation.
    std::vector<double> initial_face_mass_flux;
    // A complete same-mesh checkpoint captured at an accepted iteration may
    // legitimately contain a NONCONVERGED continuity defect and a Rhie-Chow
    // boundary pressure correction. This is not a physical convergence waiver:
    // it only changes the admission test for a verified restart initial guess.
    bool initial_face_mass_flux_from_paired_checkpoint = false;
    // Original accepted-run normalization and offset for a verifiable restart.
    // BOTH references must accompany a same-mesh U,p,phi checkpoint.
    double initial_momentum_residual_reference = 0.0;
    double initial_nonlinear_metric_reference = 0.0;
    std::size_t initial_accepted_iteration = 0;
    // Explicit opt-in for a physical free-stream backflow policy. The listed
    // patches use zeroGradient U on outward flux and the specified U on inward
    // flux, decided from the PREVIOUS accepted Rhie-Chow phi. This remains a
    // Picard-frozen BC, and the final residual is evaluated with fresh signs.
    // Limited to laminar segregated steady calculations until SST/COUPLED and
    // pressure-mixed interface contracts have independent qualifications.
    std::map<std::string, cfdx::core::Vec3> inlet_outlet_velocity_patches;
    cfdx::core::Vec3 body_force{0.0, 0.0, 0.0};
    // Optional cellwise body force per unit mass (m/s^2). Replaces the
    // uniform controls.body_force vector when supplied (Boussinesq
    // buoyancy coupling, issue #885 step 4, fiche methode 10): the
    // assembly multiplies it by the constant density exactly like the
    // uniform channel, in both the momentum predictor and the
    // independent momentum conservation audit. Supplying a nonzero
    // uniform body_force together with the cellwise channel is rejected,
    // mirroring the cellwise turbulent-viscosity contract.
    std::vector<cfdx::core::Vec3> cell_body_force;
    PressureGaugePolicy pressure_gauge_policy = PressureGaugePolicy::REFERENCE_CELL;
    std::size_t pressure_reference_cell = 0;
    double pressure_reference_value = 0.0;
    cfdx::core::LinearSolverRequest momentum_linear_solver;
    cfdx::core::LinearSolverRequest pressure_linear_solver;
    // Explicit pressure-matrix forensic capture: select at most 32 outer
    // iterations (e.g. {20,21,600} for N256). One original post-gauge CSR
    // system per selected iteration and pressure corrector; retry on the same
    // iteration may overwrite its file. No I/O unless configured.
    std::string pressure_replay_prefix;
    std::vector<std::size_t> pressure_replay_iterations;
    cfdx::core::LinearSolverRequest coupled_linear_solver;
    // PCD pressure-side solves are nested inside the coupled preconditioner.
    // They deliberately have independent solver controls so the outer coupled
    // iteration budget can never become an accidental inner Krylov budget.
    cfdx::core::LinearSolverRequest pcd_laplacian_linear_solver{
        cfdx::core::KrylovModel::GMRES,
        cfdx::core::PreconditionerModel::NativeAMG,
        40,
        false,
        cfdx::core::NullSpaceModel::None,
        cfdx::core::LinearScalingModel::None,
        0,
        false,
        {},
        {}};
    // Fp is the nested PCD pressure convection-diffusion operator. The
    // Ghia Re=100 qualification showed that FGMRES(40) repeatedly hit the
    // 100-iteration cap with true relative residuals above the unchanged
    // 1e-8 gate; an isolated FP64 experiment converged in roughly 53-56
    // iterations with restart 512. Give this inner solve enough Krylov space
    // to avoid the measured restart plateau. The iteration budget and
    // tolerance remain unchanged; the effective restart is clamped to n_rows
    // at the call site.
    cfdx::core::LinearSolverRequest pcd_convection_diffusion_linear_solver{
        cfdx::core::KrylovModel::FGMRES,
        cfdx::core::PreconditionerModel::ILU0,
        512,
        false,
        cfdx::core::NullSpaceModel::None,
        cfdx::core::LinearScalingModel::None,
        0,
        false,
        {},
        {}};
    std::size_t pcd_laplacian_max_iterations = 100;
    std::size_t pcd_convection_diffusion_max_iterations = 100;
    double pcd_laplacian_tolerance = 1e-8;
    double pcd_convection_diffusion_tolerance = 1e-8;
    bool use_bounded_convection = true;
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND;
    // Opt-in: conservative continuation from UPWIND to the requested full
    // high-order scheme; acceptance requires full order, never a ramp stage.
    ConvectionContinuationControls convection_continuation;
    DiagnosticsControls diagnostics;
    NonlinearRetryControls nonlinear_retry;
    // Opt-in nonlinear globalization/damping aids. Each is validated below and
    // disabled by default; none can relax a physical convergence threshold.
    PicardPhysicalGlobalizationControls picard_physical_globalization;
    TvdDeferredDefectMemoryControls tvd_deferred_memory;
    TvdSecantImplicitControls tvd_secant_implicit;
    // Opt-in Anderson/NGMRES predictor; the accelerated candidate is never accepted without SIMPLE reprojection.
    NonlinearMultisecantControls nonlinear_multisecant;
    std::vector<IncompressiblePointProbe> probes;
    // Optional quantities of interest gated by the convergence monitor. Each
    // gate must name a configured probe; a gate that names nothing is a
    // misconfiguration and is rejected rather than silently never satisfied.
    std::vector<cfdx::core::QoIGate> qoi_gates;
    AdaptiveRelaxationControls adaptive_relaxation;
    // Opt-in: minimise the true conservative face-flux continuity after each
    // pressure update. Defaults unchanged; requires physical A/B qualification.
    PressureContinuityLineSearchControls pressure_continuity_line_search;
    // Experimental: skew-corrected Green-Gauss pressure gradient used by
    // momentum, pressure assembly AND authoritative Rhie-Chow reconstruction.
    // Keep off until three-grid NACA and N9 physical qualification pass.
    bool pressure_skew_corrected_gradient = false;
    // Experimental: use the actual projected face position on an elongated
    // owner-neighbour pair for BOTH the predictor and pressure Schur operator.
    // Segregated SIMPLE only until the coupled Jacobian is equivalently updated.
    bool pressure_projected_face_interpolation = false;
    // V18 diagnostic ONLY: assemble the complete face-flux Jacobian by
    // pressure basis vectors. O(N*Nfaces), max 1024 cells, SIMPLE/fixed-p only.
    // Never promoted automatically and never changes physical convergence gates.
    bool v18_full_flux_pressure_oracle = false;
    // V18 scalable assembly of the same Jacobian from local Gauss stencils.
    // Still experimental; no automatic physical qualification.
    bool v18_sparse_flux_pressure_jacobian = false;
    // V19 experimental SIMPLE flux policy: compute conservative phi using
    // the FULL pressure correction, while pressure used for velocity and
    // gradients remains under-relaxed. Mirrors the SIMPLE flux/pressure
    // separation in OpenFOAM. The face flux is an independent state and
    // deliberately need not equal RC(HbyA, under-relaxed cell pressure).
    // Requires exact V18 pressure/flux Jacobian, single SIMPLE corrector.
    bool v19_full_pressure_flux_step = false;
    std::function<void(const IncompressibleProbeSample&)> probe_callback;
    // Called after each completed nonlinear iteration with the authoritative solver state.
    // Return false to stop the solve after the current iteration.
    using IterationOutputCallback = std::function<bool(
        std::size_t,
        double,
        const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>;
    IterationOutputCallback iteration_output_callback;
    // Optional PURE candidate-validation gate, evaluated BEFORE the nonlinear
    // U,p,phi transaction commits. Throw NonlinearRetryableFailure to reject
    // the whole tentative step; this is not an output/stop callback. Do not
    // perform disk I/O or irreversible side effects here; use the accepted
    // iteration callback for checkpoints, or iteration_output_callback for
    // post-commit visualization/output.
    std::function<void(std::size_t, double, const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>
        tentative_iteration_gate_callback;
    // Called only after a nonlinear U,p,phi transaction has committed.
    // Unlike diagnostic assembly hooks, this cannot expose a rejected candidate.
    // phi is the authoritative accepted Rhie-Chow face flux (kg/s), not
    // a reconstruction from U. This hook is for reproducible restart snapshots.
    using AcceptedIterationCallback = std::function<void(
        std::size_t,
        const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::FACE>&,
        double, double)>;
    AcceptedIterationCallback accepted_iteration_callback;
    // V7 opt-in checkpoint observer: only a COMMITTED state, with the exact
    // low pressure component used by velocity/RC. Allows callers to write
    // DAT6 U/p/p_low/phi without rebuilding any of these from U alone.
    using AcceptedCompensatedIterationCallback = std::function<void(
        std::size_t,
        const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::FACE>&,
        const cfdx::core::Vector&)>;
    AcceptedCompensatedIterationCallback accepted_compensated_iteration_callback;
    // Evidence hook called after the declared initialization strategy has been
    // applied and before geometry assembly or the first residual is formed.
    // References are valid only for the duration of the callback.
    using InitializedStateCallback = std::function<void(
        const cfdx::core::Mesh&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>;
    InitializedStateCallback initialized_state_callback;
    // Evidence hook at every actual momentum predictor (including retries).
    // The viscosity and flux are the exact arrays consumed by assembly; U is
    // the current carrier iterate. All references are read-only and transient.
    using MomentumPredictorCallback = std::function<void(
        std::size_t,
        const std::vector<double>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::FACE>&)>;
    MomentumPredictorCallback momentum_predictor_callback;
    // Opt-in forensic snapshot immediately before each momentum component
    // assembly. Called for every predictor attempt, including nonlinear retries.
    // The callback is observational only; references are valid during the call.
    using MomentumAssemblyForensicCallback = std::function<void(
        std::size_t,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::FACE>&,
        const std::vector<double>&)>;
    MomentumAssemblyForensicCallback momentum_assembly_forensic_callback;
    // Optional post-assembly evidence hook. The callback sees the exact scalar
    // equation passed to the momentum solve, with source fields and geometry
    // needed to compare pressure/body-force scales with the assembled operator.
    // Null by default; this callback never changes the equation.
    using MomentumEquationForensicCallback = std::function<void(
        std::size_t,
        std::size_t,
        const ScalarEquation&,
        const FvGeometry&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>;
    MomentumEquationForensicCallback momentum_equation_forensic_callback;
    // Accepted-iterate per-cell convection/diffusion/pressure/RC diagnostic.
    // The callback never changes the solver state, equations, or tolerances.
    std::function<void(std::size_t, std::size_t, const MomentumFaceBudget&,
                       const FvGeometry&)> momentum_face_budget_callback;
    std::size_t momentum_face_budget_every = 1;
    // Opt-in true residual of the RELAXED linear system and the ORIGINAL
    // frozen momentum equation at the same predictor. A small relaxed Krylov
    // residual does not imply a small original/nonlinear residual.
    std::function<void(std::size_t, std::size_t,
                       const MomentumLinearRelaxationEvidence&,
                       double, double)> momentum_linear_relaxation_callback;
    // Same H/A, pressure state and authoritative flux: only published after
    // nonlinear transaction commit, never on a rejected retry.
    std::function<void(std::size_t, const MomentumSameStateTrace&)>
        momentum_same_state_face_callback;
    // Pre-correction Picard assembly: U_old, phi_old and the exact
    // deferred-convection reconstruction. Published only after commit.
    std::function<void(std::size_t, const MomentumPicardTvdTrace&)>
        momentum_picard_tvd_face_callback;
    // Evidence-only hook: invoked at the exact Krylov call point of the
    // coupled solve - after assembly, gauge treatment, finite checks and
    // preconditioner setup - so the snapshot matrix, RHS and initial guess
    // are the very objects handed to the solver. Null by default: no
    // capture, no overhead, no behavior change. It exists so an isolated
    // benchmark can replay the REAL solved system instead of a divergent
    // replica - see tests/validation/test_pcd_schur_action_benchmark.cpp.
    std::function<void(const CoupledAssemblySnapshot&)> coupled_assembly_callback;
    // Optional per-equation diagnostic on the post-gauge coupled matrix.
    // No output and no analysis cost unless enabled. This is an observer, not
    // a new criterion for accepting a nonlinear or linear iteration.
    std::function<void(const std::string&, const cfdx::core::AnisotropyReport&)>
        anisotropy_report_callback;
};

struct IncompressibleIteration {
    std::size_t iteration = 0;
    double high_order_fraction = 1.0; // 1.0 is the only qualified final operator.
    bool full_order_eligible = true; // Held at full order for mandatory iterations.
    double momentum_residual = std::numeric_limits<double>::infinity();
    double pressure_residual = std::numeric_limits<double>::infinity();
    double continuity_l1 = std::numeric_limits<double>::infinity();
    double continuity_linf = std::numeric_limits<double>::infinity();
    double continuity_normalized = std::numeric_limits<double>::infinity();
    bool multisecant_proposed = false;
    bool multisecant_accepted = false;
    double multisecant_relative_step = 0.0;
    double multisecant_estimated_objective = 0.0;
    // Explicit normalization evidence: raw residual / max(abs(momentum RHS))
    // and raw residual / first-iteration raw residual are different metrics.
    double momentum_rhs_scale = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_reference = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_rhs_relative = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_rhs_active_component_max = std::numeric_limits<double>::infinity();
    std::array<double, 3> momentum_rhs_scale_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double, 3> momentum_equation_residual_rhs_relative_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    // Conservative (accepted Rhie-Chow) flux throughput and the globalization
    // merit/backtrack counters, filled only for iterations actually accepted.
    double conservative_face_flux_l1 = std::numeric_limits<double>::infinity();
    double globalization_merit = std::numeric_limits<double>::quiet_NaN();
    std::size_t globalization_backtracks = 0;
    double tvd_memory_source_relative = 0.0;
    // NaN unless the expensive post-SST operator is enabled. Reassembled
    // with the trial k/omega/nu_t, U, pressure gradient and authoritative phi.
    double post_sst_momentum_rhs_relative =
        std::numeric_limits<double>::quiet_NaN();
    // Returned-state certificate; unset when V8 is disabled.
    GlobalNonlinearResidualV8 global_nonlinear_v8;
    double momentum_equation_residual = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_relative = std::numeric_limits<double>::infinity();
    double velocity_change_inf = std::numeric_limits<double>::infinity();
    double pressure_change_inf = std::numeric_limits<double>::infinity();
    double effective_alpha_u = 0.0;
    double effective_alpha_p = 0.0;
    double nonlinear_convergence_metric = std::numeric_limits<double>::infinity();
    // NaN when disabled. Integrated artificial diagonal M_P (kg/s) is not a
    // physical residual and MUST NOT be included in the steady final audit.
    double v17_pseudo_cfl = std::numeric_limits<double>::quiet_NaN();
    double v17_pseudo_mass_min = std::numeric_limits<double>::quiet_NaN();
    double v17_pseudo_mass_max = std::numeric_limits<double>::quiet_NaN();
    std::size_t v17_pseudo_capped_cells = 0;
    std::size_t momentum_linear_iterations = 0;
    std::size_t pressure_linear_iterations = 0;
    std::size_t pressure_correctors_used = 0;
    // Evidence for experimental iterated PISO/PIMPLE: count predictor
    // rebuilds using corrected neighbour velocities and measure their effect.
    // Zero for the production-default frozen-predictor implementation.
    std::size_t piso_neighbor_predictor_rebuilds = 0;
    double piso_neighbor_predictor_flux_delta_linf = 0.0;
    double selected_pressure_alpha = std::numeric_limits<double>::quiet_NaN();
    double pressure_line_search_norm_before = std::numeric_limits<double>::quiet_NaN();
    double pressure_line_search_norm_after = std::numeric_limits<double>::quiet_NaN();
    double linear_tolerance_used = std::numeric_limits<double>::quiet_NaN();
    // Conservative pressure-corrected flux continuity at the end of the iteration.
    double corrected_flux_continuity_linf = std::numeric_limits<double>::infinity();
    // Continuity obtained by reconstructing the face flux from the corrected cell velocity.
    double reconstructed_velocity_continuity_linf = std::numeric_limits<double>::infinity();
    // Maximum face-flux discrepancy between the conservative flux and reconstructed U flux.
    double flux_velocity_mismatch_linf = std::numeric_limits<double>::infinity();
    // Boundary-separated evidence channel: same discrepancy restricted to
    // fixed-pressure boundary faces, where the two reconstruction paths
    // differ analytically (Rhie-Chow boundary correction plus the
    // predictor/corrector velocity response). Diagnostic only.
    double flux_velocity_mismatch_boundary_linf = std::numeric_limits<double>::infinity();
    // Normalized internal-face mismatch
    // E_phi = ||phi_RC - phi_U||_inf / max(||phi_RC||_inf, rho*max|U|*max|Sf|),
    // with an explicit physical floor so the value is comparable across
    // meshes and never divides by a vanishing flux scale.
    double flux_velocity_mismatch_normalized = std::numeric_limits<double>::infinity();
    // Operator identity: max |phi_stored - phi_rederived| where phi_rederived
    // is the same Rhie-Chow operator re-derived on the current state.
    // Expected zero up to round-off; nonzero flags a stale or mutated
    // authoritative flux, or an omitted/duplicated SIMPLEC response term.
    double rhie_chow_identity_linf = std::numeric_limits<double>::infinity();
    // Independent SIMPLEC fixed-point equivalent: the fully corrected
    // rAtU flux must match Rhie-Chow reconstructed from the same H/A and
    // physical pressure with the unmodified rAU response. This catches
    // a missing or duplicated (rAtU-rAU) correction on skewed faces.
    double simplec_plain_rau_flux_gap_linf = 0.0;
    std::size_t simplec_plain_rau_worst_face = 0;
    // SIMPLEC coefficient audit: compare the signed row sum used by the
    // production denominator with the absolute-value alternative. These are
    // diagnostics only; validation gates remain independent.
    std::array<double, 3> simplec_signed_offdiag_sum_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_abs_offdiag_sum_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_denominator_signed_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_denominator_abs_linf{0.0, 0.0, 0.0};
    std::array<double, 3> simplec_denominator_gap_linf{0.0, 0.0, 0.0};
    std::array<std::size_t, 3> simplec_denominator_worst_cell{0, 0, 0};

    // Component-wise and location-aware diagnostics for independently
    // reassembled momentum equations. These are diagnostic signals only;
    // acceptance thresholds remain owned by the validation tests.
    std::array<double, 3> momentum_equation_residual_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    double momentum_equation_residual_internal = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_boundary = std::numeric_limits<double>::infinity();
    double pressure_gradient_linf = std::numeric_limits<double>::infinity();
    double pressure_gradient_l2 = std::numeric_limits<double>::infinity();
    std::size_t momentum_residual_cell = 0;
    double momentum_residual_no_pressure = std::numeric_limits<double>::infinity();
    double momentum_pressure_contribution = std::numeric_limits<double>::infinity();
    std::string momentum_residual_patch;
    // Independent finite-volume conservation diagnostics for the authoritative
    // pressure-corrected mass flux. These are distinct from linear residuals.
    double mass_boundary_flux = std::numeric_limits<double>::infinity();
    double mass_global_cell_balance = std::numeric_limits<double>::infinity();
    double mass_local_l1 = std::numeric_limits<double>::infinity();
    double mass_local_linf = std::numeric_limits<double>::infinity();
    double mass_local_l2 = std::numeric_limits<double>::infinity();
    double mass_normalized_imbalance = std::numeric_limits<double>::infinity();
    std::size_t mass_worst_cell = 0;
    std::size_t mass_nonfinite_faces = 0;
    std::size_t boundedness_nonfinite_velocity = 0;
    std::array<double, 3> momentum_conservation_residual{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double, 3> momentum_conservation_normalized{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<std::size_t, 3> momentum_conservation_worst_cell{0,0,0};
    // True nonlinear face-flux balance: worst per-cell defect, not the
    // global signed sum (which can cancel), on the accepted physical state.
    // A transient audit includes BE/BDF2 accumulation; pseudo-time is excluded.
    std::array<double, 3> momentum_face_local_linf{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double, 3> momentum_face_local_relative{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
};

struct IncompressibleSolveResult {
    bool converged = false;
    std::size_t iterations = 0;
    std::vector<IncompressibleIteration> history;
    // Retained point-probe samples, one per configured probe per iteration.
    // This is the artifact a probe export is written from.
    std::vector<IncompressibleProbeSample> probe_samples;
    // Final conservative face mass flux used by the pressure-velocity
    // continuity operator. This is exposed as evidence so validation can
    // compare algorithms on the same authoritative face operator rather than
    // reconstructing phi from the converged cell velocity.
    cfdx::core::Field<double, cfdx::core::Location::FACE> authoritative_mass_flux;
    // V7: low part at the LAST ACCEPTED state. Empty unless opt-in V6 pressure
    // production is active. Pair with authoritative_mass_flux for DAT6 write.
    cfdx::core::Vector authoritative_pressure_low;
    double reference_momentum_residual = 0.0;
    cfdx::core::LinearSolverContextStats pressure_linear_context;
    cfdx::core::ConvergenceStatus convergence_status = cfdx::core::ConvergenceStatus::CONTINUE;
    std::string convergence_reason;
    // The plans the dispatcher actually resolved, not just the request. A
    // requested method can legitimately resolve to another one, and a report
    // cannot show that unless the resolved plan is recorded here.
    cfdx::core::LinearSolverPlan pressure_linear_plan;
    bool pressure_linear_plan_resolved = false;
    cfdx::core::LinearSolverPlan coupled_linear_plan;
    bool coupled_linear_plan_resolved = false;
    CoupledSchurModel requested_coupled_schur_model = CoupledSchurModel::BlockLocal;
    CoupledSchurModel resolved_coupled_schur_model = CoupledSchurModel::BlockLocal;
    bool coupled_schur_model_resolved = false;
};

// Applies the declared initialization strategy to the caller's fields.
//
// legacy_restart_path is the positional restart argument of
// solve_steady_incompressible. It supplies the checkpoint for a declared
// Restart mode, so existing callers that only pass the positional argument keep
// working. A caller that declares Provided or Uniform while also passing a
// legacy path gets the path ignored rather than silently applied, because that
// would make the initial state depend on an undeclared input.
inline void apply_initialization(
    const InitializationControls& init,
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    cfdx::io::DatRestartFields& restart_fields,
    const std::string& legacy_restart_path)
{
    using namespace cfdx::core;
    const std::string& restart_source =
        init.restart_path.empty() ? legacy_restart_path : init.restart_path;
    switch (init.mode) {
        case InitializationMode::Provided:
            // Caller-owned state is used as-is.
            break;
        case InitializationMode::Uniform:
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                U.component_data(0)[c] = init.uniform_velocity.x;
                U.component_data(1)[c] = init.uniform_velocity.y;
                U.component_data(2)[c] = init.uniform_velocity.z;
                p(c) = init.uniform_pressure;
            }
            break;
        case InitializationMode::Restart: {
            if (restart_source.empty())
                throw std::invalid_argument(
                    "initialization mode 'restart' requires a checkpoint path");
            (void)cfdx::io::read_dat_restart_fields(
                restart_source, mesh, U, p, restart_fields);
            break;
        }
    }
}

inline void validate_initialization_controls(const InitializationControls& init) {
    switch (init.mode) {
        case InitializationMode::Provided:
            break;
        case InitializationMode::Uniform:
            if (!std::isfinite(init.uniform_pressure) ||
                !std::isfinite(init.uniform_velocity.x) ||
                !std::isfinite(init.uniform_velocity.y) ||
                !std::isfinite(init.uniform_velocity.z))
                throw std::invalid_argument(
                    "uniform initialization state must be finite");
            break;
        case InitializationMode::Restart:
            // The checkpoint may come from init.restart_path or from the
            // positional argument; apply_initialization resolves which.
            break;
    }
}

inline void validate_incompressible_controls(
    const IncompressibleSolverControls& c, std::size_t n_cells)
{
    validate_initialization_controls(c.initialization);
    if (c.initial_face_mass_flux_from_paired_checkpoint) {
        if (c.initial_face_mass_flux.empty() ||
            !std::isfinite(c.initial_momentum_residual_reference) ||
            !(c.initial_momentum_residual_reference > 0.0) ||
            !std::isfinite(c.initial_nonlinear_metric_reference) ||
            !(c.initial_nonlinear_metric_reference > 0.0) ||
            c.initial_accepted_iteration == 0 ||
            c.initial_accepted_iteration > std::numeric_limits<std::size_t>::max() -
                c.convergence.max_iterations || c.convection_continuation.enabled)
            throw std::invalid_argument(
                "verified paired restart requires full reference/iteration contract (no continuation)");
    } else if (c.initial_momentum_residual_reference != 0.0 ||
               c.initial_nonlinear_metric_reference != 0.0 ||
               c.initial_accepted_iteration != 0)
        throw std::invalid_argument("incomplete initial reference restart contract");
    validate_coupling_controls(c.coupling);
    if (c.coupling.piso_neighbor_correction &&
        c.algorithm != PressureVelocityAlgorithm::PISO &&
        c.algorithm != PressureVelocityAlgorithm::PIMPLE)
        throw std::invalid_argument(
            "PISO neighbour correction requires PISO or PIMPLE algorithm");
    if (c.coupling.pcd_lumped_mass_commutator_experimental &&
        (c.algorithm != PressureVelocityAlgorithm::COUPLED ||
         c.coupling.schur_model != CoupledSchurModel::PCD))
        throw std::invalid_argument(
            "PCD lumped-mass commutator requires COUPLED with PCD Schur");
    if (c.coupling.pcd_reverse_lumped_mass_v39_experimental &&
        (c.algorithm != PressureVelocityAlgorithm::COUPLED ||
         c.coupling.schur_model != CoupledSchurModel::PCD ||
         c.coupling.pcd_lumped_mass_commutator_experimental))
        throw std::invalid_argument(
            "V39 reverse PCD requires COUPLED/PCD and excludes other commutator mode");
    if (c.coupling.pcd_blocklocal_recovery_v42_experimental &&
        (c.algorithm != PressureVelocityAlgorithm::COUPLED ||
         c.coupling.schur_model != CoupledSchurModel::PCD ||
         c.coupling.coupled_monolithic_ilut_fallback_v16_experimental))
        throw std::invalid_argument(
            "V42 BlockLocal recovery requires COUPLED/PCD and excludes ILUT recovery");
    if (c.coupling.coupled_schwarz_experimental &&
        (c.algorithm != PressureVelocityAlgorithm::COUPLED ||
         c.coupling.schur_model != CoupledSchurModel::BlockLocal))
        throw std::invalid_argument(
            "experimental Schwarz requires COUPLED with BlockLocal Schur");
    validate_convergence_criteria(c.convergence);
    validate_convergence_acceleration_controls(c.acceleration);
    validate_adaptive_relaxation_controls(c.adaptive_relaxation);
    validate_nonlinear_retry_controls(c.nonlinear_retry);
    validate_coupled_phase_escape_v45(
        c.coupled_period_four_damping_v45_experimental,
        c.coupled_period_four_damping_steps_v45);
    if (c.coupled_period_four_damping_v45_experimental > 0.0 &&
        !c.coupled_period_four_full_state_v44_experimental)
        throw std::invalid_argument(
            "V45 cycle-escape continuation requires V44 coherent full-state predictor");
    if (!std::isfinite(c.coupled_outer_step_v43_experimental) ||
        c.coupled_outer_step_v43_experimental < 0.0 ||
        c.coupled_outer_step_v43_experimental > 1.0)
        throw std::invalid_argument("V43 COUPLED outer step must be in [0,1]");
    if ((c.coupled_period_four_mean_v43_experimental &&
         (c.coupled_outer_step_v43_experimental > 0.0 ||
          c.coupled_period_four_full_state_v44_experimental)) ||
        (c.coupled_period_four_full_state_v44_experimental &&
         c.coupled_outer_step_v43_experimental > 0.0))
        throw std::invalid_argument(
            "V44 full-phase, V43 phase-mean and carrier damping modes are exclusive");
    // The paired pressure-high/low representation must not be mixed with
    // an unpaired gauge-pressure phase average. The V44 predictor cannot
    // yet preserve its companion low-pressure state and correction flux.
    if (c.coupled_period_four_full_state_v44_experimental &&
        c.pressure_high_low_production)
        throw std::invalid_argument(
            "V44 full-phase predictor is incompatible with paired pressure-high/low state");
    if ((c.coupled_outer_step_v43_experimental > 0.0 ||
         c.coupled_period_four_mean_v43_experimental ||
         c.coupled_period_four_full_state_v44_experimental) &&
        (c.algorithm != PressureVelocityAlgorithm::COUPLED ||
         c.transient.enabled || !c.sst_runtime ||
         !c.require_global_nonlinear_certificate_v8 ||
         !c.require_post_sst_momentum_reassembly ||
         c.initial_face_mass_flux_from_paired_checkpoint ||
         c.outer_ptc_v31_experimental || c.pseudo_momentum.enabled))
        throw std::invalid_argument(
            "V43 COUPLED outer step requires steady COUPLED/SST, V7/V8, "
            "no PTC or incomplete paired restart");
    if (c.monotone_global_merit_v14_experimental &&
        !c.joint_globalization_v13_experimental)
        throw std::invalid_argument("V14 monotone physical merit requires V13 joint globalization");
    if (c.joint_globalization_v13_experimental &&
        (!c.require_global_nonlinear_certificate_v8 ||
         !c.require_post_sst_momentum_reassembly ||
         !c.sst_runtime || c.transient.enabled ||
         c.nonlinear_retry.max_retries == 0))
        throw std::invalid_argument(
            "V13 global trial guard requires steady SST, V8 global residual, "
            "post-SST momentum reconstruction and nonlinear retries");
    if (c.outer_ptc_balance_v32_experimental &&
        (!c.outer_ptc_v31_experimental ||
         (!c.outer_ptc_balance_momentum_v32 && !c.outer_ptc_balance_sst_v32)))
        throw std::invalid_argument("V32 PTC requires V31 and at least one capped subsystem");
    if (c.outer_ptc_v31_experimental &&
        (!c.sst_runtime || c.transient.enabled || c.pseudo_momentum.enabled ||
         !c.require_global_nonlinear_certificate_v8 ||
         !c.require_post_sst_momentum_reassembly ||
         !c.joint_globalization_v13_experimental ||
         !c.monotone_global_merit_v14_experimental ||
         c.initial_face_mass_flux_from_paired_checkpoint ||
         !c.sst_runtime->transport.sst_coupled_ptc_v30_experimental ||
         !c.sst_runtime->transport.sst_jfnk_v22_experimental ||
         c.sst_runtime->allow_bounded_provisional_transport))
        throw std::invalid_argument("V31 outer PTC requires steady SST, V8/V13/V14 physical gates, V30 JFNK, no provisional transport/restart and no duplicate V17 pseudo-time");
    // V4 checkpoints contain U, p, phi and norm references, not the
    // committed deferred-TVD history or a full nonmonotone merit window.
    // Such a restart would not reproduce the nonlinear iteration map.
    if (c.initial_face_mass_flux_from_paired_checkpoint &&
        (c.tvd_deferred_memory.enabled || c.tvd_secant_implicit.enabled ||
         c.picard_physical_globalization.enabled ||
         c.nonlinear_multisecant.method != NonlinearMultisecantMethod::OFF))
        throw std::invalid_argument(
            "paired restart cannot restore TVD/Picard/multisecant nonlinear history");
    validate_picard_physical_globalization(c.picard_physical_globalization);
    validate_tvd_deferred_defect_memory(c.tvd_deferred_memory);
    validate_tvd_secant_implicit(c.tvd_secant_implicit);
    if (c.tvd_secant_implicit.enabled && c.tvd_deferred_memory.enabled)
        throw std::invalid_argument(
            "TVD secant and deferred memory require a qualified combined scheme");
    if (c.require_simplec_same_state_flux_identity &&
        c.algorithm != PressureVelocityAlgorithm::SIMPLEC)
        throw std::invalid_argument(
            "SIMPLEC same-state flux identity requires SIMPLEC algorithm");
    validate_nonlinear_multisecant(c.nonlinear_multisecant);
    if (c.nonlinear_multisecant.method != NonlinearMultisecantMethod::OFF &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE ||
         c.transient.enabled || c.sst_runtime ||
         c.pseudo_momentum.enabled || c.convection_continuation.enabled ||
         c.picard_physical_globalization.enabled ||
         c.tvd_deferred_memory.enabled || c.tvd_secant_implicit.enabled ||
         c.adaptive_relaxation.enabled ||
         c.pressure_continuity_line_search.enabled ||
         c.v19_full_pressure_flux_step ||
         (c.convection_scheme != ConvectionScheme::TVD &&
          c.convection_scheme != ConvectionScheme::TVD_VANLEER &&
          c.convection_scheme != ConvectionScheme::TVD_VANALBADA &&
          c.convection_scheme != ConvectionScheme::MUSCL_CELL_LIMITED)))
        throw std::invalid_argument(
            "nonlinear multisecant requires steady SIMPLE with full TVD, "
            "no transient/SST/ramp, and no simultaneous experimental acceleration");
    if (c.tvd_secant_implicit.enabled &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE || c.transient.enabled ||
         c.sst_runtime || c.pseudo_momentum.enabled ||
         c.convection_continuation.enabled ||
         (c.convection_scheme != ConvectionScheme::TVD &&
          c.convection_scheme != ConvectionScheme::TVD_VANLEER &&
          c.convection_scheme != ConvectionScheme::TVD_VANALBADA)))
        throw std::invalid_argument(
            "TVD implicit secant requires steady SIMPLE and the complete TVD operator");
    if (c.tvd_deferred_memory.enabled &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE || c.transient.enabled ||
         c.sst_runtime || c.pseudo_momentum.enabled ||
         c.convection_continuation.enabled ||
         (c.convection_scheme != ConvectionScheme::TVD &&
          c.convection_scheme != ConvectionScheme::TVD_VANLEER &&
          c.convection_scheme != ConvectionScheme::TVD_VANALBADA)))
        throw std::invalid_argument(
            "TVD defect memory requires steady SIMPLE, full TVD and no transient/SST/pseudo/ramp");
    if (c.picard_physical_globalization.enabled &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE || c.transient.enabled ||
         c.sst_runtime != nullptr || c.pseudo_momentum.enabled ||
         c.v19_full_pressure_flux_step || c.adaptive_relaxation.enabled))
        throw std::invalid_argument(
            "Picard physical globalization requires steady SIMPLE, no SST/pseudo, "
            "no split V19 flux and no simultaneous adaptive relaxation");
    validate_convection_continuation(c.convection_continuation,
                                     c.convergence.max_iterations);
    if (c.convection_continuation.enabled &&
        ((c.convection_scheme != ConvectionScheme::TVD &&
          c.convection_scheme != ConvectionScheme::TVD_VANLEER &&
          c.convection_scheme != ConvectionScheme::TVD_VANALBADA &&
          c.convection_scheme != ConvectionScheme::MUSCL_CELL_LIMITED &&
          c.convection_scheme != ConvectionScheme::SECOND_ORDER_UPWIND) ||
         c.algorithm != PressureVelocityAlgorithm::SIMPLE ||
         c.transient.enabled || c.sst_runtime != nullptr))
        throw std::invalid_argument(
            "convection continuation: only steady SIMPLE with SOU/TVD is qualified");
    validate_pressure_continuity_line_search(c.pressure_continuity_line_search);
    if (c.pressure_replay_prefix.empty() != c.pressure_replay_iterations.empty() ||
        c.pressure_replay_iterations.size() > 32 ||
        std::any_of(c.pressure_replay_iterations.begin(),
                    c.pressure_replay_iterations.end(),
                    [&](std::size_t i) {
                        return i == 0 || i > c.convergence.max_iterations;
                    }))
        throw std::invalid_argument(
            "pressure linear replay: require prefix and 1..32 valid iterations");
    auto selected_iterations = c.pressure_replay_iterations;
    std::sort(selected_iterations.begin(), selected_iterations.end());
    if (std::adjacent_find(selected_iterations.begin(), selected_iterations.end()) !=
        selected_iterations.end())
        throw std::invalid_argument("pressure linear replay: duplicate iteration");
    // The monolithic U-p block still assembles its analytical G operator
    // with the historical face weights. Until that block is augmented with
    // the skew-correction Jacobian, activating this option would make the
    // coupled matrix inconsistent with the reconstructed Rhie-Chow flux.
    if (c.pressure_projected_face_interpolation &&
        c.algorithm != PressureVelocityAlgorithm::SIMPLE)
        throw std::invalid_argument(
            "projected pressure face interpolation requires segregated SIMPLE");
    if (c.pressure_skew_corrected_gradient &&
        c.algorithm == PressureVelocityAlgorithm::COUPLED)
        throw std::invalid_argument(
            "skew-corrected pressure gradient is not qualified for COUPLED");
    if (c.require_post_sst_momentum_reassembly && !c.sst_runtime)
        throw std::invalid_argument(
            "post-SST momentum reassembly requires an SST production runtime");

    if (c.sst_runtime) {
        if (!c.sst_runtime->state)
            throw std::invalid_argument("SST production runtime requires an initialized state");
        if (c.sst_runtime->transport.model != TurbulenceModel::SST)
            throw std::invalid_argument("SST production runtime requires SST transport controls");
        if (c.sst_runtime->transport.density != c.density ||
            c.sst_runtime->transport.molecular_viscosity != c.kinematic_viscosity)
            throw std::invalid_argument(
                "SST transport density and molecular viscosity must match the carrier material properties");
        if (c.sst_runtime->max_transport_iterations == 0 ||
            !(c.sst_runtime->transport_tolerance > 0.0) ||
            !std::isfinite(c.sst_runtime->transport_tolerance))
            throw std::invalid_argument("SST production runtime has invalid transport controls");
        validate_sst_production_state(*c.sst_runtime->state, n_cells);
        if (!c.cell_turbulent_viscosity.empty() || c.turbulent_viscosity != 0.0)
            throw std::invalid_argument(
                "SST production runtime cannot be combined with prescribed turbulent viscosity");
    }
    if (!c.qoi_gates.empty()) {
        validate_convergence_monitor_controls(cfdx::core::ConvergenceMonitorControls{
            cfdx::core::ConvergenceGate{}, 1.0, 1, 1000, 1, 1.0, 1, 10.0, c.qoi_gates});
        for (const auto& gate : c.qoi_gates) {
            const bool declared =
                std::any_of(c.probes.begin(), c.probes.end(),
                            [&](const IncompressiblePointProbe& probe) {
                                return probe.name == gate.name;
                            });
            if (!declared)
                throw std::invalid_argument(
                    "convergence QoI gate '" + gate.name +
                    "' does not name a configured probe");
        }
    }
    if (c.v19_full_pressure_flux_step &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE ||
         (!c.v18_full_flux_pressure_oracle && !c.v18_sparse_flux_pressure_jacobian) ||
         c.pressure_continuity_line_search.enabled))
        throw std::invalid_argument(
            "V19 split pressure/flux requires SIMPLE, V18 full Jacobian, "
            "one pressure corrector and no pressure line-search");
    if (c.pressure_high_low_production &&
        std::numeric_limits<long double>::digits <=
            std::numeric_limits<double>::digits)
        throw std::invalid_argument(
            "V6 high/low pressure requires extended long double; "
            "this target needs a portable double-double face operator");
    if (c.pressure_high_low_production &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE ||
         !c.v18_sparse_flux_pressure_jacobian ||
         !c.v19_full_pressure_flux_step ||
         c.pressure_continuity_line_search.enabled ||
         c.pressure_skew_corrected_gradient ||
         c.nonlinear_multisecant.method != NonlinearMultisecantMethod::OFF))
        throw std::invalid_argument(
            "V6 pressure high/low requires SIMPLE, sparse V18 full RC Jacobian, "
            "V19 full flux, standard one-sweep gradient, no separate affine "
            "pressure search/multisecant, and DAT6 for restart");
    if (c.v18_full_flux_pressure_oracle && c.v18_sparse_flux_pressure_jacobian)
        throw std::invalid_argument("V18 reference and sparse pressure solvers are mutually exclusive");
    if (c.v18_sparse_flux_pressure_jacobian &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE ||
         c.pressure_skew_corrected_gradient))
        throw std::invalid_argument("V18 sparse pressure operator requires SIMPLE and standard gradients");
    if (c.v18_full_flux_pressure_oracle &&
        (c.algorithm != PressureVelocityAlgorithm::SIMPLE ||
         c.pressure_skew_corrected_gradient || n_cells>1024))
        throw std::invalid_argument(
            "V18 full pressure flux oracle requires SIMPLE, standard gradients, <=1024 cells");
    validate_v17_pseudo_momentum_controls(c.pseudo_momentum);
    if (c.pseudo_momentum.enabled) {
        if (c.transient.enabled)
            throw std::invalid_argument("V17 pseudo-momentum and physical transient are exclusive");
        if (c.algorithm == PressureVelocityAlgorithm::COUPLED)
            throw std::invalid_argument("V17 pseudo-momentum requires segregated pressure coupling");
    }
    if (c.transient.enabled) {
        if (!(c.transient.dt > 0.0) || !std::isfinite(c.transient.dt))
            throw std::invalid_argument("transient dual-time dt must be finite and positive");
        if (c.transient.scheme == cfdx::core::TimeScheme::BDF2 &&
            c.transient.history_valid &&
            (!(c.transient.dt_previous > 0.0) || !std::isfinite(c.transient.dt_previous)))
            throw std::invalid_argument("transient BDF2 previous dt must be finite and positive");
        if (!c.transient.previous || c.transient.previous->size() != n_cells ||
            c.transient.previous->dimension() != 3)
            throw std::invalid_argument("transient dual-time previous velocity history is missing or invalid");
        if (c.transient.scheme == cfdx::core::TimeScheme::BDF2 && c.transient.history_valid &&
            (!c.transient.previous_previous || c.transient.previous_previous->size() != n_cells ||
             c.transient.previous_previous->dimension() != 3))
            throw std::invalid_argument("transient BDF2 second velocity history is missing or invalid");
    }
    if (!(c.density > 0.0) || !std::isfinite(c.density) ||
        !(c.kinematic_viscosity >= 0.0) || !std::isfinite(c.kinematic_viscosity) ||
        !(c.turbulent_viscosity >= 0.0) || !std::isfinite(c.turbulent_viscosity))
        throw std::invalid_argument("invalid incompressible material properties");
    // Reject overflow in the effective dynamic viscosity before applying
    // initialization. Otherwise an invalid material combination can mutate
    // caller-owned U/p fields before assembly discovers non-finite coefficients.
    const auto validate_effective_viscosity = [&](double nu_t) {
        const double total_kinematic_viscosity = c.kinematic_viscosity + nu_t;
        const double effective_dynamic_viscosity =
            c.density * total_kinematic_viscosity;
        if (!std::isfinite(total_kinematic_viscosity) ||
            !std::isfinite(effective_dynamic_viscosity))
            throw std::invalid_argument(
                "effective dynamic viscosity must be finite");
    };
    validate_effective_viscosity(c.turbulent_viscosity);
    if (!c.cell_turbulent_viscosity.empty()) {
        if (c.cell_turbulent_viscosity.size() != n_cells)
            throw std::invalid_argument(
                "cellwise turbulent viscosity size does not match mesh");
        if (c.turbulent_viscosity != 0.0)
            throw std::invalid_argument(
                "specify either scalar or cellwise turbulent viscosity, not both");
        for (const double nu_t : c.cell_turbulent_viscosity) {
            if (!(nu_t >= 0.0) || !std::isfinite(nu_t))
                throw std::invalid_argument(
                    "cellwise turbulent viscosity must be finite and non-negative");
            validate_effective_viscosity(nu_t);
        }
    }
    if (!c.cell_body_force.empty()) {
        if (c.cell_body_force.size() != n_cells)
            throw std::invalid_argument(
                "cellwise body force size does not match mesh");
        if (c.body_force.x != 0.0 || c.body_force.y != 0.0 ||
            c.body_force.z != 0.0)
            throw std::invalid_argument(
                "specify either a uniform or a cellwise body force, not both");
        for (const auto& force : c.cell_body_force) {
            if (!std::isfinite(force.x) || !std::isfinite(force.y) ||
                !std::isfinite(force.z))
                throw std::invalid_argument(
                    "cellwise body force components must be finite");
        }
    }
    if (c.pressure_reference_cell >= n_cells)
        throw std::invalid_argument("pressure reference cell out of range");
    if (c.linear_max_iterations == 0 || !(c.linear_tolerance > 0.0) ||
        !std::isfinite(c.linear_tolerance) ||
        !std::isfinite(c.pressure_linear_tolerance) ||
        c.pressure_linear_tolerance < 0.0)
        throw std::invalid_argument("invalid linear solver controls");
}

// Independent, read-only pressure-correction diagnostic. A tiny backward
// error in A*x=b does NOT certify a sufficiently small local mass defect.
// Use extended precision for this diagnostic to avoid cancellation errors
// in the double-precision matrix-vector product itself. It can never alter
// the linear-solver status or the physical convergence gate.
struct PressureLinearCancellationDiagnostic {
    bool valid = false;
    double max_componentwise_backward_error = std::numeric_limits<double>::infinity();
    double max_local_mass_ratio = std::numeric_limits<double>::infinity();
};

inline PressureLinearCancellationDiagnostic diagnose_pressure_linear_cancellation(
    const cfdx::core::SparseMatrix& matrix,
    const cfdx::core::Vector& rhs,
    const cfdx::core::Vector& candidate,
    const std::vector<double>& local_incident_mass_scale,
    double required_local_mass_ratio)
{
    PressureLinearCancellationDiagnostic certificate;
    const std::size_t n = matrix.n_rows();
    if (!n || n != matrix.n_cols() || rhs.size() != n || candidate.size() != n ||
        local_incident_mass_scale.size() != n ||
        !std::isfinite(required_local_mass_ratio) ||
        !(required_local_mass_ratio > 0.0))
        return certificate;
    long double largest_backward = 0.0L, largest_local = 0.0L;
    for (std::size_t row = 0; row < n; ++row) {
        const double reference = local_incident_mass_scale[row];
        if (!(reference > 0.0) || !std::isfinite(reference) ||
            !std::isfinite(rhs(row)))
            return certificate;
        long double ax = 0.0L, absolute_terms = std::abs(static_cast<long double>(rhs(row)));
        for (std::uint32_t k = matrix.row_offsets_data()[row];
             k < matrix.row_offsets_data()[row + 1]; ++k) {
            const double a = matrix.values_data()[k];
            const double x = candidate(matrix.columns_data()[k]);
            if (!std::isfinite(a) || !std::isfinite(x))
                return certificate;
            const long double term = static_cast<long double>(a) * x;
            ax += term;
            absolute_terms += std::abs(term);
        }
        const long double defect = std::abs(ax - static_cast<long double>(rhs(row)));
        const long double eta = absolute_terms > 0.0L
            ? defect / absolute_terms : (defect == 0.0L ? 0.0L :
               std::numeric_limits<long double>::infinity());
        largest_backward = std::max(largest_backward, eta);
        largest_local = std::max(largest_local, defect / reference);
    }
    certificate.max_componentwise_backward_error =
        static_cast<double>(largest_backward);
    certificate.max_local_mass_ratio = static_cast<double>(largest_local);
    certificate.valid = largest_backward <=
            64.0L * std::numeric_limits<double>::epsilon() &&
        largest_local <= required_local_mass_ratio;
    return certificate;
}

inline void apply_pressure_gauge_policy(
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const FvGeometry& geometry,
    PressureGaugePolicy policy,
    std::size_t reference_cell,
    double reference_value,
    bool has_fixed_pressure_boundary)
{
    if (has_fixed_pressure_boundary || policy == PressureGaugePolicy::NONE) return;
    if (policy == PressureGaugePolicy::REFERENCE_CELL) {
        const double shift = reference_value - p(reference_cell);
        for (std::size_t c = 0; c < p.size(); ++c) p(c) += shift;
        return;
    }
    if (policy == PressureGaugePolicy::ZERO_MEAN) {
        double volume = 0.0, weighted_sum = 0.0;
        for (std::size_t c = 0; c < p.size(); ++c) {
            volume += geometry.cell_volumes[c];
            weighted_sum += geometry.cell_volumes[c] * p(c);
        }
        if (!(volume > 0.0) || !std::isfinite(volume) || !std::isfinite(weighted_sum))
            throw std::runtime_error("pressure gauge: invalid volume-weighted mean");
        const double mean = weighted_sum / volume;
        for (std::size_t c = 0; c < p.size(); ++c) p(c) -= mean;
    }
}

// V19: split the full pressure step used for flux correction from the
// under-relaxed physical pressure used in the momentum corrector. Both
// states have the same gauge policy. This is deliberately transactional:
// validate the full increment before modifying either output.
inline void apply_v19_pressure_flux_policy(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& pressure,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& flux_pressure,
    const cfdx::core::Vector& correction,
    double alpha,
    bool full_flux,
    const FvGeometry& geometry,
    PressureGaugePolicy gauge,
    std::size_t reference_cell,
    double reference_value,
    bool has_fixed_pressure_boundary)
{
    if (correction.size()!=pressure.size() ||
        flux_pressure.size()!=pressure.size() ||
        !(alpha>0.0 && alpha<=1.0) || !std::isfinite(alpha))
        throw std::invalid_argument("V19 pressure flux policy: invalid dimensions/relaxation");
    for (std::size_t c=0;c<pressure.size();++c)
        if (!std::isfinite(pressure(c)) || !std::isfinite(correction(c)))
            throw std::invalid_argument("V19 pressure flux policy: nonfinite pressure/increment");
    auto adjusted_pressure=pressure;
    auto adjusted_flux_pressure=pressure;
    for (std::size_t c=0;c<pressure.size();++c) {
        const double original=pressure(c);
        adjusted_pressure(c)=original+alpha*correction(c);
        if (full_flux) adjusted_flux_pressure(c)=original+correction(c);
        if (!std::isfinite(adjusted_pressure(c)) ||
            !std::isfinite(adjusted_flux_pressure(c)))
            throw std::runtime_error("V19 pressure flux policy overflow");
    }
    apply_pressure_gauge_policy(adjusted_pressure, geometry, gauge, reference_cell,
                                reference_value, has_fixed_pressure_boundary);
    if (full_flux)
        apply_pressure_gauge_policy(adjusted_flux_pressure, geometry, gauge, reference_cell,
                                    reference_value, has_fixed_pressure_boundary);
    else
        adjusted_flux_pressure=adjusted_pressure;
    pressure=std::move(adjusted_pressure);
    flux_pressure=std::move(adjusted_flux_pressure);
}

inline cfdx::core::Field<double, cfdx::core::Location::CELL>
gauss_gradient_with_boundary(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& bcs,
    bool skew_correct = false,
    bool projected_face_interpolation = false)
{
    using namespace cfdx::core;
    if (field.dimension() != 1 || field.size() != mesh.n_cells())
        throw std::invalid_argument("gauss_gradient_with_boundary: invalid field");

    Field<double, Location::CELL> grad(
        mesh.n_cells(), field.name() + "_grad_bc", field.metadata().unit + "/m", 3);
    const auto& own = mesh.ownership();

    // Opt-in fixed-point correction for skewed internal faces. Every face
    // is assembled once per cell, preserving the conservative orientation.
    // The default one-sweep Green-Gauss implementation remains unchanged.
    const std::size_t sweeps = skew_correct ? 8 : 1;
    for (std::size_t sweep = 0; sweep < sweeps; ++sweep) {
    // Avoid copying a vector field on the default (single sweep) path.
    std::optional<Field<double, Location::CELL>> previous;
    if (skew_correct && sweep > 0) previous.emplace(grad);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        Vec3 sum{0.0, 0.0, 0.0};
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = own.owner(f) == c;
            Vec3 Sf = owner ? geometry.face_area_vectors[f] : geometry.face_area_vectors[f] * -1.0;
            double vf = field(c);

            if (own.neighbour(f) >= 0) {
                // The mesh stores owner/neighbour independent of the cell being
                // visited.  When c is the neighbour, own.neighbour(f) == c;
                // using it again would therefore collapse the face interpolation
                // to field(c) and halve the Green-Gauss pressure gradient.
                const std::size_t other = owner
                    ? static_cast<std::size_t>(own.neighbour(f))
                    : static_cast<std::size_t>(own.owner(f));
                // V8: the pressure-face predictor and Rhie--Chow Schur
                // operator optionally use the longitudinal projection.  Its
                // cell gradients must use EXACTLY the same interpolation,
                // including the skew correction's reference point.  Reuse
                // the existing core geometry weight; do not invent another
                // LS/Green-Gauss operator.  Preserve the legacy mode for
                // regression comparisons and COUPLED's existing Jacobian.
                double w = 0.0;
                if (projected_face_interpolation) {
                    const double wn = cfdx::core::geometric_interpolation_weight(
                        geometry.face_centres[f], geometry.cell_centres[c],
                        &geometry.cell_centres[other]);
                    w = 1.0 - wn;
                } else {
                    const double dc = (geometry.face_centres[f] -
                                       geometry.cell_centres[c]).mag();
                    const double dn = (geometry.cell_centres[other] -
                                       geometry.face_centres[f]).mag();
                    if (!(dc + dn > 0.0))
                        throw std::runtime_error(
                            "gauss_gradient_with_boundary: degenerate face distance");
                    w = dn / (dc + dn);
                }
                vf = w * field(c) + (1.0 - w) * field(other);
                if (skew_correct && sweep > 0) {
                    const Vec3 interpolated_centre =
                        geometry.cell_centres[c] * w +
                        geometry.cell_centres[other] * (1.0 - w);
                    const Vec3 skew = geometry.face_centres[f] - interpolated_centre;
                    const Vec3 grad_face{
                        0.5 * (previous->component_data(0)[c] + previous->component_data(0)[other]),
                        0.5 * (previous->component_data(1)[c] + previous->component_data(1)[other]),
                        0.5 * (previous->component_data(2)[c] + previous->component_data(2)[other])};
                    vf += grad_face.dot(skew);
                }
            } else {
                const std::size_t p = geometry.face_patch[f];
                if (p < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(p).type == PatchType::EMPTY)
                    continue;
                if (p < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(p).name;
                    const auto it = bcs.find(name);
                    if (it != bcs.end() && it->second.type == ScalarBoundaryType::FIXED_VALUE)
                        vf = it->second.value;
                    else if (it != bcs.end() &&
                             it->second.type == ScalarBoundaryType::FIXED_GRADIENT) {
                        const Vec3 outward = Sf * (1.0 / Sf.mag());
                        const Vec3 centre_to_face =
                            geometry.face_centres[f] - geometry.cell_centres[c];
                        vf += it->second.gradient * centre_to_face.dot(outward);
                    }
                }
            }
            sum += Sf * vf;
        }
        const double invV = 1.0 / geometry.cell_volumes[c];
        grad.component_data(0)[c] = sum.x * invV;
        grad.component_data(1)[c] = sum.y * invV;
        grad.component_data(2)[c] = sum.z * invV;
        for (std::size_t d = 0; d < 3; ++d)
            if (!std::isfinite(grad.component_data(d)[c]))
                throw std::runtime_error(
                    "gauss_gradient_with_boundary: non-finite skew-corrected gradient");
    }
    }
    return grad;
}

// Matrix-free exact directional derivative of the *same finite-sweep* Gauss
// pressure-gradient reconstruction used in production. The skew correction
// is linear in pressure for fixed geometry, so homogeneous boundary data
// produces the Jacobian action without finite-difference cancellation.
// This is an observable Jacobian action, NOT permission to enable the
// incomplete monolithic COUPLED G/C assembly on skewed meshes.
inline cfdx::core::Field<double, cfdx::core::Location::CELL>
gauss_pressure_gradient_jacobian_action(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& direction,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& pressure_bcs,
    bool skew_correct = true,
    bool projected_face_interpolation = false)
{
    ScalarBoundaryConditions homogeneous_bcs = pressure_bcs;
    for (auto& [name, bc] : homogeneous_bcs) {
        (void)name;
        bc.value = 0.0;
        bc.gradient = 0.0;
    }
    return gauss_gradient_with_boundary(direction, mesh, geometry,
                                        homogeneous_bcs, skew_correct,
                                        projected_face_interpolation);
}

// The weight on the neighbour is the longitudinal projection of the actual
// face centroid. A geometric average (0.5) is not affine-exact when the two
// centre-to-face distances differ. The same weight MUST be used in HbyA,
// dAU=V/A_P, pressure gradients and the pressure-correction matrix.
inline double pressure_face_neighbour_weight(
    const FvGeometry& geometry, std::size_t face,
    std::size_t owner, std::size_t neighbour, bool projected)
{
    if (!projected) return 0.5;
    return cfdx::core::geometric_interpolation_weight(
        geometry.face_centres[face], geometry.cell_centres[owner],
        &geometry.cell_centres[neighbour]);
}

inline VelocityFaceConditions resolve_inlet_outlet_velocity_faces(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& previous_phi,
    const VelocityBoundaryConditions& velocity_bcs,
    const std::map<std::string, cfdx::core::Vec3>& inflow_values)
{
    using namespace cfdx::core;
    if (previous_phi.size() != mesh.n_faces() || previous_phi.dimension() != 1 ||
        geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument("inlet/outlet face BC: invalid flux/geometry size");
    VelocityFaceConditions resolved{ScalarBoundaryFaceConditions(mesh.n_faces()),
                                    ScalarBoundaryFaceConditions(mesh.n_faces()),
                                    ScalarBoundaryFaceConditions(mesh.n_faces())};
    for (const auto& [name, value] : inflow_values) {
        const std::size_t patch_id = mesh.boundary().find(name);
        if (patch_id >= mesh.boundary().n_patches() ||
            mesh.boundary().patch(patch_id).empty() ||
            mesh.boundary().patch(patch_id).type == PatchType::EMPTY ||
            mesh.boundary().patch(patch_id).type == PatchType::SYMMETRY)
            throw std::invalid_argument("inlet/outlet face BC: missing or impermeable patch: " + name);
        const auto bc = velocity_bcs.find(name);
        if (bc == velocity_bcs.end() ||
            bc->second.type != VelocityBoundaryCondition::Type::ZERO_GRADIENT)
            throw std::invalid_argument("inlet/outlet face BC: patch must have base U zeroGradient: " + name);
        if (!std::isfinite(value.x) || !std::isfinite(value.y) || !std::isfinite(value.z))
            throw std::invalid_argument("inlet/outlet face BC: non-finite inflow velocity");
    }
    for (std::size_t face=0; face<mesh.n_faces(); ++face) {
        if (mesh.ownership().neighbour(face) >= 0) continue;
        const auto patch_id = geometry.face_patch[face];
        if (patch_id >= mesh.boundary().n_patches()) continue;
        const auto& patch = mesh.boundary().patch(patch_id);
        if (patch.type == PatchType::EMPTY || patch.type == PatchType::SYMMETRY)
            continue;
        const auto bc = velocity_bcs.find(patch.name);
        if (bc == velocity_bcs.end())
            throw std::invalid_argument("inlet/outlet face BC: velocity BC missing: " + patch.name);
        const double f = previous_phi(face);
        if (!std::isfinite(f))
            throw std::invalid_argument("inlet/outlet face BC: non-finite authoritative flux");
        const auto inflow = inflow_values.find(patch.name);
        const bool mixed_inflow = inflow != inflow_values.end() && f < 0.0;
        const bool fixed = mixed_inflow || bc->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE;
        const Vec3 value = mixed_inflow ? inflow->second : bc->second.value;
        const std::array<double, 3> components{value.x, value.y, value.z};
        for (std::size_t d=0;d<3;++d) {
            resolved[d].specified[face] = 1;
            resolved[d].conditions[face] = {
                fixed ? ScalarBoundaryFaceCondition::Type::FIXED_VALUE
                      : ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT,
                fixed ? components[d] : 0.0};
        }
    }
    return resolved;
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {},
    bool projected_face_interpolation = false,
    const VelocityFaceConditions* face_velocity_conditions = nullptr)
{
    using namespace cfdx::core;
    (void)pressure_bcs;
    Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
    const auto& own = mesh.ownership();

    for (std::size_t c = 0; c < U.size(); ++c) {
        const double ux = U.component_data(0)[c];
        const double uy = U.component_data(1)[c];
        const double uz = U.component_data(2)[c];
        if (!std::isfinite(ux) || !std::isfinite(uy) || !std::isfinite(uz))
            throw std::runtime_error("make_mass_flux: non-finite cell velocity");
    }
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const std::size_t o=own.owner(f);
        Vec3 Uf;
        U.get(o,Uf.x,Uf.y,Uf.z);
        if (own.neighbour(f)>=0) {
            const std::size_t n=static_cast<std::size_t>(own.neighbour(f));
            Vec3 Un;
            U.get(n,Un.x,Un.y,Un.z);
            const double wn = pressure_face_neighbour_weight(
                geometry,f,o,n,projected_face_interpolation);
            Uf=Uf*(1.0-wn)+Un*wn;
        } else {
            const std::size_t p=geometry.face_patch[f];
            if (p < mesh.boundary().n_patches() &&
                (mesh.boundary().patch(p).type == PatchType::EMPTY ||
                 mesh.boundary().patch(p).type == PatchType::SYMMETRY)) {
                // Symmetry is impermeable, not a generic outflow.
                flux(f) = 0.0;
                continue;
            }
            if(p<mesh.boundary().n_patches()) {
                const auto& name=mesh.boundary().patch(p).name;
                const auto it=bcs.find(name);
                if(it!=bcs.end() && it->second.type==VelocityBoundaryCondition::Type::FIXED_VALUE)
                    Uf=it->second.value;
            }
        }
        if (own.neighbour(f) < 0 && face_velocity_conditions)
            for (std::size_t d=0;d<3;++d)
                if ((*face_velocity_conditions)[d].has(f) &&
                    (*face_velocity_conditions)[d].conditions[f].type ==
                        ScalarBoundaryFaceCondition::Type::FIXED_VALUE) {
                    const double value = (*face_velocity_conditions)[d].conditions[f].value;
                    if (d==0) Uf.x=value;
                    else if (d==1) Uf.y=value;
                    else Uf.z=value;
                }
        flux(f)=rho*Uf.dot(geometry.face_area_vectors[f]);
    }
    return flux;
}

// Orthogonal implicit coefficient shared by the Rhie-Chow face flux and the
// pressure-correction matrix. dAU_diagonal contains the volume-weighted
// directional inverse momentum diagonal (V*rAU), interpolated to the face.
// The implicit pressure difference is along the owner-to-neighbour line, so
// the matching surface measure is Sf.e, not |Sf|. They coincide only on an
// orthogonal face. The deferred non-orthogonal part is handled separately by
// the Rhie-Chow flux reconstruction.
inline double pressure_face_orthogonal_coefficient(
    const cfdx::core::Vec3& Sf,
    const cfdx::core::Vec3& dvec,
    const cfdx::core::Vec3& dAU_diagonal)
{
    const double d = dvec.mag();
    const double area = Sf.mag();
    if (!(d > 0.0) || !std::isfinite(d) || !(area > 0.0) || !std::isfinite(area))
        throw std::runtime_error("pressure_face_orthogonal_coefficient: invalid face geometry");
    const cfdx::core::Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
    const double orthogonal_area = Sf.dot(e);
    if (!(orthogonal_area > 0.0) || !std::isfinite(orthogonal_area))
        throw std::runtime_error(
            "pressure_face_orthogonal_coefficient: face orientation is inconsistent with cell-centre direction");
    const double rfn = dAU_diagonal.x * e.x * e.x +
                       dAU_diagonal.y * e.y * e.y +
                       dAU_diagonal.z * e.z * e.z;
    const double coefficient = rfn * orthogonal_area / d;
    if (!std::isfinite(coefficient))
        throw std::runtime_error("pressure_face_orthogonal_coefficient: non-finite coefficient");
    return coefficient;
}

// Deferred part of the pressure-response face flux. The implicit term uses
// only the centre-line pressure jump; this remainder makes the sum equal to
// (diag(dAU) * grad(p)) dot Sf for a linear pressure field, including skewed
// faces and anisotropic momentum response.
inline double pressure_face_nonorthogonal_flux(
    const cfdx::core::Vec3& Sf,
    const cfdx::core::Vec3& dvec,
    const cfdx::core::Vec3& dAU_diagonal,
    const cfdx::core::Vec3& grad_p)
{
    const double d = dvec.mag();
    if (!(d > 0.0) || !std::isfinite(d))
        throw std::runtime_error("pressure_face_nonorthogonal_flux: invalid centre distance");
    const cfdx::core::Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
    const double orthogonal_area = Sf.dot(e);
    const double rfn = dAU_diagonal.x * e.x * e.x +
                       dAU_diagonal.y * e.y * e.y +
                       dAU_diagonal.z * e.z * e.z;
    const double full_response_flux =
        dAU_diagonal.x * grad_p.x * Sf.x +
        dAU_diagonal.y * grad_p.y * Sf.y +
        dAU_diagonal.z * grad_p.z * Sf.z;
    const double implicit_gradient_flux =
        rfn * grad_p.dot(e) * orthogonal_area;
    const double deferred = full_response_flux - implicit_gradient_flux;
    if (!std::isfinite(deferred))
        throw std::runtime_error("pressure_face_nonorthogonal_flux: non-finite deferred flux");
    return deferred;
}

inline double rhie_chow_pressure_flux_internal(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    std::size_t face,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& grad_p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    bool projected_face_interpolation = false)
{
    using namespace cfdx::core;
    const auto nr = mesh.ownership().neighbour(face);
    if (nr < 0) return 0.0;
    const std::size_t o = mesh.ownership().owner(face);
    const std::size_t n = static_cast<std::size_t>(nr);
    const Vec3 Sf = geometry.face_area_vectors[face];
    const double area = Sf.mag();
    if (!(area > 0.0) || !std::isfinite(area))
        throw std::runtime_error("rhie_chow_pressure_flux_internal: degenerate face");
    const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
    const double d = dvec.mag();
    if (!(d > 0.0) || !std::isfinite(d))
        throw std::runtime_error("rhie_chow_pressure_flux_internal: degenerate centre distance");
    // rAU is the inverse integrated momentum diagonal (1/A_P). The volume
    // factor is introduced here, at the velocity/pressure projection stage,
    // because the pressure gradient is a volumetric momentum source. Thus
    // cx/cy/cz are dAU = V*rAU, not a redefinition of rAU itself.
    const double wn = pressure_face_neighbour_weight(
        geometry, face, o, n, projected_face_interpolation);
    const double wo = 1.0 - wn;
    const double cx = wo * geometry.cell_volumes[o] * rAU[0][o] +
                      wn * geometry.cell_volumes[n] * rAU[0][n];
    const double cy = wo * geometry.cell_volumes[o] * rAU[1][o] +
                      wn * geometry.cell_volumes[n] * rAU[1][n];
    const double cz = wo * geometry.cell_volumes[o] * rAU[2][o] +
                      wn * geometry.cell_volumes[n] * rAU[2][n];
    const Vec3 gp{
        wo * grad_p.component_data(0)[o] + wn * grad_p.component_data(0)[n],
        wo * grad_p.component_data(1)[o] + wn * grad_p.component_data(1)[n],
        wo * grad_p.component_data(2)[o] + wn * grad_p.component_data(2)[n]};
    // The face response uses the same dAU = V/A_P operator as the cell
    // velocity reconstruction. Keeping the orthogonal pressure difference
    // implicit and the non-orthogonal remainder deferred is important: the
    // continuity matrix and the conservative Rhie-Chow flux must represent
    // the same discrete operator, otherwise a residual floor can survive
    // even when the momentum equations appear converged.
    // Orthogonal pressure difference is implicit in the pressure matrix.
    // The non-orthogonal remainder is deferred explicitly using the same
    // Gauss gradient. This keeps the face flux and pressure equation on the
    // same discrete operator on arbitrary meshes.
    const double implicit_coefficient = pressure_face_orthogonal_coefficient(
        Sf, dvec, {cx, cy, cz});
    const double non_orthogonal_flux = pressure_face_nonorthogonal_flux(
        Sf, dvec, {cx, cy, cz}, gp);
    return rho * (
        implicit_coefficient * (p(n) - p(o)) +
        non_orthogonal_flux);
}

enum class RhieChowVelocityContract {
    HBYA_PREDICTOR,
    PHYSICAL_VELOCITY
};

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {},
    RhieChowVelocityContract velocity_contract =
        RhieChowVelocityContract::HBYA_PREDICTOR,
    bool skew_correct_gradient = false,
    bool projected_face_interpolation = false,
    const VelocityFaceConditions* face_velocity_conditions = nullptr)
;

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::vector<double>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs)
{
    std::array<std::vector<double>, 3> directional_rAU{rAU, rAU, rAU};
    return make_rhie_chow_mass_flux(
        mesh, geometry, U, p, directional_rAU, rho, bcs, {});
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    RhieChowVelocityContract velocity_contract,
    bool skew_correct_gradient,
    bool projected_face_interpolation,
    const VelocityFaceConditions* face_velocity_conditions)
{
    using namespace cfdx::core;
    const bool coupled_physical_velocity =
        velocity_contract == RhieChowVelocityContract::PHYSICAL_VELOCITY;
    const auto& own = mesh.ownership();
    if (p.size() != mesh.n_cells())
        throw std::invalid_argument("make_rhie_chow_mass_flux: field size mismatch");
    // Boundary pressure is passed explicitly to the gradient operator. A
    // boundary face has no neighbour coefficient: any dAU/rAU interpolation
    // must therefore use the owner-cell value rather than a ghost/zero entry.
    const auto grad_p = gauss_gradient_with_boundary(
        p, mesh, geometry, pressure_bcs, skew_correct_gradient,
        projected_face_interpolation);
    for (const auto& component : rAU) {
        if (component.size() != mesh.n_cells())
            throw std::invalid_argument("make_rhie_chow_mass_flux: inverse momentum diagonal size mismatch");
        for (const double value : component) {
            if (!std::isfinite(value) || value <= 0.0)
                throw std::invalid_argument(
                    "make_rhie_chow_mass_flux: inverse momentum diagonal must be finite and positive");
        }
    }

    // Momentum-weighted (Rhie-Chow) face flux:
    //
    //     phi_f = rho [ H/A |_f . Sf - rAU_f * grad(p)_f . Sf ]
    //
    // The first term is linearly interpolated from the cell-centred predictor.
    // The implicit pressure jump is taken along the owner-neighbour centre
    // line and multiplied by Sf.e, not |Sf|. The tangential remainder is a
    // deferred non-orthogonal correction. The same implicit coefficient is
    // assembled in the pressure-correction matrix below; it must never be
    // replaced by an unrelated least-squares fit of the cell velocity.
    //
    // For a fixed-pressure boundary the boundary pressure is known. For a
    // zero-gradient pressure boundary there is no normal pressure correction.
    Field<double, Location::FACE> flux(mesh.n_faces(), "phi", "kg/s", 1);
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const std::size_t o = own.owner(f);
        Vec3 Hf;
        U.get(o, Hf.x, Hf.y, Hf.z);
        const auto nr = own.neighbour(f);
        double correction = 0.0;
        if (nr >= 0) {
            const std::size_t n = static_cast<std::size_t>(nr);
            Vec3 Hn;
            U.get(n, Hn.x, Hn.y, Hn.z);
            const double wn = pressure_face_neighbour_weight(
                geometry, f, o, n, projected_face_interpolation);
            Hf = Hf * (1.0 - wn) + Hn * wn;
            correction = rhie_chow_pressure_flux_internal(
                mesh, geometry, f, p, grad_p, rAU, rho,
                projected_face_interpolation);
            if (coupled_physical_velocity) {
                // In the monolithic block U is the physical velocity, not H/A.
                // Add the interpolated dAU*grad(p) compensation so an affine
                // pressure field is not subtracted twice from the face flux.
                const double cx = (1.0-wn) * geometry.cell_volumes[o] * rAU[0][o] +
                                  wn * geometry.cell_volumes[n] * rAU[0][n];
                const double cy = (1.0-wn) * geometry.cell_volumes[o] * rAU[1][o] +
                                  wn * geometry.cell_volumes[n] * rAU[1][n];
                const double cz = (1.0-wn) * geometry.cell_volumes[o] * rAU[2][o] +
                                  wn * geometry.cell_volumes[n] * rAU[2][n];
                const Vec3 gp{
                    (1.0-wn)*grad_p.component_data(0)[o]+wn*grad_p.component_data(0)[n],
                    (1.0-wn)*grad_p.component_data(1)[o]+wn*grad_p.component_data(1)[n],
                    (1.0-wn)*grad_p.component_data(2)[o]+wn*grad_p.component_data(2)[n]};
                const Vec3 Sf = geometry.face_area_vectors[f];
                correction -= rho * (cx * gp.x * Sf.x +
                                     cy * gp.y * Sf.y +
                                     cz * gp.z * Sf.z);
            }
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                (mesh.boundary().patch(patch).type == PatchType::EMPTY ||
                 mesh.boundary().patch(patch).type == PatchType::SYMMETRY)) {
                // Both EMPTY and SYMMETRY faces prohibit normal mass flow.
                flux(f) = 0.0;
                continue;
            }
            if (patch < mesh.boundary().n_patches()) {
                const auto& name = mesh.boundary().patch(patch).name;
                const auto ubc = bcs.find(name);
                if (ubc != bcs.end() &&
                    ubc->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    Hf = ubc->second.value;
                if (face_velocity_conditions)
                    for (std::size_t d=0;d<3;++d)
                        if ((*face_velocity_conditions)[d].has(f) &&
                            (*face_velocity_conditions)[d].conditions[f].type ==
                                ScalarBoundaryFaceCondition::Type::FIXED_VALUE) {
                            const double value = (*face_velocity_conditions)[d].conditions[f].value;
                            if (d==0) Hf.x=value;
                            else if (d==1) Hf.y=value;
                            else Hf.z=value;
                        }

                const auto pbc = pressure_bcs.find(name);
                if (pbc != pressure_bcs.end() &&
                    pbc->second.type == ScalarBoundaryType::FIXED_VALUE) {
                    const std::size_t o = own.owner(f);
                    const Vec3 Sf = geometry.face_area_vectors[f];
                    const Vec3 dvec =
                        geometry.face_centres[f] - geometry.cell_centres[o];
                    const double d = dvec.mag();
                    if (!(d > 0.0) || !std::isfinite(d))
                        throw std::runtime_error("make_rhie_chow_mass_flux: degenerate boundary distance");
                    const Vec3 gp{
                        grad_p.component_data(0)[o],
                        grad_p.component_data(1)[o],
                        grad_p.component_data(2)[o]};
                    const Vec3 dAU_boundary{
                        geometry.cell_volumes[o] * rAU[0][o],
                        geometry.cell_volumes[o] * rAU[1][o],
                        geometry.cell_volumes[o] * rAU[2][o]};
                    const double implicit_coefficient =
                        pressure_face_orthogonal_coefficient(
                            Sf, dvec, dAU_boundary);
                    const double deferred = pressure_face_nonorthogonal_flux(
                        Sf, dvec, dAU_boundary, gp);
                    correction = rho * (
                        implicit_coefficient * (pbc->second.value - p(o)) +
                        deferred);
                    if (coupled_physical_velocity) {
                        const Vec3 gp{
                            grad_p.component_data(0)[o],
                            grad_p.component_data(1)[o],
                            grad_p.component_data(2)[o]};
                        correction -= rho * (
                            dAU_boundary.x * gp.x * Sf.x +
                            dAU_boundary.y * gp.y * Sf.y +
                            dAU_boundary.z * gp.z * Sf.z);
                    }
                }
            }
        }
        flux(f) = rho * Hf.dot(geometry.face_area_vectors[f]) - correction;
    }
    return flux;
}

// Experimental audit of the COMPLETE pressure-dependent Rhie-Chow face
// response of the same Picard-frozen SIMPLE state.  Unlike an implicit-face
// residual, this also includes the deferred/nonorthogonal Green-Gauss term.
// It is a diagnostic: the public FV fields are still FP64 and the nonlinear
// momentum/SST acceptance gates are not replaced.  Skewed multi-sweep
// gradients, SIMPLEC and COUPLED need separate linearization contracts.
//
// The derivative is assembled from the high/low *increment* using extended
// precision, so neither (high+low) nor p_old+increment is rounded before
// computing pressure gradients and differences.  Each face is written ONCE;
// owner/neighbour divergence is accumulated with opposite signs.
inline std::vector<long double> reconstruct_compensated_rhie_chow_pressure_step(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::vector<long double>& base_flux,
    const cfdx::core::Vector& high, const cfdx::core::Vector& low,
    const std::array<std::vector<double>,3>& rAU, double rho,
    const ScalarBoundaryConditions& pressure_bcs,
    bool projected_face_interpolation = false)
{
    using namespace cfdx::core;
    const std::size_t nc=mesh.n_cells(), nf=mesh.n_faces();
    if (high.size()!=nc || low.size()!=nc || base_flux.size()!=nf ||
        !(rho>0.0) || !std::isfinite(rho) ||
        std::numeric_limits<long double>::digits <=
            std::numeric_limits<double>::digits)
        throw std::invalid_argument("compensated RC: invalid dimensions/precision/density");
    for (const auto& component : rAU) {
        if (component.size()!=nc)
            throw std::invalid_argument("compensated RC: rAU dimension mismatch");
        for (double value : component)
            if (!std::isfinite(value) || !(value>0.0))
                throw std::invalid_argument("compensated RC: invalid inverse momentum diagonal");
    }
    const auto& ownership=mesh.ownership();
    const auto* offsets=mesh.cells().offsets_data();
    const auto* cell_faces=mesh.cells().faces_data();
    std::vector<std::array<long double,3>> gradient(nc);
    for (std::size_t c=0;c<nc;++c) {
        if (!std::isfinite(high(c)) || !std::isfinite(low(c)) ||
            !(geometry.cell_volumes[c]>0.0))
            throw std::invalid_argument("compensated RC: invalid pressure/volume");
        std::array<long double,3> sum{0.0L,0.0L,0.0L};
        for (Offset k=offsets[c];k<offsets[c+1];++k) {
            const auto f=cell_faces[k];
            const bool is_owner=ownership.owner(f)==c;
            const auto Sf=is_owner ? geometry.face_area_vectors[f] :
                                      geometry.face_area_vectors[f]*-1.0;
            const auto nr=ownership.neighbour(f);
            long double vf=static_cast<long double>(high(c))+low(c);
            if (nr>=0) {
                const std::size_t other=is_owner ? static_cast<std::size_t>(nr) :
                                                      ownership.owner(f);
                double w;
                if (projected_face_interpolation) {
                    const double wn=geometric_interpolation_weight(
                        geometry.face_centres[f], geometry.cell_centres[c],
                        &geometry.cell_centres[other]);
                    w=1.0-wn;
                } else {
                    const double dc=(geometry.face_centres[f]-geometry.cell_centres[c]).mag();
                    const double dn=(geometry.cell_centres[other]-geometry.face_centres[f]).mag();
                    if (!(dc+dn>0.0))
                        throw std::runtime_error("compensated RC: zero face distance");
                    w=dn/(dc+dn);
                }
                vf=static_cast<long double>(w)*
                       (static_cast<long double>(high(c))+low(c)) +
                   static_cast<long double>(1.0-w)*
                       (static_cast<long double>(high(other))+low(other));
            } else {
                const auto patch=geometry.face_patch[f];
                if (patch<mesh.boundary().n_patches() &&
                    mesh.boundary().patch(patch).type==PatchType::EMPTY)
                    continue;
                if (patch<mesh.boundary().n_patches()) {
                    const auto& name=mesh.boundary().patch(patch).name;
                    const auto bc=pressure_bcs.find(name);
                    // A correction field has homogeneous fixed-value BC and
                    // homogeneous fixed-gradient BC. Both must be zeroed.
                    if (bc!=pressure_bcs.end() &&
                        bc->second.type==ScalarBoundaryType::FIXED_VALUE)
                        vf=0.0L;
                }
            }
            sum[0]+=static_cast<long double>(Sf.x)*vf;
            sum[1]+=static_cast<long double>(Sf.y)*vf;
            sum[2]+=static_cast<long double>(Sf.z)*vf;
        }
        for (std::size_t d=0;d<3;++d) {
            gradient[c][d]=sum[d]/geometry.cell_volumes[c];
            if (!std::isfinite(gradient[c][d]))
                throw std::domain_error("compensated RC: nonfinite pressure gradient");
        }
    }
    std::vector<long double> flux(nf);
    const auto deferred = [&](const Vec3& Sf, const Vec3& dvec,
                              const Vec3& dAU,
                              const std::array<long double,3>& gp) {
        const double d=dvec.mag();
        if (!(d>0.0)) throw std::runtime_error("compensated RC: degenerate face");
        const Vec3 e=dvec*(1.0/d);
        const double ge=dAU.x*e.x*e.x+dAU.y*e.y*e.y+dAU.z*e.z*e.z;
        const double orth_area=Sf.dot(e);
        const long double projected_grad=
            gp[0]*e.x+gp[1]*e.y+gp[2]*e.z;
        return static_cast<long double>(dAU.x)*gp[0]*Sf.x+
               static_cast<long double>(dAU.y)*gp[1]*Sf.y+
               static_cast<long double>(dAU.z)*gp[2]*Sf.z-
               static_cast<long double>(ge)*projected_grad*orth_area;
    };
    for (std::size_t f=0;f<nf;++f) {
        const auto o=ownership.owner(f);
        const auto nr=ownership.neighbour(f);
        const Vec3 Sf=geometry.face_area_vectors[f];
        long double response=0.0L;
        if (nr>=0) {
            const auto n=static_cast<std::size_t>(nr);
            const auto dvec=geometry.cell_centres[n]-geometry.cell_centres[o];
            const double wn=pressure_face_neighbour_weight(
                geometry,f,o,n,projected_face_interpolation);
            const double wo=1.0-wn;
            const Vec3 dAU{
                wo*geometry.cell_volumes[o]*rAU[0][o]+wn*geometry.cell_volumes[n]*rAU[0][n],
                wo*geometry.cell_volumes[o]*rAU[1][o]+wn*geometry.cell_volumes[n]*rAU[1][n],
                wo*geometry.cell_volumes[o]*rAU[2][o]+wn*geometry.cell_volumes[n]*rAU[2][n]};
            const double orth=pressure_face_orthogonal_coefficient(Sf,dvec,dAU);
            const long double jump=
                (static_cast<long double>(high(n))-high(o))+
                (static_cast<long double>(low(n))-low(o));
            const std::array<long double,3> gp{
                wo*gradient[o][0]+wn*gradient[n][0],
                wo*gradient[o][1]+wn*gradient[n][1],
                wo*gradient[o][2]+wn*gradient[n][2]};
            response=static_cast<long double>(rho)*
                (static_cast<long double>(orth)*jump+deferred(Sf,dvec,dAU,gp));
        } else {
            const auto patch=geometry.face_patch[f];
            if (patch<mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type!=PatchType::EMPTY &&
                mesh.boundary().patch(patch).type!=PatchType::SYMMETRY) {
                const auto bc=pressure_bcs.find(mesh.boundary().patch(patch).name);
                if (bc!=pressure_bcs.end() &&
                    bc->second.type==ScalarBoundaryType::FIXED_VALUE) {
                    const Vec3 dvec=geometry.face_centres[f]-geometry.cell_centres[o];
                    const Vec3 dAU{geometry.cell_volumes[o]*rAU[0][o],
                                   geometry.cell_volumes[o]*rAU[1][o],
                                   geometry.cell_volumes[o]*rAU[2][o]};
                    const double orth=pressure_face_orthogonal_coefficient(Sf,dvec,dAU);
                    const long double jump=-static_cast<long double>(high(o))-low(o);
                    response=static_cast<long double>(rho)*
                        (static_cast<long double>(orth)*jump+
                         deferred(Sf,dvec,dAU,gradient[o]));
                }
            }
        }
        flux[f]=base_flux[f]-response;
        if (!std::isfinite(flux[f]))
            throw std::domain_error("compensated RC: nonfinite reconstructed flux");
    }
    return flux;
}

// Field overload retained for all V5 diagnostics and existing callers.
inline std::vector<long double> reconstruct_compensated_rhie_chow_pressure_step(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& base_flux,
    const cfdx::core::Vector& high, const cfdx::core::Vector& low,
    const std::array<std::vector<double>,3>& rAU, double rho,
    const ScalarBoundaryConditions& pressure_bcs,
    bool projected_face_interpolation = false)
{
    if (base_flux.dimension()!=1 || base_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("compensated RC: invalid base flux");
    std::vector<long double> base(mesh.n_faces());
    for (std::size_t f=0;f<base.size();++f) base[f]=base_flux(f);
    return reconstruct_compensated_rhie_chow_pressure_step(
        mesh,geometry,base,high,low,rAU,rho,pressure_bcs,
        projected_face_interpolation);
}

// Canonical physical face operator.  A zero CELL pressure establishes the
// actual inhomogeneous boundary flux; the split cell pressure is then applied
// as a homogeneous pressure response. No large p+dp is rounded for the flux.
inline std::vector<long double> make_rhie_chow_mass_flux_high_low(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& HbyA,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& pressure_high,
    const cfdx::core::Vector& pressure_low,
    const std::array<std::vector<double>,3>& rAU, double density,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    bool projected_face_interpolation = false,
    const VelocityFaceConditions* face_velocity_conditions = nullptr)
{
    using namespace cfdx::core;
    if (pressure_high.size()!=mesh.n_cells() ||
        pressure_low.size()!=mesh.n_cells())
        throw std::invalid_argument("RC high/low: invalid pressure state");
    Field<double,Location::CELL> zero(mesh.n_cells(),"pressure_zero","Pa",1);
    zero.fill(0.0);
    const auto zero_pressure_flux=make_rhie_chow_mass_flux(
        mesh,geometry,HbyA,zero,rAU,density,velocity_bcs,pressure_bcs,
        RhieChowVelocityContract::HBYA_PREDICTOR,false,
        projected_face_interpolation,face_velocity_conditions);
    Vector high(mesh.n_cells(),0.0);
    for (std::size_t c=0;c<mesh.n_cells();++c) high(c)=pressure_high(c);
    return reconstruct_compensated_rhie_chow_pressure_step(
        mesh,geometry,zero_pressure_flux,high,pressure_low,rAU,density,
        pressure_bcs,projected_face_interpolation);
}

// True signed FV residual: -div(phi), assembled once per oriented face and
// normalized by the ORIGINAL pressure RHS. This remains distinct from the
// preconditioner CSR residual and the accepted double transport flux.
inline std::pair<cfdx::core::Vector,double> pressure_true_face_residual(
    const cfdx::core::Mesh& mesh, const std::vector<long double>& flux,
    const cfdx::core::Vector& original_rhs)
{
    using namespace cfdx::core;
    if (flux.size()!=mesh.n_faces() || original_rhs.size()!=mesh.n_cells())
        throw std::invalid_argument("RC true face residual: dimensions mismatch");
    std::vector<long double> divergence(mesh.n_cells(),0.0L);
    for (std::size_t f=0;f<flux.size();++f) {
        if (!std::isfinite(flux[f]))
            throw std::domain_error("RC true face residual: invalid flux");
        const auto o=mesh.ownership().owner(f);
        divergence[o]+=flux[f];
        const auto nr=mesh.ownership().neighbour(f);
        if (nr>=0) divergence[static_cast<std::size_t>(nr)]-=flux[f];
    }
    Vector residual(mesh.n_cells(),0.0);
    long double residual2=0.0L, b2=0.0L;
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        const long double r=-divergence[c];
        if (!std::isfinite(r) || std::abs(r)>std::numeric_limits<double>::max() ||
            !std::isfinite(original_rhs(c)))
            throw std::domain_error("RC true face residual: invalid cell balance");
        residual(c)=static_cast<double>(r);
        residual2+=r*r;
        b2+=static_cast<long double>(original_rhs(c))*original_rhs(c);
    }
    if (!(b2>0.0L))
        throw std::invalid_argument("RC true face residual: zero reference RHS");
    const double relative=static_cast<double>(std::sqrt(residual2/b2));
    if(!std::isfinite(relative))
        throw std::domain_error("RC true face residual: invalid norm");
    return {std::move(residual),relative};
}

struct PressureTrueFaceRefinementReport {
    bool converged=false;
    std::size_t attempted=0, accepted=0;
    double initial_relative=std::numeric_limits<double>::infinity();
    double final_relative=std::numeric_limits<double>::infinity();
};

// The sparse Jacobian is a PRECONDITIONER for correction steps, never a
// physical success oracle. Every trial is validated on the full RC face flux.
// A nonconverged linear inner solve is permissible as a tentative direction,
// but the final true-FV tolerance is never relaxed or bypassed.
template <class FaceResidual, class LinearIncrement>
PressureTrueFaceRefinementReport refine_pressure_true_face_high_low(
    cfdx::core::Vector& high,cfdx::core::Vector& low,double tolerance,
    std::size_t max_passes,FaceResidual&& physical_residual,
    LinearIncrement&& linear_increment)
{
    using namespace cfdx::core;
    if (high.size()==0 || high.size()!=low.size() ||
        !(tolerance>0.0) || !std::isfinite(tolerance) || max_passes>12)
        throw std::invalid_argument("RC true face refinement: invalid arguments");
    auto [r,relative]=physical_residual(high,low);
    if(r.size()!=high.size() || !std::isfinite(relative))
        throw std::domain_error("RC true face refinement: invalid initial residual");
    PressureTrueFaceRefinementReport report;
    report.initial_relative=relative;
    report.final_relative=relative;
    for(std::size_t k=0;k<max_passes && relative>tolerance;++k) {
        ++report.attempted;
        Vector delta(high.size(),0.0);
        if (!linear_increment(r,delta)) break;
        bool accepted_step=false;
        // Backtrack the pressure defect-correction direction on the independent
        // physical full-face residual, not on FGMRES's approximate CSR norm.
        // Geometric steps require no arbitrary relaxation of the final gate.
        for (unsigned backtrack=0;backtrack<8;++backtrack) {
            const double step=std::ldexp(1.0,-static_cast<int>(backtrack));
            Vector candidate_high=high,candidate_low=low;
            bool finite=true;
            for(std::size_t c=0;c<high.size();++c)
                finite=compensated_add(candidate_high(c),candidate_low(c),
                                       step*delta(c)) && finite;
            if(!finite) continue;
            auto [trial_r,trial_relative]=physical_residual(
                candidate_high,candidate_low);
            if(!std::isfinite(trial_relative) || !(trial_relative<relative))
                continue;
            high=std::move(candidate_high);
            low=std::move(candidate_low);
            r=std::move(trial_r);
            relative=trial_relative;
            ++report.accepted;
            report.final_relative=relative;
            accepted_step=true;
            break;
        }
        if (!accepted_step) break;
    }
    report.converged=relative<=tolerance;
    return report;
}

// One-sweep Green-Gauss pressure gradient with no loss of the low component.
// This is the exact pressure-gradient contract used by the experimental full
// RC Jacobian.  Skew/multi-sweep is intentionally unsupported here.
inline cfdx::core::Field<double,cfdx::core::Location::CELL>
gauss_gradient_pressure_high_low(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& high,
    const cfdx::core::Vector& low,const ScalarBoundaryConditions& bcs,
    bool projected=false)
{
    using namespace cfdx::core;
    if (high.dimension()!=1 || high.size()!=mesh.n_cells() ||
        low.size()!=mesh.n_cells())
        throw std::invalid_argument("RC high/low gradient: invalid dimensions");
    Field<double,Location::CELL> grad(mesh.n_cells(),"pressure_hilo_gradient","Pa/m",3);
    const auto& own=mesh.ownership();
    const auto* offsets=mesh.cells().offsets_data();
    const auto* faces=mesh.cells().faces_data();
    for(std::size_t c=0;c<mesh.n_cells();++c) {
        if(!std::isfinite(high(c)) || !std::isfinite(low(c)) ||
           !(geometry.cell_volumes[c]>0.0))
            throw std::domain_error("RC high/low gradient: invalid pressure/volume");
        std::array<long double,3> sum{0.0L,0.0L,0.0L};
        for(Offset at=offsets[c];at<offsets[c+1];++at) {
            const std::size_t f=faces[at];
            const bool is_owner=own.owner(f)==c;
            const Vec3 Sf=is_owner ? geometry.face_area_vectors[f] :
                                     geometry.face_area_vectors[f]*-1.0;
            long double vf=static_cast<long double>(high(c))+low(c);
            const auto nr=own.neighbour(f);
            if(nr>=0) {
                const std::size_t other=is_owner ? static_cast<std::size_t>(nr) : own.owner(f);
                double w=0.0;
                if(projected) {
                    const double wn=geometric_interpolation_weight(
                        geometry.face_centres[f],geometry.cell_centres[c],
                        &geometry.cell_centres[other]);
                    w=1.0-wn;
                } else {
                    const double dc=(geometry.face_centres[f]-geometry.cell_centres[c]).mag();
                    const double dn=(geometry.cell_centres[other]-geometry.face_centres[f]).mag();
                    if(!(dc+dn>0.0)) throw std::invalid_argument("RC high/low gradient: degenerate face");
                    w=dn/(dc+dn);
                }
                vf=static_cast<long double>(w)*vf+
                    static_cast<long double>(1.0-w)*
                        (static_cast<long double>(high(other))+low(other));
            } else {
                const auto patch=geometry.face_patch[f];
                if(patch<mesh.boundary().n_patches() &&
                   mesh.boundary().patch(patch).type==PatchType::EMPTY)
                    continue;
                if(patch<mesh.boundary().n_patches()) {
                    const auto it=bcs.find(mesh.boundary().patch(patch).name);
                    if(it!=bcs.end() && it->second.type==ScalarBoundaryType::FIXED_VALUE)
                        vf=it->second.value;
                    else if(it!=bcs.end() && it->second.type==ScalarBoundaryType::FIXED_GRADIENT) {
                        const Vec3 unit=Sf*(1.0/Sf.mag());
                        vf+=it->second.gradient*(geometry.face_centres[f]-geometry.cell_centres[c]).dot(unit);
                    }
                }

            }
            sum[0]+=static_cast<long double>(Sf.x)*vf;
            sum[1]+=static_cast<long double>(Sf.y)*vf;
            sum[2]+=static_cast<long double>(Sf.z)*vf;
        }
        for(std::size_t d=0;d<3;++d) {
            const long double v=sum[d]/geometry.cell_volumes[c];
            if(!std::isfinite(v) || std::abs(v)>std::numeric_limits<double>::max())
                throw std::domain_error("RC high/low gradient: nonfinite component");
            grad.component_data(d)[c]=static_cast<double>(v);
        }
    }
    return grad;
}

// An under-relaxed physical pressure update MUST retain its roundoff part;
// an unrelaxed flux pressure is distinct from the momentum pressure in SIMPLE.
inline void add_pressure_high_low(
    cfdx::core::Field<double,cfdx::core::Location::CELL>& high,
    cfdx::core::Vector& low,const cfdx::core::Vector& correction_high,
    const cfdx::core::Vector& correction_low,double alpha)
{
    if (high.size()!=low.size() || low.size()!=correction_high.size() ||
        low.size()!=correction_low.size() || !(alpha>0.0 && alpha<=1.0) ||
        !std::isfinite(alpha) ||
        std::numeric_limits<long double>::digits<=std::numeric_limits<double>::digits)
        throw std::invalid_argument("RC high/low update: invalid arguments");
    for(std::size_t c=0;c<high.size();++c) {
        // Never evaluate (high+low)+increment in long double: even 80-bit
        // extended precision loses bits of the low part when |high| is large.
        // Split the alpha-scaled correction before adding using TwoSum.
        const long double increment = static_cast<long double>(alpha) *
            (static_cast<long double>(correction_high(c)) +
             static_cast<long double>(correction_low(c)));
        if (!std::isfinite(increment) ||
            std::abs(increment)>std::numeric_limits<double>::max())
            throw std::domain_error("RC high/low update: invalid increment");
        const double increment_high=static_cast<double>(increment);
        const double increment_low=static_cast<double>(
            increment-static_cast<long double>(increment_high));
        if (!cfdx::core::compensated_add(high(c),low(c),increment_high) ||
            !cfdx::core::compensated_add(high(c),low(c),increment_low))
            throw std::domain_error("RC high/low update: nonfinite pressure");
    }
}

inline double compensated_flux_continuity_relative(
    const cfdx::core::Mesh& mesh, const std::vector<long double>& flux,
    const cfdx::core::Vector& normalization_rhs)
{
    if (flux.size()!=mesh.n_faces() || normalization_rhs.size()!=mesh.n_cells())
        throw std::invalid_argument("compensated RC continuity: dimensions mismatch");
    std::vector<long double> balance(mesh.n_cells(),0.0L);
    long double norm=0.0L, rhs_norm=0.0L;
    for (std::size_t f=0;f<mesh.n_faces();++f) {
        if (!std::isfinite(flux[f]))
            throw std::domain_error("compensated RC continuity: nonfinite face flux");
        const auto o=mesh.ownership().owner(f);
        balance[o]+=flux[f];
        const auto n=mesh.ownership().neighbour(f);
        if (n>=0) balance[static_cast<std::size_t>(n)]-=flux[f];
    }
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        const long double r=normalization_rhs(c);
        rhs_norm+=r*r;
        norm+=balance[c]*balance[c];
    }
    if (!(rhs_norm>0.0L))
        throw std::invalid_argument("compensated RC continuity: zero normalization RHS");
    return static_cast<double>(std::sqrt(norm/rhs_norm));
}

// V20: decompose the INDEPENDENT face-flux discrepancy into two mechanisms.
// phi_transport is conservative and uses the full pressure step in V19.
// phi_relaxed uses the identical HbyA, rAU, BCs, and RC reconstruction but
// the pressure used by the momentum-corrected cell velocity.
// phi_velocity is deliberately built using a DIFFERENT operator (linear U
// interpolation).  Thus:
//   phi_transport - phi_velocity = (phi_transport - phi_relaxed)
//                                  + (phi_relaxed - phi_velocity).
// The first component is the pressure-relaxation split, NOT a conservation
// defect. The second is the collocated Rhie-Chow versus U interpolation
// discrepancy, which must be qualified separately near walls and outlets.
struct V20FluxVelocityDecomposition {
    double total_linf=0.0, pressure_split_linf=0.0, rhie_chow_linf=0.0;
    double total_l1=0.0, pressure_split_l1=0.0, rhie_chow_l1=0.0;
    double triangle_error_linf=0.0;
    double internal_total_linf=0.0, boundary_total_linf=0.0;
    std::size_t total_worst_face=0, split_worst_face=0, rc_worst_face=0;
    bool total_worst_boundary=false;
};

inline V20FluxVelocityDecomposition audit_v20_flux_velocity_decomposition(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& phi_transport,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& phi_relaxed,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& phi_velocity)
{
    if (phi_transport.size()!=mesh.n_faces() ||
        phi_relaxed.size()!=mesh.n_faces() ||
        phi_velocity.size()!=mesh.n_faces())
        throw std::invalid_argument("V20 face flux decomposition: inconsistent size");
    V20FluxVelocityDecomposition result;
    for (std::size_t face=0; face<mesh.n_faces(); ++face) {
        const double total=phi_transport(face)-phi_velocity(face);
        const double split=phi_transport(face)-phi_relaxed(face);
        const double rc=phi_relaxed(face)-phi_velocity(face);
        if (!std::isfinite(total) || !std::isfinite(split) || !std::isfinite(rc))
            throw std::domain_error("V20 face flux decomposition: nonfinite flux");
        const double a=std::abs(total), b=std::abs(split), c=std::abs(rc);
        result.total_l1+=a;
        result.pressure_split_l1+=b;
        result.rhie_chow_l1+=c;
        if (a>result.total_linf) {
            result.total_linf=a;
            result.total_worst_face=face;
            result.total_worst_boundary=mesh.ownership().neighbour(face)<0;
        }
        if (b>result.pressure_split_linf) {
            result.pressure_split_linf=b;
            result.split_worst_face=face;
        }
        if (c>result.rhie_chow_linf) {
            result.rhie_chow_linf=c;
            result.rc_worst_face=face;
        }
        if (mesh.ownership().neighbour(face)>=0)
            result.internal_total_linf=std::max(result.internal_total_linf,a);
        else
            result.boundary_total_linf=std::max(result.boundary_total_linf,a);
        result.triangle_error_linf=std::max(
            result.triangle_error_linf,std::abs(total-(split+rc)));
    }
    return result;
}

// V18 exact sparse Schur Jacobian of the existing SIMPLE Rhie-Chow flux.
// The Green-Gauss gradient and each face flux are linear in pressure when
// gradient skew sweeps are disabled. Construct their local derivatives, not
// a second finite-volume discretization. The resulting (generally
// nonsymmetric) matrix is div(dphi/dp), including the deferred gradient.
// This is opt-in while mesh and nonlinear qualification remain open.
inline std::vector<std::map<std::size_t,double>>
make_full_rhie_chow_pressure_rows_sparse(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::array<std::vector<double>,3>& rAU, double density,
    const ScalarBoundaryConditions& pressure_bcs,
    bool projected_face_interpolation=false)
{
    using namespace cfdx::core;
    const std::size_t nc=mesh.n_cells();
    if (!std::isfinite(density) || density<=0.0)
        throw std::invalid_argument("V18 sparse pressure Jacobian: invalid density");
    for (const auto& component:rAU) {
        if (component.size()!=nc)
            throw std::invalid_argument("V18 sparse pressure Jacobian: invalid rAU size");
        for (double value:component)
            if (!std::isfinite(value) || value<=0.0)
                throw std::invalid_argument("V18 sparse pressure Jacobian: invalid rAU");
    }
    // Green-Gauss linear map grad(p)_c = sum_j G[c][j]*p_j + constant.
    // Boundary fixed values and fixed gradients contribute only constants;
    // ZERO_GRADIENT and symmetry faces use the owner pressure.
    std::vector<std::map<std::size_t,Vec3>> G(nc);
    for(std::size_t c=0;c<nc;++c) {
        const auto begin=mesh.cells().offsets_data()[c];
        const auto end=mesh.cells().offsets_data()[c+1];
        for(auto at=begin;at<end;++at) {
            const std::size_t f=mesh.cells().faces_data()[at];
            const bool owner=mesh.ownership().owner(f)==c;
            const Vec3 Sf=(owner ? geometry.face_area_vectors[f] :
                                    geometry.face_area_vectors[f]*-1.0) *
                                    (1.0/geometry.cell_volumes[c]);
            const auto other_index=mesh.ownership().neighbour(f);
            if(other_index>=0) {
                const std::size_t other=owner
                    ? static_cast<std::size_t>(other_index)
                    : mesh.ownership().owner(f);
                double w=0.0;
                if (projected_face_interpolation) {
                    const double wn=cfdx::core::geometric_interpolation_weight(
                        geometry.face_centres[f],geometry.cell_centres[c],
                        &geometry.cell_centres[other]);
                    w=1.0-wn;
                } else {
                    const double dc=(geometry.face_centres[f]-geometry.cell_centres[c]).mag();
                    const double dn=(geometry.cell_centres[other]-geometry.face_centres[f]).mag();
                    if(!(dc+dn>0.0))
                        throw std::invalid_argument("V18 sparse pressure Jacobian: zero distance");
                    w=dn/(dc+dn);
                }
                G[c][c]+=Sf*w;
                G[c][other]+=Sf*(1.0-w);
                continue;
            }
            const auto patch=geometry.face_patch[f];
            if(patch<mesh.boundary().n_patches() &&
               mesh.boundary().patch(patch).type==PatchType::EMPTY)
                continue;
            bool fixed_value=false;
            if(patch<mesh.boundary().n_patches()) {
                const auto& name=mesh.boundary().patch(patch).name;
                const auto bc=pressure_bcs.find(name);
                fixed_value=bc!=pressure_bcs.end() &&
                            bc->second.type==ScalarBoundaryType::FIXED_VALUE;
            }
            if(!fixed_value) G[c][c]+=Sf;
        }
    }
    std::vector<std::map<std::size_t,double>> rows(nc);
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const auto o=mesh.ownership().owner(f);
        const auto nr=mesh.ownership().neighbour(f);
        const auto Sf=geometry.face_area_vectors[f];
        auto add_face_derivative=[&](std::size_t j,double dphi) {
            rows[o][j]+=dphi;
            if(nr>=0) rows[static_cast<std::size_t>(nr)][j]-=dphi;
        };
        if(nr>=0) {
            const std::size_t n=static_cast<std::size_t>(nr);
            const auto dvec=geometry.cell_centres[n]-geometry.cell_centres[o];
            const double wn=pressure_face_neighbour_weight(
                geometry,f,o,n,projected_face_interpolation);
            const double wo=1.0-wn;
            const Vec3 dAU{
                wo*geometry.cell_volumes[o]*rAU[0][o]+wn*geometry.cell_volumes[n]*rAU[0][n],
                wo*geometry.cell_volumes[o]*rAU[1][o]+wn*geometry.cell_volumes[n]*rAU[1][n],
                wo*geometry.cell_volumes[o]*rAU[2][o]+wn*geometry.cell_volumes[n]*rAU[2][n]};
            const double orth=pressure_face_orthogonal_coefficient(Sf,dvec,dAU);
            add_face_derivative(o,density*orth);
            add_face_derivative(n,-density*orth);
            const Vec3 e=dvec*(1.0/dvec.mag());
            const double ge=dAU.x*e.x*e.x+dAU.y*e.y*e.y+dAU.z*e.z*e.z;
            const double ortho_area=Sf.dot(e);
            const Vec3 Q{dAU.x*Sf.x-ge*ortho_area*e.x,
                         dAU.y*Sf.y-ge*ortho_area*e.y,
                         dAU.z*Sf.z-ge*ortho_area*e.z};
            for(const auto& [j,g]:G[o]) add_face_derivative(j,-density*wo*Q.dot(g));
            for(const auto& [j,g]:G[n]) add_face_derivative(j,-density*wn*Q.dot(g));
            continue;
        }
        const auto patch=geometry.face_patch[f];
        if(patch>=mesh.boundary().n_patches()) continue;
        const auto& patchdata=mesh.boundary().patch(patch);
        if(patchdata.type==PatchType::EMPTY || patchdata.type==PatchType::SYMMETRY)
            continue;
        const auto bc=pressure_bcs.find(patchdata.name);
        if(bc==pressure_bcs.end() || bc->second.type!=ScalarBoundaryType::FIXED_VALUE)
            continue;
        const auto dvec=geometry.face_centres[f]-geometry.cell_centres[o];
        const Vec3 dAU{geometry.cell_volumes[o]*rAU[0][o],
                       geometry.cell_volumes[o]*rAU[1][o],
                       geometry.cell_volumes[o]*rAU[2][o]};
        const double orth=pressure_face_orthogonal_coefficient(Sf,dvec,dAU);
        add_face_derivative(o,density*orth);
        const Vec3 e=dvec*(1.0/dvec.mag());
        const double ge=dAU.x*e.x*e.x+dAU.y*e.y*e.y+dAU.z*e.z*e.z;
        const double ortho_area=Sf.dot(e);
        const Vec3 Q{dAU.x*Sf.x-ge*ortho_area*e.x,
                     dAU.y*Sf.y-ge*ortho_area*e.y,
                     dAU.z*Sf.z-ge*ortho_area*e.z};
        for(const auto& [j,g]:G[o]) add_face_derivative(j,-density*Q.dot(g));
    }
    for(std::size_t c=0;c<nc;++c) {
        for(auto it=rows[c].begin();it!=rows[c].end();) {
            if(it->second==0.0) it=rows[c].erase(it);
            else {
                if (!std::isfinite(it->second))
                    throw std::runtime_error("V18 sparse pressure Jacobian non-finite entry");
                ++it;
            }
        }
        if(rows[c].find(c)==rows[c].end())
            throw std::runtime_error("V18 sparse pressure Jacobian missing diagonal");
    }
    return rows;
}

// V18 diagnostic reference Schur operator. This is deliberately O(N*Nfaces)
// and limited to small meshes: for every unit pressure perturbation it uses
// the *production* face-flux reconstruction, including deferred gradients,
// face weights, boundary conditions and directional momentum response.
// It is a numerical oracle for a future scalable analytical Jacobian, NOT
// the default production pressure solver. Fixed pressure/gradient BC values
// are zeroed so the result is the homogeneous linear pressure response.
inline std::vector<std::map<std::size_t, double>>
make_full_rhie_chow_pressure_rows_reference(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::array<std::vector<double>, 3>& rAU, double density,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    bool projected_face_interpolation = false)
{
    using namespace cfdx::core;
    const std::size_t nc=mesh.n_cells();
    if (nc==0 || nc>1024)
        throw std::invalid_argument("V18 full flux pressure oracle requires 1..1024 cells");
    Field<double,Location::CELL> zero_U(nc,"zeroH","m/s",3);
    zero_U.fill(0.0);
    Field<double,Location::CELL> delta_p(nc,"delta_p","Pa",1);
    delta_p.fill(0.0);
    auto zero_velocity_bcs=velocity_bcs;
    for (auto& [name,bc]:zero_velocity_bcs)
        if (bc.type==VelocityBoundaryCondition::Type::FIXED_VALUE)
            bc.value={0.0,0.0,0.0};
    auto homogeneous_p_bcs=pressure_bcs;
    for (auto& [name,bc]:homogeneous_p_bcs) {
        if (bc.type==ScalarBoundaryType::FIXED_VALUE) bc.value=0.0;
        if (bc.type==ScalarBoundaryType::FIXED_GRADIENT) bc.gradient=0.0;
    }
    std::vector<std::map<std::size_t,double>> rows(nc);
    std::vector<double> div(nc,0.0);
    for(std::size_t j=0;j<nc;++j) {
        delta_p(j)=1.0;
        const auto flux=make_rhie_chow_mass_flux(
            mesh,geometry,zero_U,delta_p,rAU,density,
            zero_velocity_bcs,homogeneous_p_bcs,
            RhieChowVelocityContract::HBYA_PREDICTOR,false,
            projected_face_interpolation);
        delta_p(j)=0.0;
        std::fill(div.begin(),div.end(),0.0);
        for(std::size_t f=0;f<mesh.n_faces();++f) {
            const double value=flux(f);
            if (!std::isfinite(value))
                throw std::runtime_error("V18 pressure Jacobian has non-finite face value");
            const auto o=mesh.ownership().owner(f);
            div[o]+=value;
            const auto n=mesh.ownership().neighbour(f);
            if(n>=0) div[static_cast<std::size_t>(n)]-=value;
        }
        for(std::size_t i=0;i<nc;++i)
            if(div[i]!=0.0) rows[i][j]=div[i];
    }
    for(std::size_t i=0;i<nc;++i)
        if (rows[i].find(i)==rows[i].end() || !std::isfinite(rows[i][i]))
            throw std::runtime_error("V18 pressure Jacobian missing diagonal");
    return rows;
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_flux_from_hbya(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& h_by_a,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {},
    bool skew_correct_gradient = false)
{
    return make_rhie_chow_mass_flux(
        mesh, geometry, h_by_a, p, rAU, rho, bcs, pressure_bcs,
        RhieChowVelocityContract::HBYA_PREDICTOR, skew_correct_gradient);
}

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_rhie_chow_flux_from_physical_velocity(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& physical_velocity,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& p,
    const std::array<std::vector<double>, 3>& rAU,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {},
    bool skew_correct_gradient = false)
{
    return make_rhie_chow_mass_flux(
        mesh, geometry, physical_velocity, p, rAU, rho, bcs, pressure_bcs,
        RhieChowVelocityContract::PHYSICAL_VELOCITY, skew_correct_gradient);
}

inline ScalarEquation assemble_momentum_component(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    double effective_dynamic_viscosity,
    const ScalarBoundaryConditions& velocity_bcs,
    std::size_t component,
    bool bounded,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr,
    double high_order_fraction = 1.0)
{
    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        mesh.n_cells(), "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        if (pressure_gradient_component.dimension() != 3 ||
            pressure_gradient_component.size() != mesh.n_cells())
            throw std::invalid_argument("assemble_momentum_component: pressure gradient must be a 3-component cell field");
        source(c) = body_component(c) - pressure_gradient_component.component_data(component)[c];
        sp(c) = 0.0;
    }
    if (component >= 3)
        throw std::invalid_argument("assemble_momentum_component: component index out of range");
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, effective_dynamic_viscosity,
        source, sp, velocity_bcs, bounded, nullptr, nullptr, nullptr, nullptr,
        convection_scheme, convected_field, nullptr, true, nullptr,
        high_order_fraction);
}

// SSTm has zero turbulent shear viscosity at a resolved no-slip wall.  The
// cell-centre nu_t is NOT the face value there. Reusing mu_eff(owner) at the
// Dirichlet wall produces an artificial wall force O(nu_t/nu), and makes the
// momentum equation inconsistent with the molecular-only Cf diagnostic.
// The returned field is used by predictor, acceptance residual and
// independent face-flux audit, so all three see the same wall coefficient.
inline std::vector<double> resolved_wall_face_viscosity(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const std::vector<double>& cell_mu,
    const ScalarBoundaryConditions& velocity_bcs,
    double molecular_wall_mu,
    const std::vector<double>* wall_face_mu_override = nullptr)
{
    if (wall_face_mu_override && wall_face_mu_override->size()!=mesh.n_faces())
        throw std::invalid_argument("V17 wall momentum face viscosity shape mismatch");
    if (!(molecular_wall_mu > 0.0) || !std::isfinite(molecular_wall_mu) ||
        cell_mu.size() != mesh.n_cells() ||
        geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument("SST wall diffusion: invalid molecular viscosity or geometry");
    std::vector<double> face_mu(mesh.n_faces(), 0.0);
    for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
        const auto owner = mesh.ownership().owner(face);
        if (owner >= cell_mu.size() || !std::isfinite(cell_mu[owner]) ||
            cell_mu[owner] < molecular_wall_mu)
            throw std::invalid_argument("SST wall diffusion: invalid effective cell viscosity");
        face_mu[face] = cell_mu[owner];
        if (mesh.ownership().neighbour(face) >= 0) continue;
        const auto patch_id = geometry.face_patch[face];
        if (patch_id >= mesh.boundary().n_patches()) continue;
        const auto& patch = mesh.boundary().patch(patch_id);
        const auto it = velocity_bcs.find(patch.name);
        if (patch.type == cfdx::core::PatchType::WALL &&
            it != velocity_bcs.end() &&
            it->second.type == ScalarBoundaryType::FIXED_VALUE)
            face_mu[face] = wall_face_mu_override
                ? wall_face_mu_override->at(face) : molecular_wall_mu;
        if (patch.type == cfdx::core::PatchType::WALL &&
            it != velocity_bcs.end() &&
            it->second.type == ScalarBoundaryType::FIXED_VALUE &&
            (!std::isfinite(face_mu[face]) || !(face_mu[face]>0.)))
            throw std::invalid_argument("V17 wall momentum viscosity is absent or invalid");
    }
    return face_mu;
}

// Conservative explicit contribution of div(mu_eff * grad(U)^T) to one
// momentum component. The returned vector contains integrated face fluxes;
// internal-face fluxes are oriented outward from the owner cell.
inline std::vector<double> momentum_transpose_diffusion_face_fluxes(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const VelocityBoundaryConditions& velocity_bcs,
    const std::vector<double>& cell_effective_dynamic_viscosity,
    std::size_t momentum_component)
{
    using cfdx::core::Field;
    using cfdx::core::Location;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nf = mesh.n_faces();
    if (momentum_component >= 3 || velocity.dimension() != 3 ||
        velocity.size() != nc || cell_effective_dynamic_viscosity.size() != nc ||
        geometry.cell_centres.size() != nc || geometry.face_centres.size() != nf ||
        geometry.face_area_vectors.size() != nf)
        throw std::invalid_argument(
            "momentum_transpose_diffusion_face_fluxes: field/geometry shape mismatch");
    for (const double mu : cell_effective_dynamic_viscosity) {
        if (!std::isfinite(mu) || mu < 0.0)
            throw std::invalid_argument(
                "momentum_transpose_diffusion_face_fluxes: viscosity must be finite and non-negative");
    }

    std::array<std::vector<std::array<double, 3>>, 3> gradient;
    std::array<std::vector<std::array<double, 3>>, 3> boundary_gradient;
    for (std::size_t component = 0; component < 3; ++component) {
        Field<double, Location::CELL> scalar(
            nc, "velocity_component_for_transpose_stress", "m/s", 1);
        for (std::size_t cell = 0; cell < nc; ++cell)
            scalar(cell) = velocity.component_data(component)[cell];

        ScalarBoundaryConditions scalar_bcs;
        for (const auto& [name, bc] : velocity_bcs) {
            double value = 0.0;
            if (component == 0) value = bc.value.x;
            else if (component == 1) value = bc.value.y;
            else value = bc.value.z;
            scalar_bcs[name] = {
                bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                    ? ScalarBoundaryType::FIXED_VALUE
                    : ScalarBoundaryType::ZERO_GRADIENT,
                value, 0.0};
        }
        auto grad = gauss_gradient_with_boundary(
            scalar, mesh, geometry, scalar_bcs);
        boundary_gradient[component].resize(nc);
        for (std::size_t cell = 0; cell < nc; ++cell)
            for (std::size_t direction = 0; direction < 3; ++direction)
                boundary_gradient[component][cell][direction] =
                    grad.component_data(direction)[cell];
        // The two-point Green-Gauss face value is not linear-exact on skewed
        // polyhedra. Use the rank-aware least-squares gradient on cells with at
        // least three internal neighbours; preserve boundary-aware Green-Gauss
        // values on smaller stencils (e.g. a one-cell-thick channel), where
        // prescribed boundary values complete the stencil.
        const auto grad_ls = cfdx::core::compute_gradient_least_squares(scalar, mesh, cfdx::core::BoundaryGradientPolicy::EXCLUDE_BOUNDARY);
        const auto* cell_faces = mesh.cells().faces_data();
        const auto* cell_offsets = mesh.cells().offsets_data();
        const auto& ownership = mesh.ownership();
        for (std::size_t cell = 0; cell < nc; ++cell) {
            std::size_t internal_neighbours = 0;
            for (cfdx::core::Offset k = cell_offsets[cell];
                 k < cell_offsets[cell + 1]; ++k) {
                const std::size_t face = cell_faces[k];
                if (ownership.neighbour(face) >= 0 || ownership.owner(face) != cell)
                    ++internal_neighbours;
            }
            if (internal_neighbours >= 3) {
                for (std::size_t direction = 0; direction < 3; ++direction)
                    grad.component_data(direction)[cell] =
                        grad_ls.component_data(direction)[cell];
            }
        }
        gradient[component].resize(nc);
        for (std::size_t cell = 0; cell < nc; ++cell) {
            for (std::size_t direction = 0; direction < 3; ++direction) {
                const double value = grad.component_data(direction)[cell];
                if (!std::isfinite(value))
                    throw std::runtime_error(
                        "momentum_transpose_diffusion_face_fluxes: non-finite velocity gradient");
                gradient[component][cell][direction] = value;
            }
        }
    }

    // Reconstruct the cellwise viscosity gradient once per flux field.
    // Normal-distance weighting alone is not linear-exact on a skewed face:
    // the face centroid is generally not on the owner-neighbour centre line.
    Field<double, Location::CELL> mu_cell(
        nc, "transpose_stress_mu", "Pa s", 1);
    for (std::size_t cell = 0; cell < nc; ++cell)
        mu_cell(cell) = cell_effective_dynamic_viscosity[cell];
    Field<double, Location::CELL> mu_gradient(
        nc, "transpose_stress_mu_gradient", "Pa", 3);
    mu_gradient.fill(0.0);
    const bool uniform_face_viscosity =
        cell_effective_dynamic_viscosity.size() <= 1 ||
        std::all_of(
            cell_effective_dynamic_viscosity.begin() + 1,
            cell_effective_dynamic_viscosity.end(),
            [&](double value) {
                return value == cell_effective_dynamic_viscosity.front();
            });
    if (!uniform_face_viscosity) {
        mu_gradient = cfdx::core::compute_gradient_least_squares(
            mu_cell, mesh,
            cfdx::core::BoundaryGradientPolicy::EXCLUDE_BOUNDARY);
    }

    std::vector<double> face_fluxes(nf, 0.0);
    for (std::size_t face = 0; face < nf; ++face) {
        const std::size_t owner = mesh.ownership().owner(face);
        if (owner >= nc)
            throw std::runtime_error(
                "momentum_transpose_diffusion_face_fluxes: invalid face owner");
        const auto neighbour_index = mesh.ownership().neighbour(face);
        const auto& sf = geometry.face_area_vectors[face];
        const double area = sf.mag();
        if (!(area > 0.0) || !std::isfinite(area))
            throw std::runtime_error(
                "momentum_transpose_diffusion_face_fluxes: invalid face area");
        double mu_face = cell_effective_dynamic_viscosity[owner];
        std::array<double, 3> grad_face{};
        if (neighbour_index >= 0) {
            const std::size_t neighbour = static_cast<std::size_t>(neighbour_index);
            if (neighbour >= nc || neighbour == owner)
                throw std::runtime_error(
                    "momentum_transpose_diffusion_face_fluxes: invalid face neighbour");
            const auto d_owner_vec = geometry.face_centres[face] -
                                     geometry.cell_centres[owner];
            const auto d_neighbour_vec = geometry.cell_centres[neighbour] -
                                         geometry.face_centres[face];
            const double d_owner = std::abs(d_owner_vec.dot(sf)) / area;
            const double d_neighbour = std::abs(d_neighbour_vec.dot(sf)) / area;
            const double d_sum = d_owner + d_neighbour;
            if (!(d_owner > 0.0) || !(d_neighbour > 0.0) ||
                !(d_sum > 0.0) || !std::isfinite(d_sum))
                throw std::runtime_error(
                    "momentum_transpose_diffusion_face_fluxes: invalid face interpolation distance");
            const double w_owner = d_neighbour / d_sum;
            const double w_neighbour = d_owner / d_sum;
            const auto& face_centre = geometry.face_centres[face];
            const auto& owner_centre = geometry.cell_centres[owner];
            const auto& neighbour_centre = geometry.cell_centres[neighbour];
            const double mu_owner_at_face =
                cell_effective_dynamic_viscosity[owner] +
                mu_gradient.component_data(0)[owner] * (face_centre.x - owner_centre.x) +
                mu_gradient.component_data(1)[owner] * (face_centre.y - owner_centre.y) +
                mu_gradient.component_data(2)[owner] * (face_centre.z - owner_centre.z);
            const double mu_neighbour_at_face =
                cell_effective_dynamic_viscosity[neighbour] +
                mu_gradient.component_data(0)[neighbour] * (face_centre.x - neighbour_centre.x) +
                mu_gradient.component_data(1)[neighbour] * (face_centre.y - neighbour_centre.y) +
                mu_gradient.component_data(2)[neighbour] * (face_centre.z - neighbour_centre.z);
            mu_face = w_owner * mu_owner_at_face +
                      w_neighbour * mu_neighbour_at_face;
            if (!std::isfinite(mu_face) || mu_face < 0.0)
                throw std::runtime_error(
                    "momentum_transpose_diffusion_face_fluxes: reconstructed face viscosity is invalid");
            for (std::size_t velocity_component = 0; velocity_component < 3;
                 ++velocity_component) {
                grad_face[velocity_component] =
                    w_owner * gradient[velocity_component][owner][momentum_component] +
                    w_neighbour * gradient[velocity_component][neighbour][momentum_component];
            }
            const double flux = mu_face * (
                grad_face[0] * sf.x + grad_face[1] * sf.y +
                grad_face[2] * sf.z);
            if (!std::isfinite(flux))
                throw std::runtime_error(
                    "momentum_transpose_diffusion_face_fluxes: non-finite internal stress flux");
            face_fluxes[face] = flux;
        } else {
            // Preserve prescribed boundary-value information in the boundary
            // stress flux; least-squares reconstruction is used for internal
            // faces only.
            for (std::size_t velocity_component = 0; velocity_component < 3;
                 ++velocity_component)
                grad_face[velocity_component] =
                    boundary_gradient[velocity_component][owner][momentum_component];
            const double flux = mu_face * (
                grad_face[0] * sf.x + grad_face[1] * sf.y +
                grad_face[2] * sf.z);
            if (!std::isfinite(flux))
                throw std::runtime_error(
                    "momentum_transpose_diffusion_face_fluxes: non-finite boundary stress flux");
            face_fluxes[face] = flux;
        }
    }
    return face_fluxes;
}

inline std::vector<double> momentum_transpose_diffusion_rhs(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    const VelocityBoundaryConditions& velocity_bcs,
    const std::vector<double>& cell_effective_dynamic_viscosity,
    std::size_t momentum_component)
{
    const auto face_fluxes = momentum_transpose_diffusion_face_fluxes(
        mesh, geometry, velocity, velocity_bcs,
        cell_effective_dynamic_viscosity, momentum_component);
    std::vector<double> rhs(mesh.n_cells(), 0.0);
    for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
        const std::size_t owner = mesh.ownership().owner(face);
        rhs[owner] += face_fluxes[face];
        const auto neighbour = mesh.ownership().neighbour(face);
        if (neighbour >= 0)
            rhs[static_cast<std::size_t>(neighbour)] -= face_fluxes[face];
    }
    return rhs;
}

inline ScalarEquation assemble_momentum_component(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    const std::vector<double>& cell_effective_dynamic_viscosity,
    const ScalarBoundaryConditions& velocity_bcs,
    std::size_t component,
    bool bounded,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr,
    const std::vector<double>* extra_rhs = nullptr,
    double molecular_wall_mu = -1.0,
    double high_order_fraction = 1.0,
    const ScalarBoundaryFaceConditions* face_conditions = nullptr,
    const std::vector<double>* wall_face_mu_override = nullptr)
{
    const std::size_t nc = mesh.n_cells();
    if (extra_rhs != nullptr && extra_rhs->size() != nc)
        throw std::invalid_argument(
            "assemble_momentum_component: extra RHS size does not match mesh");
    if (cell_effective_dynamic_viscosity.size() != nc)
        throw std::invalid_argument(
            "assemble_momentum_component: cell viscosity size does not match mesh");
    for (const double viscosity : cell_effective_dynamic_viscosity) {
        if (!std::isfinite(viscosity) || viscosity < 0.0)
            throw std::invalid_argument(
                "assemble_momentum_component: cell viscosity must be finite and non-negative");
    }
    if (component >= 3)
        throw std::invalid_argument("assemble_momentum_component: component index out of range");
    if (pressure_gradient_component.dimension() != 3 ||
        pressure_gradient_component.size() != nc ||
        body_component.size() != nc || body_component.dimension() != 1)
        throw std::invalid_argument(
            "assemble_momentum_component: source/pressure-gradient field shape mismatch");

    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        nc, "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        nc, "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < nc; ++c) {
        source(c) = body_component(c) -
                    pressure_gradient_component.component_data(component)[c];
        sp(c) = 0.0;
    }
    const auto face_mu = molecular_wall_mu > 0.0
        ? resolved_wall_face_viscosity(
            mesh, geometry, cell_effective_dynamic_viscosity,
            velocity_bcs, molecular_wall_mu, wall_face_mu_override)
        : std::vector<double>{};
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, 0.0, source, sp, velocity_bcs, bounded,
        nullptr, nullptr, extra_rhs, &cell_effective_dynamic_viscosity,
        convection_scheme, convected_field, face_conditions, true,
        face_mu.empty() ? nullptr : &face_mu, high_order_fraction);
}

// Production wrapper: the full velocity tensor and vector boundary conditions
// are required so the transpose-stress term cannot be silently omitted.
inline ScalarEquation assemble_momentum_component_with_transpose_diffusion(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    const std::vector<double>& cell_effective_dynamic_viscosity,
    const ScalarBoundaryConditions& component_boundary_conditions,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& full_velocity,
    const VelocityBoundaryConditions& velocity_bcs,
    std::size_t component,
    bool bounded,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr,
    double high_order_fraction = 1.0)
{
    // In the constant-viscosity incompressible limit, the transpose
    // contribution is equivalent to zero when continuity is satisfied.
    // Avoid feeding the finite discrete continuity defect back as an
    // explicit source into the established constant-viscosity solver.
    const bool uniform_viscosity =
        cell_effective_dynamic_viscosity.size() <= 1 ||
        std::all_of(
            cell_effective_dynamic_viscosity.begin() + 1,
            cell_effective_dynamic_viscosity.end(),
            [&](double mu) { return mu == cell_effective_dynamic_viscosity.front(); });
    const auto transpose_rhs = uniform_viscosity
        ? std::vector<double>(mesh.n_cells(), 0.0)
        : momentum_transpose_diffusion_rhs(
              mesh, geometry, full_velocity, velocity_bcs,
              cell_effective_dynamic_viscosity, component);
    return assemble_momentum_component(
        mesh, geometry, mass_flux, pressure_gradient_component, body_component,
        cell_effective_dynamic_viscosity, component_boundary_conditions, component,
        bounded, convection_scheme, convected_field, &transpose_rhs, -1.0,
        high_order_fraction);
}

inline ScalarEquation assemble_momentum_component(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_gradient_component,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& body_component,
    double effective_dynamic_viscosity,
    const BoundaryConstraintMap& boundary_constraints,
    std::size_t component,
    bool bounded,
    double time = 0.0,
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr,
    double high_order_fraction = 1.0)
{
    if (component >= 3)
        throw std::invalid_argument("assemble_momentum_component: component index out of range");
    if (pressure_gradient_component.dimension() != 3 ||
        pressure_gradient_component.size() != mesh.n_cells())
        throw std::invalid_argument("assemble_momentum_component: pressure gradient must be a 3-component cell field");
    cfdx::core::Field<double, cfdx::core::Location::CELL> source(
        mesh.n_cells(), "momentum_source", "N/m3", 1);
    cfdx::core::Field<double, cfdx::core::Location::CELL> sp(
        mesh.n_cells(), "momentum_sp", "kg/m3/s", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        source(c) = body_component(c) - pressure_gradient_component.component_data(component)[c];
        sp(c) = 0.0;
    }
    const std::string field = "U." + std::string(1, static_cast<char>('x' + component));
    const auto face_conditions = resolve_scalar_boundary_constraints(
        mesh, geometry.face_centres, boundary_constraints, field, time, &mass_flux);
    return assemble_scalar_equation(
        mesh, geometry, mass_flux, effective_dynamic_viscosity,
        source, sp, {}, bounded, nullptr, nullptr, nullptr, nullptr,
        convection_scheme, convected_field, &face_conditions, true, nullptr,
        high_order_fraction);
}

inline void add_dual_time_mass_term(
    ScalarEquation& equation,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& u_n,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& u_nm1,
    const FvGeometry& geometry,
    double rho,
    double dt,
    double dt_previous,
    bool bdf2,
    std::size_t component)
{
    if (!(dt > 0.0) || !std::isfinite(dt) || !(rho > 0.0) || !std::isfinite(rho))
        throw std::invalid_argument("add_dual_time_mass_term: invalid rho/dt");
    if (u_n.dimension() != 3 || u_nm1.dimension() != 3 ||
        u_n.size() != equation.rhs.size() || u_nm1.size() != equation.rhs.size())
        throw std::invalid_argument("add_dual_time_mass_term: history shape mismatch");
    if (component >= 3) throw std::invalid_argument("add_dual_time_mass_term: invalid component");

    double a0 = 1.0, a1 = 0.0, denom = 1.0, b = dt;
    if (bdf2) {
        if (!(dt_previous > 0.0) || !std::isfinite(dt_previous))
            throw std::invalid_argument("add_dual_time_mass_term: invalid previous dt for BDF2");
        cfdx::core::bdf2_coefficients(dt_previous, dt, a0, a1, denom, b);
    }
    for (std::size_t c = 0; c < equation.rhs.size(); ++c) {
        const double mass = rho * geometry.cell_volumes[c] * denom / b;
        const double source = rho * geometry.cell_volumes[c] *
            (a0 * u_n.component_data(component)[c] -
             a1 * u_nm1.component_data(component)[c]) / b;
        if (!(mass > 0.0) || !std::isfinite(mass) || !std::isfinite(source))
            throw std::runtime_error("add_dual_time_mass_term: non-finite temporal coefficient");
        auto* rows = equation.matrix.row_offsets_data();
        auto* cols = equation.matrix.columns_data();
        auto* values = equation.matrix.values_data();
        bool found = false;
        for (std::uint32_t k = rows[c]; k < rows[c + 1]; ++k) {
            if (cols[k] == c) {
                values[k] += mass;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("add_dual_time_mass_term: missing momentum diagonal");
        equation.diagonal[c] += mass;
        equation.rhs(c) += source;
    }
}


// V17: conservative face accounting uses ONLY the published RC flux. All
// internal owner/neighbour contributions cancel by topology, independently of
// a cell-velocity interpolation or a pressure-matrix residual.
struct V17MassPatchAudit {
    struct Patch { std::size_t faces=0; double signed_flux=0.0; double abs_flux=0.0; };
    std::map<std::string,Patch> patches;
    std::vector<double> cell_divergence;
    double boundary_signed=0.0, boundary_abs=0.0, cell_signed=0.0;
    double cell_l1=0.0, cell_linf=0.0, topology_error=0.0;
    std::size_t worst_cell=0;
};

inline V17MassPatchAudit audit_v17_authoritative_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& flux)
{
    if (mesh.n_cells()==0 || flux.size() != mesh.n_faces() ||
        geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument("V17 flux audit: invalid face arrays");
    V17MassPatchAudit out;
    out.cell_divergence.assign(mesh.n_cells(), 0.0);
    std::vector<long double> cell_sum(mesh.n_cells(), 0.0L);
    long double boundary_signed=0.0L, cell_signed=0.0L;
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const double value=flux(f);
        if (!std::isfinite(value))
            throw std::invalid_argument("V17 flux audit: nonfinite Rhie-Chow face mass flux");
        const std::size_t owner=mesh.ownership().owner(f);
        if (owner>=mesh.n_cells())
            throw std::invalid_argument("V17 flux audit: invalid owner");
        cell_sum[owner]+=static_cast<long double>(value);
        const auto neighbor=mesh.ownership().neighbour(f);
        if (neighbor>=0) {
            if (static_cast<std::size_t>(neighbor)>=mesh.n_cells())
                throw std::invalid_argument("V17 flux audit: invalid neighbour");
            cell_sum[static_cast<std::size_t>(neighbor)]-=static_cast<long double>(value);
        } else {
            const std::size_t pi=geometry.face_patch[f];
            if (pi>=mesh.boundary().n_patches())
                throw std::invalid_argument("V17 flux audit: unowned boundary patch");
            auto &patch=out.patches[mesh.boundary().patch(pi).name];
            ++patch.faces;
            patch.signed_flux+=value;
            patch.abs_flux+=std::abs(value);
            boundary_signed+=static_cast<long double>(value);
            out.boundary_abs+=std::abs(value);
        }
    }
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        const double val=static_cast<double>(cell_sum[c]);
        out.cell_divergence[c]=val;
        cell_signed+=cell_sum[c];
        out.cell_l1+=std::abs(val);
        if (std::abs(val)>out.cell_linf) {
            out.cell_linf=std::abs(val); out.worst_cell=c;
        }
    }
    out.cell_signed=static_cast<double>(cell_signed);
    out.boundary_signed=static_cast<double>(boundary_signed);
    out.topology_error=static_cast<double>(cell_signed-boundary_signed);
    // This is a topology closure, not a physical convergence certificate.
    if (std::abs(out.topology_error)>1e-11*std::max(1.0,out.boundary_abs))
        throw std::runtime_error("V17 flux audit: topology closure failed");
    return out;
}

// Local pseudo-time coefficient M_P = rho V / dtau, with a spectral scale
// based on the SAME authoritative mass flux and positive viscosity consumed
// by the next momentum predictor. The 2*diffusion contribution is a stability
// heuristic, not a new physical source. No face-by-face flux is modified.
inline std::vector<double> make_v17_pseudo_momentum_masses(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& flux,
    const std::vector<double>& dynamic_viscosity,
    double density, double cfl, double molecular_wall_viscosity=-1.0,
    const std::vector<double>* actual_wall_face_mu=nullptr)
{
    if (!std::isfinite(molecular_wall_viscosity) ||
        (molecular_wall_viscosity < 0.0 && molecular_wall_viscosity != -1.0))
        throw std::invalid_argument("V17 pseudo-momentum: invalid wall viscosity");
    if (!(density>0.0) || !std::isfinite(density) || !(cfl>0.0) ||
        !std::isfinite(cfl) || flux.size()!=mesh.n_faces() ||
        dynamic_viscosity.size()!=mesh.n_cells())
        throw std::invalid_argument("V17 pseudo-momentum: invalid dimensions or properties");
    // V31: the production V17 active wall operator supplies a per-face
    // viscosity. Use that SAME viscosity in the pseudo-time spectral scale.
    // Negative entries mean the face is not a V17 wall, hence fall back to
    // the existing molecular-wall policy. This does NOT change the steady FV
    // momentum equation, its Rhie-Chow flux, or its physical certificate.
    if (actual_wall_face_mu) {
        if (actual_wall_face_mu->size()!=mesh.n_faces())
            throw std::invalid_argument("V31 pseudo-momentum: wall-face viscosity size");
        for (const auto v:*actual_wall_face_mu)
            if (!std::isfinite(v) || (v<0. && v!=-1.) || v==0.)
                throw std::invalid_argument("V31 pseudo-momentum: invalid wall-face viscosity");
    }
    std::vector<double> spectral(mesh.n_cells(),0.0);
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        if (!(geometry.cell_volumes[c]>0.0) || !std::isfinite(geometry.cell_volumes[c]) ||
            !(dynamic_viscosity[c]>=0.0) || !std::isfinite(dynamic_viscosity[c]))
            throw std::invalid_argument("V17 pseudo-momentum: invalid cell material or volume");
    }
    for (std::size_t f=0;f<mesh.n_faces();++f) {
        const double phi=flux(f);
        if (!std::isfinite(phi))
            throw std::invalid_argument("V17 pseudo-momentum: nonfinite flux");
        const auto owner=mesh.ownership().owner(f);
        const auto nei=mesh.ownership().neighbour(f);
        cfdx::core::PatchType boundary_type=cfdx::core::PatchType::WALL;
        if (nei<0) {
            const std::size_t patch=geometry.face_patch[f];
            if (patch>=mesh.boundary().n_patches())
                throw std::invalid_argument("V17 pseudo-momentum: unassigned boundary face");
            boundary_type=mesh.boundary().patch(patch).type;
            // A 2-D extruded mesh has zero transport through EMPTY faces:
            // treating its tiny extrusion thickness as a viscous distance
            // added a huge artificial diagonal and blocked the steady solve.
            // Symmetry has no normal flux or normal diffusive coupling either.
            if (boundary_type==cfdx::core::PatchType::EMPTY ||
                boundary_type==cfdx::core::PatchType::SYMMETRY)
                continue;
        }
        const double d=nei>=0
            ? (geometry.cell_centres[static_cast<std::size_t>(nei)]-
               geometry.cell_centres[owner]).mag()
            : (geometry.face_centres[f]-geometry.cell_centres[owner]).mag();
        if (!(d>0.0) || !std::isfinite(d))
            throw std::invalid_argument("V17 pseudo-momentum: invalid centre-to-face geometry");
        const double mu=nei>=0
            ? 0.5*(dynamic_viscosity[owner]+dynamic_viscosity[static_cast<std::size_t>(nei)])
            : ((boundary_type==cfdx::core::PatchType::WALL &&
                actual_wall_face_mu && (*actual_wall_face_mu)[f]>=0.)
                ? (*actual_wall_face_mu)[f]
                : ((boundary_type==cfdx::core::PatchType::WALL &&
                    molecular_wall_viscosity>=0.0)
                    ? molecular_wall_viscosity : dynamic_viscosity[owner]));
        const double spectral_face=std::abs(phi)+2.0*mu*geometry.face_area_vectors[f].mag()/d;
        if (!std::isfinite(spectral_face))
            throw std::invalid_argument("V17 pseudo-momentum: overflow in local spectral mass");
        spectral[owner]+=spectral_face;
        if (nei>=0) spectral[static_cast<std::size_t>(nei)]+=spectral_face;
    }
    for (std::size_t c=0;c<mesh.n_cells();++c) {
        spectral[c]=std::max(spectral[c],1e-12*density*geometry.cell_volumes[c])/cfl;
        if (!(spectral[c]>0.0) || !std::isfinite(spectral[c]))
            throw std::runtime_error("V17 pseudo-momentum: nonfinite artificial mass");
    }
    return spectral;
}

// Optional cap for the pseudo-time diagonal. It is chosen from the three
// positive *pre-pseudo* momentum diagonals and applied identically to Ux,
// Uy, Uz, preserving isotropic artificial mass and the exact physical fixed
// point M_P*(U_new-U_old)=0. It never modifies the actual physical residual
// or any convergence threshold. A finite cap must be opt-in.
inline std::size_t cap_v17_pseudo_momentum_masses(
    std::vector<double>& masses,
    const ScalarEquation& ex, const ScalarEquation& ey,
    const ScalarEquation& ez, double maximum_fraction)
{
    if (!(maximum_fraction > 0.0) || std::isnan(maximum_fraction) ||
        masses.size() != ex.diagonal.size() ||
        masses.size() != ey.diagonal.size() ||
        masses.size() != ez.diagonal.size())
        throw std::invalid_argument("V17 pseudo-momentum: invalid diagonal cap");
    if (std::isinf(maximum_fraction)) return 0;
    std::size_t capped = 0;
    for (std::size_t i=0; i<masses.size(); ++i) {
        const double physical_diagonal = std::min({
            ex.diagonal[i], ey.diagonal[i], ez.diagonal[i]});
        const double old_mass = masses[i];
        if (!(physical_diagonal > 0.0) || !std::isfinite(physical_diagonal) ||
            !(old_mass > 0.0) || !std::isfinite(old_mass))
            throw std::invalid_argument("V17 pseudo-momentum: invalid physical diagonal or mass");
        const double maximum_mass = maximum_fraction * physical_diagonal;
        if (!(maximum_mass > 0.0) || std::isnan(maximum_mass))
            throw std::invalid_argument("V17 pseudo-momentum: invalid cap mass");
        if (old_mass > maximum_mass) {
            masses[i] = maximum_mass;
            ++capped;
        }
    }
    return capped;
}

inline void add_v17_pseudo_momentum_term(
    ScalarEquation& eq,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& old_velocity,
    const std::vector<double>& masses, std::size_t component)
{
    if (old_velocity.dimension()!=3 || component>=3 ||
        old_velocity.size()!=eq.rhs.size() || masses.size()!=eq.rhs.size())
        throw std::invalid_argument("V17 pseudo-momentum: incompatible equation");
    for (std::size_t c=0;c<masses.size();++c) {
        const double m=masses[c];
        const double source=m*old_velocity.component_data(component)[c];
        if (!(m>0.0) || !std::isfinite(m) || !std::isfinite(source))
            throw std::invalid_argument("V17 pseudo-momentum: invalid mass/source");
        bool found=false;
        for (std::uint32_t k=eq.matrix.row_offsets_data()[c];
             k<eq.matrix.row_offsets_data()[c+1];++k) {
            if (eq.matrix.columns_data()[k]==c) {
                eq.matrix.values_data()[k]+=m;
                found=true; break;
            }
        }
        if (!found) throw std::runtime_error("V17 pseudo-momentum: missing diagonal");
        eq.diagonal[c]+=m;
        eq.rhs(c)+=source;
    }
}

inline double next_v17_pseudo_cfl(
    double current, double previous_momentum, double current_momentum,
    double previous_continuity, double current_continuity,
    const V17PseudoMomentumControls& settings)
{
    validate_v17_pseudo_momentum_controls(settings);
    for (double v : {current,previous_momentum,current_momentum,
                     previous_continuity,current_continuity})
        if (!std::isfinite(v) || v<0.0)
            throw std::invalid_argument("V17 pseudo-momentum: nonfinite adaptation metric");
    double updated=current;
    if (current_momentum<=0.95*previous_momentum &&
        current_continuity<=0.95*previous_continuity)
        updated*=settings.growth;
    else if (current_momentum>1.15*previous_momentum ||
             current_continuity>1.15*previous_continuity)
        updated*=settings.shrink;
    if (!std::isfinite(updated)) updated=settings.maximum_cfl;
    return std::clamp(updated,settings.minimum_cfl,settings.maximum_cfl);
}

inline void apply_velocity_boundary_conditions(
    const cfdx::core::Mesh& mesh,
    const VelocityBoundaryConditions& bcs,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U)
{
    // Cell-centred finite-volume fields have no independent patch storage.
    // Fixed-value conditions are therefore enforced through the matrix and
    // the face flux; no cell overwrite is performed here.
    (void)mesh;
    (void)bcs;
    (void)U;
}

inline double sample_incompressible_probe(
    const IncompressiblePointProbe& probe,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p)
{
    if (mesh.n_cells() == 0)
        throw std::invalid_argument("sample_incompressible_probe: mesh has no cells");
    std::size_t nearest = 0;
    double best_distance = std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto delta = geometry.cell_centres[c] - probe.location;
        const double distance = delta.mag2();
        if (distance < best_distance) {
            best_distance = distance;
            nearest = c;
        }
    }
    switch (probe.field) {
    case IncompressibleProbeField::PRESSURE: return p(nearest);
    case IncompressibleProbeField::U_X: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.x; }
    case IncompressibleProbeField::U_Y: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.y; }
    case IncompressibleProbeField::U_Z: { cfdx::core::Vec3 v; U.get(nearest, v.x, v.y, v.z); return v.z; }
    case IncompressibleProbeField::U_MAGNITUDE: {
        cfdx::core::Vec3 velocity;
        U.get(nearest, velocity.x, velocity.y, velocity.z);
        return velocity.mag();
    }
    }
    throw std::invalid_argument("sample_incompressible_probe: unsupported field");
}

inline void relax_momentum_equation(
    ScalarEquation& equation,
    const cfdx::core::Vector& old_field,
    double alpha)
{
    if (!(alpha > 0.0 && alpha <= 1.0) ||
        old_field.size() != equation.rhs.size())
        throw std::invalid_argument("relax_momentum_equation: invalid relaxation");

    if (alpha == 1.0) return;

    const auto* rows = equation.matrix.row_offsets_data();
    const auto* cols = equation.matrix.columns_data();
    auto* values = equation.matrix.values_data();

    // Standard equation relaxation:
    //
    //   A u = b
    //
    // becomes
    //
    //   [A_off + A_P/alpha] u
    //       = b + (1-alpha)/alpha A_P u_old.
    //
    // This preserves the exact algebraic relation used later to construct
    // HbyA and rAU. A field blend after the solve would break that relation.
    for (std::size_t c = 0; c < old_field.size(); ++c) {
        const double a = equation.diagonal[c];
        if (!(a > 0.0) || !std::isfinite(a))
            throw std::runtime_error("relax_momentum_equation: invalid diagonal");

        double* diagonal_value = nullptr;
        for (std::uint32_t k = rows[c]; k < rows[c + 1]; ++k) {
            if (cols[k] == c) {
                diagonal_value = &values[k];
                break;
            }
        }
        if (!diagonal_value)
            throw std::runtime_error("relax_momentum_equation: missing diagonal entry");

        const double added = a * (1.0 / alpha - 1.0);
        *diagonal_value += added;
        equation.diagonal[c] += added;
        equation.rhs(c) += added * old_field(c);
    }
}

struct PcdPressureOperators {
    cfdx::core::SparseMatrix mass;
    cfdx::core::SparseMatrix laplacian;
    cfdx::core::SparseMatrix convection_diffusion;
    std::unique_ptr<cfdx::core::Preconditioner> laplacian_preconditioner;
    std::unique_ptr<cfdx::core::Preconditioner> convection_diffusion_preconditioner;
    std::unique_ptr<PcdFpEquilibration> convection_diffusion_equilibrated;
    std::unique_ptr<PcdFpStrongComponentAnchor> convection_diffusion_anchored;
    std::unique_ptr<PcdFpBiNullspaceV13> convection_diffusion_bisided;
};

inline cfdx::core::SparseMatrix impose_pcd_reference_row(
    const cfdx::core::SparseMatrix& input, std::size_t reference_cell)
{
    const std::size_t n = input.n_rows();
    if (n == 0 || reference_cell >= n)
        throw std::invalid_argument("PCD reference cell is outside pressure operator");
    cfdx::core::SparseMatrix result(n, n);
    const auto* ro = input.row_offsets_data();
    const auto* co = input.columns_data();
    const auto* va = input.values_data();
    for (std::size_t r = 0; r < n; ++r) {
        if (r == reference_cell) {
            result.push_back(r, r, 1.0);
            continue;
        }
        for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k)
            result.push_back(r, co[k], va[k]);
    }
    result.finalize();
    return result;
}

// PCD pressure-side operators. The convection_scheme parameter is accepted
// but deliberately NOT forwarded to assemble_scalar_equation, and this is
// not a TODO. Every ConvectionScheme except UPWIND is implemented in
// assemble_scalar_equation as a DEFERRED CORRECTION: an explicit
// right-hand-side term evaluated on the convected field, riding the SAME
// conservative first-order upwind matrix. These operators are matrix-only
// preconditioner ingredients, so the assembled matrix is identical for
// every scheme and ConvectionScheme::UPWIND below is not an approximation
// of the true pressure operator but its exact matrix for every scheme
// setting. Forwarding a non-UPWIND scheme would also throw, because a
// preconditioner assembly has no convected field to evaluate the deferred
// correction against. Audit context: issue #663 (NACA0012 convection-scheme
// sensitivity campaign).
inline PcdPressureOperators assemble_pcd_pressure_operators(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
    double rho,
    double kinematic_viscosity,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    std::size_t reference_cell,
    ConvectionScheme /*convection_scheme — unused by design, see comment above*/,
    bool pressure_space_inflow_dirichlet = false)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    if (velocity.size() != nc || velocity.dimension() != 3)
        throw std::invalid_argument("PCD pressure operator assembly: invalid velocity field");
    if (!(rho > 0.0) || !std::isfinite(rho) ||
        !(kinematic_viscosity >= 0.0) || !std::isfinite(kinematic_viscosity))
        throw std::invalid_argument("PCD pressure operator assembly: invalid fluid properties");

    PcdPressureOperators result;
    result.mass = SparseMatrix(nc, nc);
    for (std::size_t c = 0; c < nc; ++c)
        result.mass.push_back(c, c, geometry.cell_volumes[c]);
    result.mass.finalize();

    Field<double, Location::FACE> zero_flux(mesh.n_faces(), "pcd_zero_flux", "", 1);
    auto mass_flux = make_mass_flux(mesh, geometry, velocity, rho, velocity_bcs);

    // PCD pressure-space inflow boundaries are NOT necessarily the physical
    // pressure BCs. The same farfield patch can carry both inflow/outflow:
    // resolve the actual signed flux on each boundary face. For an incoming
    // face use homogeneous Dirichlet for the auxiliary Kp and Fp operators;
    // characteristic/symmetry/solid/empty faces retain their original BCs.
    // This is an opt-in Schur-preconditioner experiment only. Never rewrite
    // the coupled FV matrix or pressure solution boundary conditions here.
    ScalarBoundaryFaceConditions inflow_conditions(mesh.n_faces());
    if (pressure_space_inflow_dirichlet) {
        for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
            if (mesh.ownership().neighbour(face) >= 0 ||
                !(mass_flux(face) < 0.0)) continue;
            const auto patch_index = geometry.face_patch[face];
            if (patch_index >= mesh.boundary().n_patches()) continue;
            const auto kind = mesh.boundary().patch(patch_index).type;
            if (kind == PatchType::WALL || kind == PatchType::EMPTY ||
                kind == PatchType::SYMMETRY || kind == PatchType::PERIODIC ||

                kind == PatchType::INTERFACE) continue;
            inflow_conditions.specified[face] = 1;
            inflow_conditions.conditions[face] = {
                ScalarBoundaryFaceCondition::Type::FIXED_VALUE, 0.0};
        }
    }
    const auto* auxiliary_face_bcs = pressure_space_inflow_dirichlet
        ? &inflow_conditions : nullptr;

    Field<double, Location::CELL> zero_source(nc, "pcd_zero_source", "", 1);
    Field<double, Location::CELL> zero_implicit(nc, "pcd_zero_implicit", "", 1);

    const auto diffusion = assemble_scalar_equation(
        mesh, geometry, zero_flux, 1.0,
        zero_source, zero_implicit, pressure_bcs, true, nullptr, nullptr,
        nullptr, nullptr, ConvectionScheme::UPWIND, nullptr,
        auxiliary_face_bcs, false);

    const auto convection_diffusion = assemble_scalar_equation(
        mesh, geometry, mass_flux, rho * kinematic_viscosity,
        zero_source, zero_implicit, pressure_bcs, true, nullptr, nullptr,
        nullptr, nullptr, ConvectionScheme::UPWIND, nullptr,
        auxiliary_face_bcs, false);

    // A fixed pressure boundary already removes the constant null mode of
    // the pressure operators, so the gauge reference row is imposed only
    // when no such boundary exists. The fixed-pressure decision is made
    // here, from pressure_bcs, so every caller passes its true reference
    // cell and the caller-side gauge policy cannot drift from the
    // assembly. A gauge sentinel (>= nc) still leaves the operators
    // untouched for direct callers that pass one.
    const bool fixed_pressure_gauge =
        has_effective_fixed_pressure_boundary(mesh, pressure_bcs);
    const std::size_t gauge_cell = fixed_pressure_gauge
        ? std::numeric_limits<std::size_t>::max() : reference_cell;
    if (gauge_cell < nc) {
        result.laplacian = impose_pcd_reference_row(diffusion.matrix, gauge_cell);
        result.convection_diffusion =
            impose_pcd_reference_row(convection_diffusion.matrix, gauge_cell);
    } else {
        result.laplacian = diffusion.matrix;
        result.convection_diffusion = convection_diffusion.matrix;
    }
    return result;
}

// Coherent gauge for the nested algebraic-Schur pressure operators
// (the LSC/BFBt stabilized operator). The coupled production system
// imposes the pressure reference by eliminating the reference row AND
// column and inserting a unit diagonal (see the reduced coupled
// assembly), so the nested operator mirrors exactly that structure.
//
// The reference right-hand side entry is deliberately NOT zeroed and
// the nested correction at the reference is deliberately NOT pinned to
// zero: the outer Krylov iteration needs the gauge component of the
// preconditioned residual to correct the coupled gauge row (see the
// PR #832 note on the initial gauge row, and the measured
// frozen-residual failure documented in pcd_schur.h). A fixed pressure
// boundary already removes the constant null mode, so the caller
// passes a sentinel (>= n) and this function is a no-op in that case.
inline void apply_nested_schur_pressure_gauge(
    cfdx::core::SparseMatrix& pressure_operator, std::size_t gauge_cell)
{
    const std::size_t n = pressure_operator.n_rows();
    if (n == 0 || gauge_cell >= n)
        return;
    for (std::size_t k = pressure_operator.row_offsets_data()[gauge_cell];
         k < pressure_operator.row_offsets_data()[gauge_cell + 1]; ++k) {
        pressure_operator.values_data()[k] = 0.0;
    }
    for (std::size_t row = 0; row < n; ++row) {
        if (row == gauge_cell)
            continue;
        for (std::size_t k = pressure_operator.row_offsets_data()[row];
             k < pressure_operator.row_offsets_data()[row + 1]; ++k) {
            if (pressure_operator.columns_data()[k] == gauge_cell)
                pressure_operator.values_data()[k] = 0.0;
        }
    }
    bool ref_diag_set = false;
    for (std::size_t k = pressure_operator.row_offsets_data()[gauge_cell];
         k < pressure_operator.row_offsets_data()[gauge_cell + 1]; ++k) {
        if (pressure_operator.columns_data()[k] == gauge_cell) {
            pressure_operator.values_data()[k] = 1.0;
            ref_diag_set = true;
            break;
        }
    }
    if (!ref_diag_set) {
        cfdx::core::SparseMatrix rebuilt(n, n);
        for (std::size_t row = 0; row < n; ++row) {
            if (row == gauge_cell) {
                rebuilt.push_back(row, row, 1.0);
                continue;
            }
            for (std::size_t k = pressure_operator.row_offsets_data()[row];
                 k < pressure_operator.row_offsets_data()[row + 1]; ++k) {
                rebuilt.push_back(row, pressure_operator.columns_data()[k],
                                  pressure_operator.values_data()[k]);
            }
        }
        rebuilt.finalize();
        pressure_operator = std::move(rebuilt);
    }
}

// Best-effort contract for the nested pressure solves invoked inside a
// coupled Schur approximation action (PCD Fp/Kp, LSC/BFBt H). These
// callbacks are preconditioner applications, not standalone solves: the
// Schur action itself is only an approximation of the true Schur inverse
// (pcd_schur.h, lsc_bfbt_schur.h), and the outer coupled FGMRES owns the
// final accuracy through its independently recomputed true-residual gate.
// Requiring strict CONVERGED turns a single benign inner shortfall into a
// total outer NOT_APPLICABLE abort: on the N9 Couette campaign one Fp apply
// that could not reach the last digit of its relative tolerance killed the
// whole coupled solve after 168 successful applies, because by then the
// outer residual had shrunk enough for the nested rhs direction to hit its
// own rounding floor. A nested run that exhausted its iteration budget but
// produced a finite best iterate (the Krylov drivers restore the iterate
// with the best verified true residual) remains a valid preconditioner
// action. Only a genuinely unusable result fails the action: divergence,
// breakdown, an inapplicable inner setup, or a non-finite iterate.
inline bool nested_schur_pressure_solve_usable(
    const cfdx::core::SolverResult& result,
    const cfdx::core::Vector& x,
    std::size_t expected_size)
{
    // A solver status is not a substitute for inspecting the actual vector.
    // In particular, CONVERGED must NOT let an out-of-shape/nonfinite nested
    // preconditioner image contaminate the outer FGMRES Arnoldi basis.
    if (result.status != cfdx::core::SolverStatus::CONVERGED &&
        result.status != cfdx::core::SolverStatus::MAX_ITER_REACHED)
        return false;
    if (x.size() != expected_size || !x.is_valid())
        return false;
    for (std::size_t i = 0; i < x.size(); ++i)
        if (!std::isfinite(x(i)))
            return false;
    return std::isfinite(result.residual) &&
           std::isfinite(result.residual_relative) &&
           result.residual >= 0.0 && result.residual_relative >= 0.0;
}

inline cfdx::core::SolverResult solve_coupled_momentum_continuity(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    ScalarEquation ex,
    ScalarEquation ey,
    ScalarEquation ez,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U_old,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& p_old,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    double rho,
    double kinematic_viscosity,
    std::size_t reference_cell,
    double reference_value,
    std::size_t max_iterations,
    double tolerance,
    const DiagnosticsControls& diagnostics,
    CoupledSchurModel schur_model,
    const cfdx::core::LinearSolverRequest& solver_request,
    const IncompressibleSolverControls& controls,
    // Live momentum under-relaxation consumed only by the coupled SIMPLEC
    // Schur approximation. The coupled matrix itself stays unrelaxed;
    // Picard relaxation is applied to the field update after the solve.
    double momentum_under_relaxation,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    cfdx::core::LinearSolverPlan* resolved_linear_plan = nullptr)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nv = 3 * nc;
    const std::size_t np = nc;
    const std::size_t n = nv + np;

    SparseMatrix A(n, n);
    Vector b(n, 0.0);

    // The momentum equations received the current pressure explicitly during
    // segregated assembly. Undo that explicit contribution here: the coupled
    // matrix contains the pressure gradient implicitly instead.
    const auto grad_p_old =
        gauss_gradient_with_boundary(
            p_old, mesh, geometry, pressure_bcs,
            controls.pressure_skew_corrected_gradient,
            controls.pressure_projected_face_interpolation);
    auto add_momentum_block = [&](const ScalarEquation& eq, std::size_t component) {
        const std::size_t row_base = component * nc;
        const auto* ro = eq.matrix.row_offsets_data();
        const auto* co = eq.matrix.columns_data();
        const auto* va = eq.matrix.values_data();
        for (std::size_t r = 0; r < nc; ++r) {
            const double diagonal = eq.diagonal[r];
            if (!(diagonal > 0.0) || !std::isfinite(diagonal))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: invalid momentum diagonal at cell " +
                    std::to_string(r) + " component " + std::to_string(component));

            bool has_diagonal = false;
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t column = component * nc + co[k];
                A.push_back(row_base + r, column, va[k]);
                has_diagonal = has_diagonal || (column == row_base + r);
            }
            // The scalar equation keeps its authoritative diagonal separately
            // from the CSR entries. If a transport row reaches the coupled
            // assembler without an explicit diagonal entry, restore that
            // coefficient from the authoritative metadata instead of creating
            // an empty coupled momentum row.
            if (!has_diagonal)
                A.push_back(row_base + r, row_base + r, diagonal);
            b(row_base + r) =
                eq.rhs(r) + grad_p_old.component_data(component)[r] * geometry.cell_volumes[r];
        }
    };
    add_momentum_block(ex, 0);
    add_momentum_block(ey, 1);
    add_momentum_block(ez, 2);

    // Move the current explicit pressure term from b back into the linear
    // system. The FV Gauss gradient is linear in the pressure field, so its
    // coefficients can be assembled face-by-face.
    for (std::size_t c = 0; c < nc; ++c) {
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                // Same owner/neighbour rule as the continuity block below:
                // if c is the neighbour, nr points back to c itself.
                const std::size_t ncell = owner
                    ? static_cast<std::size_t>(nr)
                    : mesh.ownership().owner(f);
                if (ncell >= nc || ncell == c)
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face topology on face " +
                        std::to_string(f) + " cell " + std::to_string(c));
                // Match gauss_gradient_with_boundary() exactly: the
                // cell-centred pressure gradient uses a distance-weighted
                // face value, not an unconditional arithmetic average.
                // On a skewed mesh dc and dn need not be equal. Using 0.5 here
                // therefore assembled a different pressure-gradient operator
                // from the one used by the segregated momentum equation and
                // from the independent final momentum audit.
                const double dc =
                    (geometry.face_centres[f] - geometry.cell_centres[c]).mag();
                const double dn =
                    (geometry.cell_centres[ncell] - geometry.face_centres[f]).mag();
                if (!(dc > 0.0) || !(dn > 0.0) || !std::isfinite(dc) ||
                    !std::isfinite(dn))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face-centre distances on face " +
                        std::to_string(f));
                const double weight_c = dn / (dc + dn);
                const double weight_n = dc / (dc + dn);
                A.push_back(c, 3*nc + c, Sf.x * weight_c);
                A.push_back(c, 3*nc + ncell, Sf.x * weight_n);
                A.push_back(nc + c, 3*nc + c, Sf.y * weight_c);
                A.push_back(nc + c, 3*nc + ncell, Sf.y * weight_n);
                A.push_back(2*nc + c, 3*nc + c, Sf.z * weight_c);
                A.push_back(2*nc + c, 3*nc + ncell, Sf.z * weight_n);
            } else {
                const std::size_t patch = geometry.face_patch[f];
                // gauss_gradient_with_boundary() omits EMPTY faces.  The
                // coupled G block must omit them too, even on a distorted
                // nominally two-dimensional extrusion.
                if (patch < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(patch).type == PatchType::EMPTY)
                    continue;
                const ScalarBoundaryCondition* bc = nullptr;
                ScalarBoundaryCondition local_bc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto it = pressure_bcs.find(name);
                    if (it != pressure_bcs.end()) {
                        local_bc = it->second;
                        bc = &local_bc;
                    }
                }
                if (bc && bc->type == ScalarBoundaryType::FIXED_VALUE) {
                    // Known boundary pressure is moved to the momentum RHS.
                    b(c) -= Sf.x * bc->value;
                    b(nc + c) -= Sf.y * bc->value;
                    b(2*nc + c) -= Sf.z * bc->value;
                } else {
                    // The Gauss face value is p_owner + gradient *
                    // (C_face-C_cell).n for FIXED_GRADIENT, and simply
                    // p_owner for ZERO_GRADIENT.  The constant gradient
                    // contribution belongs to the RHS, not to G.
                    if (bc && bc->type == ScalarBoundaryType::FIXED_GRADIENT) {
                        const double area = Sf.mag();
                        if (!(area > 0.0) || !std::isfinite(area))
                            throw std::runtime_error(
                                "solve_coupled_momentum_continuity: invalid gradient boundary face area");
                        const Vec3 outward = Sf * (1.0 / area);
                        const Vec3 centre_to_face =
                            geometry.face_centres[f] - geometry.cell_centres[c];
                        const double offset = bc->gradient * centre_to_face.dot(outward);
                        if (!std::isfinite(offset))
                            throw std::runtime_error(
                                "solve_coupled_momentum_continuity: invalid pressure boundary gradient contribution");
                        b(c) -= Sf.x * offset;
                        b(nc + c) -= Sf.y * offset;
                        b(2*nc + c) -= Sf.z * offset;
                    }
                    A.push_back(c, 3*nc + c, Sf.x);
                    A.push_back(nc + c, 3*nc + c, Sf.y);
                    A.push_back(2*nc + c, 3*nc + c, Sf.z);
                }
            }
        }
    }

    // Continuity is assembled from the same face velocity interpolation used
    // by the pressure-based solver. The pressure Schur diagonal below is the
    // exact diagonal of the same momentum-weighted (Rhie-Chow) face-flux
    // derivative and is passed explicitly to the coupled preconditioner.
    std::vector<double> coupled_pressure_schur_diagonal(nc, 0.0);
    // Preserve the exact deferred continuity RHS contribution at the point
    // where it is assembled. This audit channel must not re-evaluate the
    // correction from a later pressure iterate.
    std::vector<double> coupled_explicit_continuity_correction(nc, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        const std::size_t row = nv + c;
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                // neighbour(f) is the opposite cell only on the owner side.
                // When c is the neighbour, neighbour(f) == c; using it directly
                // creates a self-coupled face (zero centre distance) and corrupts
                // the coupled pressure Schur coefficient.
                const std::size_t ncell = owner
                    ? static_cast<std::size_t>(nr)
                    : mesh.ownership().owner(f);
                if (ncell >= nc || ncell == c)
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face topology on face " +
                        std::to_string(f) + " cell " + std::to_string(c));
                const double half_rho = 0.5 * rho;
                A.push_back(row, c, half_rho * Sf.x);
                A.push_back(row, ncell, half_rho * Sf.x);
                A.push_back(row, nc + c, half_rho * Sf.y);
                A.push_back(row, nc + ncell, half_rho * Sf.y);
                A.push_back(row, 2*nc + c, half_rho * Sf.z);
                A.push_back(row, 2*nc + ncell, half_rho * Sf.z);

                const Vec3 rawSf = geometry.face_area_vectors[f];
                const double area = rawSf.mag();
                // Keep the geometric direction aligned with the stored
                // owner-to-neighbour face orientation even when this row is
                // the neighbour cell; the matrix coefficient must be positive
                // and identical in both continuity rows.
                const std::size_t face_owner = mesh.ownership().owner(f);
                const std::size_t face_neighbour = static_cast<std::size_t>(nr);
                const Vec3 dvec = geometry.cell_centres[face_neighbour] -
                                  geometry.cell_centres[face_owner];
                const double d = dvec.mag();
                if (!(area > 0.0) || !std::isfinite(area) ||
                    !(d > 0.0) || !std::isfinite(d))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face geometry on face " +
                        std::to_string(f));
                const Vec3 e{dvec.x / d, dvec.y / d, dvec.z / d};
                const double aox = ex.diagonal[c], anx = ex.diagonal[ncell];
                const double aoy = ey.diagonal[c], any = ey.diagonal[ncell];
                const double aoz = ez.diagonal[c], anz = ez.diagonal[ncell];
                if (!(aox > 0.0) || !(anx > 0.0) ||
                    !(aoy > 0.0) || !(any > 0.0) ||
                    !(aoz > 0.0) || !(anz > 0.0) ||
                    !std::isfinite(aox) || !std::isfinite(anx) ||
                    !std::isfinite(aoy) || !std::isfinite(any) ||
                    !std::isfinite(aoz) || !std::isfinite(anz))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid momentum diagonal on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                // The coupled Schur complement uses the same dAU = V/A_P
                // pressure response as the segregated Rhie-Chow operator.
                const double rAUx_f = 0.5 * (
                    geometry.cell_volumes[c] / aox +
                    geometry.cell_volumes[ncell] / anx);
                const double rAUy_f = 0.5 * (
                    geometry.cell_volumes[c] / aoy +
                    geometry.cell_volumes[ncell] / any);
                const double rAUz_f = 0.5 * (
                    geometry.cell_volumes[c] / aoz +
                    geometry.cell_volumes[ncell] / anz);
                const double rfn = rAUx_f*e.x*e.x + rAUy_f*e.y*e.y + rAUz_f*e.z*e.z;
                if (!(rfn > 0.0) || !std::isfinite(rfn))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face rAU on face " +
                        std::to_string(f));
                const double D = rho * pressure_face_orthogonal_coefficient(
                    rawSf, dvec, {rAUx_f, rAUy_f, rAUz_f});
                if (!(d > 0.0) || !(area > 0.0) ||
                    !std::isfinite(d) || !std::isfinite(area) ||
                    !(D > 0.0) || !std::isfinite(D)) {
                    if (diagnostics.coupled_matrix_summary) {
                        std::cerr << "COUPLED_FACE_DEBUG face=" << f
                                  << " cells=" << c << "/" << ncell
                                  << " area=" << area
                                  << " distance=" << d
                                  << " a=(" << aox << "," << anx
                                  << ";" << aoy << "," << any
                                  << ";" << aoz << "," << anz << ")"
                                  << " rAU=(" << rAUx_f << "," << rAUy_f
                                  << "," << rAUz_f << ")"
                                  << " rfn=" << rfn
                                  << " D=" << D << "\n";
                    }
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid pressure coefficient on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                }
                // phi_p = D (p_P - p_N). Only this orthogonal derivative is
                // implicit in the coupled matrix. The non-orthogonal remainder
                // is a deferred correction evaluated from p_old below, exactly
                // as in the segregated Rhie-Chow operator.
                A.push_back(row, nv + c, D);
                A.push_back(row, nv + ncell, -D);
                coupled_pressure_schur_diagonal[c] += D;

                const Vec3 gp{
                    0.5*(grad_p_old.component_data(0)[c] +
                         grad_p_old.component_data(0)[ncell]),
                    0.5*(grad_p_old.component_data(1)[c] +
                         grad_p_old.component_data(1)[ncell]),
                    0.5*(grad_p_old.component_data(2)[c] +
                         grad_p_old.component_data(2)[ncell])};
                // Reuse the same deferred tensor-response operator as SIMPLE
                // and the pressure-flux reconstruction.
                const double nonorth_flux = rho * pressure_face_nonorthogonal_flux(
                    rawSf, dvec, {rAUx_f, rAUy_f, rAUz_f}, gp);
                // The continuity equation is div(phi)=0 and the pressure
                // correction is written as phi = phi_U - D grad(p).
                // Therefore the explicit non-orthogonal pressure flux is
                // subtracted from the current continuity residual.
                // nonorth_flux is evaluated with the raw owner->neighbour
                // face orientation, so the neighbour-side visit must flip
                // its sign: like the implicit pressure part, which enters
                // every row as D*(p_self - p_other), the deferred flux
                // must enter this cell's outward divergence with the
                // face sign.
                const double compensation = rho * (
                    rAUx_f * gp.x * rawSf.x +
                    rAUy_f * gp.y * rawSf.y +
                    rAUz_f * gp.z * rawSf.z);
                const double explicit_correction = sign * (nonorth_flux - compensation);
                b(row) += explicit_correction;
                coupled_explicit_continuity_correction[c] += explicit_correction;
            } else {
                const std::size_t patch = geometry.face_patch[f];
                // The authoritative mass-flux reconstruction enforces zero
                // normal flux on both EMPTY and SYMMETRY boundaries,
                // regardless of the cell-centred velocity trial state.
                // Therefore neither face contributes a D row or RHS term.
                if (patch < mesh.boundary().n_patches() &&
                    (mesh.boundary().patch(patch).type == PatchType::EMPTY ||
                     mesh.boundary().patch(patch).type == PatchType::SYMMETRY))
                    continue;
                const VelocityBoundaryCondition* vbc = nullptr;
                VelocityBoundaryCondition local_vbc;
                const ScalarBoundaryCondition* pbc = nullptr;
                ScalarBoundaryCondition local_pbc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto vit = velocity_bcs.find(name);
                    if (vit != velocity_bcs.end()) {
                        local_vbc = vit->second;
                        vbc = &local_vbc;
                    }
                    const auto pit = pressure_bcs.find(name);
                    if (pit != pressure_bcs.end()) {
                        local_pbc = pit->second;
                        pbc = &local_pbc;
                    }
                }

                if (vbc && vbc->type == VelocityBoundaryCondition::Type::FIXED_VALUE) {
                    b(row) -= rho * vbc->value.dot(Sf);
                } else {
                    A.push_back(row, c, rho * Sf.x);
                    A.push_back(row, nc + c, rho * Sf.y);
                    A.push_back(row, 2*nc + c, rho * Sf.z);
                }

                if (pbc && pbc->type == ScalarBoundaryType::FIXED_VALUE) {
                    const Vec3 rawSf = geometry.face_area_vectors[f];
                    const Vec3 dvec = geometry.face_centres[f] -
                                      geometry.cell_centres[c];
                    const double a_x = ex.diagonal[c];
                    const double a_y = ey.diagonal[c];
                    const double a_z = ez.diagonal[c];
                    if (!(a_x > 0.0) || !(a_y > 0.0) || !(a_z > 0.0) ||
                        !std::isfinite(a_x) || !std::isfinite(a_y) || !std::isfinite(a_z))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary momentum diagonal on face " +
                            std::to_string(f) + " cell " + std::to_string(c));
                    const double rAUx = geometry.cell_volumes[c] / a_x;
                    const double rAUy = geometry.cell_volumes[c] / a_y;
                    const double rAUz = geometry.cell_volumes[c] / a_z;
                    const double D = rho * pressure_face_orthogonal_coefficient(
                        rawSf, dvec, {rAUx, rAUy, rAUz});
                    if (!(D > 0.0) || !std::isfinite(D))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary pressure coefficient on face " +
                            std::to_string(f));
                    A.push_back(row, nv + c, D);
                    coupled_pressure_schur_diagonal[c] += D;
                    b(row) += D * pbc->value;
                    const Vec3 gp_boundary{
                        grad_p_old.component_data(0)[c],
                        grad_p_old.component_data(1)[c],
                        grad_p_old.component_data(2)[c]};
                    const double deferred_boundary = pressure_face_nonorthogonal_flux(
                        rawSf, dvec, {rAUx, rAUy, rAUz}, gp_boundary);
                    // Match the boundary's deferred pressure response to the
                    // internal-face treatment; otherwise COUPLED drops the
                    // tangential pressure flux at fixed-pressure boundaries.
                    const double compensation_boundary = rho * (
                        rAUx * gp_boundary.x * rawSf.x +
                        rAUy * gp_boundary.y * rawSf.y +
                        rAUz * gp_boundary.z * rawSf.z);
                    const double explicit_boundary_correction =
                        rho * deferred_boundary - compensation_boundary;
                    b(row) += explicit_boundary_correction;
                    coupled_explicit_continuity_correction[c] += explicit_boundary_correction;
                }
            }
        }
    }

    // A pure-Neumann pressure field has one gauge null mode. Eliminate the
    // reference pressure column and replace the reference continuity row by
    // p_ref = reference_value. This keeps the coupled system deterministic.
    const bool has_fixed_pressure =
        has_effective_fixed_pressure_boundary(mesh, pressure_bcs);

    if (!has_fixed_pressure) {
        // Finalize the assembled COO matrix before reading its CSR arrays.
        // Without this explicit finalize(), row_offsets_/values_/columns_ still
        // describe an empty CSR shell, so gauge elimination would silently
        // discard every assembled momentum and continuity coefficient.
        A.finalize();

        // Rebuild the matrix while eliminating the known gauge column.
        SparseMatrix reduced(n, n);
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const auto* va = A.values_data();
        for (std::size_t r = 0; r < n; ++r) {
            if (r == nv + reference_cell) {
                reduced.push_back(r, nv + reference_cell, 1.0);
                b(r) = reference_value;
                continue;
            }
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t col = co[k];
                if (col == nv + reference_cell) {
                    b(r) -= va[k] * reference_value;
                } else {
                    reduced.push_back(r, col, va[k]);
                }
            }
        }
        reduced.finalize();
        A = std::move(reduced);
    } else {
        A.finalize();
    }

    Vector x(n, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        x(c) = U_old.component_data(0)[c];
        x(nc + c) = U_old.component_data(1)[c];
        x(2*nc + c) = U_old.component_data(2)[c];
        x(nv + c) = p_old(c);
    }

    // A reference row exists only for pure-Neumann pressure.  When a fixed
    // pressure boundary is present the row is physical continuity and must
    // retain its genuine warm-start value.  Use the argument that also pins
    // the assembled matrix, not an unrelated control field.
    if (!has_fixed_pressure)
        x(nv + reference_cell) = reference_value;

    // Reject non-finite coefficients before entering GMRES. A bad local
    // coefficient must be diagnosed at assembly time rather than surfacing as
    // an opaque Krylov failure.
    {
        const auto* values = A.values_data();
        for (std::size_t k = 0; k < A.nnz(); ++k) {
            if (!std::isfinite(values[k]))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite matrix coefficient at nnz " +
                    std::to_string(k));
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(b(i)))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite RHS at row " +
                    std::to_string(i));
        }
    }

    // Evidence-only: preserve the raw (pre-pin) Rhie-Chow diagonal for the
    // coupled_assembly_callback when one is installed. The pinned form
    // stored in coupled_pressure_schur_diagonal is diagnostic-only and is
    // never consumed by production preconditioning.
    std::vector<double> pressure_schur_diagonal_raw_copy;
    if (controls.coupled_assembly_callback)
        pressure_schur_diagonal_raw_copy = coupled_pressure_schur_diagonal;

    // Optional algebraic microscope. It is deliberately computed from
    // the exact post-gauge matrix passed to GMRES, so it exposes the actual
    // M/G/D/C blocks rather than reconstructed proxy operators.
    if (diagnostics.coupled_matrix_summary) {
        double m2 = 0.0, g2 = 0.0, d2 = 0.0, c2 = 0.0;
        std::size_t zero_rows = 0;
        std::vector<std::size_t> row_nnz(n, 0);
        std::vector<std::size_t> col_nnz(n, 0);
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const auto* va = A.values_data();

        for (std::size_t r = 0; r < n; ++r) {
            row_nnz[r] = ro[r + 1] - ro[r];
            if (row_nnz[r] == 0) ++zero_rows;
            for (std::size_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t col = co[k];
                const double value = va[k];
                col_nnz[col]++;
                if (r < nv && col < nv)
                    m2 += value * value;
                else if (r < nv && col >= nv)
                    g2 += value * value;
                else if (r >= nv && col < nv)
                    d2 += value * value;
                else
                    c2 += value * value;
            }
        }
        std::size_t zero_cols = 0;
        for (const auto count : col_nnz)
            if (count == 0) ++zero_cols;

        double schur_diag_min = std::numeric_limits<double>::infinity();
        double schur_diag_max = 0.0;
        double schur_mismatch_max = 0.0;
        std::size_t schur_near_zero = 0;
        for (std::size_t c = 0; c < nc; ++c) {
            const std::size_t pr = nv + c;
            double effective = 0.0;
            for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k) {
                const std::size_t u = co[k];
                if (u >= nv) continue;
                const double d = va[k];
                double inv_m = 0.0;
                bool found = false;
                for (std::size_t mk = ro[u]; mk < ro[u + 1]; ++mk) {
                    if (co[mk] == u) {
                        const double diag = va[mk];
                        if (std::isfinite(diag) && std::abs(diag) > 1e-30) {
                            inv_m = 1.0 / diag;
                            found = true;
                        }
                        break;
                    }
                }
                if (!found) continue;
                for (std::size_t gk = ro[u]; gk < ro[u + 1]; ++gk) {
                    const std::size_t q = co[gk];
                    if (q == pr)
                        effective -= d * inv_m * va[gk];
                }
            }
            for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k)
                if (co[k] == pr) effective += va[k];

            schur_diag_min = std::min(schur_diag_min, effective);
            schur_diag_max = std::max(schur_diag_max, std::abs(effective));
            if (std::abs(effective) <=
                1e-12 * std::max(1.0, schur_diag_max))
                ++schur_near_zero;
            if (c < coupled_pressure_schur_diagonal.size() &&
                std::isfinite(coupled_pressure_schur_diagonal[c]))
                schur_mismatch_max = std::max(
                    schur_mismatch_max,
                    std::abs(effective - coupled_pressure_schur_diagonal[c]));
        }

        std::cerr << "COUPLED_MATRIX_SUMMARY n=" << n
                  << " nnz=" << A.nnz()
                  << " M_l2=" << std::sqrt(m2)
                  << " G_l2=" << std::sqrt(g2)
                  << " D_l2=" << std::sqrt(d2)
                  << " C_l2=" << std::sqrt(c2)
                  << " zero_rows=" << zero_rows
                  << " zero_cols=" << zero_cols
                  << " schur_diag_min=" << schur_diag_min
                  << " schur_diag_absmax=" << schur_diag_max
                  << " schur_near_zero=" << schur_near_zero
                  << " rhie_chow_schur_vs_algebraic_diag_delta=" << schur_mismatch_max
                  << " gauge_row=" << (nv + reference_cell)
                  << "\n";
    }

    // The coupled matrix is a saddle-point system. A cell-local 4x4
    // inverse is generally singular because continuity rows have no
    // pressure-pressure diagonal. Select the N8 linear plan before constructing
    // its preconditioner so the production path uses exactly the resolved policy.
    // Only a pure-Neumann system has the reference continuity row replaced by
    // the pressure gauge equation. With a fixed-pressure boundary all continuity
    // rows remain physical, so retain the assembled Rhie-Chow Schur diagonal.
    if (!has_fixed_pressure)
        coupled_pressure_schur_diagonal[reference_cell] = 1.0;
    // The production Couette system is only 512 unknowns. A restart of
    // 128 can stagnate on the nonsymmetric saddle-point spectrum even when
    // the assembled system is nonsingular. Use one full Krylov space for
    // small/medium coupled systems so a restart boundary cannot manufacture
    // an artificial MAX_ITER failure.
    // Feed the dispatcher the characteristics measured from the assembled
    // production CSR matrix. The convenience overload only knows the problem
    // kind and equation count; it cannot see coefficient scaling, symmetry,
    // diagonal dominance or the actual sparsity pattern.
    const auto matrix_characteristics =
        cfdx::core::measure_matrix_characteristics(A, true);
    const auto solver_plan = select_linear_solver(
        LinearProblemKind::CoupledPressureVelocity,
        matrix_characteristics,
        solver_request);
    if (solver_plan.krylov != KrylovModel::GMRES &&
        solver_plan.krylov != KrylovModel::FGMRES)
        throw std::invalid_argument("coupled solver currently requires GMRES or FGMRES");
    if (solver_plan.preconditioner != PreconditionerModel::CoupledBlockSchur &&
        solver_plan.preconditioner != PreconditionerModel::PCD &&
        solver_plan.preconditioner != PreconditionerModel::MGR &&
        solver_plan.preconditioner != PreconditionerModel::BraessSarazin)
        throw std::invalid_argument(
            "coupled solver requires coupled_block_schur, pcd, mgr, or braess_sarazin");

    // N8 owns the linear-solver/preconditioner policy. The coupled FV
    // assembler remains N9-owned, but an explicit BlockSchur request must use
    // the N8 AMG-backed Schur implementation rather than the legacy scalar
    // preconditioner. MGR remains a separate N8-qualified path until its
    // coupled dispatch metadata is integrated.
    //
    // Selection reads the *resolved* plan, not the request, so an automatically
    // resolved policy takes the same branch an explicit one would.
    const bool use_n8_block_schur =
        solver_plan.preconditioner == PreconditionerModel::CoupledBlockSchur ||
        solver_plan.preconditioner == PreconditionerModel::BraessSarazin;
    if (solver_plan.preconditioner == PreconditionerModel::BraessSarazin &&
        schur_model != CoupledSchurModel::BlockLocal)
        throw std::invalid_argument(
            "Braess-Sarazin requires block_local Schur; PCD/LSC/BFBt use their dedicated operators");
    // The Schur family is the production method selector. The outer
    // preconditioner family selects the surrounding block infrastructure;
    // PCD therefore remains reachable even when the outer plan is the generic
    // coupled-block Schur container.
    const bool use_pcd =
        schur_model == CoupledSchurModel::PCD;
    const bool use_lsc_bfbt =
        schur_model == CoupledSchurModel::LSC ||
        schur_model == CoupledSchurModel::BFBT;
    const bool use_simple_schur =
        schur_model == CoupledSchurModel::SIMPLE ||
        schur_model == CoupledSchurModel::SIMPLEC;
    std::unique_ptr<Preconditioner> coupled_preconditioner;
    // Experimental opt-in (#1183): cell-wise overlapping Schwarz for the COUPLED
    // BlockLocal plan. Fail-closed: no fallback preconditioner. Setup runs
    // before the dispatch chain so the experimental path can replace the N8
    // BlockSchur AMG assignment below without altering the chain opening.
    std::unique_ptr<Preconditioner> coupled_schwarz_preconditioner;
    if (controls.coupling.coupled_schwarz_experimental) {
        if (!use_n8_block_schur || schur_model != CoupledSchurModel::BlockLocal)
            throw std::invalid_argument(
                "experimental Schwarz requires the coupled_block_schur linear plan and BlockLocal Schur");
        CellSchwarzConfig config;
        config.cells_per_subdomain = controls.coupling.coupled_schwarz_cells_per_subdomain;
        config.overlap_hops = controls.coupling.coupled_schwarz_overlap_hops;
        config.max_local_unknowns = controls.coupling.coupled_schwarz_max_local_unknowns;
        config.restriction = controls.coupling.coupled_schwarz_weighted
            ? CellSchwarzRestriction::WeightedAdditive
            : CellSchwarzRestriction::Restricted;
        auto sch = std::make_unique<CellOverlappingSchwarzPreconditioner>(nc, config);
        if (!sch->setup(A))
            throw std::runtime_error("COUPLED Schwarz setup failed: local block singular or oversized; no fallback");
        coupled_schwarz_preconditioner = std::move(sch);
    }
    // Evidence-only alias of the production PCD pressure operators, kept
    // alive for the coupled_assembly_callback (the branch-local pcd_ops
    // shared_ptr goes out of scope below). Assigned only when a callback is
    // installed and the plan resolves to PCD; null otherwise.
    std::shared_ptr<PcdPressureOperators> evidence_pcd_ops;
    const auto extract_coupled_block = [&](std::size_t row_block,
                                           std::size_t col_block,
                                           std::size_t row_size,
                                           std::size_t col_size) {
        SparseMatrix block(row_size, col_size);
        const std::size_t velocity_size = 3 * nc;
        const std::size_t row_offset = row_block == 0 ? 0 : velocity_size;
        const std::size_t col_offset = col_block == 0 ? 0 : velocity_size;
        for (std::size_t r = 0; r < row_size; ++r) {
            const std::size_t gr = row_offset + r;
            for (std::size_t k = A.row_offsets_data()[gr];
                 k < A.row_offsets_data()[gr + 1]; ++k) {
                const std::size_t gc = A.columns_data()[k];
                if (gc >= col_offset && gc < col_offset + col_size)
                    block.push_back(r, gc - col_offset, A.values_data()[k]);
            }
        }
        block.finalize();
        return block;
    };
    if (controls.anisotropy_report_callback) {
        AnisotropyGeometry cell_geometry;
        cell_geometry.cell_centers.reserve(nc);
        for (const auto& c : geometry.cell_centres)
            cell_geometry.cell_centers.push_back({c.x, c.y, c.z});
        // Face-quality metrics for actual internal faces, not guessed mesh
        // aspect ratios. A pressure Schur block is not equivalent to C alone.
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            const auto neighbor = mesh.ownership().neighbour(f);
            if (neighbor < 0) continue;
            const auto owner = static_cast<std::size_t>(mesh.ownership().owner(f));
            const auto nb = static_cast<std::size_t>(neighbor);
            const auto d = geometry.cell_centres[nb] - geometry.cell_centres[owner];
            const auto sf = geometry.face_area_vectors[f];
            const double dm = d.mag(), sm = sf.mag();
            if (!(dm > 0.0) || !(sm > 0.0))
                throw std::runtime_error("anisotropy diagnostic: invalid internal face geometry");
            const double cosine = std::clamp(sf.dot(d) / (dm * sm), -1.0, 1.0);
            cell_geometry.face_nonorthogonality_degrees.push_back(
                std::acos(cosine) * (180.0 / 3.14159265358979323846));
            const auto at_face = geometry.face_centres[f] - geometry.cell_centres[owner];
            const auto offset = at_face - d * (at_face.dot(d) / (dm * dm));
            cell_geometry.face_skewness.push_back(offset.mag() / dm);
        }
        // Inspect the component momentum operators, whose geometric
        // connection graph has one entry per physical cell. Never feed cell
        // centers into the 4N saddle matrix: those are different unknowns.
        for (std::size_t comp = 0; comp < 3; ++comp) {
            SparseMatrix component(nc, nc);
            for (std::size_t i = 0; i < nc; ++i) {
                const auto r = comp * nc + i;
                for (std::size_t k = A.row_offsets_data()[r];
                     k < A.row_offsets_data()[r+1]; ++k) {
                    const auto c = A.columns_data()[k];
                    if (c >= comp * nc && c < (comp + 1) * nc)
                        component.push_back(i, c - comp * nc, A.values_data()[k]);
                }
            }
            component.finalize();
            const auto result = analyze_anisotropy(component, &cell_geometry);
            controls.anisotropy_report_callback(
                comp == 0 ? "momentum_x" : comp == 1 ? "momentum_y" : "momentum_z",
                result);
        }
        controls.anisotropy_report_callback("coupled_full_4N", analyze_anisotropy(A));
    }
    if (use_simple_schur) {
        if (!use_n8_block_schur)
            throw std::invalid_argument(
                "SIMPLE/SIMPLEC coupled Schur models require the coupled_block_schur linear plan");

        const auto mode =
            schur_model == CoupledSchurModel::SIMPLE
                ? SimplerSchurMode::SIMPLE
                : SimplerSchurMode::SIMPLEC;
        // The coupled matrix is deliberately unrelaxed, but the SIMPLEC
        // consistent denominator is defined on the under-relaxed momentum
        // row A_P/alpha - sum(a_PN), exactly as the segregated SIMPLEC
        // assembles it from relax_momentum_equation() output. On the
        // unrelaxed steady row the plain row sum A_P - sum(a_PN) is
        // numerically zero for interior cells (only boundary/source terms
        // survive), which made SimplerSchurApproximation::setup() reject
        // the matrix and abort the N9 campaign at the first COUPLED
        // simplec case. Feed the live momentum under-relaxation to the
        // SIMPLEC Schur approximation; SIMPLE keeps its plain diag(Auu)
        // denominator unchanged.
        auto simpler = std::make_unique<SimplerSchurApproximation>(
            mode, std::nullopt,
            schur_model == CoupledSchurModel::SIMPLEC
                ? momentum_under_relaxation
                : 1.0);
        CoupledBlockSchurOptions options;
        options.max_action_amplification =
            controls.coupling.schur_action_growth_guard_experimental ? 1e10 : 0.0;
        options.factorization = CoupledSchurFactorization::Full;
        options.velocity_approximation = CoupledSchurVelocityApproximation::Block;
        options.schur_approximation =
            schur_model == CoupledSchurModel::SIMPLE
                ? CoupledSchurApproximationModel::SIMPLE
                : CoupledSchurApproximationModel::SIMPLEC;

        auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        schur->set_simpler_schur(std::move(simpler));
        if (!schur->setup(A)) {
            throw std::runtime_error(
                std::string("N9 SIMPLE/SIMPLEC Schur setup failed: ") +
                schur->last_error());
        }
        coupled_preconditioner = std::move(schur);
    } else if (controls.coupling.coupled_schwarz_experimental) {
        // Experimental Schwarz (#1183): set up before the dispatch chain;
        // fail-closed there (misconfiguration rejects before any solve).
        coupled_preconditioner = std::move(coupled_schwarz_preconditioner);
    } else if (use_lsc_bfbt) {
        if (!use_n8_block_schur)
            throw std::invalid_argument(
                "LSC/BFBt coupled Schur models require the coupled_block_schur linear plan");

        const SparseMatrix Auu = extract_coupled_block(0, 0, 3 * nc, 3 * nc);
        const SparseMatrix G = extract_coupled_block(0, 1, 3 * nc, nc);
        const SparseMatrix D = extract_coupled_block(1, 0, nc, 3 * nc);
        const SparseMatrix C = extract_coupled_block(1, 1, nc, nc);

        std::vector<double> q_diagonal(3 * nc, 1.0);
        if (schur_model == CoupledSchurModel::LSC) {
            for (std::size_t i = 0; i < 3 * nc; ++i) {
                double diagonal = 0.0;
                for (std::size_t k = Auu.row_offsets_data()[i];
                     k < Auu.row_offsets_data()[i + 1]; ++k) {
                    if (Auu.columns_data()[k] == i)
                        diagonal += Auu.values_data()[k];
                }
                if (!(diagonal > 0.0) || !std::isfinite(diagonal))
                    throw std::runtime_error(
                        "LSC Schur setup requires positive finite momentum diagonal");
                q_diagonal[i] = diagonal;
            }
        }

        const BlockOperator lsc_blocks(Auu, G, D, C);
        SparseMatrix pressure_operator =
            LscBfbtSchurApproximation::assemble_stabilized_pressure_operator(
                lsc_blocks,
                schur_model == CoupledSchurModel::LSC
                    ? LscBfbtSchurApproximation::Mode::LSC
                    : LscBfbtSchurApproximation::Mode::BFBT,
                q_diagonal);

        // Gauge the nested LSC/BFBt pressure operator ONLY when no fixed
        // pressure boundary already removes the constant null mode (same
        // sentinel policy as the reduced coupled system and the outer
        // preconditioner), and mirror that reduced system exactly:
        // reference row AND reference column eliminated with a unit
        // diagonal. The nested solve keeps the raw right-hand side, so
        // the gauge component passes through to the outer Krylov
        // iteration (pinning it freezes the coupled gauge row residual;
        // see apply_nested_schur_pressure_gauge).
        const std::size_t lsc_gauge_cell = has_fixed_pressure
            ? std::numeric_limits<std::size_t>::max() : reference_cell;
        apply_nested_schur_pressure_gauge(pressure_operator, lsc_gauge_cell);

        // The nested LSC/BFBt operator H = D Q^-1 G - C is nonsymmetric
        // by construction and is only consumed by the nested GMRES solve
        // below, whose true-residual checks remain authoritative. Request
        // the general nonsymmetric preconditioner policy: structural setup
        // checks stay fail-closed, but the SPD/CG smoother certificate
        // must not reject the hierarchy for its non-positive diagonals.
        auto pressure_preconditioner =
            std::shared_ptr<cfdx::core::Preconditioner>(
                cfdx::core::make_scalar_preconditioner(
                    cfdx::core::PreconditionerModel::NativeAMG,
                    /*constant_null_space=*/false,
                    /*require_smoother_certificate=*/false).release());
        if (!pressure_preconditioner ||
            !pressure_preconditioner->setup(pressure_operator)) {
            throw std::runtime_error(
                "N8 LSC/BFBt pressure AMG setup failed");
        }

        auto solve_pressure = [pressure_operator,
                               pressure_preconditioner,
                               max_iterations, tolerance](
            const Vector& rhs, Vector& x) mutable {
            const int restart = static_cast<int>(std::min<std::size_t>(
                64, pressure_operator.n_rows()));
            const auto result = solve_gmres(
                pressure_operator, rhs, x, restart,
                max_iterations, tolerance,
                pressure_preconditioner.get());
            return nested_schur_pressure_solve_usable(
                result, x, pressure_operator.n_rows());
        };

        auto algebraic_schur = std::make_unique<LscBfbtSchurApproximation>(
            schur_model == CoupledSchurModel::LSC
                ? LscBfbtSchurApproximation::Mode::LSC
                : LscBfbtSchurApproximation::Mode::BFBT,
            std::move(solve_pressure),
            std::move(q_diagonal));

        CoupledBlockSchurOptions options;
        options.max_action_amplification =
            controls.coupling.schur_action_growth_guard_experimental ? 1e10 : 0.0;
        options.factorization = CoupledSchurFactorization::Full;
        options.velocity_approximation = CoupledSchurVelocityApproximation::Block;
        options.schur_approximation =
            schur_model == CoupledSchurModel::LSC
                ? CoupledSchurApproximationModel::LSC
                : CoupledSchurApproximationModel::BFBT;

        auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        schur->set_algebraic_schur(std::move(algebraic_schur));
        if (!schur->setup(A)) {
            throw std::runtime_error(
                std::string("N9 LSC/BFBt Schur setup failed: ") +
                schur->last_error());
        }
        coupled_preconditioner = std::move(schur);
    } else if (use_pcd) {
        // Trace the actual nested PCD Kp/Fp solves independently of the global
        // (expensive) coupled assembly summary. This is read-only evidence.
        const bool pcd_inner_trace = [] {
            const char* v = std::getenv("CFDX_BFS_PCD_INNER_TRACE");
            return v != nullptr && std::string(v) == "1";
        }();
        auto pcd_ops = std::make_shared<PcdPressureOperators>(
            assemble_pcd_pressure_operators(
                mesh, geometry, U_old, rho, kinematic_viscosity, velocity_bcs, pressure_bcs,
                reference_cell, ConvectionScheme::UPWIND/* exact matrix for every scheme, see assemble_pcd_pressure_operators */,
                controls.coupling.pcd_pressure_inflow_dirichlet_experimental));
        if (controls.coupling.pcd_consistent_gauge_experimental && !has_fixed_pressure) {
            apply_nested_schur_pressure_gauge(pcd_ops->laplacian, reference_cell);
            apply_nested_schur_pressure_gauge(pcd_ops->convection_diffusion, reference_cell);
        }

        // V40: optional sparse comparison of the actual PCD auxiliary rows
        // against C-D diag(Auu)^-1 G on this SAME frozen coupled operator.
        // In particular, a large raw row norm ratio is not, by itself,
        // a BC/stencil defect; the reported fit removes the best scalar.
        // No physical matrix, SST field or preconditioner is changed.
        if (const char* csv_path=std::getenv("CFDX_PCD_SCHUR_ROW_CSV_V40")) {
            if (*csv_path) {
                const auto fit=audit_pcd_schur_rows_v40(
                    A, pcd_ops->convection_diffusion,pcd_ops->laplacian,nc);
                const std::string temp=std::string(csv_path)+".partial";
                {
                    std::ofstream csv(temp, std::ios::trunc);
                    if (!csv) throw std::runtime_error("V40 PCD Schur audit: CSV open failed");
                    csv << "cell,sp_l2,fp_l2,kp_l2,fp_best_multiplier,kp_best_multiplier,"
                           "fp_scaled_shape_defect,kp_scaled_shape_defect,sp_nonzeros,fp_nonzeros,gauge_like\n";
                    csv << std::setprecision(17);
                    for (const auto& r:fit)
                        csv << r.cell << ',' << r.sp_l2 << ',' << r.fp_l2 << ','
                            << r.kp_l2 << ',' << r.fp_best_multiplier << ','
                            << r.kp_best_multiplier << ',' << r.fp_scaled_shape_defect
                            << ',' << r.kp_scaled_shape_defect << ',' << r.sp_nonzeros
                            << ',' << r.fp_nonzeros << ',' << (r.gauge_like?1:0) << '\n';
                    csv.flush();
                    if (!csv) throw std::runtime_error("V40 PCD Schur audit: CSV write failed");
                }
                if (std::rename(temp.c_str(),csv_path)!=0)
                    throw std::runtime_error("V40 PCD Schur audit: CSV rename failed");
                const auto worst=std::max_element(fit.begin(),fit.end(),
                    [](const auto& a,const auto& b) {
                        return a.fp_scaled_shape_defect<b.fp_scaled_shape_defect;
                    });
                std::cerr << "COUPLED_PCD_SCHUR_FIT_V40 cells=" << nc
                          << " worst_cell=" << worst->cell
                          << " worst_scaled_fp_shape=" << worst->fp_scaled_shape_defect
                          << " exact_schur=0 physical_converged=0 original_A_unchanged=1\n";
            }
        }

        // Opt-in read-only capture of Mp/Kp/Fp at the exact production state.
        // External diagnostics may compare their commutator with the ORIGINAL
        // coupled saddle-point A. No residual or preconditioner is changed.
        if (const char* prefix = std::getenv("CFDX_PCD_OPERATORS_CAPTURE_PREFIX_V38")) {
            if (*prefix) {
                const auto slot = capture_pcd_pressure_operators_v38(
                    prefix, pcd_ops->mass, pcd_ops->laplacian,
                    pcd_ops->convection_diffusion);
                if (slot)
                    std::cerr << "COUPLED_PCD_PRESSURE_OPERATORS_V38 slot=" << slot
                              << " rows=" << pcd_ops->mass.n_rows()
                              << " prefix=" << prefix
                              << " placeholder_rhs=1 original_coupled_operator_unchanged=1\n";
            }
        }

        const auto prepare_pcd_solver =
            [](const cfdx::core::SparseMatrix& matrix,
               const cfdx::core::LinearSolverRequest& request,
               bool is_fp_operator) {
                if (request.krylov == cfdx::core::KrylovModel::Auto ||
                    request.preconditioner == cfdx::core::PreconditionerModel::Auto)
                    throw std::invalid_argument(
                        "PCD nested solvers require explicit Krylov and preconditioner models");
                // SPD/CG nested solves keep the fail-closed smoother
                // certificate; nested nonsymmetric Krylov solves use the
                // relaxed policy (outer true-residual checks remain
                // authoritative).
                auto preconditioner =
                    cfdx::core::make_scalar_preconditioner(
                        request.preconditioner,
                        /*constant_null_space=*/false,
                        /*require_smoother_certificate=*/
                        request.krylov == cfdx::core::KrylovModel::CG);
                // Apply the anti-stagnation coarsening policy ONLY to the
                // nonsymmetric nested Fp AMG. Kp/CG and all other AMG users
                // retain their historical multigrid hierarchy unchanged.
                if (is_fp_operator &&
                    request.preconditioner == cfdx::core::PreconditionerModel::NativeAMG &&
                    request.krylov == cfdx::core::KrylovModel::FGMRES) {
                    auto* native = dynamic_cast<cfdx::core::NativeAMGPreconditioner*>(
                        preconditioner.get());
                    if (!native) throw std::logic_error("PCD Fp AMG type mismatch");
                    native->configure_max_coarsening_fraction(0.90);
                    native->configure_bounded_nonsymmetric_action(1e6);
                }
                if (!preconditioner || !preconditioner->setup(matrix))
                    throw std::runtime_error(
                        "PCD nested pressure preconditioner setup failed");
                return preconditioner;
            };

        pcd_ops->laplacian_preconditioner =
            prepare_pcd_solver(pcd_ops->laplacian, controls.pcd_laplacian_linear_solver, false);
        if (controls.coupling.pcd_fp_recirculation_anchor_experimental &&
            controls.coupling.pcd_fp_row_equilibration_experimental)
            throw std::invalid_argument(
                "PCD Fp recirculation anchor and row equilibration cannot be combined");
        if (controls.coupling.pcd_fp_recirculation_anchor_experimental) {
            pcd_ops->convection_diffusion_anchored =
                std::make_unique<PcdFpStrongComponentAnchor>(
                    anchor_pcd_fp_strong_components(pcd_ops->convection_diffusion));
            std::cerr << "COUPLED_PCD_FP_AUX_ANCHORS count="
                      << pcd_ops->convection_diffusion_anchored->anchored_rows.size()
                      << " rows=";
            for (const auto row : pcd_ops->convection_diffusion_anchored->anchored_rows)
                std::cerr << row << ",";
            std::cerr << " original_operator_unchanged=1\n";
        }
        const bool basis_by_directory =
            !controls.coupling.pcd_fp_bisided_basis_directory_v15_experimental.empty();
        const bool basis_by_file =
            !controls.coupling.pcd_fp_bisided_basis_file_experimental.empty();
        if (basis_by_directory && basis_by_file)
            throw std::invalid_argument("V15 select Fp basis directory OR single file, not both");
        if (controls.coupling.pcd_fp_refine_projected_v16_experimental &&
            !controls.coupling.pcd_fp_full_rank_action_v15_experimental)
            throw std::invalid_argument("V16 projected Fp refinement requires V15 full-rank action");
        if (controls.coupling.pcd_fp_full_rank_action_v15_experimental &&
            !basis_by_directory && !basis_by_file)
            throw std::invalid_argument("V15 full-rank Fp action needs a certified basis source");
        if (basis_by_directory || basis_by_file) {
            if (!pcd_ops->convection_diffusion_anchored ||
                controls.coupling.pcd_fp_row_equilibration_experimental ||
                controls.coupling.pcd_fp_second_pass_experimental ||
                controls.coupling.pcd_lumped_mass_commutator_experimental ||
                controls.coupling.pcd_reverse_lumped_mass_v39_experimental ||
                controls.pcd_convection_diffusion_linear_solver.krylov !=
                    cfdx::core::KrylovModel::FGMRES)
                throw std::invalid_argument(
                    "V14 bilateral Fp deflation requires the anchor and FGMRES; no scaling, second pass or commutator");
            pcd_ops->convection_diffusion_bisided =
                std::make_unique<PcdFpBiNullspaceV13>(
                    basis_by_directory
                        ? load_pcd_fp_basis_registry_v15(

                            controls.coupling.pcd_fp_bisided_basis_directory_v15_experimental,
                            pcd_ops->convection_diffusion)
                        : load_pcd_fp_bisided_basis_v13(
                            controls.coupling.pcd_fp_bisided_basis_file_experimental,
                            pcd_ops->convection_diffusion));
            std::cerr << "COUPLED_PCD_BILATERAL_V14 modes="
                      << pcd_ops->convection_diffusion_bisided->left.size()
                      << " left_defect="
                      << pcd_ops->convection_diffusion_bisided->max_left_defect
                      << " right_defect="
                      << pcd_ops->convection_diffusion_bisided->max_right_defect
                      << " basis_registry=" << (basis_by_directory ? 1 : 0)
                      << " full_rank_action=" << (controls.coupling.pcd_fp_full_rank_action_v15_experimental ? 1 : 0)
                      << " original_fp_unchanged=1\n";
        }
        if (controls.coupling.pcd_fp_row_equilibration_experimental)
            pcd_ops->convection_diffusion_equilibrated =
                std::make_unique<PcdFpEquilibration>(
                    equilibrate_pcd_fp_rows(pcd_ops->convection_diffusion));
        const auto& fp_solve_matrix = pcd_ops->convection_diffusion_anchored
            ? pcd_ops->convection_diffusion_anchored->matrix
            : pcd_ops->convection_diffusion_equilibrated
                ? pcd_ops->convection_diffusion_equilibrated->matrix
                : pcd_ops->convection_diffusion;
        pcd_ops->convection_diffusion_preconditioner =
            prepare_pcd_solver(
                fp_solve_matrix, controls.pcd_convection_diffusion_linear_solver, true);
        if (controls.diagnostics.coupled_matrix_summary) {
            const auto stats = equilibrate_pcd_fp_rows(pcd_ops->convection_diffusion);
            std::cerr << "COUPLED_PCD_FP_MATRIX"
                      << " rows=" << stats.matrix.n_rows()
                      << " nnz=" << stats.matrix.nnz()
                      << " smallest_row_max=" << stats.smallest_row_max
                      << " largest_row_max=" << stats.largest_row_max
                      << " smallest_diagonal_ratio=" << stats.smallest_diagonal_ratio
                      << " largest_diagonal_ratio=" << stats.largest_diagonal_ratio
                      << " equilibrated=" << (pcd_ops->convection_diffusion_equilibrated ? 1 : 0)
                      << " pressure_inflow_dirichlet=" << (controls.coupling.pcd_pressure_inflow_dirichlet_experimental ? 1 : 0)
                      << " fp_inner_budget=" << controls.pcd_convection_diffusion_max_iterations
                      << " fp_inner_tolerance=" << controls.pcd_convection_diffusion_tolerance
                      << "\n";
        }

        const auto solve_kp = [pcd_ops, controls, pcd_inner_trace](
            const Vector& rhs, Vector& x) {
            const auto& matrix = pcd_ops->laplacian;
            // A complete PCD commutator validates Kp^{-1} against a separate
            // true-residual certificate at 1e-9 (+1e-12 absolute). A nested
            // solve configured at 1e-8 can legitimately return CONVERGED yet
            // fail that certificate, making the Schur action NOT_APPLICABLE.
            // Request a strictly smaller inner target for this mode only;
            // never loosen either the inverse certificate or the physical gate.
            const double kp_effective_tolerance =
                (controls.coupling.pcd_lumped_mass_commutator_experimental ||
                 controls.coupling.pcd_reverse_lumped_mass_v39_experimental)
                    ? std::min(controls.pcd_laplacian_tolerance, 1e-10)
                    : controls.pcd_laplacian_tolerance;
            auto* preconditioner = pcd_ops->laplacian_preconditioner.get();
            const int restart = std::min<int>(
                controls.pcd_laplacian_linear_solver.gmres_restart,
                static_cast<int>(matrix.n_rows()));
            const Vector kp_initial = x;
            cfdx::core::SolverResult result;
            if (controls.pcd_laplacian_linear_solver.krylov ==
                cfdx::core::KrylovModel::CG) {
                result = solve_cg_prepared(
                    matrix, rhs, x, *preconditioner,
                    controls.pcd_laplacian_max_iterations,
                    kp_effective_tolerance);
            } else if (controls.pcd_laplacian_linear_solver.krylov ==
                       cfdx::core::KrylovModel::GMRES) {
                result = solve_gmres(
                    matrix, rhs, x, restart,
                    controls.pcd_laplacian_max_iterations,
                    kp_effective_tolerance,
                    preconditioner, {}, {}, false);
            } else if (controls.pcd_laplacian_linear_solver.krylov ==
                       cfdx::core::KrylovModel::FGMRES) {
                result = solve_fgmres(
                    matrix, rhs, x, restart,
                    controls.pcd_laplacian_max_iterations,
                    kp_effective_tolerance,
                    preconditioner, {}, {}, false);
            } else {
                throw std::invalid_argument(
                    "unsupported PCD Laplacian Krylov model");
            }
            const auto kp_quality = pcd_fp_true_residual(
                matrix, rhs, kp_initial, x);
            const bool kp_usable = nested_schur_pressure_solve_usable(
                result, x, matrix.n_rows());
            if (controls.diagnostics.coupled_matrix_summary || pcd_inner_trace) {
                std::cerr << "COUPLED_NESTED_PCD_KP"
                          << " requested_tolerance=" << controls.pcd_laplacian_tolerance
                          << " effective_tolerance=" << kp_effective_tolerance
                          << " status=" << static_cast<int>(result.status)
                          << " iterations=" << result.iterations
                          << " reported_residual=" << result.residual
                          << " reported_relative=" << result.residual_relative
                          << " independent_true_l2=" << kp_quality.final_l2
                          << " independent_true_relative_rhs=" << kp_quality.relative_rhs
                          << " independent_true_linf=" << kp_quality.final_linf
                          << " accepted=" << (kp_usable ? 1 : 0) << "\n";
            }
            return kp_usable;
        };

        // Exact nested-system capture: the production CSR, RHS and warm start,
        // BEFORE equilibration, for identical offline ILU0/AMG comparisons.
        // Never alter the Schur action or the outer nonlinear acceptance gate.
        const auto solve_fp = [pcd_ops, controls, pcd_inner_trace](
            const Vector& rhs, Vector& x) {
            const auto& original_matrix = pcd_ops->convection_diffusion;
            if (const char* prefix = std::getenv("CFDX_PCD_FP_CAPTURE_PREFIX")) {
                if (*prefix) {
                    const auto slot = reserve_pcd_fp_capture_slot(false);
                    if (slot) write_pcd_fp_capture(
                        std::string(prefix) + "_rhs" +
                            std::to_string(slot) + ".cfp",
                        original_matrix, rhs, x);
                }
            }
            const auto& solving_matrix = pcd_ops->convection_diffusion_anchored
                ? pcd_ops->convection_diffusion_anchored->matrix
                : pcd_ops->convection_diffusion_equilibrated
                    ? pcd_ops->convection_diffusion_equilibrated->matrix
                    : original_matrix;
            const bool full_rank_deflation =
                controls.coupling.pcd_fp_full_rank_action_v15_experimental &&
                static_cast<bool>(pcd_ops->convection_diffusion_bisided);
            const Vector compatible_rhs = pcd_ops->convection_diffusion_bisided
                ? pcd_ops->convection_diffusion_bisided->compatible_rhs(rhs) : rhs;
            // V14 discarded W^T rhs and the matching Z component.  The
            // resulting PCD map annihilated all left-null directions. V15
            // optional full-rank action uses the ORIGINAL rhs with the
            // regularised *auxiliary* operator (Fp + W Z^T), never with the
            // singular original Fp. Outer coupled A,b are untouched.
            const Vector solving_rhs = full_rank_deflation ? rhs
                : pcd_ops->convection_diffusion_bisided
                ? compatible_rhs : pcd_ops->convection_diffusion_anchored
                ? pcd_ops->convection_diffusion_anchored->rhs(rhs)
                : pcd_ops->convection_diffusion_equilibrated
                    ? pcd_fp_scaled_rhs(rhs,
                        pcd_ops->convection_diffusion_equilibrated->row_multipliers)
                    : rhs;
            if (pcd_ops->convection_diffusion_bisided && !full_rank_deflation)
                x=pcd_ops->convection_diffusion_bisided->canonical_solution(x);
            else if (!pcd_ops->convection_diffusion_bisided &&
                     pcd_ops->convection_diffusion_anchored)
                pcd_ops->convection_diffusion_anchored->pin_initial(x);
            auto* preconditioner = pcd_ops->convection_diffusion_preconditioner.get();
            // Native AMG can SET UP on this nonsymmetric Fp, then return an
            // image A*M^-1*r amplified by 1e14: stop explicitly before a
            // catastrophic Arnoldi loss, without silent ILU0 substitution.
            if (controls.pcd_convection_diffusion_linear_solver.preconditioner ==
                cfdx::core::PreconditionerModel::NativeAMG) {
                Vector initial_residual(rhs.size(), 0.0);
                const auto* ro = solving_matrix.row_offsets_data();
                const auto* co = solving_matrix.columns_data();
                const auto* va = solving_matrix.values_data();
                for (std::size_t row = 0; row < rhs.size(); ++row) {
                    long double ax = 0.0L;
                    for (std::size_t k = ro[row]; k < ro[row+1]; ++k)
                        ax += static_cast<long double>(va[k]) * x(co[k]);
                    initial_residual(row) = solving_rhs(row) - static_cast<double>(ax);
                }
                const auto probe = pcd_fp_action_probe(
                    solving_matrix, *preconditioner, initial_residual);
                const bool safe = probe.applied && probe.finite &&
                    probe.amplification <= 1e6;
                // Do not produce tens of thousands of probe log lines on a
                // successful physical N9 campaign unless explicitly requested.
                const char* trace_env = std::getenv("CFDX_PCD_FP_AMG_TRACE");
                if (!safe || (trace_env && std::string(trace_env) == "1"))
                    std::cerr << "COUPLED_PCD_FP_AMG_ACTION"
                              << " applied=" << (probe.applied ? 1 : 0)
                              << " correction_l2=" << probe.correction_l2
                              << " image_l2=" << probe.image_l2
                              << " amplification=" << probe.amplification
                              << " safe=" << (safe ? 1 : 0) << "\n";
                if (!safe) return false;
            }
            const int restart = std::min<int>(
                controls.pcd_convection_diffusion_linear_solver.gmres_restart,
                static_cast<int>(solving_matrix.n_rows()));
            const Vector x_initial = x;
            const auto solve_once = [&](Vector& iterate, double requested_tolerance) {
                if (pcd_ops->convection_diffusion_bisided) {
                    const auto& projector=*pcd_ops->convection_diffusion_bisided;
                    cfdx::core::LinearOperator op;
                    op.size=original_matrix.n_rows();
                    op.apply=[&](const Vector& in, Vector& out) {
                        out=projector.deflated_matvec(original_matrix,in);
                    };
                    op.apply_fp64=op.apply;
                    // FGMRES handles a variable anchored ILU/AMG action; its
                    // operator is the certified full-rank low-rank update.
                    return solve_fgmres(op,solving_rhs,iterate,restart,
                        controls.pcd_convection_diffusion_max_iterations,
                        requested_tolerance,
                        preconditioner);
                }
                if (controls.pcd_convection_diffusion_linear_solver.krylov ==
                    cfdx::core::KrylovModel::GMRES)
                    return solve_gmres(solving_matrix, solving_rhs, iterate, restart,
                        controls.pcd_convection_diffusion_max_iterations,
                        requested_tolerance,
                        preconditioner, {}, {}, false);
                if (controls.pcd_convection_diffusion_linear_solver.krylov ==
                    cfdx::core::KrylovModel::FGMRES)
                    return solve_fgmres(solving_matrix, solving_rhs, iterate, restart,
                        controls.pcd_convection_diffusion_max_iterations,
                        requested_tolerance,
                        preconditioner, {}, {}, false);
                throw std::invalid_argument("PCD convection-diffusion requires GMRES or FGMRES");
            };
            auto result = solve_once(x, controls.pcd_convection_diffusion_tolerance);
            if (pcd_ops->convection_diffusion_bisided && !full_rank_deflation)
                x=pcd_ops->convection_diffusion_bisided->canonical_solution(x);
            auto quality = pcd_fp_true_residual(original_matrix, rhs, x_initial, x);
            auto auxiliary_quality = pcd_ops->convection_diffusion_bisided
                ? pcd_fp_deflated_true_residual_v14(
                      original_matrix,*pcd_ops->convection_diffusion_bisided,
                      solving_rhs,x_initial,x)
                : pcd_fp_true_residual(
                      solving_matrix,solving_rhs,x_initial,x);
            // Also preserve the first three troublesome production RHS. The
            // first Krylov RHS can be easy even when later Schur actions fail.
            if (const char* prefix = std::getenv("CFDX_PCD_FP_CAPTURE_DIFFICULT_PREFIX")) {
                if (*prefix && (!std::isfinite(quality.relative_rhs) ||
                                quality.relative_rhs > 1e-3)) {
                    const auto slot = reserve_pcd_fp_capture_slot(true);
                    if (slot) write_pcd_fp_capture(
                        std::string(prefix) + "_rhs" +
                            std::to_string(slot) + ".cfp",
                        original_matrix, rhs, x_initial);
                }
            }
            const auto usable = [&](const cfdx::core::SolverResult& candidate,
                                    const PcdFpTrueResidual& q,
                                    const PcdFpTrueResidual& aq) {
                const bool vector_ok = nested_schur_pressure_solve_usable(
                    candidate, x, original_matrix.n_rows());
                if (pcd_ops->convection_diffusion_bisided) {
                    // The original incompatible component is NOT discarded
                    // from evidence or the outer coupled physical system.
                    const auto& projector=*pcd_ops->convection_diffusion_bisided;
                    const double projected=projector.projected_residual_relative(
                        original_matrix,rhs,x);
                    return vector_ok && candidate.status==SolverStatus::CONVERGED &&
                        std::isfinite(aq.relative_rhs) &&
                        aq.relative_rhs <= controls.pcd_convection_diffusion_tolerance &&
                        std::isfinite(projected) &&
                        projected<=controls.pcd_convection_diffusion_tolerance &&
                        std::isfinite(q.relative_rhs) &&
                        q.final_l2 <= q.initial_l2;
                }
                if (pcd_ops->convection_diffusion_anchored) {
                    // This is explicitly a PCD preconditioner approximation.
                    // The auxiliary residual must reach the ORIGINAL inner
                    // tolerance; original Fp is singular and cannot meet it.
                    // Reject any action that increases the original true
                    // residual, and never report the auxiliary one as original.
                    return vector_ok && candidate.status == SolverStatus::CONVERGED &&
                        std::isfinite(aq.relative_rhs) &&
                        aq.relative_rhs <= controls.pcd_convection_diffusion_tolerance &&
                        std::isfinite(q.initial_l2) && q.initial_l2 > 0.0 &&
                        std::isfinite(q.final_l2) &&
                        q.final_l2 < 0.999 * q.initial_l2;
                }
                return vector_ok && std::isfinite(q.relative_rhs) &&
                    (!controls.coupling.pcd_fp_quality_guard_experimental ||
                     q.relative_rhs <= 1e-3);
            };
            if (controls.coupling.pcd_fp_second_pass_experimental &&
                result.status == SolverStatus::MAX_ITER_REACHED &&
                quality.relative_rhs > controls.pcd_convection_diffusion_tolerance &&
                x.is_valid()) {
                // Restart from the best first-pass iterate with an additional
                // independent Krylov cycle; keep only an improvement measured
                // on the UNSCALED, original Fp matrix.
                const Vector best_x = x;
                const auto best_quality = quality;
                const auto best_result = result;
                auto retry_result = solve_once(x, controls.pcd_convection_diffusion_tolerance);
                const auto retry_quality = pcd_fp_true_residual(
                    original_matrix, rhs, x_initial, x);
                const bool improved = std::isfinite(retry_quality.final_l2) &&
                    retry_quality.final_l2 < best_quality.final_l2 &&
                    nested_schur_pressure_solve_usable(
                        retry_result, x, original_matrix.n_rows());
                if (improved) {
                    quality = retry_quality;
                    auxiliary_quality = pcd_fp_true_residual(
                        solving_matrix, solving_rhs, x_initial, x);
                    result = retry_result;
                } else {
                    x = best_x;
                    quality = best_quality;
                    auxiliary_quality = pcd_fp_true_residual(
                        solving_matrix, solving_rhs, x_initial, x);
                    result = best_result;
                }
                if (controls.diagnostics.coupled_matrix_summary)
                    std::cerr << "COUPLED_PCD_FP_SECOND_PASS"
                              << " improved=" << (improved ? 1 : 0)
                              << " initial_true_relative_rhs=" << best_quality.relative_rhs
                              << " returned_true_relative_rhs=" << quality.relative_rhs
                              << " retry_iterations=" << retry_result.iterations << "\n";
            }
            bool accepted = usable(result, quality, auxiliary_quality);
            if (!accepted && full_rank_deflation &&
                controls.coupling.pcd_fp_refine_projected_v16_experimental) {
                // Do NOT loosen either certificate. A restarted Krylov call
                // at the SAME tolerance is a no-op once the auxiliary norm
                // passes. Request a strictly tighter auxiliary solve, then
                // recompute BOTH residuals against their original operators.
                const auto report = refine_pcd_fp_v16(
                    x, result, controls.pcd_convection_diffusion_tolerance, 2,
                    [&](Vector& iterate, double stricter) {
                        return solve_once(iterate, stricter);
                    },
                    [&](const cfdx::core::SolverResult& candidate, const Vector& iterate) {
                        (void)iterate;
                        quality = pcd_fp_true_residual(original_matrix, rhs, x_initial, x);
                        auxiliary_quality = pcd_fp_deflated_true_residual_v14(
                            original_matrix, *pcd_ops->convection_diffusion_bisided,
                            solving_rhs, x_initial, x);
                        return usable(candidate, quality, auxiliary_quality);
                    });
                accepted = report.certified;
                std::cerr << "COUPLED_PCD_FP_REFINE_V16 extra_solves="
                          << report.extra_solves
                          << " requested_tolerance=" << report.last_requested_tolerance
                          << " deflated_true_relative=" << auxiliary_quality.relative_rhs
                          << " projected_original_true_relative="
                          << pcd_ops->convection_diffusion_bisided->projected_residual_relative(
                                 original_matrix, rhs, x)
                          << " accepted=" << (accepted ? 1 : 0)
                          << " original_operator_unchanged=1\n";
            }
            if (controls.diagnostics.coupled_matrix_summary || pcd_inner_trace)
                std::cerr << "COUPLED_NESTED_PCD_FP"
                          << " status=" << static_cast<int>(result.status)
                          << " iterations=" << result.iterations
                          << " reported_residual=" << result.residual
                          << " reported_relative=" << result.residual_relative
                          << " original_initial_l2=" << quality.initial_l2
                          << " original_true_l2=" << quality.final_l2
                          << " original_true_relative_rhs=" << quality.relative_rhs
                          << " original_true_relative_initial=" << quality.relative_initial
                          << " original_true_linf=" << quality.final_linf
                          << " bilateral_deflation=" << (pcd_ops->convection_diffusion_bisided ? 1 : 0)
                          << " full_rank_action=" << (full_rank_deflation ? 1 : 0)
                          << " left_incompatible_floor=" << (pcd_ops->convection_diffusion_bisided
                              ? pcd_ops->convection_diffusion_bisided->incompatibility_relative(rhs) : 0.0)
                          << " projected_original_true_relative=" << (pcd_ops->convection_diffusion_bisided
                              ? pcd_ops->convection_diffusion_bisided->projected_residual_relative(
                                  original_matrix,rhs,x) : -1.0)
                          << " auxiliary_operator=" << (pcd_ops->convection_diffusion_bisided ? "deflated"
                              : pcd_ops->convection_diffusion_anchored ? "anchored"
                              : pcd_ops->convection_diffusion_equilibrated ? "equilibrated" : "original")
                          << " auxiliary_anchored=" << (pcd_ops->convection_diffusion_anchored &&
                               !pcd_ops->convection_diffusion_bisided ? 1 : 0)
                          << " preconditioner_anchored=" << (pcd_ops->convection_diffusion_anchored ? 1 : 0)
                          << " auxiliary_true_relative_rhs=" << auxiliary_quality.relative_rhs
                          << " auxiliary_anchors=" << (pcd_ops->convection_diffusion_anchored
                              ? pcd_ops->convection_diffusion_anchored->anchored_rows.size() : 0)
                          << " accepted=" << (accepted ? 1 : 0) << "\n";
            return accepted;
        };

        auto pcd = std::make_unique<PcdSchurApproximation>(
            pcd_ops->mass, pcd_ops->laplacian, pcd_ops->convection_diffusion,
            std::move(solve_kp), std::move(solve_fp),
            std::nullopt,
            has_fixed_pressure ? std::optional<std::size_t>()
                               : std::optional<std::size_t>(reference_cell),
            controls.coupling.pcd_reverse_lumped_mass_v39_experimental
                ? PcdSchurActionMode::ReverseLumpedMassCommutatorExperimental
                : controls.coupling.pcd_lumped_mass_commutator_experimental
                    ? PcdSchurActionMode::LumpedMassCommutator
                    : PcdSchurActionMode::LegacySingleConvectionDiffusionInverse);
        // Evidence-only: expose the production operators to the
        // coupled_assembly_callback when one is installed (kept alive via
        // evidence_pcd_ops past this branch scope).
        if (controls.coupled_assembly_callback)
            evidence_pcd_ops = pcd_ops;

        CoupledBlockSchurOptions options;
        options.max_action_amplification =
            controls.coupling.schur_action_growth_guard_experimental ? 1e10 : 0.0;
        options.factorization = CoupledSchurFactorization::Full;
        options.velocity_approximation = CoupledSchurVelocityApproximation::Block;
        switch (schur_model) {
            case CoupledSchurModel::BlockLocal:
                options.schur_approximation = CoupledSchurApproximationModel::BlockLocal;
                break;
            case CoupledSchurModel::PCD:
                options.schur_approximation = CoupledSchurApproximationModel::PCD;
                break;
            case CoupledSchurModel::LSC:
            case CoupledSchurModel::BFBT:
            case CoupledSchurModel::SIMPLE:
            case CoupledSchurModel::SIMPLEC:
                throw std::invalid_argument(
                    "selected algebraic Schur model is handled by its dedicated production path");
            default:
                throw std::invalid_argument("unsupported coupled Schur model");
        }

        auto schur = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        if (schur_model == CoupledSchurModel::PCD)
            schur->set_pcd_schur(std::move(pcd));
        // Set the preconditioner up here for the same reason the BlockSchur
        // branch below does. solve_fgmres also calls setup() on its
        // preconditioner before the first apply, so this is not what makes PCD
        // reachable; performing it explicitly is what makes a setup failure
        // report its own cause. Without it, a rejected PCD setup surfaces only
        // as SolverStatus::NOT_APPLICABLE from the Krylov driver and the
        // PcdSchurApproximation last_error() is lost. reset() preserves
        // pcd_schur_, so the driver's second setup is clean.
        if (!schur->setup(A)) {
            throw std::runtime_error(
                std::string("N8 coupled PCD Schur setup failed: ") +
                schur->last_error());
        }
        coupled_preconditioner = std::move(schur);
    } else if (use_n8_block_schur) {
        CoupledBlockSchurOptions options;
        options.max_action_amplification =
            controls.coupling.schur_action_growth_guard_experimental ? 1e10 : 0.0;
        if (solver_plan.preconditioner == PreconditionerModel::BraessSarazin) {
            options.velocity_approximation = CoupledSchurVelocityApproximation::Diagonal;
            options.velocity_inverse_scale = 1.0 / solver_request.braess_sarazin_alpha;
        }
        auto schur_amg = std::make_unique<CoupledBlockSchurAMGPreconditioner>(
            nc, options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        if (!schur_amg->setup(A)) {
            throw std::runtime_error(
                std::string("N8 coupled BlockSchur AMG setup failed: ") +
                schur_amg->last_error());
        }
        coupled_preconditioner = std::move(schur_amg);
    } else {
        std::vector<std::size_t> fine_variables(nv);
        std::vector<std::size_t> coarse_variables(np);
        std::iota(fine_variables.begin(), fine_variables.end(), 0);
        std::iota(coarse_variables.begin(), coarse_variables.end(), nv);
        coupled_preconditioner =
            std::make_unique<NativeMGRPreconditioner>(
                std::move(fine_variables), std::move(coarse_variables));
    }

    const bool automatic_coupled = solver_request.krylov == KrylovModel::Auto &&
        solver_request.preconditioner == PreconditionerModel::Auto;
    if (resolved_linear_plan != nullptr)
        *resolved_linear_plan = solver_plan;
    const int gmres_restart = automatic_coupled
        ? static_cast<int>(A.n_rows() <= 4096
                               ? std::min<std::size_t>(A.n_rows(),
                                                       std::max<std::size_t>(1024, max_iterations))
                               : std::min<std::size_t>(A.n_rows(), 512))
        : std::min<int>(solver_request.gmres_restart, static_cast<int>(A.n_rows()));
    // The coupled acceptance problem is small enough for a full Krylov space.
    // Do not let the global adaptive-restart policy silently clamp the restart
    // back to its generic [10,40] range: that was the direct cause of the
    // observed residual plateau around 1.33e-6. The restart must also be the
    // FULL space for small/medium coupled systems (<= 1024 unknowns): the
    // 12x16 channel system has 768 unknowns, so a 512 cap truncated the
    // Krylov space below the system dimension and manufactured a plateau
    // near 1e-7 on the Poiseuille Picard system, converting a solvable
    // 1e-12 gate into an artificial MAX_ITER failure. Medium coupled systems
    // (1024 < n <= 4096) get a restart of max(1024, max_iterations), capped
    // by the system dimension: the 16x32 cavity Picard system has 4096
    // unknowns, and a 512 cap truncated the Krylov space mid-budget
    // (1000 iterations = 512 + 488), stalling the true relative residual at
    // 3.66e-12 against the 1e-12 gate with three identical relaxation
    // retries. A 1024 restart is never truncated within the iteration
    // budget, so convergence is decided by the budget, not by a restart
    // boundary. Systems above 4096 unknowns keep the 512 cap for memory.
    // Keep the restart fixed for this block solve; adaptive restart remains
    // available to generic GMRES callers.
    cfdx::core::KrylovControls coupled_gmres_controls;
    coupled_gmres_controls.restart_min = gmres_restart;
    coupled_gmres_controls.restart_max = gmres_restart;
    coupled_gmres_controls.adaptive_restart = false;
    if (diagnostics.coupled_matrix_summary) {
        std::cerr << "COUPLED_PRECONDITIONER name=" << coupled_preconditioner->name()
                  << " restart=" << gmres_restart
                  << " adaptive_restart=0\n";
    }

    // Evidence-only capture hook, at the exact Krylov call point (null by
    // default; see IncompressibleSolverControls::coupled_assembly_callback).
    // Everything the solver consumes is now final: assembly, gauge
    // treatment, finite checks and preconditioner setup are done, and no
    // mutation happens between this hook and the solve calls below. The
    // snapshot therefore aliases the very A, b and x objects handed to the
    // Krylov driver - the captured system IS the solved system, by
    // construction rather than by convention. The consumer must still copy
    // anything it wants to keep: the references die with this scope, and x
    // becomes the solution below.
    if (controls.coupled_assembly_callback) {
        const CoupledAssemblySnapshot snapshot{
            &A, &b, &x,
            &coupled_pressure_schur_diagonal,
            &pressure_schur_diagonal_raw_copy,
            evidence_pcd_ops ? &evidence_pcd_ops->mass : nullptr,
            evidence_pcd_ops ? &evidence_pcd_ops->laplacian : nullptr,
            evidence_pcd_ops ? &evidence_pcd_ops->convection_diffusion : nullptr,
            &U_old, reference_cell, has_fixed_pressure,
            {{&ex, &ey, &ez}}, &p_old};
        controls.coupled_assembly_callback(snapshot);
    }

    // Retain the exact initial guess only when explicit small-system ILUT
    // recovery is enabled.  Ordinary six-Schur production behavior is unchanged.
    const bool allow_v16_ilut_recovery =
        controls.coupling.coupled_monolithic_ilut_fallback_v16_experimental;
    if (allow_v16_ilut_recovery && A.n_rows() > 4096)
        throw std::invalid_argument(
            "V16 monolithic ILUT recovery is restricted to <=4096 unknowns; "
            "use native Schur/AMG/MGR for large CPU/GPU meshes");
    const bool allow_v42_blocklocal_recovery =
        controls.coupling.pcd_blocklocal_recovery_v42_experimental;
    if (allow_v42_blocklocal_recovery && (schur_model != CoupledSchurModel::PCD ||
        allow_v16_ilut_recovery))
        throw std::invalid_argument(
            "V42 BlockLocal recovery requires PCD and cannot combine with ILUT");
    const Vector original_guess_for_recovery =
        (allow_v16_ilut_recovery || allow_v42_blocklocal_recovery) ? x : Vector{};
    const Vector& initial_guess_v16 = original_guess_for_recovery;

    const std::size_t primary_budget = allow_v42_blocklocal_recovery &&
        controls.coupling.pcd_blocklocal_preflight_iterations_v42 != 0
        ? std::min(max_iterations,
            controls.coupling.pcd_blocklocal_preflight_iterations_v42)
        : max_iterations;
    SolverResult result;
    if (solver_plan.krylov == KrylovModel::FGMRES) {
        result = solve_fgmres(
            A, b, x, gmres_restart, primary_budget, tolerance,
            coupled_preconditioner.get(), coupled_gmres_controls);
    } else {
        result = solve_gmres(
            A, b, x, gmres_restart, primary_budget, tolerance,
            coupled_preconditioner.get(), coupled_gmres_controls);
    }

    // V42: after a failed PCD solve, perform a separately named BlockLocal
    // retry on IDENTICAL A,b and ORIGINAL x0. This is not a better PCD action
    // and must not be counted as PCD convergence in six-method comparisons.
    // No dense factorization or small-mesh cap; uses existing AMG BlockSchur.
    // Explicit opt-in only, never used by any of the other five variants.
    if (allow_v42_blocklocal_recovery && result.status != SolverStatus::CONVERGED) {
        const auto first_status=result.status;
        const auto first_iterations=result.iterations;
        CoupledBlockSchurOptions retry_options;
        retry_options.max_action_amplification =
            controls.coupling.schur_action_growth_guard_experimental ? 1e10 : 0.0;
        CoupledBlockSchurAMGPreconditioner blocklocal(
            nc, retry_options, has_fixed_pressure
                ? std::numeric_limits<std::size_t>::max() : reference_cell);
        if (!blocklocal.setup(A)) {
            std::cerr << "COUPLED_PCD_BLOCKLOCAL_RECOVERY_V42 attempted=1 adopted=0"
                      << " reason=setup_failed"
                      << " primary_status=" << static_cast<int>(first_status)
                      << " original_A_unchanged=1 physical_qualified=0\n";
        } else {
            Vector candidate=original_guess_for_recovery;
            const auto alternative=solve_fgmres(
                A,b,candidate,gmres_restart,max_iterations,tolerance,
                &blocklocal,coupled_gmres_controls);
            const auto quality=certify_pcd_blocklocal_recovery_v42(
                A,b,candidate,nv,tolerance,
                alternative.status==SolverStatus::CONVERGED);
            std::cerr << std::scientific << std::setprecision(12)
                      << "COUPLED_PCD_BLOCKLOCAL_RECOVERY_V42 attempted=1"
                      << " adopted=" << (quality.adopt?1:0)
                      << " primary_model=PCD recovery_model=BlockLocal"
                      << " primary_status=" << static_cast<int>(first_status)
                      << " primary_iterations=" << first_iterations
                      << " primary_budget=" << primary_budget
                      << " recovery_status=" << static_cast<int>(alternative.status)
                      << " recovery_iterations=" << alternative.iterations
                      << " true_relative=" << quality.true_relative
                      << " momentum_relative="
                      << quality.blocks.momentum.relative_to_block_rhs()
                      << " continuity_and_gauge_relative="
                      << quality.blocks.continuity_and_gauge.relative_to_block_rhs()
                      << " unchanged_tolerance=" << tolerance
                      << " original_A_unchanged=1 physical_qualified=0\n";
            if (quality.adopt) {
                x=std::move(candidate);
                result=alternative;
                result.iterations+=first_iterations;
                if (resolved_linear_plan) {
                    resolved_linear_plan->krylov=KrylovModel::FGMRES;
                    resolved_linear_plan->preconditioner=PreconditionerModel::CoupledBlockSchur;
                    resolved_linear_plan->automatic_preconditioner=false;
                    resolved_linear_plan->reason+=
                        "; explicit V42 BlockLocal recovery after failed PCD";
                }
            }
        }
    }

    // V16: a separately selected, explicit experimental recovery of a
    // NONCONVERGED exterior Schur action, never an implicit solver-selection
    // change. The fallback approximates the full saddle matrix directly and
    // is accepted ONLY by the independent full original-system A*x-b gate.
    if (allow_v16_ilut_recovery && result.status != SolverStatus::CONVERGED) {
        const auto original_status = result.status;
        const auto original_iterations = result.iterations;
        CoupledMonolithicIlutV16 ilut(1.e-5,20);
        if (!ilut.setup(A)) {
            std::cerr << "COUPLED_MONOLITHIC_ILUT_V16 attempted=1 adopted=0 "
                         "reason=setup_failed original_schur="
                      << to_string(schur_model) << "\n";
        } else {
            Vector candidate=initial_guess_v16;
            cfdx::core::KrylovControls recovery_controls;
            const int recovery_restart=static_cast<int>(
                std::min<std::size_t>(160,A.n_rows()));
            recovery_controls.restart_min=recovery_restart;
            recovery_controls.restart_max=recovery_restart;
            recovery_controls.adaptive_restart=false;
            auto recovery=solve_fgmres(A,b,candidate,recovery_restart,
                                       max_iterations,tolerance,&ilut,recovery_controls);
            const auto certificate=independent_coupled_block_residual_v16(
                A,b,candidate,nv);
            const auto true_relative=static_cast<double>(
                certificate.relative_to_global_rhs());
            const bool adopted = recovery.status==SolverStatus::CONVERGED &&
                std::isfinite(true_relative) && true_relative<=tolerance;
            std::cerr << "COUPLED_MONOLITHIC_ILUT_V16 attempted=1 adopted="
                      << (adopted?1:0)
                      << " original_schur=" << to_string(schur_model)
                      << " primary_status=" << static_cast<int>(original_status)
                      << " primary_iterations=" << original_iterations
                      << " recovery_status=" << static_cast<int>(recovery.status)
                      << " recovery_iterations=" << recovery.iterations
                      << " recovery_true_relative=" << true_relative
                      << " momentum_relative_block="
                      << certificate.momentum.relative_to_block_rhs()
                      << " continuity_relative_block="
                      << certificate.continuity_and_gauge.relative_to_block_rhs()
                      << " regularized_preconditioner_pivots="
                      << ilut.regularized_pivots()
                      << " target_unchanged=" << tolerance
                      << " physical_qualified=0\n";
            if (adopted) {
                result=recovery;
                result.iterations+=original_iterations;
                x=std::move(candidate);
                if (resolved_linear_plan) {
                    resolved_linear_plan->krylov=KrylovModel::FGMRES;
                    resolved_linear_plan->preconditioner=PreconditionerModel::ILUT;
                    resolved_linear_plan->automatic_preconditioner=false;
                    resolved_linear_plan->reason+=
                        "; explicit V16 monolithic ILUT recovery after failed Schur";
                }
            }
        }
    }

    // Evidence-only: inspect the pressure/continuity rows on the exact
    // monolithic candidate returned by Krylov, before the caller applies
    // separate U/p relaxation. Exclude the pinned reference row for pure
    // Neumann pressure, because that row is a gauge equation, not continuity.
    // This is deliberately diagnostic only: it does not alter SolverResult
    // or any physical acceptance criterion.
    if (diagnostics.coupled_matrix_summary) {
        double continuity_row_linf = 0.0;
        double continuity_row_signed = 0.0;
        std::size_t continuity_row_cell = 0;
        for (std::size_t c = 0; c < nc; ++c) {
            if (!has_fixed_pressure && c == reference_cell)
                continue;
            const std::size_t row = nv + c;
            double ax = 0.0;
            for (std::size_t k = A.row_offsets_data()[row];
                 k < A.row_offsets_data()[row + 1]; ++k)
                ax += A.values_data()[k] * x(A.columns_data()[k]);
            const double residual = ax - b(row);
            if (std::abs(residual) > continuity_row_linf) {
                continuity_row_linf = std::abs(residual);
                continuity_row_signed = residual;
                continuity_row_cell = c;
            }
        }
        std::cerr << "COUPLED_PRE_RELAXATION_CONTINUITY_ROW"
                  << " linf=" << std::scientific << std::setprecision(8)
                  << continuity_row_linf
                  << " signed=" << continuity_row_signed
                  << " cell=" << continuity_row_cell
                  << " status=" << static_cast<int>(result.status)
                  << " iterations=" << result.iterations
                  << " reported_relative=" << result.residual_relative
                  << "\n";
        {
            const std::size_t row = nv + continuity_row_cell;
            double velocity_block = 0.0;
            double pressure_block = 0.0;
            for (std::size_t k = A.row_offsets_data()[row];
                 k < A.row_offsets_data()[row + 1]; ++k) {
                const std::size_t column = A.columns_data()[k];
                const double term = A.values_data()[k] * x(column);
                if (column < nv) velocity_block += term;
                else pressure_block += term;
            }
            std::cerr << "COUPLED_CONTINUITY_ROW_BREAKDOWN"
                      << " cell=" << continuity_row_cell
                      << " velocity_block=" << velocity_block
                      << " pressure_block=" << pressure_block
                      << " rhs_including_explicit_correction=" << b(row)
                      << " lhs_minus_rhs=" << (velocity_block + pressure_block - b(row))
                      << "\n";
        }
    }

    // Recompute the independent true residual even for a failed Krylov
    // attempt. An FGMRES status is not an independent measurement of A*x-b.
    const auto krylov_status_before_certificate = result.status;
    const double krylov_reported_before_certificate = result.residual_relative;

    // Independently recompute ||b-Ax|| after the linear solve.  The Krylov
    // result already carries a true residual, but keeping an explicit
    // post-solve matrix check here makes the coupled algebraic gate independent
    // from the nonlinear FV diagnostics and exposes any stale/inconsistent
    // solver bookkeeping.
    // V16: attribute failure to momentum or continuity rows *on the exact
    // same candidate and CSR*, without modifying the solver convergence gate.
    const auto block_certificate_v16 = independent_coupled_block_residual_v16(
        A, b, x, nv);
    const long double coupled_matrix_residual = block_certificate_v16.global_l2;
    const double coupled_matrix_relative = static_cast<double>(
        block_certificate_v16.relative_to_global_rhs());
    if (krylov_status_before_certificate != SolverStatus::CONVERGED ||
        diagnostics.coupled_matrix_summary) {
        std::cerr << "COUPLED_BLOCK_TRUE_RESIDUAL_V16"
                  << " momentum_l2=" << block_certificate_v16.momentum.absolute_l2
                  << " momentum_rhs_l2=" << block_certificate_v16.momentum.rhs_l2
                  << " momentum_relative_block="
                  << block_certificate_v16.momentum.relative_to_block_rhs()
                  << " momentum_worst_row=" << block_certificate_v16.momentum.worst_row
                  << " momentum_linf=" << block_certificate_v16.momentum.linf
                  << " continuity_l2="
                  << block_certificate_v16.continuity_and_gauge.absolute_l2
                  << " continuity_rhs_l2="
                  << block_certificate_v16.continuity_and_gauge.rhs_l2
                  << " continuity_relative_block="
                  << block_certificate_v16.continuity_and_gauge.relative_to_block_rhs()
                  << " continuity_worst_row="
                  << block_certificate_v16.continuity_and_gauge.worst_row
                  << " continuity_linf=" << block_certificate_v16.continuity_and_gauge.linf
                  << " momentum_squared_defect_fraction="
                  << block_certificate_v16.momentum_fraction_of_squared_defect()
                  << " full_relative=" << coupled_matrix_relative
                  << " krylov_status=" << static_cast<int>(krylov_status_before_certificate)
                  << " iterations=" << result.iterations << "\n";
    }
    if (!std::isfinite(coupled_matrix_residual) ||
        !std::isfinite(coupled_matrix_relative)) {
        result.status = SolverStatus::DIVERGED;
        result.residual = std::numeric_limits<double>::infinity();
        result.residual_relative = std::numeric_limits<double>::infinity();
        return result;
    }
    // The independently recomputed residual is the authoritative algebraic
    // gate. Never accept a Krylov-reported convergence if b-Ax is still above
    // the requested coupled linear tolerance.
    result.residual = static_cast<double>(coupled_matrix_residual);
    result.residual_relative = coupled_matrix_relative;
    if (diagnostics.coupled_matrix_summary) {
        std::cerr << "COUPLED_LINEAR_INDEPENDENT_V15"
                  << " krylov_status=" << static_cast<int>(krylov_status_before_certificate)
                  << " krylov_reported=" << krylov_reported_before_certificate
                  << " long_double_true_relative=" << coupled_matrix_relative
                  << " target=" << tolerance << "\n";
    }
    if (krylov_status_before_certificate != SolverStatus::CONVERGED)
        return result;
    if (coupled_matrix_relative > tolerance) {
        result.status = SolverStatus::DIVERGED;
        return result;
    }
    if (diagnostics.coupled_matrix_summary) {
        std::cerr << "COUPLED_LINEAR_TRUE_RESIDUAL norm="
                  << coupled_matrix_residual
                  << " relative=" << coupled_matrix_relative
                  << " gmres_reported=" << krylov_reported_before_certificate
                  << "\n";
    }

    for (std::size_t c = 0; c < nc; ++c) {
        U.component_data(0)[c] = x(c);
        U.component_data(1)[c] = x(nc + c);
        U.component_data(2)[c] = x(2*nc + c);
        p(c) = x(nv + c);
    }

    // Cellwise operator accounting for the exact monolithic candidate. This
    // is diagnostic-only: decompose A*x-b and the conservative Rhie-Chow
    // flux on the same U,p state, without changing acceptance or solver data.
    if (diagnostics.coupled_matrix_summary) {
        std::size_t row_cell = 0;
        double row_linf = -1.0;
        for (std::size_t cell = 0; cell < nc; ++cell) {
            if (!has_fixed_pressure && cell == reference_cell)
                continue;
            const std::size_t row = nv + cell;
            double ax = 0.0;
            for (std::size_t k = A.row_offsets_data()[row];
                 k < A.row_offsets_data()[row + 1]; ++k)
                ax += A.values_data()[k] * x(A.columns_data()[k]);
            const double residual = ax - b(row);
            if (std::abs(residual) > row_linf) {
                row_linf = std::abs(residual);
                row_cell = cell;
            }
        }

        std::array<std::vector<double>, 3> rAU;
        for (auto& component : rAU)
            component.resize(nc);
        for (std::size_t cell = 0; cell < nc; ++cell) {
            rAU[0][cell] = 1.0 / ex.diagonal[cell];
            rAU[1][cell] = 1.0 / ey.diagonal[cell];
            rAU[2][cell] = 1.0 / ez.diagonal[cell];
        }
        const auto grad_p = gauss_gradient_with_boundary(
            p, mesh, geometry, pressure_bcs,
            controls.pressure_skew_corrected_gradient,
            controls.pressure_projected_face_interpolation);
        const auto flux = make_rhie_chow_flux_from_physical_velocity(
            mesh, geometry, U, p, rAU, rho, velocity_bcs, pressure_bcs,
            controls.pressure_skew_corrected_gradient);
        const auto* cell_faces = mesh.cells().faces_data();
        const auto* cell_offsets = mesh.cells().offsets_data();

        auto report_cell = [&](std::size_t cell, const char* role) {
            const std::size_t row = nv + cell;
            double matrix_velocity = 0.0;
            double matrix_pressure = 0.0;
            for (std::size_t k = A.row_offsets_data()[row];
                 k < A.row_offsets_data()[row + 1]; ++k) {
                const std::size_t column = A.columns_data()[k];
                const double contribution = A.values_data()[k] * x(column);
                if (column < nv)
                    matrix_velocity += contribution;
                else
                    matrix_pressure += contribution;
            }
            const double matrix_residual = matrix_velocity + matrix_pressure - b(row);
            double flux_velocity = 0.0;
            double flux_pressure_orthogonal = 0.0;
            double flux_pressure_nonorthogonal = 0.0;
            double flux_pressure_compensation = 0.0;
            double flux_total = 0.0;

            for (Offset k = cell_offsets[cell]; k < cell_offsets[cell + 1]; ++k) {
                const std::size_t face = cell_faces[k];
                const double face_sign = mesh.ownership().owner(face) == cell ? 1.0 : -1.0;
                const std::size_t owner = mesh.ownership().owner(face);
                const auto neighbour = mesh.ownership().neighbour(face);
                const Vec3 sf = geometry.face_area_vectors[face];
                Vec3 uf{};
                U.get(owner, uf.x, uf.y, uf.z);
                if (neighbour >= 0) {
                    const std::size_t n = static_cast<std::size_t>(neighbour);
                    Vec3 un{};
                    U.get(n, un.x, un.y, un.z);
                    uf = (uf + un) * 0.5;
                    const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[owner];
                    const Vec3 dAU{
                        0.5 * (geometry.cell_volumes[owner] * rAU[0][owner] +
                               geometry.cell_volumes[n] * rAU[0][n]),
                        0.5 * (geometry.cell_volumes[owner] * rAU[1][owner] +
                               geometry.cell_volumes[n] * rAU[1][n]),
                        0.5 * (geometry.cell_volumes[owner] * rAU[2][owner] +
                               geometry.cell_volumes[n] * rAU[2][n])};
                    const Vec3 gp{
                        0.5 * (grad_p.component_data(0)[owner] + grad_p.component_data(0)[n]),
                        0.5 * (grad_p.component_data(1)[owner] + grad_p.component_data(1)[n]),
                        0.5 * (grad_p.component_data(2)[owner] + grad_p.component_data(2)[n])};
                    const double coefficient = pressure_face_orthogonal_coefficient(sf, dvec, dAU);
                    const double deferred = pressure_face_nonorthogonal_flux(sf, dvec, dAU, gp);
                    flux_velocity += face_sign * rho * uf.dot(sf);
                    flux_pressure_orthogonal -= face_sign * rho * coefficient * (p(n) - p(owner));
                    flux_pressure_nonorthogonal -= face_sign * rho * deferred;
                    flux_pressure_compensation += face_sign * rho *
                        (dAU.x * gp.x * sf.x + dAU.y * gp.y * sf.y + dAU.z * gp.z * sf.z);
                } else {
                    const std::size_t patch = geometry.face_patch[face];
                    if (patch < mesh.boundary().n_patches() &&
                        mesh.boundary().patch(patch).type == PatchType::EMPTY)
                        continue;
                    if (patch < mesh.boundary().n_patches()) {
                        const auto& name = mesh.boundary().patch(patch).name;
                        const auto ubc = velocity_bcs.find(name);
                        if (ubc != velocity_bcs.end() &&
                            ubc->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                            uf = ubc->second.value;
                        const auto pbc = pressure_bcs.find(name);
                        if (pbc != pressure_bcs.end() &&
                            pbc->second.type == ScalarBoundaryType::FIXED_VALUE) {
                            const Vec3 dvec = geometry.face_centres[face] -
                                              geometry.cell_centres[owner];
                            const Vec3 dAU{
                                geometry.cell_volumes[owner] * rAU[0][owner],
                                geometry.cell_volumes[owner] * rAU[1][owner],
                                geometry.cell_volumes[owner] * rAU[2][owner]};
                            const Vec3 gp{
                                grad_p.component_data(0)[owner],
                                grad_p.component_data(1)[owner],
                                grad_p.component_data(2)[owner]};
                            const double coefficient =
                                pressure_face_orthogonal_coefficient(sf, dvec, dAU);
                            const double deferred = pressure_face_nonorthogonal_flux(sf, dvec, dAU, gp);
                            flux_pressure_orthogonal -= face_sign * rho * coefficient *
                                (pbc->second.value - p(owner));
                            flux_pressure_nonorthogonal -= face_sign * rho * deferred;
                            flux_pressure_compensation += face_sign * rho *
                                (dAU.x * gp.x * sf.x + dAU.y * gp.y * sf.y + dAU.z * gp.z * sf.z);
                        }
                    }
                    flux_velocity += face_sign * rho * uf.dot(sf);
                }
                flux_total += face_sign * flux(face);
            }
            std::cerr << "COUPLED_CELL_OPERATOR_AUDIT role=" << role
                      << " cell=" << cell
                      << " matrix_velocity=" << matrix_velocity
                      << " matrix_pressure_orthogonal=" << matrix_pressure
                      << " matrix_rhs=" << b(row)
                      << " explicit_deferred_rhs=" << coupled_explicit_continuity_correction[cell]
                      << " rhs_without_deferred=" << (b(row) - coupled_explicit_continuity_correction[cell])
                      << " matrix_residual=" << matrix_residual
                      << " flux_velocity=" << flux_velocity
                      << " flux_pressure_orthogonal=" << flux_pressure_orthogonal
                      << " flux_pressure_nonorthogonal=" << flux_pressure_nonorthogonal
                      << " flux_pressure_compensation=" << flux_pressure_compensation
                      << " flux_sum_components=" << (flux_velocity + flux_pressure_orthogonal +
                            flux_pressure_nonorthogonal + flux_pressure_compensation)
                      << " flux_authoritative=" << flux_total << "\n";
        };

        double flux_linf = -1.0;
        std::size_t flux_cell = 0;
        for (std::size_t cell = 0; cell < nc; ++cell) {
            double divergence = 0.0;
            for (Offset k = cell_offsets[cell]; k < cell_offsets[cell + 1]; ++k) {
                const std::size_t face = cell_faces[k];
                divergence += mesh.ownership().owner(face) == cell ? flux(face) : -flux(face);
            }
            if (std::abs(divergence) > flux_linf) {
                flux_linf = std::abs(divergence);
                flux_cell = cell;
            }
        }
        report_cell(row_cell, "continuity_row_argmax");
        if (flux_cell != row_cell)
            report_cell(flux_cell, "flux_divergence_argmax");
        std::cerr << "COUPLED_CELL_OPERATOR_AUDIT_MAXIMA"
                  << " continuity_row_cell=" << row_cell
                  << " continuity_row_linf=" << row_linf
                  << " flux_cell=" << flux_cell
                  << " flux_divergence_linf=" << flux_linf << "\n";
    }
    return result;
}

inline IncompressibleSolveResult solve_steady_incompressible(
    const cfdx::core::Mesh& mesh,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    cfdx::core::Field<double, cfdx::core::Location::CELL>& p,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    const IncompressibleSolverControls& controls = {},
    const std::string& restart_path = {},
    cfdx::io::DatRestartFields restart_fields = {})
{
    using namespace cfdx::core;
    if (U.dimension() != 3 || U.size() != mesh.n_cells() ||
        p.dimension() != 1 || p.size() != mesh.n_cells())
        throw std::invalid_argument("solve_steady_incompressible: invalid fields");

    // Validate all controls and derived material properties before initialization
    // can mutate the caller-provided solution fields.
    validate_incompressible_controls(controls, mesh.n_cells());
    // Do not let a phantom fixed-pressure BC suppress the pressure gauge.
    // Validate this before Uniform/Restart initialization changes U or p.
    const bool has_fixed_pressure_boundary =
        has_effective_fixed_pressure_boundary(mesh, pressure_bcs);
    // Validate the caller-owned SST bundle before initialization can mutate
    // the carrier fields. Derived SST quantities must be refreshed only after
    // initialization/restart has established the actual velocity field.
    const FvGeometry geometry = build_fv_geometry(mesh);
    if (!controls.inlet_outlet_velocity_patches.empty()) {
        if (controls.tvd_deferred_memory.enabled ||
            controls.tvd_secant_implicit.enabled ||
            controls.momentum_same_state_face_callback ||
            controls.momentum_picard_tvd_face_callback)
            throw std::invalid_argument(
                "mixed inlet/outlet is not qualified with legacy TVD defect or same-state face audit");

        if (controls.algorithm == PressureVelocityAlgorithm::COUPLED ||
            controls.sst_runtime || controls.transient.enabled ||
            controls.pressure_continuity_line_search.enabled ||
            controls.v18_sparse_flux_pressure_jacobian ||
            controls.v18_full_flux_pressure_oracle ||
            controls.pressure_projected_face_interpolation ||
            controls.convection_continuation.enabled)
            throw std::invalid_argument(
                "inlet/outlet mixed velocity is currently segregated steady laminar SIMPLE/SIMPLEC only");
        for (const auto& [name, velocity] : controls.inlet_outlet_velocity_patches) {
            const auto it = pressure_bcs.find(name);
            if (it == pressure_bcs.end() || it->second.type != ScalarBoundaryType::FIXED_VALUE)
                throw std::invalid_argument(
                    "inlet/outlet mixed velocity requires fixed pressure at outlet: " + name);
        }
    }
    auto mixed_velocity_faces = [&](const Field<double,Location::FACE>& flux)
        -> std::optional<VelocityFaceConditions> {
        if (controls.inlet_outlet_velocity_patches.empty()) return std::nullopt;
        return resolve_inlet_outlet_velocity_faces(
            mesh, geometry, flux, velocity_bcs,
            controls.inlet_outlet_velocity_patches);
    };
    if (controls.sst_runtime) {
        validate_sst_production_state(*controls.sst_runtime->state, mesh.n_cells());
        const auto& runtime = *controls.sst_runtime;
        const SSTRuntimeBoundaryContract boundary_contract{
            runtime.k_boundary_conditions, runtime.omega_boundary_conditions,
            runtime.wall_omega_sst_patches};
        require_supported_sst_runtime_boundaries(
            boundary_contract, mesh, geometry, &runtime.omega_face_values);
    }

    // In an explicitly opted-in SST restart, the DAT checkpoint is the
    // source of truth for k and omega just as it is for U and p. Wire those
    // named fields into the same atomic reader instead of requiring callers
    // to remember a separate restore helper. Preserve unrelated optional
    // fields (e.g. temperature) supplied by the caller.
    auto effective_restart_fields = restart_fields;
    // DAT6 must restore a single atomic U,p,phi,p_low snapshot. The reader
    // validates the ordered face/cell graph, finite values, required fields
    // and rejects partial/truncated checkpoints without touching caller state.
    cfdx::core::Field<double,cfdx::core::Location::CELL> restored_pressure_low;
    cfdx::core::Field<double,cfdx::core::Location::FACE> restored_face_flux;
    if (controls.pressure_high_low_production &&
        controls.initialization.mode == InitializationMode::Restart) {
        if (effective_restart_fields.pressure_low ||
            effective_restart_fields.face_mass_flux)
            throw std::invalid_argument(
                "V7 split pressure restart owns p_low and phi restore targets");
        restored_pressure_low=cfdx::core::Field<double,cfdx::core::Location::CELL>(
            mesh.n_cells(),"p_low","Pa",1);
        restored_face_flux=cfdx::core::Field<double,cfdx::core::Location::FACE>(
            mesh.n_faces(),"mdot","kg/s",1);
        effective_restart_fields.pressure_low=&restored_pressure_low;
        effective_restart_fields.face_mass_flux=&restored_face_flux;
    }
    if (controls.sst_runtime &&
        controls.initialization.mode == InitializationMode::Restart) {
        auto& sst_state = *controls.sst_runtime->state;
        if (effective_restart_fields.k && effective_restart_fields.k != &sst_state.k)
            throw std::invalid_argument(
                "SST DAT restart k target must be the active runtime state");
        if (effective_restart_fields.omega &&
            effective_restart_fields.omega != &sst_state.omega)
            throw std::invalid_argument(
                "SST DAT restart omega target must be the active runtime state");
        effective_restart_fields.k = &sst_state.k;
        effective_restart_fields.omega = &sst_state.omega;
    }

    apply_initialization(
        controls.initialization, mesh, U, p, effective_restart_fields, restart_path);
    if (controls.initialized_state_callback)
        controls.initialized_state_callback(mesh, U, p);

    // Authoritative dynamic viscosity used by predictor assembly, nonlinear
    // acceptance residuals, and the independent momentum conservation audit.
    std::vector<double> mu_eff(mesh.n_cells());
    // The legacy uniform-viscosity incompressible path uses mu*Laplacian(U).
    // Cellwise/SST viscosity requires the full symmetric stress; its transpose
    // term is conservative and explicit in each nonlinear predictor.
    const bool include_transpose_stress = controls.sst_runtime != nullptr ||
        !controls.cell_turbulent_viscosity.empty();
    if (controls.sst_runtime) {
        std::map<std::string, std::array<double, 3>> fixed_velocity_values;
        for (const auto& [patch, bc] : velocity_bcs) {
            if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                fixed_velocity_values.emplace(
                    patch, std::array<double, 3>{bc.value.x, bc.value.y, bc.value.z});
        }
        // Initialization may replace U (uniform, profile, or DAT restart).
        // Refresh wall distance, strain and nu_t against that final field so
        // the very first momentum predictor cannot use stale SST viscosity.
        refresh_sst_production_derived_state_from_mesh(
            *controls.sst_runtime->state, mesh, controls.sst_runtime->transport,
            U, fixed_velocity_values, controls.sst_runtime->wall_patch_names);
    }
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const double nu_t = controls.sst_runtime
            ? controls.sst_runtime->state->nu_t(cell)
            : (controls.cell_turbulent_viscosity.empty()
                ? controls.turbulent_viscosity
                : controls.cell_turbulent_viscosity[cell]);
        mu_eff[cell] = controls.density * (controls.kinematic_viscosity + nu_t);
        if (!(mu_eff[cell] >= 0.0) || !std::isfinite(mu_eff[cell]))
            throw std::invalid_argument(
                "effective cellwise dynamic viscosity must be finite and non-negative");
    }

    cfdx::core::ConvergenceMonitorControls monitor_controls;
    // This monitor receives a dimensionless maximum of physical, RHS-normalized
    // residuals and normalized field changes. During an operator homotopy the
    // initial reference is not the same equation as the final one, so its
    // relative reduction is only a diagnostic. For the opt-in continuation,
    // use the ORIGINAL per-criterion relative threshold as the absolute bound
    // for that already-normalized synthetic metric. The independent momentum
    // component, pressure, continuity and field-change gates remain mandatory.
    monitor_controls.residual.absolute_tolerance = convection_monitor_residual_threshold(
        controls.convection_continuation, controls.convergence.absolute_tolerance,
        controls.convergence.relative_tolerance);
    monitor_controls.residual.relative_tolerance = controls.convergence.relative_tolerance;
    monitor_controls.conservation_tolerance = controls.convergence.continuity_tolerance;
    monitor_controls.minimum_iterations = std::min<std::size_t>(2, controls.convergence.max_iterations);
    monitor_controls.max_iterations = controls.convergence.max_iterations;
    // Conservative defaults: stagnation is diagnostic, not a shortcut around
    // the configured residual/conservation gates. Divergence is detected only
    // after repeated growth and never by relaxing a tolerance.
    monitor_controls.stagnation_window = 10;
    monitor_controls.stagnation_relative_improvement =
        controls.convergence.stagnation_relative_improvement;
    monitor_controls.divergence_window = 3;
    monitor_controls.divergence_growth_factor = 10.0;
    monitor_controls.qoi_gates = controls.qoi_gates;
    cfdx::core::ConvergenceMonitor convergence_monitor(monitor_controls);
    if (controls.initial_face_mass_flux_from_paired_checkpoint)
        convergence_monitor.restore_initial_reference(
            controls.initial_nonlinear_metric_reference);

    const std::size_t configured_pcorr =
        (controls.algorithm == PressureVelocityAlgorithm::PISO ||
         controls.algorithm == PressureVelocityAlgorithm::PIMPLE)
            ? static_cast<std::size_t>(controls.coupling.n_pressure_correctors)
            : (controls.algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP
                ? static_cast<std::size_t>(controls.coupling.n_fractional_steps)
                : 1u);

    const std::size_t minimum_outer_correctors =
        controls.algorithm == PressureVelocityAlgorithm::PIMPLE
            ? static_cast<std::size_t>(controls.coupling.n_outer_correctors)
            : 1u;

    IncompressibleSolveResult result;
    result.reference_momentum_residual = controls.initial_momentum_residual_reference;
    double v17_pseudo_cfl=controls.pseudo_momentum.initial_cfl;
    result.requested_coupled_schur_model = controls.coupling.schur_model;
    result.resolved_coupled_schur_model = controls.coupling.schur_model;
    result.coupled_schur_model_resolved = true;
    LinearSolverPlan pressure_plan;
    std::unique_ptr<Preconditioner> pressure_preconditioner;
    std::unique_ptr<ReusableCgContext> pressure_context;
    std::unique_ptr<NullSpaceProjector> pressure_null_space;
    const bool use_projected_pressure_null_space =
        !has_fixed_pressure_boundary &&
        controls.pressure_linear_solver.null_space ==
            NullSpaceModel::Constant;
    const bool reduce_pressure_gauge =
        !has_fixed_pressure_boundary &&
        !use_projected_pressure_null_space &&
        mesh.n_cells() > 1;
    if (controls.algorithm != PressureVelocityAlgorithm::COUPLED) {
        pressure_plan = select_linear_solver(
            LinearProblemKind::PressurePoisson,
            reduce_pressure_gauge ? mesh.n_cells() - 1 : mesh.n_cells(),
            controls.pressure_linear_solver);
        result.pressure_linear_plan = pressure_plan;
        result.pressure_linear_plan_resolved = true;
        // The reusable context is intentionally dimension-agnostic until its
        // first prepare(matrix). For reference-cell gauge reduction that first
        // matrix has N-1 rows, so the preconditioner must not be configured to
        // project a constant null space: the gauge has already removed it.
        // Project only on the full-size, explicitly projected-null-space path.
        // Reuse AMG across pressure correctors even when symmetric algebraic
        // scaling is active; the dispatcher below transforms the matrix and
        // verifies the physical residual after unscaling the pressure field.
        if (pressure_plan.krylov == KrylovModel::CG) {
            const bool project_pressure_null_space =
                pressure_plan.null_space == NullSpaceModel::Constant &&
                !reduce_pressure_gauge;
            pressure_preconditioner =
                make_scalar_preconditioner(
                    pressure_plan.preconditioner,
                    project_pressure_null_space);
            pressure_context = std::make_unique<ReusableCgContext>(
                *pressure_preconditioner);
            if (project_pressure_null_space) {
                pressure_null_space = std::make_unique<NullSpaceProjector>(
                    NullSpaceProjector::constant(mesh.n_cells()));
            }
        }
    }

    auto build_hbya = [&](const ScalarEquation& eq,
                          const Vector& field,
                          const Field<double, Location::CELL>& gradp,
                          std::size_t component,
                          const std::vector<double>& rAU) {
        std::vector<double> hbya(mesh.n_cells(), 0.0);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double h = eq.rhs(c) + gradp.component_data(component)[c] * geometry.cell_volumes[c];
            const auto begin = eq.matrix.row_offsets_data()[c];
            const auto end = eq.matrix.row_offsets_data()[c + 1];
            for (std::uint32_t k = begin; k < end; ++k) {
                const auto col = eq.matrix.columns_data()[k];
                if (col != c)
                    h -= eq.matrix.values_data()[k] * field(col);
            }
            // HbyA is always the plain momentum-inverse combination
            // rAU*H, for SIMPLEC included. The SIMPLEC consistent
            // operator (rAtU) belongs to the pressure equation and its
            // face fluxes only: the over-relaxation terms added to
            // phiHbyA and to the corrected flux are sized so that
            // (dAtU - dAU) cancels the rAtU response of the pressure
            // operator, which collapses the SIMPLEC fixed point back
            // onto the plain Rhie-Chow state built with rAU. Folding
            // (rAtU - rAU)*grad(p) into the cell HbyA instead leaves
            // the converged state on a displaced continuity operator
            // carrying the spurious source
            // div(rho*(dAtU - dAU)*snGrad(p)*|Sf|) of O(Laplacian(p))
            // scale: invisible on Couette (constant p) and Poiseuille
            // (linear p) but far outside the algorithm-equivalence
            // gate on the lid-driven cavity.
            hbya[c] = rAU[c] * h;
            if (!std::isfinite(hbya[c]))
                throw std::runtime_error("solve_steady_incompressible: non-finite HbyA");
        }
        return hbya;
    };

    auto make_hbya_field = [&](const std::array<std::vector<double>,3>& hbya) {
        Field<double, Location::CELL> field(
            mesh.n_cells(), "HbyA", "m/s", 3);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            field.component_data(0)[c] = hbya[0][c];
            field.component_data(1)[c] = hbya[1][c];
            field.component_data(2)[c] = hbya[2][c];
        }
        return field;
    };

    // Face pressure projection coefficient.
    //
    // Momentum diagonals are integrated finite-volume coefficients A_P [N s/m].
    // rAU = 1/A_P therefore has units of m/(N s), while the velocity response
    // to a pressure gradient is dAU = V*rAU.  The same dAU operator must be
    // used by cell velocity reconstruction, Rhie-Chow face fluxes and the
    // pressure-correction matrix.
    auto directional_face_coefficient =
        [&](const std::array<std::vector<double>,3>& coeff,
            std::size_t f, std::size_t o, std::size_t n,
            const Vec3& Sf, const Vec3& dvec) {
            const double wn = pressure_face_neighbour_weight(
                geometry, f, o, n,
                controls.pressure_projected_face_interpolation);
            const double wo = 1.0 - wn;
            const double cx = wo * geometry.cell_volumes[o] * coeff[0][o] +
                              wn * geometry.cell_volumes[n] * coeff[0][n];
            const double cy = wo * geometry.cell_volumes[o] * coeff[1][o] +
                              wn * geometry.cell_volumes[n] * coeff[1][n];
            const double cz = wo * geometry.cell_volumes[o] * coeff[2][o] +
                              wn * geometry.cell_volumes[n] * coeff[2][n];
            return pressure_face_orthogonal_coefficient(Sf, dvec, {cx, cy, cz});
        };

    // Sidecar pressure low part: committed only with the nonlinear physical
    // state. A failed iteration MUST roll this back alongside U/p/phi.
    Vector pressure_low_state(
        controls.pressure_high_low_production ? mesh.n_cells() : 0,0.0);
    if (controls.pressure_high_low_production &&
        controls.initialization.mode == InitializationMode::Restart) {
        for (std::size_t cell=0; cell<mesh.n_cells(); ++cell)
            pressure_low_state(cell)=restored_pressure_low(cell);
    }
    Field<double, Location::FACE> mass_flux =
        make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs,
                       {}, controls.pressure_projected_face_interpolation);
    if (!controls.initial_face_mass_flux.empty()) {
        if (controls.initialization.mode == InitializationMode::Restart ||
            controls.initial_face_mass_flux.size() != mesh.n_faces())
            throw std::invalid_argument(
                "hybrid initial face mass flux requires current-mesh, non-restart initialization");
        // Retain the original open-flow boundary normalization; only closed
        // cavities fall back to the internal face throughput (N7).
        double prescribed_boundary_throughput = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            const double mdot = controls.initial_face_mass_flux[f];
            if (!std::isfinite(mdot))
                throw std::invalid_argument(
                    "hybrid initial face mass flux contains nonfinite values");
            if (mesh.ownership().neighbour(f) >= 0) continue;
            const auto patch_index = geometry.face_patch[f];
            if (patch_index >= mesh.boundary().n_patches())
                throw std::invalid_argument(
                    "hybrid flux has an unlabelled boundary face");
            const auto& name = mesh.boundary().patch(patch_index).name;
            const auto it = velocity_bcs.find(name);
            if (it != velocity_bcs.end() &&
                it->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE) {
                const double reference = mass_flux(f);
                const double allowed = 1e-11 * std::max(1.0, std::abs(reference));
                // For a pressure-Dirichlet face the production RC operator
                // adds an implicit + nonorthogonal pressure contribution even
                // with fixed U. The accepted flux is not equal to U_bc dot Sf.
                // Only a verified paired restart may replay that difference;
                // legacy hybrid flux input keeps the stronger existing rule.
                const auto pbc = pressure_bcs.find(name);
                const bool pressure_flux_boundary =
                    pbc != pressure_bcs.end() &&
                    pbc->second.type == ScalarBoundaryType::FIXED_VALUE;
                if (std::abs(mdot - reference) > allowed &&
                    !(controls.initial_face_mass_flux_from_paired_checkpoint &&
                      pressure_flux_boundary))
                    throw std::invalid_argument(
                        "hybrid flux violates a fixed-velocity boundary face");
                prescribed_boundary_throughput += std::abs(reference);
            }
        }
        const auto continuity = initial_face_flux_continuity_report(
            mesh, controls.initial_face_mass_flux,
            prescribed_boundary_throughput);
        // A checkpoint represents an accepted iterate, not necessarily a
        // converged one. Preserve its exact phi and separately audit continuity
        // during the solve; reject pathological values, never relax the final
        // nonlinear convergence or independent face-momentum gates.
        const double admissible =
            controls.initial_face_mass_flux_from_paired_checkpoint ? 1e-1 : 1e-7;
        if (!(continuity.normalized <= admissible))
            throw std::invalid_argument(
                "hybrid initial mass flux fails local continuity admission gate");
        for (std::size_t f=0; f<mesh.n_faces(); ++f)
            mass_flux(f)=controls.initial_face_mass_flux[f];
        std::cout << "HYBRID_INITIAL_FACE_FLUX continuity_l1=" << continuity.local_l1
                  << " normalized=" << continuity.normalized << " (NOT a final convergence gate)\n";
    }

    // A restart must replay the continuous iteration sequence. Iteration n+1
    // assembles momentum from the Rhie-Chow flux reconstructed with
    // iteration n's rAU; that flux is not re-derivable from the restored U and
    // p, so a checkpoint that carries it (DatRestartFields::face_mass_flux)
    // replaces the simple-interpolation initial flux. Without this, a
    // restart produces a different second iteration than the continuous run.
    if (controls.initialization.mode == InitializationMode::Restart &&
        effective_restart_fields.face_mass_flux) {
        const auto& restored_flux = *effective_restart_fields.face_mass_flux;
        if (restored_flux.dimension() != 1 || restored_flux.size() != mesh.n_faces())
            throw std::invalid_argument(
                "solve_steady_incompressible: restart face mass flux does not match mesh");
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            if (!std::isfinite(restored_flux(f)))
                throw std::invalid_argument(
                    "solve_steady_incompressible: restart face mass flux has non-finite value");
        mass_flux = restored_flux;
    }

    const bool freeze_state = controls.diagnostics.freeze_state_probe;
    Field<double, Location::FACE> phi_used_in_momentum;
    std::array<std::vector<double>, 3> frozen_rAU;
    std::array<std::vector<double>, 3> frozen_hbyA;
    Field<double, Location::CELL> frozen_grad_p;
    bool frozen_state_valid = false;

    // All multi-secant samples are committed only after an accepted full
    // SIMPLE iteration, never from a rejected nonlinear trial. The predictor
    // vector has [Ux,Uy,Uz,p,phi] in that order; phi is an INDEPENDENT,
    // conservative state variable and is not silently replaced by interp(U).
    NonlinearMultisecantHistory multisecant_history(controls.nonlinear_multisecant);
    const std::size_t nc_multisecant = mesh.n_cells();
    const std::size_t nf_multisecant = mesh.n_faces();
    const auto pack_nonlinear_state = [&](const auto& velocity, const auto& pressure,
                                           const auto& face_flux) {
        std::vector<double> flat(4 * nc_multisecant + nf_multisecant);
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t cell = 0; cell < nc_multisecant; ++cell)
                flat[d*nc_multisecant+cell]=velocity.component_data(d)[cell];
        for (std::size_t cell = 0; cell < nc_multisecant; ++cell)
            flat[3*nc_multisecant+cell]=pressure(cell);
        for (std::size_t face = 0; face < nf_multisecant; ++face)
            flat[4*nc_multisecant+face]=face_flux(face);
        return flat;
    };
    const auto restore_nonlinear_state = [&](const std::vector<double>& flat) {
        if (flat.size() != 4 * nc_multisecant + nf_multisecant)
            throw std::invalid_argument("multisecant: incompatible state size");
        for (const double x : flat)
            if (!std::isfinite(x))
                throw std::invalid_argument("multisecant: nonfinite proposed field");
        for (std::size_t d = 0; d < 3; ++d)
            for (std::size_t cell = 0; cell < nc_multisecant; ++cell)
                U.component_data(d)[cell]=flat[d*nc_multisecant+cell];
        for (std::size_t cell = 0; cell < nc_multisecant; ++cell)
            p(cell)=flat[3*nc_multisecant+cell];
        for (std::size_t face = 0; face < nf_multisecant; ++face)
            mass_flux(face)=flat[4*nc_multisecant+face];
    };
    // Freeze *physical* normalization once, so NGMRES cannot obtain a false
    // reduction by changing the RHS scale from one iteration to the next.
    std::array<double,3> multisecant_rhs_reference{1.0,1.0,1.0};
    double multisecant_flux_reference=1.0;
    bool multisecant_scales_valid=false;

    GlobalNonlinearContractionLedgerV14 committed_global_merit_v14;
    CoupledPhaseEscapeLedgerV45 coupled_phase_escape_v45;
    const bool v43_field_probe_enabled = controls.sst_runtime &&
        (controls.coupled_period_four_mean_v43_experimental ||
         controls.coupled_period_four_full_state_v44_experimental ||
         (std::getenv("CFDX_OUTER_FIELDS_V43") &&
          std::string(std::getenv("CFDX_OUTER_FIELDS_V43"))=="1"));
    std::unique_ptr<OuterFieldPeriodProbeV43> v43_field_probe;
    if (v43_field_probe_enabled)
        v43_field_probe=std::make_unique<OuterFieldPeriodProbeV43>(
            mesh.n_cells(),mesh.n_faces());
    double effective_alpha_u = controls.coupling.alpha_u;
    double effective_alpha_p = controls.coupling.alpha_p;
    std::size_t adaptive_u_improvement_streak = 0;
    std::size_t adaptive_u_degradation_streak = 0;
    std::size_t adaptive_p_improvement_streak = 0;
    std::size_t adaptive_p_degradation_streak = 0;
    std::size_t adaptive_u_recovery_cooldown = 0;
    std::size_t adaptive_p_recovery_cooldown = 0;

    // Only the last COMMITTED effective TVD source is cached. This is an
    // exponentially filtered nonlinear fixed-point variable, not an extra
    // source added to the final accepted-state physics. Retrying the same nonlinear
    // iteration must never reuse a deferred flux from a rejected attempt.
    std::array<std::vector<double>,3> committed_tvd_deferred;
    std::array<std::vector<double>,3> committed_tvd_raw;
    std::array<std::vector<double>,3> committed_tvd_two_steps_ago;
    bool committed_tvd_two_steps_valid = false;
    std::array<std::vector<double>,3> committed_tvd_input_velocity;
    bool committed_tvd_deferred_valid = false;

    // Sample only committed physical states, never rejected Newton/Picard trials.
    const bool cell_trace_enabled = !controls.diagnostics.cell_convergence_trace_path.empty();
    std::unique_ptr<CellConvergenceTracker> cell_tracker;
    std::vector<std::array<double, 3>> cell_trace_momentum;
    std::vector<double> cell_trace_continuity;
    std::vector<std::array<double, 3>> cell_trace_centers;
    if (cell_trace_enabled) {
        if (mesh.n_cells() > 5'000'000)
            throw std::invalid_argument("cell convergence diagnostics limited to 5 million cells");
        cell_trace_momentum.resize(mesh.n_cells());
        cell_trace_continuity.resize(mesh.n_cells());
        cell_trace_centers.resize(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            cell_trace_centers[c] = {geometry.cell_centres[c].x,
                geometry.cell_centres[c].y, geometry.cell_centres[c].z};
        cell_tracker = std::make_unique<CellConvergenceTracker>(
            controls.diagnostics.cell_convergence_trace_path,
            controls.diagnostics.cell_convergence_mesh_source,
            mesh.n_cells(), mesh.n_faces(),
            controls.diagnostics.cell_convergence_every, 24,
            controls.convergence.relative_tolerance,
            controls.diagnostics.cell_convergence_full_map_path);
    }
    std::array<double, 3> pending_cell_rhs_scales{1.0, 1.0, 1.0};
    std::size_t last_committed_cell_iteration = 0;

    for (std::size_t iter = 1; iter <= controls.convergence.max_iterations; ++iter) {
        const double outer_high_order_fraction =
            convection_high_order_fraction(controls.convection_continuation, iter);
        const bool full_order_eligible =
            convection_full_order_eligible(controls.convection_continuation, iter);
        const bool same_convection_operator_as_previous =
            convection_residuals_share_operator(
                controls.convection_continuation, iter);
        bool stop_after_iteration = false;
        NonlinearRetryController retry_controller(controls.nonlinear_retry);
        NonlinearStateRollback transaction(U, p);
        const auto mass_flux_snapshot = mass_flux;
        const Vector pressure_low_snapshot=pressure_low_state;
        bool iteration_completed = false;
        // A rejected multisecant predictor gets exactly one ordinary Picard
        // retry at UNCHANGED alpha_u/alpha_p. Reset the multisecant history;
        // a rejected prediction must not seed the next nonlinear proposal.
        bool picard_fallback_after_multisecant = false;
        while (!iteration_completed) {
            // A retry must restore not only physical fields but also any
            // convergence/history state tentatively updated before the outer
            // transaction commits.
            stop_after_iteration = false;
            transaction.begin();
            std::optional<MomentumSameStateTrace> pending_same_state_trace;
            std::optional<MomentumPicardTvdTrace> pending_picard_tvd_trace;
            std::array<std::vector<double>,3> pending_effective_tvd_deferred;
            std::array<std::vector<double>,3> pending_raw_tvd_deferred;
            std::array<std::vector<double>,3> pending_tvd_input_velocity;
            std::array<double,3> pending_memory_source_linf{0.0,0.0,0.0};
            std::array<std::optional<MomentumLinearRelaxationEvidence>,3>
                pending_linear_relaxation;
            const double pending_linear_tolerance = controls.linear_tolerance;
            double pending_alpha_u = controls.coupling.alpha_u;
            const auto convergence_monitor_snapshot = convergence_monitor;
            const auto result_history_size = result.history.size();
            const auto probe_samples_size = result.probe_samples.size();
            const auto previous_result_converged = result.converged;
            const auto previous_result_status = result.convergence_status;
            const auto previous_result_reason = result.convergence_reason;
            const auto previous_result_iterations = result.iterations;
            const auto v30_ptc_attempt_snapshot = controls.sst_runtime
                ? std::optional<SSTV30ContinuationState>(
                    controls.sst_runtime->coupled_ptc_v30) : std::nullopt;
            const auto sst_attempt_snapshot = controls.sst_runtime
                ? std::optional<SSTProductionState>(
                    snapshot_sst_production_state(*controls.sst_runtime->state))
                : std::nullopt;
            // Visible to the exception handlers: pressure/linear failures may
            // occur before a candidate reaches the physical-residual gate.
            bool multisecant_trial_used = false;
            bool coupled_linear_failed_v37 = false;
            try {
        // Reversible V43 trial initializer: a full four-field cycle can have
        // very small ||x_n-x_{n-4}|| while failing every physical tolerance.
        // A mean of the four committed phases is a candidate for breaking
        // that fixed-point orbit, not a solution or a new residual oracle.
        // The complete COUPLED/SST solve follows, so its final Rhie-Chow
        // face flux and returned-state FV defects are independently tested.
        if ((controls.coupled_period_four_mean_v43_experimental ||
             controls.coupled_period_four_full_state_v44_experimental) &&
            v43_field_probe && v43_field_probe->period_four_candidate()) {
            const double gauge=p(controls.pressure_reference_cell);
            for (std::size_t cell=0;cell<mesh.n_cells();++cell) {
                for (std::size_t component=0;component<3;++component)
                    U.component_data(component)[cell]=
                        v43_field_probe->average_last_four(component,cell);
                p(cell)=gauge+v43_field_probe->average_last_four(3,cell);
            }
            for (std::size_t face=0;face<mesh.n_faces();++face)
                mass_flux(face)=v43_field_probe->average_last_four(4,face);
            if (controls.coupled_period_four_full_state_v44_experimental) {
                // Do not average nu_t: nu_t(k,omega,|S|,F2) is nonlinear.
                // Averaging U/p/phi alone (V43) also leaves the frozen
                // momentum viscosity inconsistent with its trial SST fields.
                // Set k/omega, then reconstruct ALL derived SST fields and
                // the viscosity actually consumed by COUPLED assembly.
                auto& state=*controls.sst_runtime->state;
                for (std::size_t cell=0;cell<mesh.n_cells();++cell) {
                    const auto phase=make_coupled_sst_phase_candidate_v44(
                        v43_field_probe->average_last_four(5,cell),
                        v43_field_probe->average_last_four(6,cell));
                    state.k(cell)=phase.k;
                    state.omega(cell)=phase.omega;
                }
                std::map<std::string,std::array<double,3>> fixed_velocity_values;
                for (const auto& pair:velocity_bcs)
                    if (pair.second.type==VelocityBoundaryCondition::Type::FIXED_VALUE)
                        fixed_velocity_values.emplace(pair.first,
                            std::array<double,3>{pair.second.value.x,
                                pair.second.value.y,pair.second.value.z});
                refresh_sst_production_derived_state_from_mesh(
                    state,mesh,controls.sst_runtime->transport,
                    U,fixed_velocity_values,controls.sst_runtime->wall_patch_names);
                for (std::size_t cell=0;cell<mesh.n_cells();++cell) {
                    mu_eff[cell]=controls.density*(
                        controls.kinematic_viscosity+state.nu_t(cell));
                    if (!(mu_eff[cell]>0.0) || !std::isfinite(mu_eff[cell]))
                        throw std::runtime_error("V44 full-phase invalid viscosity");
                }
                if (controls.diagnostics.iteration_trace)
                    std::cerr<<"COUPLED_PERIOD4_FULL_V44 iteration="<<iter
                             <<" retry="<<retry_controller.retries()
                             <<" k_omega_averaged=1 nut_recomputed=1"
                             <<" strain_recomputed=1 pressure_gauge_preserved=1"
                             <<" committed=0 physical_qualified=0\n";
            }
            if (controls.diagnostics.iteration_trace)
                std::cerr<<"COUPLED_PERIOD4_PREDICTOR_V43 iteration="<<iter
                         <<" retry="<<retry_controller.retries()
                         <<" mean_phases=4 predicted_accepted=0"
                         <<" true_residual_must_be_reassembled=1"
                         <<" final_flux_rederived_by_coupled=1\n";
        }
        MultisecantProposal multisecant_proposal;
        if (controls.nonlinear_multisecant.method != NonlinearMultisecantMethod::OFF &&
            iter >= controls.nonlinear_multisecant.start_iteration &&
            !picard_fallback_after_multisecant &&
            retry_controller.retries() == 0 && full_order_eligible &&
            same_convection_operator_as_previous) {
            multisecant_proposal = multisecant_history.propose();
            if (multisecant_proposal.available) {
                // The extrapolated state is a *predictor*, not an accepted
                // Rhie-Chow state. Every attempted SIMPLE step completely
                // rebuilds H/A, pressure and its conservative face flux.
                // Check the discrete divergence before using any blended
                // flux; each previously committed member was conservative.
                std::vector<double> blended_divergence(nc_multisecant,0.0);
                for (std::size_t f=0;f<nf_multisecant;++f) {
                    const double phi = multisecant_proposal.state[4*nc_multisecant+f];
                    const auto o=mesh.ownership().owner(f);
                    blended_divergence[o]+=phi;
                    const auto n=mesh.ownership().neighbour(f);
                    if (n>=0) blended_divergence[static_cast<std::size_t>(n)]-=phi;
                }
                double div_l1=0.0, phi_l1=0.0;
                for (const double v:blended_divergence) div_l1+=std::abs(v);
                for (std::size_t f=0;f<nf_multisecant;++f)
                    phi_l1+=std::abs(multisecant_proposal.state[4*nc_multisecant+f]);
                const double prev_cont = result.history.empty() ? 0.0 :
                    result.history.back().continuity_normalized;
                if (std::isfinite(div_l1) && std::isfinite(phi_l1) &&
                    phi_l1>0.0 && div_l1/phi_l1 <=
                       std::max(10.0*controls.convergence.continuity_tolerance,
                                2.0*prev_cont)) {
                    restore_nonlinear_state(multisecant_proposal.state);
                    multisecant_trial_used=true;
                } else {
                    multisecant_proposal.reason =
                        "candidate conservative flux divergence exceeds previous state";
                }
            }
        }
        const std::vector<double> multisecant_input =
            controls.nonlinear_multisecant.method == NonlinearMultisecantMethod::OFF
            ? std::vector<double>{} : pack_nonlinear_state(U,p,mass_flux);
        double linear_tolerance = controls.linear_tolerance;
        if (controls.acceleration.adaptive_linear_tolerance) {
            const auto& forcing = controls.acceleration.linear_forcing;
            if (result.history.empty()) {
                linear_tolerance = forcing.eta_max;
            } else {
                const auto nonlinear_residual = std::max(
                    result.history.back().momentum_residual,
                    result.history.back().continuity_normalized);
                const auto previous_residual = result.history.size() > 1
                    ? std::max(result.history[result.history.size()-2].momentum_residual,
                               result.history[result.history.size()-2].continuity_normalized)
                    : std::max(nonlinear_residual, forcing.residual_floor);
                linear_tolerance = nonlinear_forcing_tolerance(
                    nonlinear_residual,
                    std::max(previous_residual, forcing.residual_floor), forcing);
            }
        }

        if (controls.adaptive_relaxation.enabled) {
            // Adapt the momentum and pressure relaxation independently. SIMPLE
            // couples the two equations, but their residuals measure different
            // numerical mechanisms: momentum relaxation is controlled from the
            // momentum-equation residual, while pressure relaxation is
            // controlled from the normalized continuity imbalance. A single
            // aggregate metric must not force both factors in the same direction.
            const auto window = controls.adaptive_relaxation.adaptation_window;
            if (window > 0 && result.history.size() > window &&
                (result.history.size() - 1) % window == 0) {
                const auto& current = result.history.back();
                const auto& reference = result.history[result.history.size() - 1 - window];

                if (adaptive_u_recovery_cooldown > 0) {
                    --adaptive_u_recovery_cooldown;
                    adaptive_u_improvement_streak = 0;
                    adaptive_u_degradation_streak = 0;
                } else if (std::isfinite(current.momentum_equation_residual_relative) &&
                           std::isfinite(reference.momentum_equation_residual_relative) &&
                           reference.momentum_equation_residual_relative > 0.0) {
                    effective_alpha_u = adapt_relaxation_factor_windowed(
                        effective_alpha_u,
                        reference.momentum_equation_residual_relative,
                        current.momentum_equation_residual_relative,
                        controls.adaptive_relaxation.min_alpha_u,
                        controls.adaptive_relaxation.max_alpha_u,
                        controls.adaptive_relaxation,
                        adaptive_u_improvement_streak,
                        adaptive_u_degradation_streak);
                } else {
                    adaptive_u_improvement_streak = 0;
                    adaptive_u_degradation_streak = 0;
                }

                if (adaptive_p_recovery_cooldown > 0) {
                    --adaptive_p_recovery_cooldown;
                    adaptive_p_improvement_streak = 0;
                    adaptive_p_degradation_streak = 0;
                } else {
                    // V19's exact flux projection makes continuity an
                    // arithmetic-noise metric. Adapt the PRESSURE split from
                    // the nonlinear momentum/RC coupling merit instead.
                    const double current_merit = v21_pressure_relaxation_merit(
                        controls.v19_full_pressure_flux_step,
                        current.continuity_normalized,
                        current.momentum_equation_residual_rhs_relative,
                        current.flux_velocity_mismatch_normalized);
                    const double reference_merit = v21_pressure_relaxation_merit(
                        controls.v19_full_pressure_flux_step,
                        reference.continuity_normalized,
                        reference.momentum_equation_residual_rhs_relative,
                        reference.flux_velocity_mismatch_normalized);
                    if (reference_merit > 0.0) {
                    effective_alpha_p = adapt_relaxation_factor_windowed(
                        effective_alpha_p,
                        reference_merit,
                        current_merit,
                        controls.adaptive_relaxation.min_alpha_p,
                        controls.adaptive_relaxation.max_alpha_p,
                        controls.adaptive_relaxation,
                        adaptive_p_improvement_streak,
                        adaptive_p_degradation_streak);
                    } else {
                        adaptive_p_improvement_streak = 0;
                        adaptive_p_degradation_streak = 0;
                    }
                }
            }
        }

        double last_selected_pressure_alpha = std::numeric_limits<double>::quiet_NaN();
        double last_line_search_before = std::numeric_limits<double>::quiet_NaN();
        double last_line_search_after = std::numeric_limits<double>::quiet_NaN();
        std::size_t pressure_correctors = configured_pcorr;
        std::size_t piso_neighbor_predictor_rebuilds = 0;
        double piso_neighbor_predictor_flux_delta_linf = 0.0;
        const bool supports_multiple_pressure_corrections =
            controls.algorithm == PressureVelocityAlgorithm::PISO ||
            controls.algorithm == PressureVelocityAlgorithm::PIMPLE ||
            controls.algorithm == PressureVelocityAlgorithm::FRACTIONAL_STEP;
        if (supports_multiple_pressure_corrections &&
            controls.acceleration.adaptive_pressure_correctors) {
            pressure_correctors = result.history.empty()
                ? controls.acceleration.pressure_correctors.max_correctors
                : choose_pressure_correctors(
                      result.history.back().continuity_normalized,
                      controls.acceleration.pressure_correctors);
        }
        // Evaluate CFL from two *accepted* preceding states only; never from
        // a rejected nonlinear retry. Fixed for every attempt of this iter.
        if (controls.pseudo_momentum.enabled && result.history.size()>=2) {
            const auto& previous=result.history[result.history.size()-2];
            const auto& latest=result.history.back();
            v17_pseudo_cfl=next_v17_pseudo_cfl(
                v17_pseudo_cfl,
                previous.momentum_equation_residual_rhs_relative,
                latest.momentum_equation_residual_rhs_relative,
                previous.continuity_normalized,latest.continuity_normalized,
                controls.pseudo_momentum);
        }
        double v32_trial_fraction=1.0;
        if (controls.outer_ptc_v31_experimental) {
            auto& controller=controls.sst_runtime->outer_ptc_v31;
            controller.validate();
            v17_pseudo_cfl=controller.cfl;
            auto& sst_controller=controls.sst_runtime->coupled_ptc_v30;
            sst_controller.cfl=controller.cfl;
            // V31 froze the SST inner CFL, removing V30's safeguard during
            // the nonlinear k-omega solve. V32 freezes only the *external*
            // CFL; inner adaptation occurs in an attempt-local copy and is
            // rolled back unless the entire outer transaction commits.
            sst_controller.frozen_by_outer_v31=
                !controls.outer_ptc_balance_v32_experimental;
            // A ratio, not a dimensional mass: the same cap applies to the
            // pre-pseudo diagonals of momentum and k/omega. The CFL remains
            // fixed during a nonlinear retry and is published only at commit.
            if (controls.outer_ptc_balance_v32_experimental) {
                v32_trial_fraction=controls.sst_runtime->outer_ptc_balance_v32
                    .trial_fraction(retry_controller.retries());
                sst_controller.v32_max_diagonal_fraction=
                    controls.outer_ptc_balance_sst_v32
                        ? v32_trial_fraction
                        : std::numeric_limits<double>::infinity();
                if (controls.diagnostics.iteration_trace)
                    std::cerr<<"OUTER_PTC_BALANCE_V32 iteration="<<iter
                             <<" retry="<<retry_controller.retries()
                             <<" cfl="<<v17_pseudo_cfl
                             <<" fraction="<<v32_trial_fraction
                             <<" momentum_cap="<<controls.outer_ptc_balance_momentum_v32
                             <<" sst_cap="<<controls.outer_ptc_balance_sst_v32
                             <<" qualified=0\n";
            } else {
                sst_controller.v32_max_diagonal_fraction=
                    std::numeric_limits<double>::infinity();
            }
        }
        const auto U_old = U;
        const auto p_old = p;

        // Exact conservative flux consumed by the momentum assembly.
        if (freeze_state)
            phi_used_in_momentum = mass_flux;

        auto grad_p = controls.pressure_high_low_production
            ? gauss_gradient_pressure_high_low(
                mesh,geometry,p,pressure_low_state,pressure_bcs,
                controls.pressure_projected_face_interpolation)
            : gauss_gradient_with_boundary(
                p,mesh,geometry,pressure_bcs,
                controls.pressure_skew_corrected_gradient,
                controls.pressure_projected_face_interpolation);

        Field<double, Location::CELL> body_x(mesh.n_cells(), "body_x", "N/m3", 1);
        Field<double, Location::CELL> body_y(mesh.n_cells(), "body_y", "N/m3", 1);
        Field<double, Location::CELL> body_z(mesh.n_cells(), "body_z", "N/m3", 1);
        const bool cellwise_body_force = !controls.cell_body_force.empty();
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const auto& cell_force = cellwise_body_force
                ? controls.cell_body_force[c]
                : controls.body_force;
            body_x(c) = controls.density * cell_force.x;
            body_y(c) = controls.density * cell_force.y;
            body_z(c) = controls.density * cell_force.z;
        }

        ScalarBoundaryConditions ubc_x, ubc_y, ubc_z;
        Field<double, Location::CELL> u_x_field(mesh.n_cells(), "Ux_reconstruction", "m/s", 1);
        Field<double, Location::CELL> u_y_field(mesh.n_cells(), "Uy_reconstruction", "m/s", 1);
        Field<double, Location::CELL> u_z_field(mesh.n_cells(), "Uz_reconstruction", "m/s", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            u_x_field(c) = U.component_data(0)[c];
            u_y_field(c) = U.component_data(1)[c];
            u_z_field(c) = U.component_data(2)[c];
        }
        for (const auto& [name, bc] : velocity_bcs) {
            ubc_x[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.x, 0.0};
            ubc_y[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.y, 0.0};
            ubc_z[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.z, 0.0};
        }

        // Picard-freeze the flux sign from the preceding accepted RC face
        // flux. Every momentum matrix and H/A flux uses these same signs.
        const auto face_velocity_conditions = mixed_velocity_faces(mass_flux);
        // Momentum predictor uses the conservative corrected flux from the
        // previous outer iteration. This preserves the SIMPLE phi state.

        // A Picard step uses the previous accepted wall state to assemble
        // momentum, then post-SST residuals use the final candidate wall state.
        std::optional<SSTWallOperatorV17> active_wall_operator_v17;
        if (controls.sst_runtime && controls.sst_runtime->active_openfoam_wall_operator_v17) {
            std::map<std::string,std::array<double,3>> wall_velocity;
            for (const auto& [name,bc] : velocity_bcs)
                if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    wall_velocity.emplace(name,std::array<double,3>{bc.value.x,bc.value.y,bc.value.z});
            const auto ref=prepare_openfoam_wall_triplet_v15(
                mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                U,*controls.sst_runtime->state,wall_velocity,
                controls.sst_runtime->transport);
            active_wall_operator_v17=make_sst_wall_operator_v17(
                mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                ref,controls.sst_runtime->transport,
                controls.sst_runtime->omega_wall_owner_constraints_v18);
        }
        const auto* active_wall_face_mu_v17 = active_wall_operator_v17
            ? &active_wall_operator_v17->momentum_face_mu : nullptr;
        Field<double, Location::CELL> transpose_force;
        if (include_transpose_stress) {
            const auto stress_flux = transpose_viscous_stress_flux(
                mesh, geometry, U, mu_eff, {ubc_x, ubc_y, ubc_z},
                controls.sst_runtime
                    ? controls.density * controls.kinematic_viscosity : -1.0,
                    active_wall_face_mu_v17);
            transpose_force = viscous_face_flux_divergence(mesh, geometry, stress_flux);
            for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
                body_x(cell) += transpose_force(cell,0);
                body_y(cell) += transpose_force(cell,1);
                body_z(cell) += transpose_force(cell,2);
            }
        }

        if (controls.momentum_assembly_forensic_callback)
            controls.momentum_assembly_forensic_callback(
                iter, p, grad_p, U, mass_flux, mu_eff);
        if (controls.momentum_predictor_callback)
            controls.momentum_predictor_callback(iter, mu_eff, U, mass_flux);

        const double sst_molecular_wall_mu = controls.sst_runtime
            ? controls.density * controls.kinematic_viscosity : -1.0;
        auto ex = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme,
            &u_x_field, nullptr, sst_molecular_wall_mu,
            outer_high_order_fraction,
            face_velocity_conditions ? &(*face_velocity_conditions)[0] : nullptr,
            active_wall_face_mu_v17);
        auto ey = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme,
            &u_y_field, nullptr, sst_molecular_wall_mu,
            outer_high_order_fraction,
            face_velocity_conditions ? &(*face_velocity_conditions)[1] : nullptr,
            active_wall_face_mu_v17);
        auto ez = assemble_momentum_component(
            mesh, geometry, mass_flux, grad_p, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme,
            &u_z_field, nullptr, sst_molecular_wall_mu,
            outer_high_order_fraction,
            face_velocity_conditions ? &(*face_velocity_conditions)[2] : nullptr,
            active_wall_face_mu_v17);
        if (controls.tvd_deferred_memory.enabled ||
            controls.tvd_secant_implicit.enabled) {
            const std::array<
                const Field<double, Location::CELL>*,3> components{
                    &u_x_field, &u_y_field, &u_z_field};
            const std::array<const ScalarBoundaryConditions*,3> boundaries{
                &ubc_x, &ubc_y, &ubc_z};
            const std::array<ScalarEquation*,3> equations{&ex,&ey,&ez};
            for (std::size_t d = 0; d < 3; ++d) {
                auto raw_trace = audit_momentum_picard_tvd_faces(
                    mesh, geometry, mass_flux, *components[d], *boundaries[d],
                    controls.convection_scheme, outer_high_order_fraction);
                std::vector<double> raw(mesh.n_faces(), 0.0);
                for (const auto& record : raw_trace.faces)
                    raw[record.face] = record.deferred_tvd_flux;
                std::vector<double> u_input(mesh.n_cells(),0.0);
                for (std::size_t cell=0;cell<mesh.n_cells();++cell)
                    u_input[cell]=(*components[d])(cell);
                const auto applied = apply_tvd_deferred_defect_memory(
                    mesh, *equations[d], raw,
                    controls.tvd_deferred_memory.enabled && committed_tvd_deferred_valid
                        ? &committed_tvd_deferred[d] : nullptr,
                    controls.tvd_deferred_memory.enabled
                        ? controls.tvd_deferred_memory.previous_weight : 0.0,
                    controls.tvd_deferred_memory.oscillation_only &&
                        committed_tvd_two_steps_valid
                        ? &committed_tvd_two_steps_ago[d] : nullptr,
                    controls.tvd_deferred_memory.oscillation_only &&
                        committed_tvd_two_steps_valid);
                if (controls.tvd_secant_implicit.enabled && committed_tvd_deferred_valid) {
                    const auto regularization = add_tvd_secant_implicit_diagonal(
                        mesh, *equations[d], raw, committed_tvd_raw[d],
                        u_input, committed_tvd_input_velocity[d], mass_flux,
                        controls.tvd_secant_implicit);
                    if (controls.diagnostics.iteration_trace && d==1)
                        std::cerr << "TVD_SECANT_IMPLICIT iteration=" << iter
                                  << " regularized_cells=" << regularization.regularized_cells
                                  << " max_diagonal=" << regularization.max_diagonal_added
                                  << " max_spectral_fraction="
                                  << regularization.max_ratio_to_face_flux << '\n';
                }
                pending_effective_tvd_deferred[d] = applied.used;
                pending_raw_tvd_deferred[d] = raw;
                pending_tvd_input_velocity[d] = std::move(u_input);
                pending_memory_source_linf[d] = applied.max_cell_source_difference;
                if (d == 1 && controls.momentum_picard_tvd_face_callback) {
                    // Report both physical raw and effective lag-filtered flux
                    // rather than falsely labelling the raw value as assembled.
                    for (auto& record : raw_trace.faces) {
                        record.effective_deferred_flux = applied.used[record.face];
                        record.memory_filter_delta =
                            record.effective_deferred_flux - record.deferred_tvd_flux;
                    }
                    pending_picard_tvd_trace = std::move(raw_trace);
                }
            }
        } else if (controls.momentum_picard_tvd_face_callback) {
            // This is the same U_y and phi passed to assemble_momentum_component
            // above, including the iteration's continuation fraction. It is
            // intentionally distinct from the accepted-state reconstruction.
            pending_picard_tvd_trace = audit_momentum_picard_tvd_faces(
                mesh,geometry,mass_flux,u_y_field,ubc_y,
                controls.convection_scheme,outer_high_order_fraction);
        }


        if (controls.momentum_equation_forensic_callback) {
            controls.momentum_equation_forensic_callback(
                iter, 0, ex, geometry, grad_p, body_x);
            controls.momentum_equation_forensic_callback(
                iter, 1, ey, geometry, grad_p, body_y);
            controls.momentum_equation_forensic_callback(
                iter, 2, ez, geometry, grad_p, body_z);
        }

        if (controls.transient.enabled) {
            const auto& u_n = *controls.transient.previous;
            const auto& u_nm1 = controls.transient.history_valid
                ? *controls.transient.previous_previous : u_n;
            const bool use_bdf2 = controls.transient.scheme == cfdx::core::TimeScheme::BDF2 &&
                                  controls.transient.history_valid;
            add_dual_time_mass_term(ex, u_n, u_nm1, geometry, controls.density,
                                    controls.transient.dt, controls.transient.dt_previous,
                                    use_bdf2, 0);
            add_dual_time_mass_term(ey, u_n, u_nm1, geometry, controls.density,
                                    controls.transient.dt, controls.transient.dt_previous,
                                    use_bdf2, 1);
            add_dual_time_mass_term(ez, u_n, u_nm1, geometry, controls.density,
                                    controls.transient.dt, controls.transient.dt_previous,
                                    use_bdf2, 2);
        }
        double v17_mass_min=std::numeric_limits<double>::quiet_NaN();
        double v17_mass_max=std::numeric_limits<double>::quiet_NaN();
        std::size_t v17_capped_cells=0;
        if (controls.pseudo_momentum.enabled || controls.outer_ptc_v31_experimental) {
            if (active_wall_face_mu_v17 && !controls.outer_ptc_v31_experimental)
                throw std::invalid_argument(
                    "V17 active wall functions do not yet support legacy pseudo-momentum");
            auto masses=make_v17_pseudo_momentum_masses(
                mesh,geometry,mass_flux,mu_eff,controls.density,v17_pseudo_cfl,
                controls.sst_runtime
                    ? controls.density*controls.kinematic_viscosity : -1.0,
                controls.outer_ptc_v31_experimental
                    ? active_wall_face_mu_v17 : nullptr);
            v17_capped_cells = cap_v17_pseudo_momentum_masses(
                masses, ex, ey, ez,
                controls.outer_ptc_v31_experimental
                    ? (controls.outer_ptc_balance_v32_experimental &&
                       controls.outer_ptc_balance_momentum_v32
                           ? v32_trial_fraction : 1.0)
                    : controls.pseudo_momentum.maximum_diagonal_fraction);
            add_v17_pseudo_momentum_term(ex,U_old,masses,0);
            add_v17_pseudo_momentum_term(ey,U_old,masses,1);
            add_v17_pseudo_momentum_term(ez,U_old,masses,2);
            const auto bounds=std::minmax_element(masses.begin(),masses.end());
            v17_mass_min=*bounds.first;
            v17_mass_max=*bounds.second;
        }

        // Preserve the pre-solve velocity for equation relaxation. HbyA is
        // reconstructed after the momentum solve from the solved neighbour
        // values, which is the algebraic H/A split of the assembled equation.
        // The pre-solve state is used only for the under-relaxation source term.
        Vector ux(mesh.n_cells()), uy(mesh.n_cells()), uz(mesh.n_cells());
        Vector ux_old(mesh.n_cells()), uy_old(mesh.n_cells()), uz_old(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            ux(c) = U.component_data(0)[c];
            uy(c) = U.component_data(1)[c];
            uz(c) = U.component_data(2)[c];
            ux_old(c) = ux(c);
            uy_old(c) = uy(c);
            uz_old(c) = uz(c);
        }

        cfdx::core::SolverResult rx{}, ry{}, rz{};
        double pressure_residual = std::numeric_limits<double>::infinity();
        std::size_t pressure_iterations = 0;
        std::array<double, 3> simplec_signed_offdiag_sum_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_abs_offdiag_sum_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_denominator_signed_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_denominator_abs_audit{0.0, 0.0, 0.0};
        std::array<double, 3> simplec_denominator_gap_audit{0.0, 0.0, 0.0};
        std::array<std::size_t, 3> simplec_denominator_worst_cell_audit{0, 0, 0};
        // Rhie-Chow operator identity audit (common scope): each algorithm
        // branch evaluates it right after its authoritative flux rebuild.
        // Evidence channel only; no acceptance criterion consumes it.
        double coupled_outer_step_applied_v43 = 1.0;
        const auto phase_escape_trial_v45 = coupled_phase_escape_v45.trial(
            iter, v43_field_probe && v43_field_probe->period_four_candidate(),
            controls.coupled_period_four_damping_v45_experimental,
            controls.coupled_period_four_damping_steps_v45,
            retry_controller.alpha_u(1.0));
        bool rhie_chow_identity_evaluated = false;
        double rhie_chow_identity_audit_linf = 0.0;
        double simplec_plain_rau_gap = 0.0;
        std::size_t simplec_plain_rau_worst = 0;
        // V20 diagnostic-only reference for the V19 split pressure state.
        std::optional<Field<double,Location::FACE>> v20_relaxed_pressure_flux;

        if (controls.algorithm == PressureVelocityAlgorithm::COUPLED) {
            // V10: opt-in, phase-only evidence.  This is deliberately outside
            // the solver's residual/acceptance logic and cannot certify it.
            if (std::getenv("CFDX_BFS_STAGE_TRACE"))
                std::cerr << "BFS_PHASE_V10 iteration=" << iter
                          << " phase=COUPLED_LINEAR_BEGIN\n";
            const auto coupled_result = solve_coupled_momentum_continuity(
                mesh, geometry, ex, ey, ez, U_old, p_old,
                velocity_bcs, pressure_bcs, controls.density,
                controls.kinematic_viscosity,
                controls.pressure_reference_cell,
                controls.pressure_reference_value,
                controls.coupling.coupled_max_iterations,
                controls.coupling.coupled_linear_tolerance,
                controls.diagnostics,
                controls.coupling.schur_model,
                controls.coupled_linear_solver,
                controls,
                retry_controller.alpha_u(effective_alpha_u),
                U, p, &result.coupled_linear_plan);
            result.coupled_linear_plan_resolved = true;
            if (std::getenv("CFDX_BFS_STAGE_TRACE"))
                std::cerr << "BFS_PHASE_V10 iteration=" << iter
                          << " phase=COUPLED_LINEAR_END status="
                          << static_cast<int>(coupled_result.status)
                          << " iterations=" << coupled_result.iterations
                          << " true_relative=" << coupled_result.residual_relative
                          << " actual_krylov=" << cfdx::core::to_string(result.coupled_linear_plan.krylov)
                          << " actual_preconditioner=" << cfdx::core::to_string(result.coupled_linear_plan.preconditioner)
                          << " resolved_schur=" << cfdx::physics::to_string(result.resolved_coupled_schur_model)
                          << "\n";
            if (coupled_result.status != cfdx::core::SolverStatus::CONVERGED) {
                // Report the stalled linear state in scientific notation.
                // std::to_string truncates to six fixed decimals, so a
                // MAX_ITER stall printed as "residual=0.000000" hid the
                // plateau level and made the failure impossible to calibrate
                // against the coupled linear gate or the nested PCD budget.
                const auto coupled_linear_failure_value = [](double value) {
                    std::ostringstream out;
                    out << std::scientific << std::setprecision(6) << value;
                    return out.str();
                };
                coupled_linear_failed_v37 = true;
                throw NonlinearRetryableFailure(
                    "solve_steady_incompressible: coupled momentum-continuity solve did not converge "
                    "(status=" + std::to_string(static_cast<int>(coupled_result.status)) +
                    ", iterations=" + std::to_string(coupled_result.iterations) +
                    ", residual=" + coupled_linear_failure_value(coupled_result.residual) +
                    ", relative=" + coupled_linear_failure_value(coupled_result.residual_relative) + ")");
            }

            // Diagnose the full monolithic candidate before independent U/p
            // under-relaxation. This separates a block/flux inconsistency from
            // continuity loss introduced by the nonlinear update itself.
            if (controls.diagnostics.coupled_matrix_summary) {
                const auto* cell_faces = mesh.cells().faces_data();
                const auto* cell_offsets = mesh.cells().offsets_data();
                std::array<std::vector<double>, 3> candidate_rAU;
                for (auto& component : candidate_rAU)
                    component.resize(mesh.n_cells());
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    candidate_rAU[0][c] = 1.0 / ex.diagonal[c];
                    candidate_rAU[1][c] = 1.0 / ey.diagonal[c];
                    candidate_rAU[2][c] = 1.0 / ez.diagonal[c];
                }
                const auto candidate_flux = make_rhie_chow_flux_from_physical_velocity(
                    mesh, geometry, U, p, candidate_rAU,
                    controls.density, velocity_bcs, pressure_bcs,
                    controls.pressure_skew_corrected_gradient);
                double candidate_continuity_linf = 0.0;
                double candidate_continuity_signed = 0.0;
                std::size_t candidate_continuity_cell = 0;
                double candidate_global_divergence_sum = 0.0;
                double candidate_boundary_flux_sum = 0.0;
                double candidate_boundary_flux_abs_sum = 0.0;
                double candidate_gauge_cell_divergence = 0.0;
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    double divergence = 0.0;
                    for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
                        const std::size_t face = cell_faces[k];
                        divergence += mesh.ownership().owner(face) == c
                            ? candidate_flux(face) : -candidate_flux(face);
                    }
                    candidate_global_divergence_sum += divergence;
                    if (c == controls.pressure_reference_cell)
                        candidate_gauge_cell_divergence = divergence;
                    if (std::abs(divergence) > candidate_continuity_linf) {
                        candidate_continuity_linf = std::abs(divergence);
                        candidate_continuity_signed = divergence;
                        candidate_continuity_cell = c;
                    }
                }
                // Every internal face cancels exactly in the global signed
                // divergence sum. The remaining sum must equal the signed
                // boundary flux; report both independently to audit topology,
                // face orientation, and the reference-pressure cell.
                for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
                    if (mesh.ownership().neighbour(face) >= 0) continue;
                    // Include every topological boundary face. Empty/extruded
                    // faces should contribute zero geometrically; excluding
                    // them here would hide a mesh/flux inconsistency.
                    candidate_boundary_flux_sum += candidate_flux(face);
                    candidate_boundary_flux_abs_sum += std::abs(candidate_flux(face));
                }
                std::cerr << "COUPLED_GLOBAL_CONSERVATION_AUDIT"
                          << " divergence_sum=" << candidate_global_divergence_sum
                          << " signed_boundary_flux_sum=" << candidate_boundary_flux_sum
                          << " boundary_flux_abs_sum=" << candidate_boundary_flux_abs_sum
                          << " divergence_minus_boundary_flux="
                          << (candidate_global_divergence_sum - candidate_boundary_flux_sum)
                          << " gauge_cell=" << controls.pressure_reference_cell
                          << " gauge_cell_divergence=" << candidate_gauge_cell_divergence
                          << " max_cell=" << candidate_continuity_cell
                          << " max_cell_divergence_signed=" << candidate_continuity_signed
                          << " max_cell_divergence_linf=" << candidate_continuity_linf
                          << "\n";

                const auto candidate_grad_p =
                    gauss_gradient_with_boundary(
            p, mesh, geometry, pressure_bcs,
            controls.pressure_skew_corrected_gradient,
            controls.pressure_projected_face_interpolation);
                double velocity_flux_div = 0.0;
                double pressure_orth_div = 0.0;
                double pressure_nonorth_div = 0.0;
                double physical_velocity_compensation_div = 0.0;
                for (Offset k = cell_offsets[candidate_continuity_cell];
                     k < cell_offsets[candidate_continuity_cell + 1]; ++k) {
                    const std::size_t face = cell_faces[k];
                    const bool owner = mesh.ownership().owner(face) == candidate_continuity_cell;
                    const double sign = owner ? 1.0 : -1.0;
                    const auto nr = mesh.ownership().neighbour(face);
                    double orth = 0.0;
                    double nonorth = 0.0;
                    double compensation = 0.0;
                    if (nr >= 0) {
                        const std::size_t face_owner = mesh.ownership().owner(face);
                        const std::size_t face_neighbour = static_cast<std::size_t>(nr);
                        const Vec3 sf = geometry.face_area_vectors[face];
                        const Vec3 dvec = geometry.cell_centres[face_neighbour] -
                                          geometry.cell_centres[face_owner];
                        const Vec3 cx{
                            0.5 * (geometry.cell_volumes[face_owner] * candidate_rAU[0][face_owner] +
                                   geometry.cell_volumes[face_neighbour] * candidate_rAU[0][face_neighbour]),
                            0.5 * (geometry.cell_volumes[face_owner] * candidate_rAU[1][face_owner] +
                                   geometry.cell_volumes[face_neighbour] * candidate_rAU[1][face_neighbour]),
                            0.5 * (geometry.cell_volumes[face_owner] * candidate_rAU[2][face_owner] +
                                   geometry.cell_volumes[face_neighbour] * candidate_rAU[2][face_neighbour])};
                        const Vec3 gp{
                            0.5 * (candidate_grad_p.component_data(0)[face_owner] +
                                   candidate_grad_p.component_data(0)[face_neighbour]),
                            0.5 * (candidate_grad_p.component_data(1)[face_owner] +
                                   candidate_grad_p.component_data(1)[face_neighbour]),
                            0.5 * (candidate_grad_p.component_data(2)[face_owner] +
                                   candidate_grad_p.component_data(2)[face_neighbour])};
                        // In the PHYSICAL_VELOCITY contract the authoritative
                        // flux is U_f.Sf minus the pressure jump and deferred
                        // pressure flux, plus the dAU*grad(p) compensation.
                        // Keep these signs identical to make_rhie_chow_mass_flux.
                        orth = -controls.density *
                            pressure_face_orthogonal_coefficient(sf, dvec, cx) *
                            (p(face_neighbour) - p(face_owner));
                        nonorth = -controls.density *
                            pressure_face_nonorthogonal_flux(sf, dvec, cx, gp);
                        compensation = controls.density *
                            (cx.x * gp.x * sf.x + cx.y * gp.y * sf.y +
                             cx.z * gp.z * sf.z);
                    }
                    pressure_orth_div += sign * orth;
                    pressure_nonorth_div += sign * nonorth;
                    physical_velocity_compensation_div += sign * compensation;
                    velocity_flux_div += sign *
                        (candidate_flux(face) - orth - nonorth - compensation);
                }
                std::cerr << "COUPLED_PRE_RELAXATION_FLUX_BREAKDOWN"
                          << " cell=" << candidate_continuity_cell

                          << " velocity_interpolation_div=" << velocity_flux_div
                          << " pressure_orthogonal_div=" << pressure_orth_div
                          << " pressure_nonorthogonal_div=" << pressure_nonorth_div
                          << " physical_velocity_compensation_div=" << physical_velocity_compensation_div
                          << " sum=" << (velocity_flux_div + pressure_orth_div +
                                         pressure_nonorth_div + physical_velocity_compensation_div)
                          << " authoritative_div=" << candidate_continuity_signed
                          << "\n";
                std::cerr << "COUPLED_PRE_RELAXATION_FLUX_CONTINUITY linf="
                          << std::scientific << std::setprecision(8)
                          << candidate_continuity_linf
                          << " signed=" << candidate_continuity_signed
                          << " cell=" << candidate_continuity_cell
                          << " alpha_u=" << retry_controller.alpha_u(effective_alpha_u)
                          << " alpha_p=" << retry_controller.alpha_p(effective_alpha_p)
                          << "\n";
            }

            // V43: a guarded nonlinear field update AFTER a certified
            // monolithic linear solve. This is Picard globalization, NOT a
            // second way of certifying A*x=b. Retries actually shorten the
            // carrier step, rather than changing only the segregated alpha_u
            // and alpha_p that COUPLED does not apply. We do not interpolate
            // old and new face flux: those may correspond to different rAU.
            // The authoritative Rhie-Chow flux is reconstructed below from
            // this very U,p candidate and the current frozen momentum
            // diagonals. It is audited by the existing physical gates.
            if (controls.coupled_outer_step_v43_experimental > 0.0 ||
                phase_escape_trial_v45.active) {
                coupled_outer_step_applied_v43=phase_escape_trial_v45.active ?
                    phase_escape_trial_v45.beta : coupled_outer_beta_v43(
                        controls.coupled_outer_step_v43_experimental,
                        retry_controller.alpha_u(1.0));
                for (std::size_t cell=0;cell<mesh.n_cells();++cell) {
                    for (std::size_t component=0;component<3;++component) {
                        const double previous=U_old.component_data(component)[cell];
                        const double candidate=U.component_data(component)[cell];
                        U.component_data(component)[cell]=coupled_outer_blend_v43(
                            previous,candidate,coupled_outer_step_applied_v43);
                    }
                    p(cell)=coupled_outer_blend_v43(
                        p_old(cell),p(cell),coupled_outer_step_applied_v43);
                }
                if (controls.diagnostics.iteration_trace && phase_escape_trial_v45.active)
                    std::cerr<<"COUPLED_PHASE_ESCAPE_V45 iteration="<<iter
                             <<" retry="<<retry_controller.retries()
                             <<" triggered="<<(phase_escape_trial_v45.triggered?1:0)
                             <<" remaining="<<phase_escape_trial_v45.remaining_including_this
                             <<" beta="<<coupled_outer_step_applied_v43
                             <<" linear_candidate_certified=1"
                             <<" updated_linear_solution_certified=0"
                             <<" physical_qualified=0\n";
                if (controls.diagnostics.iteration_trace && !phase_escape_trial_v45.active)
                    std::cerr<<"COUPLED_CARRIER_STEP_V43 iteration="<<iter
                             <<" retry="<<retry_controller.retries()
                             <<" beta="<<coupled_outer_step_applied_v43
                             <<" linear_candidate_certified=1"
                             <<" updated_linear_solution_certified=0"
                             <<" flux_rederived_from_trial_state=1"
                             <<" physical_qualified=0\n";
            }
            apply_pressure_gauge_policy(
                p, geometry, controls.pressure_gauge_policy,
                controls.pressure_reference_cell, controls.pressure_reference_value,
                has_fixed_pressure_boundary);

            // The coupled solve has already enforced the discrete continuity
            // equation. Rebuild the authoritative flux from the same RC
            // interpolation used by the block matrix for diagnostics and for
            // the next nonlinear iteration.
            std::array<std::vector<double>,3> coupled_rAU;
            coupled_rAU[0].resize(mesh.n_cells());
            coupled_rAU[1].resize(mesh.n_cells());
            coupled_rAU[2].resize(mesh.n_cells());
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                if (!(ex.diagonal[c] > 0.0) || !std::isfinite(ex.diagonal[c]) ||
                    !(ey.diagonal[c] > 0.0) || !std::isfinite(ey.diagonal[c]) ||
                    !(ez.diagonal[c] > 0.0) || !std::isfinite(ez.diagonal[c]))
                    throw std::runtime_error(
                        "solve_steady_incompressible: invalid coupled momentum diagonal after solve");
                // make_rhie_chow_mass_flux() and rhie_chow_pressure_flux_internal()
                // take the INVERSE momentum diagonal rAU = 1/A_P and apply the
                // cell volume internally (cx = 0.5*(V*rAU) = dAU). Passing V/A_P
                // here double-counted the volume and shrank the Rhie-Chow
                // pressure correction of the authoritative coupled flux by a
                // factor V, silently removing the checkerboard suppression.
                // On the constant/linear-pressure channel cases the correction
                // vanishes and the defect was invisible (Couette/Poiseuille
                // passed); on the lid-driven cavity the assembled continuity
                // rows and the authoritative flux disagreed by div(interp(U)),
                // a permanent continuity floor (~8.5e-3 against the 1e-7 gate)
                // that no iteration budget or linear tolerance could remove:
                // the velocity froze at the Picard fixed point and the
                // stagnation detector terminated the case.
                coupled_rAU[0][c] = 1.0 / ex.diagonal[c];
                coupled_rAU[1][c] = 1.0 / ey.diagonal[c];
                coupled_rAU[2][c] = 1.0 / ez.diagonal[c];
            }
            mass_flux = make_rhie_chow_flux_from_physical_velocity(
                mesh, geometry, U, p, coupled_rAU,
                controls.density, velocity_bcs, pressure_bcs,
                    controls.pressure_skew_corrected_gradient);

            // Operator identity audit: re-derive the authoritative flux from
            // the same state through the same production call and compare
            // against the stored field. Zero unless the stored flux was
            // mutated after the rebuild.
            {
                const auto rederived_flux = make_rhie_chow_flux_from_physical_velocity(
                    mesh, geometry, U, p, coupled_rAU,
                    controls.density, velocity_bcs, pressure_bcs,
                    controls.pressure_skew_corrected_gradient);
                for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                    rhie_chow_identity_audit_linf = std::max(
                        rhie_chow_identity_audit_linf,
                        std::abs(mass_flux(f) - rederived_flux(f)));
                rhie_chow_identity_evaluated = true;
            }

            // A coupled iteration has no segregated momentum/pressure Krylov
            // sub-solves. The block residual is the authoritative linear metric.
            // The independent nonlinear momentum residual below remains the
            // final physical acceptance metric.
        } else {
        relax_momentum_equation(ex, ux_old, retry_controller.alpha_u(effective_alpha_u));
        relax_momentum_equation(ey, uy_old, retry_controller.alpha_u(effective_alpha_u));
        relax_momentum_equation(ez, uz_old, retry_controller.alpha_u(effective_alpha_u));

        // Equation relaxation is already part of the matrix. The Krylov solve
        // therefore returns the solution of the relaxed equation directly;
        // applying alpha_u again here would double-relax the predictor.
        ScalarSolveControls momentum_solve_controls;
        momentum_solve_controls.max_iterations = controls.linear_max_iterations;
        momentum_solve_controls.tolerance = linear_tolerance;
        momentum_solve_controls.relaxation = 1.0;
        momentum_solve_controls.problem = LinearProblemKind::Momentum;
        momentum_solve_controls.linear_solver = controls.momentum_linear_solver;
        rx = solve_scalar_equation(ex, ux, momentum_solve_controls);
        ry = solve_scalar_equation(ey, uy, momentum_solve_controls);
        rz = solve_scalar_equation(ez, uz, momentum_solve_controls);

        auto require_linear_convergence = [](const char* component, const auto& solve) {
            if (solve.status != cfdx::core::SolverStatus::CONVERGED) {
                std::ostringstream message;
                message << "solve_steady_incompressible: " << component
                        << " momentum solve did not converge (status="
                        << static_cast<int>(solve.status)
                        << ", iterations=" << solve.iterations
                        << ", residual=" << std::scientific << std::setprecision(3)
                        << solve.residual
                        << ", relative=" << solve.residual_relative << ")";
                throw NonlinearRetryableFailure(message.str());
            }
        };
        require_linear_convergence("Ux", rx);
        require_linear_convergence("Uy", ry);
        require_linear_convergence("Uz", rz);
        if (controls.momentum_linear_relaxation_callback) {
            pending_alpha_u = retry_controller.alpha_u(effective_alpha_u);
            pending_linear_relaxation[0] = audit_momentum_linear_relaxation(
                ex, ux, ux_old, pending_alpha_u);
            pending_linear_relaxation[1] = audit_momentum_linear_relaxation(
                ey, uy, uy_old, pending_alpha_u);
            pending_linear_relaxation[2] = audit_momentum_linear_relaxation(
                ez, uz, uz_old, pending_alpha_u);
        }

        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            U.component_data(0)[c] = ux(c);
            U.component_data(1)[c] = uy(c);
            U.component_data(2)[c] = uz(c);
        }

        std::array<std::vector<double>,3> rAU, rAtU, hbya;
        for (std::size_t d = 0; d < 3; ++d) {
            rAU[d].resize(mesh.n_cells());
            rAtU[d].resize(mesh.n_cells());
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                const ScalarEquation* eqs[3] = {&ex, &ey, &ez};
                const double a = eqs[d]->diagonal[c];
                if (!(a > 0.0) || !std::isfinite(a))
                    throw std::runtime_error(
                        "solve_steady_incompressible: invalid momentum diagonal");
                rAU[d][c] = 1.0 / a;
                rAtU[d][c] = rAU[d][c];
                if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
                    double signed_offdiag_sum = 0.0;
                    double abs_offdiag_sum = 0.0;
                    const auto begin = eqs[d]->matrix.row_offsets_data()[c];
                    const auto end = eqs[d]->matrix.row_offsets_data()[c + 1];
                    for (std::uint32_t k = begin; k < end; ++k) {
                        const auto col = eqs[d]->matrix.columns_data()[k];
                        // Boundary contributions are assembled into the owner
                        // diagonal/source; off-diagonal entries are therefore
                        // genuine internal neighbour coefficients.
                        if (col != c) {
                            const double value = eqs[d]->matrix.values_data()[k];
                            signed_offdiag_sum += value;
                            abs_offdiag_sum += std::abs(value);
                        }
                    }
                    const double diagonal = 1.0 / rAU[d][c];
                    const double denom =
                        simplec_consistent_diagonal(diagonal, signed_offdiag_sum);
                    // Independent audit counterpart. With the expected non-positive
                    // stored off-diagonals this equals the signed-contract
                    // denominator; a nonzero gap flags a matrix-sign mismatch.
                    const double abs_denom =
                        diagonal - abs_offdiag_sum;
                    simplec_signed_offdiag_sum_audit[d] =
                        std::max(simplec_signed_offdiag_sum_audit[d],
                                 std::abs(signed_offdiag_sum));
                    simplec_abs_offdiag_sum_audit[d] =
                        std::max(simplec_abs_offdiag_sum_audit[d], abs_offdiag_sum);
                    simplec_denominator_signed_audit[d] =
                        std::max(simplec_denominator_signed_audit[d], std::abs(denom));
                    simplec_denominator_abs_audit[d] =
                        std::max(simplec_denominator_abs_audit[d], std::abs(abs_denom));
                    const double gap = std::abs(denom - abs_denom);
                    if (gap > simplec_denominator_gap_audit[d])
                        simplec_denominator_worst_cell_audit[d] = c;
                    simplec_denominator_gap_audit[d] =
                        std::max(simplec_denominator_gap_audit[d], gap);
                    rAtU[d][c] = 1.0 / denom;
                }
            }
        }

        hbya[0] = build_hbya(ex, ux, grad_p, 0, rAU[0]);
        hbya[1] = build_hbya(ey, uy, grad_p, 1, rAU[1]);
        hbya[2] = build_hbya(ez, uz, grad_p, 2, rAU[2]);
        auto HbyA = make_hbya_field(hbya);

        if (freeze_state) {
            frozen_rAU = rAU;
            frozen_hbyA = hbya;
            frozen_grad_p = grad_p;
            frozen_state_valid = true;
        }

        // rAU is the inverse integrated momentum diagonal (1/A_P).
        // Pressure-gradient terms in the integrated momentum equation therefore
        // carry the cell volume explicitly when reconstructing cell velocity.
        // phiHbyA is the authoritative predictor flux. The Rhie-Chow term
        // uses the same directional inverse momentum coefficient as the
        // pressure equation, avoiding a second, incompatible face operator.
        // V19: last conservative flux may correspond to the unrelaxed
        // pressure step, not the under-relaxed physical pressure.
        auto v19_flux_pressure_state = p;
        Vector v19_flux_pressure_low(
            controls.pressure_high_low_production ? mesh.n_cells() : 0,0.0);
        auto phiHbyA = make_mass_flux(
            mesh, geometry, HbyA, controls.density, velocity_bcs,
            {}, controls.pressure_projected_face_interpolation,
            face_velocity_conditions ? &*face_velocity_conditions : nullptr);

        if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nr = mesh.ownership().neighbour(f);
                const std::size_t o = mesh.ownership().owner(f);
                if (nr < 0) {
                    const std::size_t patch = geometry.face_patch[f];
                    if (patch >= mesh.boundary().n_patches() ||
                        mesh.boundary().patch(patch).type == PatchType::EMPTY)
                        continue;
                    const auto pbc = pressure_bcs.find(
                        mesh.boundary().patch(patch).name);
                    if (pbc == pressure_bcs.end() ||
                        pbc->second.type != ScalarBoundaryType::FIXED_VALUE)
                        continue;
                    const Vec3 Sf = geometry.face_area_vectors[f];
                    const Vec3 dvec = geometry.face_centres[f] -
                                      geometry.cell_centres[o];
                    const Vec3 delta_dAU{
                        geometry.cell_volumes[o] * (rAtU[0][o] - rAU[0][o]),
                        geometry.cell_volumes[o] * (rAtU[1][o] - rAU[1][o]),
                        geometry.cell_volumes[o] * (rAtU[2][o] - rAU[2][o])};
                    const Vec3 grad_p_cell{
                        grad_p.component_data(0)[o],
                        grad_p.component_data(1)[o],
                        grad_p.component_data(2)[o]};
                    const double correction_coeff = pressure_face_orthogonal_coefficient(
                        Sf, dvec, delta_dAU);
                    const double correction_deferred = pressure_face_nonorthogonal_flux(
                        Sf, dvec, delta_dAU, grad_p_cell);
                    phiHbyA(f) += controls.density * (
                        correction_coeff * (pbc->second.value - p(o)) +
                        correction_deferred);
                    continue;
                }
                const std::size_t n = static_cast<std::size_t>(nr);
                const auto Sf = geometry.face_area_vectors[f];
                const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
                const double dx = 0.5*(geometry.cell_volumes[o]*rAtU[0][o] +
                                       geometry.cell_volumes[n]*rAtU[0][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[0][o] +
                                       geometry.cell_volumes[n]*rAU[0][n]);
                const double dy = 0.5*(geometry.cell_volumes[o]*rAtU[1][o] +
                                       geometry.cell_volumes[n]*rAtU[1][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[1][o] +
                                       geometry.cell_volumes[n]*rAU[1][n]);
                const double dz = 0.5*(geometry.cell_volumes[o]*rAtU[2][o] +
                                       geometry.cell_volumes[n]*rAtU[2][n]) -
                                  0.5*(geometry.cell_volumes[o]*rAU[2][o] +
                                       geometry.cell_volumes[n]*rAU[2][n]);
                const Vec3 delta_dAU{dx, dy, dz};
                const double correction_coeff = pressure_face_orthogonal_coefficient(
                    Sf, dvec, delta_dAU);
                const Vec3 grad_p_face{
                    0.5 * (grad_p.component_data(0)[o] + grad_p.component_data(0)[n]),
                    0.5 * (grad_p.component_data(1)[o] + grad_p.component_data(1)[n]),
                    0.5 * (grad_p.component_data(2)[o] + grad_p.component_data(2)[n])};
                // SIMPLEC changes the full face response, not only its
                // centre-line coefficient. Include the matching deferred
                // difference or skewed-face continuity retains an O(grad p)
                // mismatch between rAtU and rAU.
                const double correction_deferred = pressure_face_nonorthogonal_flux(
                    Sf, dvec, delta_dAU, grad_p_face);
                phiHbyA(f) += controls.density * (
                    correction_coeff * (p(n)-p(o)) + correction_deferred);
            }
        }

    // SIMPLEC's pressure equation uses the consistent response
            // rAtU, while the corrected cell velocity uses rAU. Convert the
            // reconstructed rAtU Rhie-Chow flux back to the rAU response by
            // adding the full directional difference (rAtU-rAU), not only
            // its orthogonal pressure-jump part. Omitting the deferred
            // non-orthogonal remainder leaves a continuity floor on skewed
            // meshes. Fixed-pressure boundary faces require the same response
            // conversion; EMPTY and zero-gradient-pressure faces contribute
            // no pressure correction.
    const auto add_simplec_consistent_response =
        [&](cfdx::core::Field<double, cfdx::core::Location::FACE>& flux) {
            const auto response_grad_p =
                gauss_gradient_with_boundary(
            p, mesh, geometry, pressure_bcs,
            controls.pressure_skew_corrected_gradient,
            controls.pressure_projected_face_interpolation);
            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                    const auto nr = mesh.ownership().neighbour(f);
                    const std::size_t o = mesh.ownership().owner(f);
                    const Vec3 Sf = geometry.face_area_vectors[f];
                    if (nr >= 0) {
                        const std::size_t n = static_cast<std::size_t>(nr);
                        const Vec3 dvec =
                            geometry.cell_centres[n] - geometry.cell_centres[o];
                        const Vec3 delta_dAU{
                            0.5*(geometry.cell_volumes[o]*rAtU[0][o] +
                                 geometry.cell_volumes[n]*rAtU[0][n]) -
                            0.5*(geometry.cell_volumes[o]*rAU[0][o] +
                                 geometry.cell_volumes[n]*rAU[0][n]),
                            0.5*(geometry.cell_volumes[o]*rAtU[1][o] +
                                 geometry.cell_volumes[n]*rAtU[1][n]) -
                            0.5*(geometry.cell_volumes[o]*rAU[1][o] +
                                 geometry.cell_volumes[n]*rAU[1][n]),
                            0.5*(geometry.cell_volumes[o]*rAtU[2][o] +
                                 geometry.cell_volumes[n]*rAtU[2][n]) -
                            0.5*(geometry.cell_volumes[o]*rAU[2][o] +
                                 geometry.cell_volumes[n]*rAU[2][n])};
                        const double orthogonal = pressure_face_orthogonal_coefficient(
                            Sf, dvec, delta_dAU);
                        const Vec3 grad_p_face{
                            0.5*(response_grad_p.component_data(0)[o] +
                                 response_grad_p.component_data(0)[n]),
                            0.5*(response_grad_p.component_data(1)[o] +
                                 response_grad_p.component_data(1)[n]),
                            0.5*(response_grad_p.component_data(2)[o] +
                                 response_grad_p.component_data(2)[n])};
                        const double deferred = pressure_face_nonorthogonal_flux(
                            Sf, dvec, delta_dAU, grad_p_face);
                        flux(f) += controls.density *
                            (orthogonal * (p(n)-p(o)) + deferred);
                        continue;
                    }

                    const std::size_t patch = geometry.face_patch[f];
                    if (patch >= mesh.boundary().n_patches() ||
                        mesh.boundary().patch(patch).type == PatchType::EMPTY)
                        continue;
                    const auto pbc = pressure_bcs.find(
                        mesh.boundary().patch(patch).name);
                    if (pbc == pressure_bcs.end() ||
                        pbc->second.type != ScalarBoundaryType::FIXED_VALUE)
                        continue;
                    const Vec3 dvec =
                        geometry.face_centres[f] - geometry.cell_centres[o];
                    const Vec3 delta_dAU{
                        geometry.cell_volumes[o] * (rAtU[0][o] - rAU[0][o]),
                        geometry.cell_volumes[o] * (rAtU[1][o] - rAU[1][o]),
                        geometry.cell_volumes[o] * (rAtU[2][o] - rAU[2][o])};
                    const double orthogonal = pressure_face_orthogonal_coefficient(
                        Sf, dvec, delta_dAU);
                    const Vec3 grad_p_boundary{
                        response_grad_p.component_data(0)[o],
                        response_grad_p.component_data(1)[o],
                        response_grad_p.component_data(2)[o]};
                    const double deferred = pressure_face_nonorthogonal_flux(
                        Sf, dvec, delta_dAU, grad_p_boundary);
                    flux(f) += controls.density *
                        (orthogonal * (pbc->second.value - p(o)) + deferred);
                }
            };


        for (std::size_t corr = 0; corr < pressure_correctors; ++corr) {
            const std::size_t nc = mesh.n_cells();
            SparseMatrix A(nc, nc);
            Vector b(nc, 0.0);
            std::vector<std::map<std::size_t, double>> rows(nc);
            std::vector<double> diag(nc, 0.0);
            std::vector<double> continuity(nc, 0.0);
            // The legacy path uses exactly its original double reductions;
            // only an explicit opt-in constructs extended accumulators.
            const bool compensated_assembly = controls.pressure_compensated_assembly;
            // Optional exact face connectivity for the diagnostic high/low
            // pressure Schur operator. Stored only on explicitly requested
            // validation runs; all ordinary production paths are unchanged.
            std::vector<double> pressure_face_coefficients;
            if (controls.pressure_compensated_refinement_diagnostic)
                pressure_face_coefficients.assign(mesh.n_faces(), 0.0);
            std::unique_ptr<PressureFaceAssemblyAccumulator> pressure_accumulator;
            if (compensated_assembly)
                pressure_accumulator = std::make_unique<PressureFaceAssemblyAccumulator>(nc);
            const auto add_continuity = [&](std::size_t cell, double value) {
                if (compensated_assembly)
                    pressure_accumulator->add_continuity(cell, value);
                else
                    continuity[cell] += value;
            };
            const auto current_grad_p =
                gauss_gradient_with_boundary(
            p, mesh, geometry, pressure_bcs,
            controls.pressure_skew_corrected_gradient,
            controls.pressure_projected_face_interpolation);

            const auto* cell_faces = mesh.cells().faces_data();
            const auto* cell_offsets = mesh.cells().offsets_data();
            for (std::size_t c = 0; c < nc; ++c) {
                const Offset off = cell_offsets[c];
                const Offset count = cell_offsets[c + 1] - off;
                for (Offset k = 0; k < count; ++k) {
                    const std::size_t f = cell_faces[off+k];
                    add_continuity(c, mesh.ownership().owner(f) == c
                        ? phiHbyA(f) : -phiHbyA(f));
                    if (mesh.ownership().neighbour(f) >= 0) {
                        const double phi_nonorth = rhie_chow_pressure_flux_internal(
                            mesh, geometry, f, p, current_grad_p,
                            controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                            controls.density,
                            controls.pressure_projected_face_interpolation);
                        add_continuity(c, mesh.ownership().owner(f) == c
                            ? -phi_nonorth : phi_nonorth);
                    }
                }
            }

            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const auto nr = mesh.ownership().neighbour(f);
                if (nr < 0) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const std::size_t n = static_cast<std::size_t>(nr);
                const Vec3 Sf = geometry.face_area_vectors[f];
                const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
                const auto& pressure_coeff =
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;
                const double coeff = controls.density * directional_face_coefficient(
                    pressure_coeff, f, o, n, Sf, dvec);
                if (controls.pressure_compensated_refinement_diagnostic)
                    pressure_face_coefficients[f] = coeff;
                if (compensated_assembly)
                    pressure_accumulator->add_internal_pair(o, n, coeff);
                else {
                    diag[o] += coeff;
                    diag[n] += coeff;
                }
                rows[o][n] -= coeff;
                rows[n][o] -= coeff;
            }

            for (std::size_t c = 0; c < nc; ++c) {
                if (!compensated_assembly) {
                    rows[c][c] += diag[c];
                    b(c) = -continuity[c];
                }
            }

            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                if (mesh.ownership().neighbour(f) >= 0) continue;
                const std::size_t patch = geometry.face_patch[f];
                if (patch >= mesh.boundary().n_patches()) continue;
                const auto& name = mesh.boundary().patch(patch).name;
                const auto it = pressure_bcs.find(name);
                if (it == pressure_bcs.end() ||
                    it->second.type != ScalarBoundaryType::FIXED_VALUE) continue;
                const std::size_t o = mesh.ownership().owner(f);
                const Vec3 dvec = geometry.face_centres[f] - geometry.cell_centres[o];
                const Vec3 Sf = geometry.face_area_vectors[f];
                const auto& pressure_coeff =
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;
                const Vec3 dAU_boundary{
                    geometry.cell_volumes[o] * pressure_coeff[0][o],
                    geometry.cell_volumes[o] * pressure_coeff[1][o],
                    geometry.cell_volumes[o] * pressure_coeff[2][o]};
                const double boundary_coeff = pressure_face_orthogonal_coefficient(
                    Sf, dvec, dAU_boundary);
                if (controls.pressure_compensated_refinement_diagnostic)
                    pressure_face_coefficients[f] = controls.density * boundary_coeff;
                if (!compensated_assembly)
                    rows[o][o] += controls.density * boundary_coeff;
                // Explicit pressure flux at fixed-pressure boundaries (#894):
                // the half-Dirichlet diagonal term above is only consistent
                // when the boundary mass flux is reconstructed as
                //   phi_b = phiHbyA_b - rho*boundary_coeff*(p_bc - p_P).
                // boundary_coeff uses the same centre-line projection as the
                // boundary flux reconstruction below. Without the matching
                // source term the discretised equation forces the boundary
                // flux to zero and continuity stalls on a nonzero floor.
                const Vec3 gp_boundary{
                    current_grad_p.component_data(0)[o],
                    current_grad_p.component_data(1)[o],
                    current_grad_p.component_data(2)[o]};
                const double deferred_boundary = pressure_face_nonorthogonal_flux(
                    Sf, dvec, dAU_boundary, gp_boundary);
                const double boundary_source = controls.density * (
                    boundary_coeff * (it->second.value - p(o)) +
                    deferred_boundary);
                if (compensated_assembly)
                    pressure_accumulator->add_pressure_boundary(
                        o, controls.density * boundary_coeff, boundary_source);
                else
                    b(o) += boundary_source;
            }
            if (compensated_assembly) {
                for (std::size_t c = 0; c < nc; ++c) {
                    rows[c][c] += pressure_accumulator->diagonal(c);
                    b(c) = pressure_accumulator->rhs(c);
                }
            }

            // Read-only decomposition of the exact pressure RHS. The
            // production continuity assembly above is deliberately left
            // untouched; this second pass attributes its terms per cell/face.
            if (controls.momentum_assembly_forensic_callback) {
                std::vector<double> rhs_predictor(nc, 0.0);
                std::vector<double> rhs_pressure_implicit(nc, 0.0);
                std::vector<double> rhs_pressure_deferred(nc, 0.0);
                std::vector<double> rhs_boundary_source(nc, 0.0);
                double predictor_internal_linf = 0.0;
                double predictor_boundary_linf = 0.0;
                double pressure_implicit_linf = 0.0;
                double pressure_deferred_linf = 0.0;
                double boundary_source_linf = 0.0;
                double pressure_coefficient_min = std::numeric_limits<double>::infinity();
                double pressure_coefficient_max = 0.0;
                std::size_t predictor_internal_face = 0;
                std::size_t predictor_boundary_face = 0;
                std::size_t pressure_implicit_face = 0;
                std::size_t pressure_deferred_face = 0;
                std::size_t boundary_source_face = 0;
                const auto& pressure_coeff =
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;

                for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                    const std::size_t o = mesh.ownership().owner(f);
                    const auto nr = mesh.ownership().neighbour(f);
                    const double predictor = phiHbyA(f);
                    rhs_predictor[o] -= predictor;
                    if (nr >= 0) {
                        const std::size_t n = static_cast<std::size_t>(nr);
                        rhs_predictor[n] += predictor;
                        predictor_internal_linf = std::max(
                            predictor_internal_linf, std::abs(predictor));
                        if (std::abs(predictor) == predictor_internal_linf)
                            predictor_internal_face = f;

                        const Vec3 Sf = geometry.face_area_vectors[f];
                        const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
                        const double coeff = controls.density *
                            directional_face_coefficient(pressure_coeff, f, o, n, Sf, dvec);
                        // The forensic RHS audit must use the same projected
                        // dAU/gradient as the production Schur complement,
                        // or a diagnostic callback can report a false mismatch.
                        const double wn = pressure_face_neighbour_weight(
                            geometry, f, o, n,
                            controls.pressure_projected_face_interpolation);
                        const double wo = 1.0 - wn;
                        const Vec3 dAU{
                            wo * geometry.cell_volumes[o] * pressure_coeff[0][o] +
                            wn * geometry.cell_volumes[n] * pressure_coeff[0][n],
                            wo * geometry.cell_volumes[o] * pressure_coeff[1][o] +
                            wn * geometry.cell_volumes[n] * pressure_coeff[1][n],
                            wo * geometry.cell_volumes[o] * pressure_coeff[2][o] +
                            wn * geometry.cell_volumes[n] * pressure_coeff[2][n]};
                        const Vec3 grad_face{
                            wo * current_grad_p.component_data(0)[o] +
                            wn * current_grad_p.component_data(0)[n],
                            wo * current_grad_p.component_data(1)[o] +
                            wn * current_grad_p.component_data(1)[n],
                            wo * current_grad_p.component_data(2)[o] +
                            wn * current_grad_p.component_data(2)[n]};
                        const double implicit_flux = coeff * (p(n) - p(o));
                        const double deferred_flux = controls.density *
                            pressure_face_nonorthogonal_flux(Sf, dvec, dAU, grad_face);
                        rhs_pressure_implicit[o] += implicit_flux;
                        rhs_pressure_implicit[n] -= implicit_flux;
                        rhs_pressure_deferred[o] += deferred_flux;
                        rhs_pressure_deferred[n] -= deferred_flux;
                        if (std::abs(implicit_flux) > pressure_implicit_linf) {
                            pressure_implicit_linf = std::abs(implicit_flux);
                            pressure_implicit_face = f;
                        }
                        if (std::abs(deferred_flux) > pressure_deferred_linf) {
                            pressure_deferred_linf = std::abs(deferred_flux);
                            pressure_deferred_face = f;
                        }
                        pressure_coefficient_min = std::min(
                            pressure_coefficient_min, coeff);
                        pressure_coefficient_max = std::max(
                            pressure_coefficient_max, std::abs(coeff));
                    } else {
                        predictor_boundary_linf = std::max(
                            predictor_boundary_linf, std::abs(predictor));
                        if (std::abs(predictor) == predictor_boundary_linf)
                            predictor_boundary_face = f;
                        const std::size_t patch = geometry.face_patch[f];
                        if (patch >= mesh.boundary().n_patches()) continue;
                        const auto& name = mesh.boundary().patch(patch).name;
                        const auto pbc = pressure_bcs.find(name);
                        if (pbc == pressure_bcs.end() ||
                            pbc->second.type != ScalarBoundaryType::FIXED_VALUE)
                            continue;
                        const Vec3 dvec =
                            geometry.face_centres[f] - geometry.cell_centres[o];
                        const Vec3 Sf = geometry.face_area_vectors[f];
                        const Vec3 dAU{
                            geometry.cell_volumes[o] * pressure_coeff[0][o],
                            geometry.cell_volumes[o] * pressure_coeff[1][o],
                            geometry.cell_volumes[o] * pressure_coeff[2][o]};
                        const double coeff = controls.density *
                            pressure_face_orthogonal_coefficient(Sf, dvec, dAU);
                        const Vec3 grad_boundary{
                            current_grad_p.component_data(0)[o],
                            current_grad_p.component_data(1)[o],
                            current_grad_p.component_data(2)[o]};
                        const double deferred = controls.density *
                            pressure_face_nonorthogonal_flux(
                                Sf, dvec, dAU, grad_boundary);
                        const double source = coeff * (pbc->second.value - p(o)) +
                                              deferred;
                        rhs_boundary_source[o] += source;
                        if (std::abs(source) > boundary_source_linf) {
                            boundary_source_linf = std::abs(source);
                            boundary_source_face = f;
                        }
                    }
                }

                std::size_t rhs_cell = 0;
                double rhs_linf = 0.0;
                double continuity_linf = 0.0;
                double decomposition_error = 0.0;
                for (std::size_t c = 0; c < nc; ++c) {
                    if (std::abs(b(c)) > rhs_linf) {
                        rhs_linf = std::abs(b(c));
                        rhs_cell = c;
                    }
                    continuity_linf = std::max(continuity_linf, std::abs(continuity[c]));
                    const double reconstructed = rhs_predictor[c] +
                        rhs_pressure_implicit[c] + rhs_pressure_deferred[c] +
                        rhs_boundary_source[c];
                    decomposition_error = std::max(
                        decomposition_error, std::abs(b(c) - reconstructed));
                }
                const auto& rhs_point = geometry.cell_centres[rhs_cell];
                std::cerr << std::setprecision(17)
                          << "NACA_PRESSURE_RHS"
                          << " iter=" << iter
                          << " corrector=" << corr
                          << " rhs_linf=" << rhs_linf
                          << " rhs_cell=" << rhs_cell
                          << " rhs_xyz=" << rhs_point.x << ":" << rhs_point.y << ":" << rhs_point.z
                          << " continuity_linf=" << continuity_linf
                          << " predictor_rhs_cell=" << rhs_predictor[rhs_cell]
                          << " implicit_pressure_rhs_cell=" << rhs_pressure_implicit[rhs_cell]
                          << " deferred_pressure_rhs_cell=" << rhs_pressure_deferred[rhs_cell]
                          << " boundary_source_rhs_cell=" << rhs_boundary_source[rhs_cell]
                          << " decomposition_error_linf=" << decomposition_error
                          << " predictor_internal_linf=" << predictor_internal_linf
                          << " predictor_internal_face=" << predictor_internal_face
                          << " predictor_boundary_linf=" << predictor_boundary_linf
                          << " predictor_boundary_face=" << predictor_boundary_face
                          << " pressure_implicit_linf=" << pressure_implicit_linf
                          << " pressure_implicit_face=" << pressure_implicit_face
                          << " pressure_deferred_linf=" << pressure_deferred_linf
                          << " pressure_deferred_face=" << pressure_deferred_face
                          << " boundary_source_linf=" << boundary_source_linf
                          << " boundary_source_face=" << boundary_source_face
                          << " pressure_coeff_min=" << pressure_coefficient_min
                          << " pressure_coeff_abs_max=" << pressure_coefficient_max
                          << '\n';

                if (iter == 1 || iter % 10 == 0) {
                    const auto* cell_faces = mesh.cells().faces_data();
                    const auto* cell_offsets = mesh.cells().offsets_data();
                    const Offset off = cell_offsets[rhs_cell];
                    const Offset count = cell_offsets[rhs_cell + 1] - off;
                    for (Offset k = 0; k < count; ++k) {
                        const std::size_t f = cell_faces[off + k];
                        const std::size_t o = mesh.ownership().owner(f);
                        const auto nr = mesh.ownership().neighbour(f);
                        const bool is_owner = o == rhs_cell;
                        const double orientation = is_owner ? 1.0 : -1.0;
                        const double rhs_phi = -orientation * phiHbyA(f);
                        double implicit_flux = 0.0;
                        double deferred_flux = 0.0;
                        double boundary_source = 0.0;
                        double coeff = 0.0;
                        std::string patch_name = "internal";
                        if (nr >= 0) {
                            const std::size_t n = static_cast<std::size_t>(nr);
                            const Vec3 Sf = geometry.face_area_vectors[f];
                            const Vec3 dvec = geometry.cell_centres[n] - geometry.cell_centres[o];
                            const Vec3 dAU{
                                0.5 * (geometry.cell_volumes[o] * pressure_coeff[0][o] +
                                       geometry.cell_volumes[n] * pressure_coeff[0][n]),
                                0.5 * (geometry.cell_volumes[o] * pressure_coeff[1][o] +
                                       geometry.cell_volumes[n] * pressure_coeff[1][n]),
                                0.5 * (geometry.cell_volumes[o] * pressure_coeff[2][o] +
                                       geometry.cell_volumes[n] * pressure_coeff[2][n])};
                            const Vec3 grad_face{
                                0.5 * (current_grad_p.component_data(0)[o] +
                                       current_grad_p.component_data(0)[n]),
                                0.5 * (current_grad_p.component_data(1)[o] +
                                       current_grad_p.component_data(1)[n]),
                                0.5 * (current_grad_p.component_data(2)[o] +
                                       current_grad_p.component_data(2)[n])};
                            coeff = controls.density *
                                pressure_face_orthogonal_coefficient(Sf, dvec, dAU);
                            implicit_flux = coeff * (p(n) - p(o));
                            deferred_flux = controls.density *
                                pressure_face_nonorthogonal_flux(Sf, dvec, dAU, grad_face);
                        } else {
                            const std::size_t patch = geometry.face_patch[f];
                            if (patch < mesh.boundary().n_patches()) {
                                patch_name = mesh.boundary().patch(patch).name;
                                const auto pbc = pressure_bcs.find(patch_name);
                                if (pbc != pressure_bcs.end() &&
                                    pbc->second.type == ScalarBoundaryType::FIXED_VALUE) {
                                    const Vec3 dvec =
                                        geometry.face_centres[f] - geometry.cell_centres[o];
                                    const Vec3 Sf = geometry.face_area_vectors[f];
                                    const Vec3 dAU{
                                        geometry.cell_volumes[o] * pressure_coeff[0][o],
                                        geometry.cell_volumes[o] * pressure_coeff[1][o],
                                        geometry.cell_volumes[o] * pressure_coeff[2][o]};
                                    coeff = controls.density *
                                        pressure_face_orthogonal_coefficient(Sf, dvec, dAU);
                                    const Vec3 grad_boundary{
                                        current_grad_p.component_data(0)[o],
                                        current_grad_p.component_data(1)[o],
                                        current_grad_p.component_data(2)[o]};
                                    deferred_flux = controls.density *
                                        pressure_face_nonorthogonal_flux(
                                            Sf, dvec, dAU, grad_boundary);
                                    boundary_source = coeff * (pbc->second.value - p(o)) +
                                                      deferred_flux;
                                }
                            }
                        }
                        std::cerr << std::setprecision(17)
                                  << "NACA_PRESSURE_FACE"
                                  << " iter=" << iter
                                  << " corrector=" << corr
                                  << " cell=" << rhs_cell
                                  << " face=" << f
                                  << " owner=" << o
                                  << " neighbour=" << (nr >= 0 ? std::to_string(nr) : "boundary")
                                  << " patch=" << patch_name
                                  << " phiHbyA=" << phiHbyA(f)
                                  << " rhs_phi=" << rhs_phi
                                  << " implicit_pressure_flux=" << implicit_flux
                                  << " rhs_implicit=" << orientation * implicit_flux
                                  << " deferred_pressure_flux=" << deferred_flux
                                  << " rhs_deferred=" << orientation * deferred_flux
                                  << " boundary_pressure_source=" << boundary_source
                                  << " pressure_coeff=" << coeff
                                  << '\n';
                    }
                }
            }

            // The V6 canonical pressure state includes the LOW sidecar of the
            // previously accepted pressure. Keep its full face flux as the
            // single residual operator for the new pressure increment.
            std::vector<long double> canonical_base_flux;
            if (controls.pressure_high_low_production) {
                canonical_base_flux=make_rhie_chow_mass_flux_high_low(
                    mesh,geometry,HbyA,p,pressure_low_state,rAU,controls.density,
                    velocity_bcs,pressure_bcs,
                    controls.pressure_projected_face_interpolation,
                    face_velocity_conditions ? &*face_velocity_conditions : nullptr);
            }
            // Experimental reference matrix includes the derivative of every
            // pressure-dependent deferred Rhie-Chow gradient. Retain the same
            // physical RHS: minus divergence of the authoritative flux.
            if (controls.v18_full_flux_pressure_oracle ||
                controls.v18_sparse_flux_pressure_jacobian) {
                if (!has_fixed_pressure_boundary)
                    throw std::invalid_argument(
                        "V18 pressure flux oracle currently requires fixed-pressure BC");
                if (controls.v18_full_flux_pressure_oracle)
                    rows=make_full_rhie_chow_pressure_rows_reference(
                        mesh,geometry,rAU,controls.density,velocity_bcs,
                        pressure_bcs,controls.pressure_projected_face_interpolation);
                else
                    rows=make_full_rhie_chow_pressure_rows_sparse(
                        mesh,geometry,rAU,controls.density,pressure_bcs,
                        controls.pressure_projected_face_interpolation);
                auto base_flux=make_rhie_chow_mass_flux(
                    mesh,geometry,HbyA,p,rAU,controls.density,
                    velocity_bcs,pressure_bcs,
                    RhieChowVelocityContract::HBYA_PREDICTOR,false,
                    controls.pressure_projected_face_interpolation);
                const auto audit=audit_v17_authoritative_flux(mesh,geometry,base_flux);
                double rhs_error=0.0;
                for(std::size_t c=0;c<nc;++c) {
                    rhs_error+=std::abs(b(c)+audit.cell_divergence[c]);
                    b(c)=-audit.cell_divergence[c];
                }
                if (rhs_error>1e-9*std::max(1.0,audit.cell_l1))
                    throw std::runtime_error("V18 full-flux pressure RHS inconsistent");
                if (controls.pressure_high_low_production) {
                    // Replace the FP64 RHS by the canonical split-pressure
                    // minus divergence. The initial p_low=0 agrees with V18,
                    // subsequent iterations must not discard its history.
                    std::vector<long double> divergence(nc,0.0L);
                    for(std::size_t f=0;f<mesh.n_faces();++f) {
                        const auto o=mesh.ownership().owner(f);
                        divergence[o]+=canonical_base_flux[f];
                        const auto n=mesh.ownership().neighbour(f);
                        if(n>=0) divergence[static_cast<std::size_t>(n)]-=canonical_base_flux[f];
                    }
                    for(std::size_t c=0;c<nc;++c) {
                        const long double r=-divergence[c];
                        if(!std::isfinite(r) || std::abs(r)>std::numeric_limits<double>::max())
                            throw std::domain_error("RC high/low: invalid base divergence");
                        b(c)=static_cast<double>(r);
                    }
                }
                std::cerr << std::scientific << std::setprecision(12)
                          << "V18_FULL_FLUX_OPERATOR mode="
                          << (controls.v18_full_flux_pressure_oracle ? "reference" : "sparse")
                          << " iteration=" << iter
                          << " corrector=" << corr
                          << " cells=" << nc
                          << " rhs_l1=" << audit.cell_l1
                          << " rhs_contract_error=" << rhs_error
                          << " qualified=0\n";
            }
            Vector p_corr(nc, 0.0);
            cfdx::core::SolverResult rp;
            const auto solve_pressure_system =
                [&](const SparseMatrix& matrix,
                    const Vector& rhs,
                    Vector& correction) {
                    // Capture selected difficult *actual* post-gauge pressure
                    // systems, including their original RHS and initial guess.
                    // This hook also covers reusable CG, which bypasses dispatch.
                    if (!controls.pressure_replay_prefix.empty() &&
                        std::find(controls.pressure_replay_iterations.begin(),
                                  controls.pressure_replay_iterations.end(), iter) !=
                        controls.pressure_replay_iterations.end()) {
                        const auto path = controls.pressure_replay_prefix +
                            "_i" + std::to_string(iter) +
                            "_c" + std::to_string(corr) + ".cfdxlin";
                        cfdx::core::write_linear_system_capture(
                            path, matrix, rhs, correction,
                            cfdx::core::LinearProblemKind::PressurePoisson,
                            controls.pressure_linear_solver.null_space,
                            controls.pressure_linear_solver.scaling,
                            controls.pressure_linear_max_iterations != 0
                                ? controls.pressure_linear_max_iterations
                                : controls.linear_max_iterations,
                            controls.pressure_linear_tolerance > 0.0
                                ? controls.pressure_linear_tolerance
                                : linear_tolerance);
                    }
                    // The reusable pressure-CG context bypasses generic dispatch.
                    // Capture the first ACTUAL pressure matrix before the solve
                    // (including the post-gauge reduction and true initial guess).
                    // No overhead and no file creation without explicit opt-in.
                    if (iter == 1 && corr == 0 &&
                        !controls.pressure_linear_solver.linear_capture_path.empty()) {
                        cfdx::core::write_linear_system_capture(
                            controls.pressure_linear_solver.linear_capture_path,
                            matrix, rhs, correction,
                            cfdx::core::LinearProblemKind::PressurePoisson,
                            controls.pressure_linear_solver.null_space,
                            controls.pressure_linear_solver.scaling,
                            controls.pressure_linear_max_iterations != 0
                                ? controls.pressure_linear_max_iterations
                                : controls.linear_max_iterations,
                            controls.pressure_linear_tolerance > 0.0
                                ? controls.pressure_linear_tolerance
                                : linear_tolerance);
                    }
                    if (controls.v18_full_flux_pressure_oracle ||
                        controls.v18_sparse_flux_pressure_jacobian) {
                        auto request=controls.pressure_linear_solver;
                        request.krylov=cfdx::core::KrylovModel::FGMRES;
                        if (request.preconditioner != cfdx::core::PreconditionerModel::AnisotropicLineAMG)
                            request.preconditioner=cfdx::core::PreconditionerModel::ILU0;
                        request.gmres_restart=std::min<int>(512,static_cast<int>(nc));
                        request.null_space=cfdx::core::NullSpaceModel::None;
                        return solve_linear_system(
                            matrix,rhs,correction,LinearProblemKind::General,
                            request,
                            controls.pressure_linear_max_iterations!=0
                                ? controls.pressure_linear_max_iterations
                                : controls.linear_max_iterations,
                            controls.pressure_linear_tolerance>0.0
                                ? controls.pressure_linear_tolerance
                                : linear_tolerance).result;
                    }
                    if (pressure_context) {
                        return solve_reusable_scaled_cg(
                            matrix, rhs, correction, *pressure_context,
                            controls.pressure_linear_solver.scaling,
                            controls.pressure_linear_max_iterations != 0
                                ? controls.pressure_linear_max_iterations
                                : controls.linear_max_iterations,
                            controls.pressure_linear_tolerance > 0.0
                                ? controls.pressure_linear_tolerance
                                : linear_tolerance,
                            pressure_null_space.get());
                    }
                    return solve_linear_system(
                        matrix, rhs, correction,
                        LinearProblemKind::PressurePoisson,
                        controls.pressure_linear_solver,
                        controls.pressure_linear_max_iterations != 0
                            ? controls.pressure_linear_max_iterations
                            : controls.linear_max_iterations,
                        controls.pressure_linear_tolerance > 0.0
                            ? controls.pressure_linear_tolerance
                            : linear_tolerance).result;
                };
            if (std::getenv("CFDX_BFS_STAGE_TRACE"))
                std::cerr << "BFS_PHASE_V10 iteration=" << iter
                          << " corrector=" << corr
                          << " phase=PRESSURE_LINEAR_BEGIN\n";
            std::unique_ptr<SparseMatrix> reduced_pressure_matrix;
            std::unique_ptr<Vector> reduced_pressure_rhs;
            const SparseMatrix* solved_pressure_matrix = &A;
            const Vector* solved_pressure_rhs = &b;
            if (has_fixed_pressure_boundary ||
                use_projected_pressure_null_space) {
                // With a prescribed pressure boundary the pressure-correction
                // operator is SPD. For pure Neumann pressure, an explicitly
                // requested constant null-space model instead keeps the
                // singular operator and projects the CG/GMRES Krylov space.
                for (std::size_t row = 0; row < nc; ++row)
                    for (const auto& [col, value] : rows[row])
                        A.push_back(row, col, value);
                A.finalize();
                rp = solve_pressure_system(A, b, p_corr);
            } else if (nc == 1) {
                // Eliminating the only pressure unknown would produce an empty
                // linear system. Keep the trivial gauge equation instead.
                A.push_back(0, 0, 1.0);
                A.finalize();
                b(0) = 0.0;
                rp = solve_pressure_system(A, b, p_corr);
            } else {
                // Pure-Neumann pressure has one gauge null mode. Replacing one
                // row by p_ref=0 breaks symmetry. Eliminate the reference DOF
                // instead: the remaining principal submatrix is SPD for a
                // connected cavity and remains eligible for the pressure
                // solver policy and its true-residual convergence checks.
                const std::size_t ref = controls.pressure_reference_cell;
                reduced_pressure_matrix =
                    std::make_unique<SparseMatrix>(nc - 1, nc - 1);
                reduced_pressure_rhs =
                    std::make_unique<Vector>(nc - 1, 0.0);
                const auto reduced_index = [ref](std::size_t cell) {
                    return cell < ref ? cell : cell - 1;
                };
                for (std::size_t row = 0; row < nc; ++row) {
                    if (row == ref) continue;
                    const std::size_t rr = reduced_index(row);
                    (*reduced_pressure_rhs)(rr) = b(row);
                    for (const auto& [col, value] : rows[row]) {
                        if (col == ref) continue;
                        reduced_pressure_matrix->push_back(
                            rr, reduced_index(col), value);
                    }
                }
                reduced_pressure_matrix->finalize();
                solved_pressure_matrix = reduced_pressure_matrix.get();
                solved_pressure_rhs = reduced_pressure_rhs.get();

                Vector reduced_corr(nc - 1, 0.0);
                rp = solve_pressure_system(
                    *solved_pressure_matrix, *solved_pressure_rhs, reduced_corr);
                if (rp.status == cfdx::core::SolverStatus::CONVERGED) {
                    for (std::size_t cell = 0; cell < nc; ++cell)
                        if (cell != ref)
                            p_corr(cell) = reduced_corr(reduced_index(cell));
                    p_corr(ref) = 0.0;
                }
            }
            if (std::getenv("CFDX_BFS_STAGE_TRACE"))
                std::cerr << "BFS_PHASE_V10 iteration=" << iter
                          << " corrector=" << corr
                          << " phase=PRESSURE_LINEAR_END status="
                          << static_cast<int>(rp.status)
                          << " iterations=" << rp.iterations
                          << " true_relative=" << rp.residual_relative
                          << "\n";
            // The original V18 sparse matrix is an inexact Jacobian for the
            // REAL full Rhie-Chow face residual.  Neither CG/FGMRES status nor
            // an algebraic CSR residual can promote a pressure correction.
            Vector p_corr_low(nc,0.0);
            if (controls.pressure_high_low_production) {
                const double csr_relative=rp.residual_relative;
                const auto physical_residual=[&](const Vector& hi,const Vector& lo) {
                    const auto flux=reconstruct_compensated_rhie_chow_pressure_step(
                        mesh,geometry,canonical_base_flux,hi,lo,rAU,
                        controls.density,pressure_bcs,
                        controls.pressure_projected_face_interpolation);
                    return pressure_true_face_residual(mesh,flux,b);
                };
                const auto inner_direction=[&](const Vector& defect,Vector& delta) {
                    const auto inner=solve_pressure_system(A,defect,delta);
                    // Even a Krylov MAX_ITER direction must be checked on the
                    // true FV faces; DIVERGED is never eligible.
                    if(inner.status==SolverStatus::DIVERGED ||
                       inner.status==SolverStatus::NOT_APPLICABLE) return false;
                    for(std::size_t c=0;c<delta.size();++c)
                        if(!std::isfinite(delta(c))) return false;
                    return true;
                };
                const double physical_pressure_tolerance=
                    controls.pressure_linear_tolerance>0.0
                    ? controls.pressure_linear_tolerance : linear_tolerance;
                const auto true_report=refine_pressure_true_face_high_low(
                    p_corr,p_corr_low,physical_pressure_tolerance,12,
                    physical_residual,inner_direction);
                std::cerr << std::scientific << std::setprecision(16)
                          << "PRESSURE_HILO_PRODUCTION iteration=" << iter
                          << " corrector=" << corr
                          << " original_csr_relative=" << csr_relative
                          << " initial_face_relative=" << true_report.initial_relative
                          << " final_face_relative=" << true_report.final_relative
                          << " refinements_attempted=" << true_report.attempted
                          << " refinements_accepted=" << true_report.accepted
                          << " true_face_converged=" << true_report.converged
                          << " requested_tolerance=" << physical_pressure_tolerance
                          << " qualified_bfs=0\n";
                // A MAX_ITER algebraic solve may supply a usable direction,
                // but only a certified TRUE full face residual can promote it.
                rp.residual_relative=true_report.final_relative;
                long double b2=0.0L;
                for(std::size_t c=0;c<nc;++c)
                    b2+=static_cast<long double>(b(c))*b(c);
                rp.residual=static_cast<double>(
                    std::sqrt(b2)*true_report.final_relative);
                rp.status=true_report.converged
                    ? SolverStatus::CONVERGED : SolverStatus::MAX_ITER_REACHED;
            }
            // Optional snapshot of the ORIGINAL pressure matrix and right-hand
            // side. Independent of expensive momentum forensic logging.
            // Does not modify the operator, initial state or convergence gate.
            if (iter == 1 && corr == 0) {
                if (const char* prefix = std::getenv("CFDX_PRESSURE_DUMP_PREFIX");
                    prefix && *prefix) {
                    const auto& mat = *solved_pressure_matrix;
                    const auto& rhs = *solved_pressure_rhs;
                    std::ofstream mm(std::string(prefix) + ".mtx");
                    std::ofstream bb(std::string(prefix) + ".rhs");
                    if (!mm || !bb)
                        throw std::runtime_error("CFDX pressure dump: failed to open output");
                    mm << "%%MatrixMarket matrix coordinate real general\n";
                    mm << mat.n_rows() << " " << mat.n_cols() << " " << mat.nnz() << "\n";
                    mm << std::setprecision(17);
                    bb << std::setprecision(17);

                    const auto* offsets = mat.row_offsets_data();
                    const auto* columns = mat.columns_data();
                    const auto* values = mat.values_data();
                    for (std::size_t row = 0; row < mat.n_rows(); ++row) {
                        bb << rhs(row) << "\n";
                        for (std::uint32_t k = offsets[row]; k < offsets[row+1]; ++k)
                            mm << row + 1 << " " << columns[k] + 1
                               << " " << values[k] << "\n";
                    }
                }
            }
            const auto pressure_row_to_cell = [&](std::size_t row) {
                if (solved_pressure_matrix->n_rows() == nc) return row;
                const std::size_t ref = controls.pressure_reference_cell;
                return row < ref ? row : row + 1;
            };
            if (controls.momentum_assembly_forensic_callback) {
                const auto& matrix = *solved_pressure_matrix;
                const auto& rhs = *solved_pressure_rhs;
                double correction_linf = 0.0;
                double pressure_linf = 0.0;
                double rhs_linf = 0.0;
                double row_abs_sum_linf = 0.0;
                std::size_t correction_cell = 0;
                std::size_t pressure_cell = 0;
                std::size_t rhs_cell = 0;
                std::size_t diagonal_cell = 0;
                std::size_t minimum_diagonal_cell = 0;
                std::size_t rhs_diagonal_ratio_cell = 0;
                std::size_t minimum_dominance_cell = 0;
                double diagonal_min = std::numeric_limits<double>::infinity();
                double diagonal_max = 0.0;
                double minimum_diagonal_rhs = 0.0;
                double minimum_diagonal_correction = 0.0;
                double rhs_diagonal_ratio_linf = 0.0;
                double diagonal_dominance_min = std::numeric_limits<double>::infinity();
                double row_sum_linf = 0.0;
                const auto* offsets = matrix.row_offsets_data();
                const auto* columns = matrix.columns_data();
                const auto* values = matrix.values_data();
                for (std::size_t cell = 0; cell < nc; ++cell) {
                    if (std::abs(p_corr(cell)) > correction_linf) {
                        correction_linf = std::abs(p_corr(cell));
                        correction_cell = cell;
                    }
                    if (std::abs(p(cell)) > pressure_linf) {
                        pressure_linf = std::abs(p(cell));
                        pressure_cell = cell;
                    }
                }
                for (std::size_t row = 0; row < rhs.size(); ++row) {
                    if (std::abs(rhs(row)) > rhs_linf) {
                        rhs_linf = std::abs(rhs(row));
                        rhs_cell = pressure_row_to_cell(row);
                    }
                    double diagonal = 0.0;
                    double row_sum = 0.0;
                    double row_abs_sum = 0.0;
                    double off_diagonal_abs_sum = 0.0;
                    for (std::uint32_t k = offsets[row]; k < offsets[row + 1]; ++k) {
                        row_sum += values[k];
                        row_abs_sum += std::abs(values[k]);
                        if (columns[k] == row) {
                            diagonal = values[k];
                        } else {
                            off_diagonal_abs_sum += std::abs(values[k]);
                        }
                    }
                    if (diagonal < diagonal_min) {
                        diagonal_min = diagonal;
                        minimum_diagonal_cell = pressure_row_to_cell(row);
                        minimum_diagonal_rhs = rhs(row);
                        minimum_diagonal_correction =
                            p_corr(minimum_diagonal_cell);
                    }
                    if (std::abs(diagonal) > diagonal_max) {
                        diagonal_max = std::abs(diagonal);
                        diagonal_cell = pressure_row_to_cell(row);
                    }
                    if (std::abs(diagonal) > 0.0 &&
                        std::isfinite(rhs(row))) {
                        const double ratio = std::abs(rhs(row)) / std::abs(diagonal);
                        if (ratio > rhs_diagonal_ratio_linf) {
                            rhs_diagonal_ratio_linf = ratio;
                            rhs_diagonal_ratio_cell = pressure_row_to_cell(row);
                        }
                    }
                    const double dominance = std::abs(diagonal) /
                        std::max(off_diagonal_abs_sum,
                                 std::numeric_limits<double>::min());
                    if (dominance < diagonal_dominance_min) {
                        diagonal_dominance_min = dominance;
                        minimum_dominance_cell = pressure_row_to_cell(row);
                    }
                    row_sum_linf = std::max(row_sum_linf, std::abs(row_sum));
                    row_abs_sum_linf = std::max(row_abs_sum_linf, row_abs_sum);
                }
                std::cerr << std::setprecision(17)
                          << "NACA_PRESSURE_CORRECTION"
                          << " iter=" << iter
                          << " corrector=" << corr
                          << " status=" << static_cast<int>(rp.status)
                          << " iterations=" << rp.iterations
                          << " residual=" << rp.residual
                          << " relative=" << rp.residual_relative
                          << " rows=" << matrix.n_rows()
                          << " nnz=" << matrix.nnz()
                          << " rhs_linf=" << rhs_linf
                          << " rhs_cell=" << rhs_cell
                          << " correction_linf=" << correction_linf
                          << " correction_cell=" << correction_cell
                          << " pressure_before_linf=" << pressure_linf
                          << " pressure_cell=" << pressure_cell
                          << " diagonal_min=" << diagonal_min
                          << " diagonal_abs_max=" << diagonal_max
                          << " diagonal_cell=" << diagonal_cell
                          << " minimum_diagonal_cell=" << minimum_diagonal_cell
                          << " minimum_diagonal_rhs=" << minimum_diagonal_rhs
                          << " minimum_diagonal_correction=" << minimum_diagonal_correction
                          << " rhs_diagonal_ratio_linf=" << rhs_diagonal_ratio_linf
                          << " rhs_diagonal_ratio_cell=" << rhs_diagonal_ratio_cell
                          << " diagonal_dominance_min=" << diagonal_dominance_min
                          << " minimum_dominance_cell=" << minimum_dominance_cell
                          << " row_sum_linf=" << row_sum_linf
                          << " row_abs_sum_linf=" << row_abs_sum_linf
                          << '\n';

                // Optional Matrix-Market snapshot of the exact matrix and
                // RHS passed to Krylov. No dump is written unless requested.
                const char* dump_prefix = std::getenv("CFDX_NACA_MATRIX_DUMP_PREFIX");
                if (dump_prefix && *dump_prefix && iter == 1 && corr == 0) {
                    const std::string dump_stem =
                        std::string(dump_prefix) + ".cells" + std::to_string(nc);
                    std::ofstream matrix_file(dump_stem + ".mtx");
                    std::ofstream rhs_file(dump_stem + ".rhs.csv");
                    if (matrix_file && rhs_file) {
                        matrix_file << "%%MatrixMarket matrix coordinate real general\n";
                        matrix_file << matrix.n_rows() << " " << matrix.n_rows()
                                    << " " << matrix.nnz() << "\n";
                        for (std::size_t row = 0; row < matrix.n_rows(); ++row) {
                            for (std::uint32_t k = offsets[row]; k < offsets[row + 1]; ++k)
                                matrix_file << row + 1 << " " << columns[k] + 1
                                            << " " << std::setprecision(17) << values[k] << "\n";
                        }
                        rhs_file << "row,cell,rhs\n";
                        for (std::size_t row = 0; row < rhs.size(); ++row)
                            rhs_file << row << "," << pressure_row_to_cell(row) << ","
                                     << std::setprecision(17) << rhs(row) << "\n";
                        std::cerr << "NACA_PRESSURE_MATRIX_DUMP prefix="
                                  << dump_stem << " rows=" << matrix.n_rows()
                                  << " nnz=" << matrix.nnz() << '\n';
                    } else {
                        std::cerr << "NACA_PRESSURE_MATRIX_DUMP_ERROR prefix="
                                  << dump_prefix << '\n';
                    }
                }
            }
            pressure_residual = rp.residual_relative;
            pressure_iterations = rp.iterations;
            // Forensic only: an SPD pressure operator can be resolved to
            // backward-error O(eps) yet still violate cellwise continuity.
            // Measure this explicitly; NEVER bypass rp.status with it.
            if (std::getenv("CFDX_PRESSURE_ROUNDOFF_DIAGNOSTICS") &&
                rp.status != cfdx::core::SolverStatus::CONVERGED &&
                has_fixed_pressure_boundary &&
                solved_pressure_matrix->n_rows() == mesh.n_cells()) {
                std::vector<double> local_scale(mesh.n_cells(), 0.0);
                const auto* faces = mesh.cells().faces_data();
                const auto* offsets = mesh.cells().offsets_data();
                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                    for (cfdx::core::Offset k = offsets[cell];
                         k < offsets[cell + 1]; ++k)
                        local_scale[cell] += std::abs(phiHbyA(faces[k]));
                const auto diagnostic = diagnose_pressure_linear_cancellation(
                    *solved_pressure_matrix, *solved_pressure_rhs, p_corr,
                    local_scale, controls.convergence.continuity_tolerance);
                const auto leakage = cfdx::core::diagnose_pressure_roundoff_leakage(
                    *solved_pressure_matrix, *solved_pressure_rhs, p_corr);
                std::cerr << std::scientific << std::setprecision(9)
                          << "PRESSURE_GHOST_ROW_SUM_DIAGNOSTIC"
                          << " iteration=" << iter
                          << " corrector=" << corr
                          << " valid=" << leakage.valid
                          << " anchor_rows=" << leakage.pressure_anchor_rows
                          << " roundoff_rows=" << leakage.negligible_roundoff_rows
                          << " worst_row=" << leakage.worst_row
                          << " spurious_relative_l2=" << leakage.spurious_leakage_relative_l2
                          << " physical_face_certified=0 qualified=0\n";
                std::cerr << std::scientific << std::setprecision(9)
                          << "PRESSURE_ROUNDOFF_DIAGNOSTIC iteration=" << iter
                          << " corrector=" << corr
                          << " original_true_relative=" << rp.residual_relative
                          << " backward_error="
                          << diagnostic.max_componentwise_backward_error
                          << " max_local_mass_ratio="
                          << diagnostic.max_local_mass_ratio
                          << " local_gate_pass=" << diagnostic.valid
                          << " qualified=0\n";
            }
            // Opt-in high/low refinement on the EXACT assembled production
            // matrix. This is an observation only: a reduced split CSR
            // residual cannot certify the actual FP64 pressure field, its
            // Rhie-Chow face flux or any nonlinear convergence gate.
            if (controls.pressure_compensated_refinement_diagnostic &&
                rp.status != cfdx::core::SolverStatus::CONVERGED &&
                has_fixed_pressure_boundary &&
                solved_pressure_matrix->n_rows() == nc &&
                iter == 1 && corr == 0) {
                const auto& original_matrix = *solved_pressure_matrix;
                const auto& original_rhs = *solved_pressure_rhs;
                try {
                    cfdx::core::Vector high = p_corr;
                    cfdx::core::Vector low(nc, 0.0);
                    cfdx::core::NativeAMGPreconditioner amg(
                        cfdx::core::NativeAMGMethod::BoomerStyle);
                    if (!amg.setup(original_matrix)) {
                        std::cerr << "PRESSURE_HILO_LIVE status=SETUP_REJECTED"
                                  << " physical_fv_qualified=0\n";
                    } else {
                        const double final_tolerance =
                            controls.pressure_linear_tolerance > 0.0
                                ? controls.pressure_linear_tolerance
                                : linear_tolerance;
                        const auto split = cfdx::core::refine_pressure_high_low(
                            original_matrix, original_rhs, high, low,
                            final_tolerance, 2,
                            [&](const cfdx::core::SparseMatrix& mat,
                                const cfdx::core::Vector& residual,
                                cfdx::core::Vector& increment) {
                                const auto inner = cfdx::core::solve_cg_prepared(
                                    mat, residual, increment, amg, 350, 1e-7);
                                return inner.status !=
                                           cfdx::core::SolverStatus::NOT_APPLICABLE &&
                                       inner.status !=
                                           cfdx::core::SolverStatus::DIVERGED;
                            });
                        // Solve the IMPLICIT face pressure operator as a
                        // distinct object, using the exact live owner/neighbour
                        // and fixed-pressure boundary coefficients. The CSR's
                        // tiny floating row sums must not become fake face
                        // leakage on a Neumann cell. Only a complete FV
                        // Rhie-Chow/momentum update may qualify the CFD run.
                        std::vector<cfdx::core::ConservativePressureFaceOperator::Face>
                            face_entries;
                        face_entries.reserve(mesh.n_faces());
                        for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
                            const double coefficient = pressure_face_coefficients[face];
                            // Preserve EVERY native face slot, even zero-pressure
                            // BC faces, to compare against the authoritative
                            // Rhie-Chow face flux without any reindexing.
                            face_entries.push_back({
                                mesh.ownership().owner(face),
                                mesh.ownership().neighbour(face), coefficient});
                        }
                        const cfdx::core::ConservativePressureFaceOperator
                            physical_implicit_faces(nc, std::move(face_entries));
                        const auto face_report =
                            cfdx::core::refine_pressure_on_implicit_faces(
                                physical_implicit_faces, original_matrix,
                                original_rhs, high, low, final_tolerance, 2,
                                [&](const cfdx::core::SparseMatrix& mat,
                                    const cfdx::core::Vector& residual,
                                    cfdx::core::Vector& increment) {
                                    const auto inner = cfdx::core::solve_cg_prepared(
                                        mat, residual, increment, amg, 350, 1e-7);
                                    return inner.status !=
                                               cfdx::core::SolverStatus::NOT_APPLICABLE &&
                                           inner.status !=
                                               cfdx::core::SolverStatus::DIVERGED;
                                });
                        // Independent bridge: the pressure RHS must be the
                        // negative divergence of the SAME frozen Rhie-Chow
                        // flux used before the correction.  Reconstruct both
                        // the implicit FACE step and the full (nonorthogonal)
                        // pressure response from high/low.  A tiny matrix
                        // residual alone never promotes the solver status.
                        if (controls.algorithm == PressureVelocityAlgorithm::SIMPLE &&
                            !controls.pressure_skew_corrected_gradient) {
                            const auto base_rhie_flux=make_rhie_chow_mass_flux(
                                mesh,geometry,HbyA,p,rAU,controls.density,
                                velocity_bcs,pressure_bcs,
                                RhieChowVelocityContract::HBYA_PREDICTOR,false,
                                controls.pressure_projected_face_interpolation,
                                face_velocity_conditions ? &*face_velocity_conditions : nullptr);
                            std::vector<double> base_face_flux(mesh.n_faces());
                            for (std::size_t face=0;face<mesh.n_faces();++face)
                                base_face_flux[face]=base_rhie_flux(face);
                            const auto bridge=physical_implicit_faces.audit_implicit_flux_bridge(
                                original_rhs,high,low,base_face_flux,final_tolerance);
                            const auto full_compensated=reconstruct_compensated_rhie_chow_pressure_step(
                                mesh,geometry,base_rhie_flux,high,low,rAU,
                                controls.density,pressure_bcs,
                                controls.pressure_projected_face_interpolation);
                            const double full_rc_relative=compensated_flux_continuity_relative(
                                mesh,full_compensated,original_rhs);
                            const auto flux_gap=physical_implicit_faces.diagnose_full_flux_gap(
                                original_rhs,high,low,base_face_flux,full_compensated);
                            auto rounded_p=p;
                            for (std::size_t c=0;c<nc;++c) {
                                rounded_p(c)+=high(c);
                                rounded_p(c)+=low(c);
                            }
                            const auto rounded_flux=make_rhie_chow_mass_flux(
                                mesh,geometry,HbyA,rounded_p,rAU,controls.density,
                                velocity_bcs,pressure_bcs,
                                RhieChowVelocityContract::HBYA_PREDICTOR,false,
                                controls.pressure_projected_face_interpolation,
                                face_velocity_conditions ? &*face_velocity_conditions : nullptr);
                            std::vector<long double> rounded_face_flux(mesh.n_faces());
                            long double mismatch2=0.0L, rhs2=0.0L;
                            for (std::size_t f=0;f<mesh.n_faces();++f) {
                                rounded_face_flux[f]=rounded_flux(f);
                                const long double delta=full_compensated[f]-rounded_flux(f);
                                mismatch2+=delta*delta;
                            }
                            for (std::size_t c=0;c<nc;++c)
                                rhs2+=static_cast<long double>(original_rhs(c))*original_rhs(c);
                            const double rounded_rc_relative=compensated_flux_continuity_relative(
                                mesh,rounded_face_flux,original_rhs);
                            std::cerr << std::scientific << std::setprecision(16)
                                      << "PRESSURE_HILO_RHIE_CHOW_BRIDGE"
                                      << " iter=" << iter << " corrector=" << corr
                                      << " rhs_base_flux_relative="
                                      << bridge.rhs_base_flux_relative
                                      << " implicit_linearized_continuity_relative="
                                      << bridge.implicit_flux_continuity_relative
                                      << " implicit_operator_relative="
                                      << bridge.implicit_operator_residual_relative
                                      << " full_rc_compensated_continuity_relative="
                                      << full_rc_relative
                                      << " full_rc_fp64_continuity_relative="
                                      << rounded_rc_relative
                                      << " deferred_divergence_relative="
                                      << flux_gap.deferred_divergence_relative
                                      << " deferred_face_l2_over_rhs="
                                      << flux_gap.deferred_face_l2_relative
                                      << " worst_deferred_face="
                                      << flux_gap.worst_deferred_face
                                      << " worst_deferred_cell="
                                      << flux_gap.worst_deferred_cell
                                      << " worst_full_divergence_cell="
                                      << flux_gap.worst_divergence_cell
                                      << " worst_full_divergence_x="
                                      << geometry.cell_centres[flux_gap.worst_divergence_cell].x
                                      << " worst_full_divergence_y="
                                      << geometry.cell_centres[flux_gap.worst_divergence_cell].y
                                      << " worst_full_divergence_z="
                                      << geometry.cell_centres[flux_gap.worst_divergence_cell].z
                                      << " full_rc_face_rounding_l2_over_rhs="
                                      << static_cast<double>(std::sqrt(mismatch2/rhs2))
                                      << " rhs_base_contract_pass="
                                      << bridge.rhs_matches_base_flux
                                      << " full_rc_qualified=0"
                                      << " pressure_state_unchanged=1\n";
                        } else {
                            std::cerr << "PRESSURE_HILO_RHIE_CHOW_BRIDGE"
                                      << " iter=" << iter << " corrector=" << corr
                                      << " status=UNSUPPORTED_ALGORITHM_OR_SKEW"
                                      << " full_rc_qualified=0 pressure_state_unchanged=1\n";
                        }
                        std::cerr << std::scientific << std::setprecision(16)
                                  << "PRESSURE_HILO_FACES_LIVE"
                                  << " iter=" << iter << " corrector=" << corr
                                  << " active_faces="
                                  << physical_implicit_faces.active_face_count()
                                  << " before_relative="
                                  << face_report.initial_face_relative
                                  << " after_relative="
                                  << face_report.final_face_relative
                                  << " original_csr_after_relative="
                                  << face_report.final_csr_relative
                                  << " accepted_refinements="
                                  << face_report.accepted_refinements
                                  << " implicit_face_converged="
                                  << face_report.implicit_face_converged
                                  << " original_csr_converged="
                                  << face_report.original_csr_converged
                                  << " full_rhie_chow_qualified=0"
                                  << " original_status_unchanged=1\n";
                        std::cerr << std::scientific << std::setprecision(16)
                                  << "PRESSURE_HILO_LIVE"
                                  << " iter=" << iter << " corrector=" << corr
                                  << " n=" << nc
                                  << " original_double_relative=" << rp.residual_relative
                                  << " initial_split_relative="
                                  << split.initial_relative_residual
                                  << " final_split_relative="
                                  << split.final_relative_residual
                                  << " final_fp64_relative="
                                  << split.final_double_only_residual
                                  << " accepted_refinements="
                                  << split.passes_accepted
                                  << " algebraic_converged="
                                  << split.algebraic_converged
                                  << " physical_fv_qualified=0"
                                  << " original_status_unchanged=1\n";
                    }
                } catch (const std::exception& e) {
                    std::cerr << "PRESSURE_HILO_LIVE status=DIAGNOSTIC_REJECTED"
                              << " detail=" << e.what()
                              << " physical_fv_qualified=0\n";
                }
            }

            if (rp.status != cfdx::core::SolverStatus::CONVERGED) {
                if (controls.momentum_assembly_forensic_callback) {
                    const auto& matrix = *solved_pressure_matrix;
                    const auto& rhs = *solved_pressure_rhs;
                    double rhs_linf = 0.0;
                    double diagonal_min = std::numeric_limits<double>::infinity();
                    double diagonal_max = 0.0;
                    double row_sum_linf = 0.0;
                    double row_abs_sum_linf = 0.0;
                    std::size_t rhs_cell = 0;
                    std::size_t diagonal_cell = 0;
                    std::size_t nonfinite_matrix = 0;
                    std::size_t nonfinite_rhs = 0;
                    const auto* offsets = matrix.row_offsets_data();
                    const auto* columns = matrix.columns_data();
                    const auto* values = matrix.values_data();
                    for (std::size_t row = 0; row < matrix.n_rows(); ++row) {
                        const double rhs_value = rhs(row);
                        if (std::isfinite(rhs_value)) {
                            if (std::abs(rhs_value) > rhs_linf) {
                                rhs_linf = std::abs(rhs_value);
                                rhs_cell = pressure_row_to_cell(row);
                            }
                        } else {
                            ++nonfinite_rhs;
                        }
                        double diagonal = 0.0;
                        double row_sum = 0.0;
                        double row_abs_sum = 0.0;
                        for (std::uint32_t k = offsets[row]; k < offsets[row + 1]; ++k) {
                            const double value = values[k];
                            if (!std::isfinite(value)) {
                                ++nonfinite_matrix;
                                continue;
                            }
                            row_sum += value;
                            row_abs_sum += std::abs(value);
                            if (columns[k] == row) diagonal = value;
                        }
                        if (std::isfinite(diagonal)) {
                            diagonal_min = std::min(diagonal_min, diagonal);
                            if (std::abs(diagonal) > diagonal_max) {
                                diagonal_max = std::abs(diagonal);
                                diagonal_cell = pressure_row_to_cell(row);
                            }
                        } else {
                            ++nonfinite_matrix;
                        }
                        row_sum_linf = std::max(row_sum_linf, std::abs(row_sum));
                        row_abs_sum_linf = std::max(row_abs_sum_linf, row_abs_sum);
                    }
                    const auto& rhs_point = geometry.cell_centres[rhs_cell];
                    double pressure_linf = 0.0;
                    double grad_pressure_linf = 0.0;
                    double predictor_flux_linf = 0.0;
                    for (std::size_t cell = 0; cell < nc; ++cell) {
                        pressure_linf = std::max(pressure_linf, std::abs(p(cell)));
                        for (std::size_t d = 0; d < 3; ++d)
                            grad_pressure_linf = std::max(
                                grad_pressure_linf,
                                std::abs(current_grad_p.component_data(d)[cell]));
                    }
                    for (std::size_t face = 0; face < mesh.n_faces(); ++face)
                        predictor_flux_linf = std::max(
                            predictor_flux_linf, std::abs(phiHbyA(face)));
                    std::cerr << std::setprecision(17)
                              << "NACA_PRESSURE_CORRECTION_FAILURE"
                              << " iter=" << iter
                              << " corrector=" << corr
                              << " status=" << static_cast<int>(rp.status)
                              << " iterations=" << rp.iterations
                              << " residual=" << rp.residual
                              << " relative=" << rp.residual_relative
                              << " rows=" << matrix.n_rows()
                              << " nnz=" << matrix.nnz()
                              << " rhs_linf=" << rhs_linf
                              << " rhs_cell=" << rhs_cell
                              << " rhs_xyz=" << rhs_point.x << ":" << rhs_point.y << ":" << rhs_point.z
                              << " diag_min=" << diagonal_min
                              << " diag_abs_max=" << diagonal_max
                              << " diag_cell=" << diagonal_cell
                              << " row_sum_linf=" << row_sum_linf
                              << " row_abs_sum_linf=" << row_abs_sum_linf
                              << " max_pressure=" << pressure_linf
                              << " max_grad_pressure=" << grad_pressure_linf
                              << " max_predictor_flux=" << predictor_flux_linf
                              << " nonfinite_matrix=" << nonfinite_matrix
                              << " nonfinite_rhs=" << nonfinite_rhs
                              << '\n';
                }
                std::ostringstream message;
                message << "solve_steady_incompressible: pressure-correction solve did not converge "
                        << "(status=" << static_cast<int>(rp.status)
                        << ", iterations=" << rp.iterations
                        << ", residual=" << std::scientific
                        << std::setprecision(std::numeric_limits<double>::max_digits10)
                        << rp.residual
                        << ", relative=" << rp.residual_relative << ")";
                throw NonlinearRetryableFailure(message.str());
            }

            // V17: audit the CURRENT nonlinear flux before applying the
            // pressure correction. Do not substitute the pressure-matrix
            // residual for the post-correction continuity imbalance.
            std::optional<V17MassPatchAudit> v17_pressure_before;
            if (controls.diagnostics.flux_patch_audit) {
                auto before_flux=make_rhie_chow_mass_flux(
                    mesh,geometry,HbyA,p,
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                    controls.density,velocity_bcs,pressure_bcs,
                    RhieChowVelocityContract::HBYA_PREDICTOR,
                    controls.pressure_skew_corrected_gradient,
                    controls.pressure_projected_face_interpolation);
                if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC)
                    add_simplec_consistent_response(before_flux);
                v17_pressure_before=audit_v17_authoritative_flux(mesh,geometry,before_flux);
            }

            // Relax the physical pressure exactly once. The reference is a
            // gauge: enforce it by a uniform shift, never by overwriting one
            // cell and creating an artificial local pressure jump.
            const double nominal_alpha_p = retry_controller.alpha_p(effective_alpha_p);
            double alpha_p = nominal_alpha_p;
            std::optional<PressureAffineFluxTrials> affine_line_search;
            if (controls.pressure_continuity_line_search.enabled) {
                std::vector<double> base_pressure(nc);
                for (std::size_t c = 0; c < nc; ++c) base_pressure[c] = p(c);
                const auto restore_base_pressure = [&]() {
                    for (std::size_t c = 0; c < nc; ++c) p(c) = base_pressure[c];
                };
                const auto trial_flux_at = [&](double candidate_alpha) {
                    for (std::size_t c = 0; c < nc; ++c)
                        p(c) = base_pressure[c] + candidate_alpha * p_corr(c);
                    apply_pressure_gauge_policy(
                        p, geometry, controls.pressure_gauge_policy,
                        controls.pressure_reference_cell, controls.pressure_reference_value,
                        has_fixed_pressure_boundary);
                    auto trial_flux = make_rhie_chow_mass_flux(
                        mesh, geometry, HbyA, p,
                        controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                        controls.density, velocity_bcs, pressure_bcs,
                        RhieChowVelocityContract::HBYA_PREDICTOR,
                        controls.pressure_skew_corrected_gradient,
                        controls.pressure_projected_face_interpolation);
                    if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC)
                        add_simplec_consistent_response(trial_flux);
                    return trial_flux;
                };
                const auto continuity_score = [&](const auto& flux) {
                    std::vector<double> balance(nc, 0.0);
                    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                        const double value = flux(f);
                        if (!std::isfinite(value))
                            return std::numeric_limits<double>::infinity();
                        balance[mesh.ownership().owner(f)] += value;
                        const auto neighbour = mesh.ownership().neighbour(f);
                        if (neighbour >= 0)
                            balance[static_cast<std::size_t>(neighbour)] -= value;
                    }
                    double norm = 0.0;
                    for (double v : balance) {
                        if (!std::isfinite(v)) return std::numeric_limits<double>::infinity();
                        norm += std::abs(v);
                    }
                    return norm;
                };
                // All trial pressure evaluations are transactional. On ANY
                // rejected selection or unexpected exception, p is restored.
                try {
                    if (controls.diagnostics.flux_patch_audit) {
                        constexpr double epsilon=1e-5;
                        const auto f0=audit_v17_authoritative_flux(
                            mesh,geometry,trial_flux_at(0.0));
                        const auto f1=audit_v17_authoritative_flux(
                            mesh,geometry,trial_flux_at(epsilon));
                        // Independent operator contract: (a) RHS must equal
                        // minus the divergence of the *same* pressure/flux
                        // state, (b) A*p_corr=b, (c) the measured face-flux
                        // Jacobian must equal A*p_corr. This distinguishes a
                        // wrong source from a wrong face operator.
                        double rhs_flux_mismatch_l1=0.0;
                        double matrix_solve_mismatch_l1=0.0;
                        double flux_jacobian_mismatch_l1=0.0;
                        double sum_rhs_l1=0.0, sum_derivative_l1=0.0;
                        std::size_t worst_rhs_cell=0, worst_jac_cell=0;
                        double rhs_mismatch_linf=0.0, jac_mismatch_linf=0.0;
                        for(std::size_t c=0;c<nc;++c) {
                            double matrix_times_correction=0.0;
                            for(const auto& [col,val]:rows[c])
                                matrix_times_correction+=val*p_corr(col);
                            const double rhs_err=b(c)+f0.cell_divergence[c];
                            const double solve_err=matrix_times_correction-b(c);
                            const double derivative=(f1.cell_divergence[c]-f0.cell_divergence[c])/epsilon;
                            const double jac_err=derivative-matrix_times_correction;
                            rhs_flux_mismatch_l1+=std::abs(rhs_err);
                            matrix_solve_mismatch_l1+=std::abs(solve_err);
                            flux_jacobian_mismatch_l1+=std::abs(jac_err);
                            sum_rhs_l1+=std::abs(b(c));
                            sum_derivative_l1+=std::abs(derivative);
                            if(std::abs(rhs_err)>rhs_mismatch_linf) {
                                rhs_mismatch_linf=std::abs(rhs_err);
                                worst_rhs_cell=c;
                            }
                            if(std::abs(jac_err)>jac_mismatch_linf) {
                                jac_mismatch_linf=std::abs(jac_err);
                                worst_jac_cell=c;
                            }
                        }
                        std::cerr << std::scientific << std::setprecision(13)
                                  << "V17_PRESSURE_OPERATOR iteration=" << iter
                                  << " corrector=" << corr
                                  << " rhs_flux_mismatch_l1=" << rhs_flux_mismatch_l1
                                  << " matrix_solve_mismatch_l1=" << matrix_solve_mismatch_l1
                                  << " flux_jacobian_mismatch_l1=" << flux_jacobian_mismatch_l1
                                  << " rhs_l1=" << sum_rhs_l1
                                  << " derivative_l1=" << sum_derivative_l1
                                  << " rhs_mismatch_linf=" << rhs_mismatch_linf
                                  << " worst_rhs_cell=" << worst_rhs_cell
                                  << " jac_mismatch_linf=" << jac_mismatch_linf
                                  << " worst_jac_cell=" << worst_jac_cell
                                  << " qualified=0\n";
                        double mismatch_l1=0.0, expected_l1=0.0;
                        double mismatch_linf=0.0;
                        std::size_t worst=0;
                        double worst_before=0.0, worst_jvp=0.0;
                        for (std::size_t c=0;c<nc;++c) {
                            const double before=f0.cell_divergence[c];
                            const double derivative=(f1.cell_divergence[c]-before)/epsilon;
                            const double mismatch=derivative+before;
                            mismatch_l1+=std::abs(mismatch);
                            expected_l1+=std::abs(before);
                            if(std::abs(mismatch)>mismatch_linf) {
                                mismatch_linf=std::abs(mismatch);
                                worst=c;worst_before=before;worst_jvp=derivative;
                            }
                        }
                        std::cerr << std::scientific << std::setprecision(13)
                                  << "V17_PRESSURE_JVP iteration=" << iter
                                  << " corrector=" << corr
                                  << " epsilon=" << epsilon
                                  << " mismatch_l1=" << mismatch_l1
                                  << " baseline_l1=" << expected_l1
                                  << " ratio=" << mismatch_l1/std::max(expected_l1,1e-300)
                                  << " mismatch_linf=" << mismatch_linf
                                  << " worst_cell=" << worst
                                  << " before=" << worst_before
                                  << " actual_derivative=" << worst_jvp
                                  << " qualified=0\n";
                        restore_base_pressure();
                    }
                    if (controls.pressure_continuity_line_search.affine_flux_trials) {
                        const auto f0 = trial_flux_at(0.0);
                        const auto f1 = trial_flux_at(nominal_alpha_p);
                        std::vector<double> v0, v1;
                        std::vector<std::size_t> owner;
                        std::vector<std::int64_t> neighbour;
                        v0.reserve(mesh.n_faces()); v1.reserve(mesh.n_faces());
                        owner.reserve(mesh.n_faces()); neighbour.reserve(mesh.n_faces());
                        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                            v0.push_back(f0(f)); v1.push_back(f1(f));
                            owner.push_back(mesh.ownership().owner(f));
                            neighbour.push_back(mesh.ownership().neighbour(f));
                        }
                        affine_line_search.emplace(v0, v1, owner, neighbour, nc,
                                                   nominal_alpha_p);
                    }
                    const auto score = [&](double candidate_alpha) {
                        const double value=affine_line_search
                            ? affine_line_search->continuity_norm(candidate_alpha)
                            : continuity_score(trial_flux_at(candidate_alpha));
                        if (controls.diagnostics.flux_patch_audit)
                            std::cerr << std::scientific << std::setprecision(13)
                                      << "V17_PRESSURE_TRIAL iteration=" << iter
                                      << " corrector=" << corr
                                      << " alpha=" << candidate_alpha
                                      << " continuity_l1=" << value
                                      << " qualified=0\n";
                        return value;
                    };
                    const auto selection = select_pressure_continuity_step(
                        nominal_alpha_p, controls.pressure_continuity_line_search, score);
                    alpha_p = selection.alpha;
                    last_line_search_before = selection.baseline_norm;
                    last_line_search_after = selection.accepted_norm;
                } catch (...) {
                    restore_base_pressure();
                    throw;
                }
                restore_base_pressure();
            }
            last_selected_pressure_alpha = alpha_p;
            if (controls.pressure_high_low_production) {
                // Transactional high/low update of BOTH SIMPLE pressures:
                // alpha*p_corr for velocity and full p_corr for transport.
                auto physical_high=p;
                auto physical_low=pressure_low_state;
                auto flux_high=p;
                auto flux_low=pressure_low_state;
                add_pressure_high_low(physical_high,physical_low,
                                      p_corr,p_corr_low,alpha_p);
                add_pressure_high_low(flux_high,flux_low,
                                      p_corr,p_corr_low,1.0);
                p=std::move(physical_high);
                pressure_low_state=std::move(physical_low);
                v19_flux_pressure_state=std::move(flux_high);
                v19_flux_pressure_low=std::move(flux_low);
            } else {
                apply_v19_pressure_flux_policy(
                    p, v19_flux_pressure_state, p_corr, alpha_p,
                    controls.v19_full_pressure_flux_step,
                    geometry, controls.pressure_gauge_policy,
                    controls.pressure_reference_cell, controls.pressure_reference_value,
                    has_fixed_pressure_boundary);
            }

            auto corrected_grad_p = controls.pressure_high_low_production
                ? gauss_gradient_pressure_high_low(
                    mesh,geometry,p,pressure_low_state,pressure_bcs,
                    controls.pressure_projected_face_interpolation)
                : gauss_gradient_with_boundary(
                    p,mesh,geometry,pressure_bcs,
                    controls.pressure_skew_corrected_gradient,
                    controls.pressure_projected_face_interpolation);

            // Rebuild the corrected velocity from the same rAU operator
            // that built the cell HbyA, so the converged velocity
            // satisfies the momentum equation. The SIMPLEC consistent
            // operator (rAtU) belongs to the pressure equation and its
            // face fluxes only. This is the momentum correction; no
            // least-squares flux fitting is permitted.
            for (std::size_t c = 0; c < nc; ++c) {
                U.component_data(0)[c] =
                    HbyA.component_data(0)[c] - rAU[0][c] * corrected_grad_p.component_data(0)[c] * geometry.cell_volumes[c];
                U.component_data(1)[c] =
                    HbyA.component_data(1)[c] - rAU[1][c] * corrected_grad_p.component_data(1)[c] * geometry.cell_volumes[c];
                U.component_data(2)[c] =
                    HbyA.component_data(2)[c] - rAU[2][c] * corrected_grad_p.component_data(2)[c] * geometry.cell_volumes[c];
            }

            // The conservative face flux is rebuilt from the same
            // momentum-weighted interpolation used by the pressure equation.
            // This is the single authoritative phi state for continuity.
            if (controls.pressure_high_low_production) {
                const auto exact_flux=make_rhie_chow_mass_flux_high_low(
                    mesh,geometry,HbyA,v19_flux_pressure_state,
                    v19_flux_pressure_low,rAU,controls.density,
                    velocity_bcs,pressure_bcs,
                    controls.pressure_projected_face_interpolation,
                    face_velocity_conditions ? &*face_velocity_conditions : nullptr);
                const auto exact_balance=pressure_true_face_residual(
                    mesh,exact_flux,b);
                for(std::size_t f=0;f<mesh.n_faces();++f)
                    mass_flux(f)=static_cast<double>(exact_flux[f]);
                std::vector<long double> rounded_flux(mesh.n_faces());
                for(std::size_t f=0;f<mesh.n_faces();++f)
                    rounded_flux[f]=mass_flux(f);
                const double double_relative=pressure_true_face_residual(
                    mesh,rounded_flux,b).second;
                std::cerr << std::scientific << std::setprecision(16)
                          << "PRESSURE_HILO_PHYSICAL_FLUX iteration=" << iter
                          << " corrector=" << corr
                          << " exact_relative=" << exact_balance.second
                          << " stored_fp64_relative=" << double_relative
                          << " alpha_p=" << alpha_p
                          << " physical_pressure_low_max=" << [&] {
                            double maximum=0.0;
                            for(std::size_t c=0;c<nc;++c)
                                maximum=std::max(maximum,std::abs(pressure_low_state(c)));
                            return maximum;
                          }()
                          << " accepted_iteration=0\n";
                const double physical_pressure_tolerance=
                    controls.pressure_linear_tolerance>0.0
                    ? controls.pressure_linear_tolerance : linear_tolerance;
                if (exact_balance.second>physical_pressure_tolerance ||
                    double_relative>physical_pressure_tolerance)
                    throw NonlinearRetryableFailure(
                        "V6 compensated pressure: full Rhie-Chow transport flux "
                        "does not satisfy the original pressure tolerance");
            } else {
                mass_flux = make_rhie_chow_mass_flux(
                    mesh, geometry, HbyA, v19_flux_pressure_state,
                    controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU,
                    controls.density, velocity_bcs, pressure_bcs,
                    RhieChowVelocityContract::HBYA_PREDICTOR,
                    controls.pressure_skew_corrected_gradient,
                    controls.pressure_projected_face_interpolation,
                    face_velocity_conditions ? &*face_velocity_conditions : nullptr);
            }

            if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC)
                add_simplec_consistent_response(mass_flux);

            if (controls.v19_full_pressure_flux_step && controls.diagnostics.flux_patch_audit) {
                // Both fluxes use identical predictor, rAU, face BCs and
                // geometry. The only difference is full vs relaxed pressure.
                const auto relaxed_flux=make_rhie_chow_mass_flux(
                    mesh,geometry,HbyA,p,rAU,controls.density,
                    velocity_bcs,pressure_bcs,
                    RhieChowVelocityContract::HBYA_PREDICTOR,
                    controls.pressure_skew_corrected_gradient,
                    controls.pressure_projected_face_interpolation,
                    face_velocity_conditions ? &*face_velocity_conditions : nullptr);
                double l1_delta=0.0, linf_delta=0.0;
                for(std::size_t f=0;f<mesh.n_faces();++f) {
                    const double delta=std::abs(mass_flux(f)-relaxed_flux(f));
                    l1_delta+=delta;
                    linf_delta=std::max(linf_delta,delta);
                }
                v20_relaxed_pressure_flux=relaxed_flux;
                const auto a=audit_v17_authoritative_flux(mesh,geometry,mass_flux);
                const auto b=audit_v17_authoritative_flux(mesh,geometry,relaxed_flux);
                std::cerr << std::scientific << std::setprecision(13)
                          << "V19_SPLIT_PRESSURE_FLUX iteration=" << iter
                          << " corrector=" << corr
                          << " alpha_p=" << alpha_p
                          << " full_flux_continuity_l1=" << a.cell_l1
                          << " relaxed_flux_continuity_l1=" << b.cell_l1
                          << " face_delta_l1=" << l1_delta
                          << " face_delta_linf=" << linf_delta
                          << " qualified=0\n";
            }
            if (v17_pressure_before) {
                const auto after=audit_v17_authoritative_flux(mesh,geometry,mass_flux);
                const double ratio=after.cell_l1/std::max(v17_pressure_before->cell_l1,1e-300);
                std::cerr << std::scientific << std::setprecision(13)
                          << "V17_PRESSURE_CORRECTION iteration=" << iter
                          << " corrector=" << corr
                          << " alpha_p=" << alpha_p
                          << " before_l1=" << v17_pressure_before->cell_l1
                          << " after_l1=" << after.cell_l1
                          << " before_linf=" << v17_pressure_before->cell_linf
                          << " after_linf=" << after.cell_linf
                          << " reduction_ratio=" << ratio
                          << " pressure_linear_relative=" << rp.residual_relative
                          << " qualified=0\n";
            }

            // Final authoritative reconstruction. This independent check
            // refuses a nonlinear boundary/response path rather than quietly
            // accepting an invalid affine trial on an unsupported contract.
            if (affine_line_search) {
                std::vector<double> exact_flux;
                exact_flux.reserve(mesh.n_faces());
                for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                    exact_flux.push_back(mass_flux(f));
                if (!affine_line_search->agrees_with_flux(alpha_p, exact_flux))
                    throw NonlinearRetryableFailure(
                        "pressure affine line-search: authoritative flux mismatch");
            }

            // Issa corrector form: every pressure corrector is driven by the
            // ORIGINAL predictor flux together with the current pressure,
            // div(phiHbyA - RC(p)); the assembly at the top of this loop
            // already subtracts the Rhie-Chow pressure flux of the current
            // p. Overwriting phiHbyA with the corrected flux here would
            // subtract RC(p) twice in the next corrector's continuity source
            // (div(phiHbyA) - 2*div(RC(p))), which displaces the fixed point
            // and leaves a permanent continuity floor on every
            // multi-corrector algorithm (PISO, PIMPLE, FRACTIONAL_STEP).
            // The momentum matrix stays frozen inside the inner loop; the
            // relaxed pressure update is what carries the first correction
            // into the next one. In the optional iterated-PISO variant,
            // refresh the diagonal-inverse predictor H/A using the just-
            // corrected cell velocities. This restores the neighbour terms
            // of the frozen momentum operator before the next pressure solve.
            // The next pressure RHS and its face flux must BOTH use this new
            // predictor. Do not feed the already pressure-corrected flux back
            // as phiHbyA: the pressure gradient would then be applied twice.
            if (controls.coupling.piso_neighbor_correction &&
                corr + 1 < pressure_correctors) {
                for (std::size_t c = 0; c < nc; ++c) {
                    ux(c) = U.component_data(0)[c];
                    uy(c) = U.component_data(1)[c];
                    uz(c) = U.component_data(2)[c];
                }
                hbya[0] = build_hbya(ex, ux, grad_p, 0, rAU[0]);
                hbya[1] = build_hbya(ey, uy, grad_p, 1, rAU[1]);
                hbya[2] = build_hbya(ez, uz, grad_p, 2, rAU[2]);
                HbyA = make_hbya_field(hbya);
                const auto old_predictor_flux = phiHbyA;
                phiHbyA = make_mass_flux(
                    mesh, geometry, HbyA, controls.density, velocity_bcs,
                    {}, controls.pressure_projected_face_interpolation,
                    face_velocity_conditions ? &*face_velocity_conditions : nullptr);
                double predictor_delta_linf = 0.0;
                for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                    predictor_delta_linf = std::max(
                        predictor_delta_linf,
                        std::abs(phiHbyA(f) - old_predictor_flux(f)));
                ++piso_neighbor_predictor_rebuilds;
                piso_neighbor_predictor_flux_delta_linf = std::max(
                    piso_neighbor_predictor_flux_delta_linf, predictor_delta_linf);
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "PISO_NEIGHBOUR_UPDATE iteration=" << iter
                              << " corrector=" << corr
                              << " predictor_flux_delta_linf=" << predictor_delta_linf
                              << " (experimental; qualification=0)\n";
            }
        }

        // Operator identity audit: re-derive the authoritative flux from the
        // CURRENT HbyA/flux-pressure/rAtU(rAU) state through the same operator primitives
        // the pressure equation balanced (including the SIMPLEC
        // consistent-response term) and compare with the stored mass_flux.
        // Expected zero up to round-off; a nonzero value flags a stale or
        // mutated flux, or an omitted/duplicated SIMPLEC response term.
        {
            auto rederived_flux = make_rhie_chow_mass_flux(
                mesh, geometry, HbyA, v19_flux_pressure_state,
                controls.algorithm == PressureVelocityAlgorithm::SIMPLEC
                    ? rAtU : rAU,
                controls.density, velocity_bcs, pressure_bcs,
                RhieChowVelocityContract::HBYA_PREDICTOR,
                controls.pressure_skew_corrected_gradient,
                controls.pressure_projected_face_interpolation,
                face_velocity_conditions ? &*face_velocity_conditions : nullptr);
            if (controls.pressure_high_low_production) {
                const auto split_flux=make_rhie_chow_mass_flux_high_low(
                    mesh,geometry,HbyA,v19_flux_pressure_state,
                    v19_flux_pressure_low,rAU,controls.density,
                    velocity_bcs,pressure_bcs,
                    controls.pressure_projected_face_interpolation,
                    face_velocity_conditions ? &*face_velocity_conditions : nullptr);
                for(std::size_t f=0;f<mesh.n_faces();++f)
                    rederived_flux(f)=static_cast<double>(split_flux[f]);
            }
            if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC)
                add_simplec_consistent_response(rederived_flux);
            for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                rhie_chow_identity_audit_linf = std::max(
                    rhie_chow_identity_audit_linf,
                    std::abs(mass_flux(f) - rederived_flux(f)));
            rhie_chow_identity_evaluated = true;
        }
        if (controls.momentum_same_state_face_callback) {
            // Capture the EXACT live predictor and pressure states before
            // the nonlinear step can advance or replace the flux. Keep it
            // pending until the transaction is committed.
            const auto predictor_base = make_mass_flux(
                mesh, geometry, HbyA, controls.density, velocity_bcs,
                {}, controls.pressure_projected_face_interpolation);
            const auto& response_inverse =
                controls.algorithm == PressureVelocityAlgorithm::SIMPLEC ? rAtU : rAU;
            const auto pressure_relaxed_flux = make_rhie_chow_mass_flux(
                mesh, geometry, HbyA, p, response_inverse, controls.density,
                velocity_bcs, pressure_bcs,
                RhieChowVelocityContract::HBYA_PREDICTOR,
                controls.pressure_skew_corrected_gradient,
                controls.pressure_projected_face_interpolation);
            const auto pressure_full_flux = make_rhie_chow_mass_flux(
                mesh, geometry, HbyA, v19_flux_pressure_state, response_inverse,
                controls.density, velocity_bcs, pressure_bcs,
                RhieChowVelocityContract::HBYA_PREDICTOR,
                controls.pressure_skew_corrected_gradient,
                controls.pressure_projected_face_interpolation);
            // Independent internal-face pressure Schur action at the same
            // flux-pressure state and rAU/rAtU. Boundaries are not compared
            // to the internal-face formula (they have different BC terms).
            const auto full_grad = gauss_gradient_with_boundary(
                v19_flux_pressure_state, mesh, geometry, pressure_bcs,
                controls.pressure_skew_corrected_gradient,
                controls.pressure_projected_face_interpolation);
            cfdx::core::Field<double,cfdx::core::Location::FACE> independent_pressure(
                mesh.n_faces(), "rc_pressure_face", "kg/s", 1);
            independent_pressure.fill(0.0);
            for (std::size_t f=0; f<mesh.n_faces(); ++f)
                if (mesh.ownership().neighbour(f)>=0)
                    independent_pressure(f) = rhie_chow_pressure_flux_internal(
                        mesh, geometry, f, v19_flux_pressure_state, full_grad,
                        response_inverse, controls.density,
                        controls.pressure_projected_face_interpolation);
            cfdx::core::Field<double,cfdx::core::Location::CELL> uy_trace(
                mesh.n_cells(), "Uy_same_state", "m/s", 1);
            for (std::size_t cell=0; cell<mesh.n_cells(); ++cell)
                uy_trace(cell)=U.component_data(1)[cell];
            pending_same_state_trace = audit_momentum_same_state_faces(
                mesh, geometry, mass_flux, uy_trace, ubc_y,
                controls.convection_scheme, outer_high_order_fraction,
                predictor_base, pressure_relaxed_flux, pressure_full_flux,
                &independent_pressure);
            const double op_scale = [&] {
                double scale = 1.0;
                for (std::size_t f=0; f<mesh.n_faces(); ++f)
                    scale = std::max(scale, std::abs(mass_flux(f)));
                return scale;
            }();
            if (pending_same_state_trace->max_pressure_operator_error >
                5e-10*op_scale)
                throw std::runtime_error("same-state RC pressure response/operator mismatch");
        }

        if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC) {
            // DIFFERENT operator path from the assembled rAtU correction:
            // use the plain momentum inverse rAU on the very same accepted
            // HbyA/pressure/face BCs. No stale pressure, nor relaxed-U flux.
            const auto plain_flux = make_rhie_chow_mass_flux(
                mesh, geometry, HbyA, v19_flux_pressure_state,
                rAU, controls.density, velocity_bcs, pressure_bcs,
                RhieChowVelocityContract::HBYA_PREDICTOR,
                controls.pressure_skew_corrected_gradient,
                controls.pressure_projected_face_interpolation,
                face_velocity_conditions ? &*face_velocity_conditions : nullptr);
            double gap = 0.0, phi_scale = 1.0;
            std::size_t worst_face = 0;
            for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                const double error = std::abs(mass_flux(f)-plain_flux(f));
                if (error > gap) { gap = error; worst_face = f; }
                phi_scale = std::max(phi_scale, std::abs(mass_flux(f)));
            }
            simplec_plain_rau_gap = gap;
            simplec_plain_rau_worst = worst_face;
            if (controls.require_simplec_same_state_flux_identity &&
                gap > 1e-10 * phi_scale)
                throw NonlinearRetryableFailure(
                    "SIMPLEC same-state rAtU/rAU flux identity mismatch");
        }

        }
        
        // Independent final momentum residual: evaluate the actual corrected
        // state against freshly assembled equations. Do not use a U reconstructed
        // from the conservative flux as the residual state.
        Field<double, Location::CELL> final_ux_field(mesh.n_cells(), "Ux_final", "m/s", 1);
        Field<double, Location::CELL> final_uy_field(mesh.n_cells(), "Uy_final", "m/s", 1);

        Field<double, Location::CELL> final_uz_field(mesh.n_cells(), "Uz_final", "m/s", 1);
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            final_ux_field(c) = U.component_data(0)[c];
            final_uy_field(c) = U.component_data(1)[c];
            final_uz_field(c) = U.component_data(2)[c];
        }
        auto final_grad_p = controls.pressure_high_low_production
            ? gauss_gradient_pressure_high_low(
                mesh,geometry,p,pressure_low_state,pressure_bcs,
                controls.pressure_projected_face_interpolation)
            : gauss_gradient_with_boundary(
                p,mesh,geometry,pressure_bcs,
                controls.pressure_skew_corrected_gradient,
                controls.pressure_projected_face_interpolation);
        if (include_transpose_stress) {
            const auto final_stress_flux = transpose_viscous_stress_flux(
                mesh, geometry, U, mu_eff, {ubc_x, ubc_y, ubc_z},
                controls.sst_runtime
                    ? controls.density * controls.kinematic_viscosity : -1.0,
                    active_wall_face_mu_v17);
            const auto final_transpose_force = viscous_face_flux_divergence(
                mesh, geometry, final_stress_flux);
            for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
                body_x(cell) += final_transpose_force(cell,0)-transpose_force(cell,0);
                body_y(cell) += final_transpose_force(cell,1)-transpose_force(cell,1);
                body_z(cell) += final_transpose_force(cell,2)-transpose_force(cell,2);
            }
        }
        const auto final_face_velocity_conditions = mixed_velocity_faces(mass_flux);
        bool mixed_modes_stable = true;
        if (face_velocity_conditions && final_face_velocity_conditions) {
            for (const auto& [name, value] : controls.inlet_outlet_velocity_patches) {
                const auto& patch = mesh.boundary().patch(mesh.boundary().find(name));
                for (const auto face : patch.face_ids)
                    if ((*face_velocity_conditions)[0].conditions[face].type !=
                        (*final_face_velocity_conditions)[0].conditions[face].type)
                        mixed_modes_stable = false;
            }
        }
        auto final_ex = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme,
            &final_ux_field, nullptr, sst_molecular_wall_mu,
            outer_high_order_fraction,
            final_face_velocity_conditions ? &(*final_face_velocity_conditions)[0] : nullptr);
        auto final_ey = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme,
            &final_uy_field, nullptr, sst_molecular_wall_mu,
            outer_high_order_fraction,
            final_face_velocity_conditions ? &(*final_face_velocity_conditions)[1] : nullptr);
        auto final_ez = assemble_momentum_component(
            mesh, geometry, mass_flux, final_grad_p, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme,
            &final_uz_field, nullptr, sst_molecular_wall_mu,
            outer_high_order_fraction,
            final_face_velocity_conditions ? &(*final_face_velocity_conditions)[2] : nullptr);

        // The final nonlinear acceptance residual must represent the same
        // production equation that was solved. For transient runs the
        // physical-time derivative is part of that equation; omitting it here
        // would measure only the steady spatial operator and reject every
        // legitimate transient state whose temporal acceleration balances the
        // spatial momentum residual. Reuse the exact production mass-term
        // assembly, including variable-step BDF2 coefficients and the accepted
        // history, so the acceptance metric is algebraically consistent with
        // the coupled solve above.
        if (controls.transient.enabled) {
            const auto& u_n = *controls.transient.previous;
            const auto& u_nm1 = controls.transient.history_valid
                ? *controls.transient.previous_previous : u_n;
            const bool use_bdf2 =
                controls.transient.scheme == cfdx::core::TimeScheme::BDF2 &&
                controls.transient.history_valid;
            add_dual_time_mass_term(
                final_ex, u_n, u_nm1, geometry, controls.density,
                controls.transient.dt, controls.transient.dt_previous,
                use_bdf2, 0);
            add_dual_time_mass_term(
                final_ey, u_n, u_nm1, geometry, controls.density,
                controls.transient.dt, controls.transient.dt_previous,
                use_bdf2, 1);
            add_dual_time_mass_term(
                final_ez, u_n, u_nm1, geometry, controls.density,
                controls.transient.dt, controls.transient.dt_previous,
                use_bdf2, 2);
        }

        Vector final_ux(mesh.n_cells()), final_uy(mesh.n_cells()), final_uz(mesh.n_cells());
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            final_ux(c) = U.component_data(0)[c];
            final_uy(c) = U.component_data(1)[c];
            final_uz(c) = U.component_data(2)[c];
        }
        const double final_momentum_residual = std::max({
            scalar_equation_residual_inf(final_ex, final_ux),
            scalar_equation_residual_inf(final_ey, final_uy),
            scalar_equation_residual_inf(final_ez, final_uz)
        });

        // Recompute the cell-wise residuals explicitly so a large independent
        // residual can be localized to a component and to an interior/boundary
        // cell. This is intentionally redundant with scalar_equation_residual_inf:
        // the purpose is forensic diagnostics, not a second numerical path.
        const auto cell_has_boundary_face = [&](std::size_t cell) {
            const Offset off = mesh.cells().offsets_data()[cell];
            const Offset count = mesh.cells().offsets_data()[cell + 1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t face = mesh.cells().faces_data()[off + k];
                if (mesh.ownership().neighbour(face) < 0) {
                    const std::size_t patch = geometry.face_patch[face];
                    if (patch < mesh.boundary().n_patches() &&
                        mesh.boundary().patch(patch).type == PatchType::EMPTY)
                        continue;
                    return true;
                }
            }
            return false;
        };
        struct ResidualDiagnostic {
            double global = 0.0;
            double interior = 0.0;
            double boundary = 0.0;
            std::size_t max_cell = 0;
            double no_pressure = 0.0;
            double pressure_contribution = 0.0;
            std::string patch;
        };
        const auto residual_diagnostics =
            [&](const ScalarEquation& equation,
                const Vector& solution,
                std::size_t component) {
                ResidualDiagnostic diagnostic;
                for (std::size_t row = 0; row < mesh.n_cells(); ++row) {
                    double ax = 0.0;
                    const auto begin = equation.matrix.row_offsets_data()[row];
                    const auto end = equation.matrix.row_offsets_data()[row + 1];
                    for (std::uint32_t k = begin; k < end; ++k)
                        ax += equation.matrix.values_data()[k] *
                              solution(equation.matrix.columns_data()[k]);
                    const double full_signed = ax - equation.rhs(row);
                    const double residual = std::abs(full_signed);
                    if (cell_trace_enabled)
                        cell_trace_momentum[row][component] = residual;
                    if (residual > diagnostic.global) {
                        diagnostic.global = residual;
                        diagnostic.max_cell = row;
                    }
                    if (cell_has_boundary_face(row))
                        diagnostic.boundary = std::max(diagnostic.boundary, residual);
                    else
                        diagnostic.interior = std::max(diagnostic.interior, residual);

                    // final_ex is assembled with source = body - grad(p).
                    // Remove the pressure source only for diagnosis, never for
                    // the acceptance metric, to distinguish pressure/transport
                    // inconsistency from boundary transport inconsistency.
                    const double pressure = final_grad_p.component_data(component)[row] *
                                            geometry.cell_volumes[row];
                    const double transport_signed = full_signed - pressure;
                    diagnostic.no_pressure =
                        std::max(diagnostic.no_pressure, std::abs(transport_signed));
                    diagnostic.pressure_contribution =
                        std::max(diagnostic.pressure_contribution, std::abs(pressure));
                }
                const std::size_t row = diagnostic.max_cell;
                const Offset off = mesh.cells().offsets_data()[row];
                const Offset count = mesh.cells().offsets_data()[row + 1] - off;
                for (Offset k = 0; k < count; ++k) {
                    const std::size_t face = mesh.cells().faces_data()[off + k];
                    if (mesh.ownership().neighbour(face) < 0) {
                        const std::size_t patch = geometry.face_patch[face];
                        if (patch >= mesh.boundary().n_patches())
                            continue;
                        // EMPTY faces are geometric extrusion faces, not
                        // physical boundary patches. Do not let the first
                        // front/back face mask the physical patch adjacent to
                        // the worst momentum-residual cell.
                        if (mesh.boundary().patch(patch).type == PatchType::EMPTY)
                            continue;
                        diagnostic.patch = mesh.boundary().patch(patch).name;
                        break;
                    }
                }
                return diagnostic;
            };
        const auto rx_diag = residual_diagnostics(final_ex, final_ux, 0);
        const auto ry_diag = residual_diagnostics(final_ey, final_uy, 1);
        const auto rz_diag = residual_diagnostics(final_ez, final_uz, 2);

        double pressure_gradient_linf = 0.0;
        double pressure_gradient_l2_sum = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            const double gx = final_grad_p.component_data(0)[c];
            const double gy = final_grad_p.component_data(1)[c];
            const double gz = final_grad_p.component_data(2)[c];
            const double gp2 = gx*gx + gy*gy + gz*gz;
            pressure_gradient_linf = std::max(pressure_gradient_linf, std::sqrt(gp2));
            pressure_gradient_l2_sum += gp2;
        }
        const double pressure_gradient_l2 =
            std::sqrt(pressure_gradient_l2_sum /
                      static_cast<double>(std::max<std::size_t>(mesh.n_cells(), 1)));

        double l1 = 0.0, linf = 0.0;
        const auto* cell_faces = mesh.cells().faces_data();
        const auto* cell_offsets = mesh.cells().offsets_data();
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double div = 0.0;
            const Offset off = cell_offsets[c];
            const Offset count = cell_offsets[c+1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t f = cell_faces[off+k];
                div += mesh.ownership().owner(f) == c ? mass_flux(f) : -mass_flux(f);
            }
            if (cell_trace_enabled)
                cell_trace_continuity[c] = std::abs(div);
            l1 += std::abs(div);
            linf = std::max(linf, std::abs(div));
        }

        constexpr double tiny = 1e-300;
        double velocity_change_inf = 0.0, pressure_change_inf = 0.0;
        double velocity_scale = 0.0;
        double pressure_min = std::numeric_limits<double>::infinity();
        double pressure_max = -std::numeric_limits<double>::infinity();
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            for (std::size_t d = 0; d < 3; ++d) {
                velocity_change_inf = std::max(
                    velocity_change_inf,
                    std::abs(U.component_data(d)[c]-U_old.component_data(d)[c]));
                velocity_scale = std::max(
                    velocity_scale, std::abs(U.component_data(d)[c]));
            }
            const double pressure_delta = has_fixed_pressure_boundary
                ? std::abs(p(c) - p_old(c))
                : std::abs((p(c) - p(controls.pressure_reference_cell)) -
                            (p_old(c) - p_old(controls.pressure_reference_cell)));
            pressure_change_inf = std::max(pressure_change_inf, pressure_delta);
            pressure_min = std::min(pressure_min, p(c));
            pressure_max = std::max(pressure_max, p(c));
        }
        velocity_change_inf /= std::max(velocity_scale, tiny);
        const double pressure_range = pressure_max - pressure_min;
        // Use a physical pressure scale when the converged solution is nearly uniform.
        // Otherwise round-off-level pressure changes are divided by a vanishing range
        // and can prevent convergence forever on constant-pressure cases.
        const double pressure_scale = std::max(
            {pressure_range, controls.density * velocity_scale * velocity_scale, 1.0e-30});
        pressure_change_inf /= pressure_scale;

        const std::array<const ScalarEquation*, 3> momentum_equations{
            &final_ex, &final_ey, &final_ez};
        std::array<double, 3> momentum_rhs_scale_components{0.0, 0.0, 0.0};
        for (std::size_t d = 0; d < momentum_equations.size(); ++d) {
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                momentum_rhs_scale_components[d] = std::max(
                    momentum_rhs_scale_components[d],
                    std::abs(momentum_equations[d]->rhs(c)));
            }
            momentum_rhs_scale_components[d] =
                std::max(momentum_rhs_scale_components[d], tiny);
        }
        const double momentum_rhs_scale = std::max({
            momentum_rhs_scale_components[0],
            momentum_rhs_scale_components[1],
            momentum_rhs_scale_components[2]});

        if (cell_trace_enabled)
            pending_cell_rhs_scales = momentum_rhs_scale_components;
        IncompressibleIteration h;
        if (controls.pseudo_momentum.enabled || controls.outer_ptc_v31_experimental) {
            h.v17_pseudo_cfl=v17_pseudo_cfl;
            h.v17_pseudo_mass_min=v17_mass_min;
            h.v17_pseudo_mass_max=v17_mass_max;
            h.v17_pseudo_capped_cells=v17_capped_cells;
        }
        h.selected_pressure_alpha = last_selected_pressure_alpha;
        h.pressure_line_search_norm_before = last_line_search_before;
        h.pressure_line_search_norm_after = last_line_search_after;
        h.iteration = iter;
        h.high_order_fraction = outer_high_order_fraction;
        h.full_order_eligible = full_order_eligible;
        h.momentum_residual = final_momentum_residual / momentum_rhs_scale;
        h.pressure_residual = pressure_residual;
        h.continuity_l1 = l1;
        h.continuity_linf = linf;
        h.momentum_equation_residual = final_momentum_residual;
        h.momentum_rhs_scale = momentum_rhs_scale;
        h.momentum_rhs_scale_components = momentum_rhs_scale_components;
        h.momentum_equation_residual_rhs_relative =
            final_momentum_residual / momentum_rhs_scale;
        h.momentum_equation_residual_rhs_relative_components = {
            rx_diag.global / momentum_rhs_scale_components[0],
            ry_diag.global / momentum_rhs_scale_components[1],
            rz_diag.global / momentum_rhs_scale_components[2]};
        const auto momentum_component_certificate =
            certify_momentum_residual_components(
                {rx_diag.global, ry_diag.global, rz_diag.global},
                momentum_rhs_scale_components,
                controls.convergence.relative_tolerance);
        h.momentum_equation_residual_rhs_active_component_max =
            momentum_component_certificate.max_active_component_relative;
        if (controls.tvd_deferred_memory.enabled) {
            for (std::size_t d=0; d<3; ++d)
                h.tvd_memory_source_relative = std::max(
                    h.tvd_memory_source_relative,
                    pending_memory_source_linf[d] /
                        std::max(momentum_rhs_scale_components[d], 1e-30));
            if (controls.diagnostics.iteration_trace)
                std::cerr << "TVD_DEFECT_MEMORY iteration=" << iter
                          << " previous_valid=" << committed_tvd_deferred_valid
                          << " weight=" << controls.tvd_deferred_memory.previous_weight
                          << " source_relative=" << h.tvd_memory_source_relative << '\n';
        }
        double flux_reference = 0.0;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
            flux_reference += std::abs(mass_flux(f));
        h.continuity_normalized = l1 / std::max(flux_reference, tiny);
        h.conservative_face_flux_l1 = flux_reference;
        if (controls.algorithm == PressureVelocityAlgorithm::COUPLED)
            h.pressure_residual = h.continuity_normalized;
        if (iter == 1 && !controls.initial_face_mass_flux_from_paired_checkpoint)
            result.reference_momentum_residual =
                std::max(final_momentum_residual, tiny);
        h.momentum_equation_residual_reference =
            std::max(result.reference_momentum_residual, tiny);
        h.momentum_equation_residual_relative =
            final_momentum_residual / h.momentum_equation_residual_reference;
        h.velocity_change_inf = velocity_change_inf;
        h.pressure_change_inf = pressure_change_inf;
        // Report the relaxation factors actually used on this attempt, including
        // any retry-controller scaling. These are the values needed to replay
        // an accepted iteration deterministically from its history record.
        // The monolithic COUPLED candidate is applied without post-solve
        // damping. Report its actual update factor, not the segregated
        // alpha_u/alpha_p controls, which do not scale this block solution.
        h.effective_alpha_u = controls.algorithm == PressureVelocityAlgorithm::COUPLED
            ? coupled_outer_step_applied_v43 : retry_controller.alpha_u(effective_alpha_u);
        h.effective_alpha_p = controls.algorithm == PressureVelocityAlgorithm::COUPLED
            ? coupled_outer_step_applied_v43 : retry_controller.alpha_p(effective_alpha_p);
        h.simplec_signed_offdiag_sum_linf = simplec_signed_offdiag_sum_audit;
        h.simplec_abs_offdiag_sum_linf = simplec_abs_offdiag_sum_audit;
        h.simplec_denominator_signed_linf = simplec_denominator_signed_audit;
        h.simplec_denominator_abs_linf = simplec_denominator_abs_audit;
        h.simplec_denominator_gap_linf = simplec_denominator_gap_audit;
        h.simplec_denominator_worst_cell = simplec_denominator_worst_cell_audit;
        h.momentum_linear_iterations =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED
                ? 0 : std::max({rx.iterations, ry.iterations, rz.iterations});
        h.pressure_linear_iterations =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED ? 0 : pressure_iterations;
        h.pressure_correctors_used =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED ? 0 : pressure_correctors;
        h.piso_neighbor_predictor_rebuilds = piso_neighbor_predictor_rebuilds;
        h.piso_neighbor_predictor_flux_delta_linf =
            piso_neighbor_predictor_flux_delta_linf;
        h.linear_tolerance_used = controls.algorithm == PressureVelocityAlgorithm::COUPLED
            ? controls.coupling.coupled_linear_tolerance : linear_tolerance;
        h.corrected_flux_continuity_linf = linf;
        // Independent flux/velocity audit (restores the measurement channel
        // that commit 7c7dadf5 made tautological with
        // reconstructed_flux = mass_flux). Two reconstruction paths are
        // compared:
        //   channel RC (authoritative): mass_flux, the Rhie-Chow flux built
        //     from HbyA/p plus, for SIMPLEC, the consistent-response term;
        //   channel U (independent): the flux rebuilt from the corrected
        //     cell velocity with the plain linear interpolation of
        //     make_mass_flux, which shares no operator with Rhie-Chow.
        // The mismatch is split internal / fixed-pressure boundary faces:
        // on fixed-pressure boundary faces the two channels differ
        // analytically by the Rhie-Chow boundary correction plus the
        // predictor/corrector velocity response, so a persistent value there
        // is an expected artifact of the boundary contract, not a
        // conservation error. A persistent internal-face mismatch at
        // convergence is an evidence channel about the interpolation
        // identity. Neither channel is an acceptance metric: no criterion
        // consumes them (see the N9/S9 test suite for the sensitivity
        // proofs).
        auto reconstructed_flux = make_mass_flux(
            mesh, geometry, U, controls.density, velocity_bcs,
            {}, controls.pressure_projected_face_interpolation);
        double reconstructed_linf = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
            double div_reconstructed = 0.0;
            const Offset off = cell_offsets[c];
            const Offset count = cell_offsets[c+1] - off;
            for (Offset k = 0; k < count; ++k) {
                const std::size_t f = cell_faces[off+k];
                div_reconstructed +=
                    mesh.ownership().owner(f) == c ? reconstructed_flux(f) : -reconstructed_flux(f);
            }
            reconstructed_linf = std::max(reconstructed_linf, std::abs(div_reconstructed));
        }
        double flux_mismatch_linf = 0.0;
        double flux_mismatch_boundary_linf = 0.0;
        double flux_mismatch_scale = 0.0;
        double max_velocity_abs = 0.0;
        double max_face_area = 0.0;
        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
            for (std::size_t d = 0; d < 3; ++d)
                max_velocity_abs = std::max(
                    max_velocity_abs, std::abs(U.component_data(d)[c]));
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            const double mismatch =
                std::abs(mass_flux(f) - reconstructed_flux(f));
            flux_mismatch_scale = std::max(flux_mismatch_scale, std::abs(mass_flux(f)));
            max_face_area = std::max(
                max_face_area, geometry.face_area_vectors[f].mag());
            if (mesh.ownership().neighbour(f) >= 0) {
                flux_mismatch_linf = std::max(flux_mismatch_linf, mismatch);
                continue;
            }
            const std::size_t patch = geometry.face_patch[f];
            if (patch >= mesh.boundary().n_patches())
                continue;
            if (mesh.boundary().patch(patch).type == PatchType::EMPTY)
                continue;
            const auto& patch_name = mesh.boundary().patch(patch).name;
            const auto pbc = pressure_bcs.find(patch_name);
            if (pbc != pressure_bcs.end() &&
                pbc->second.type == ScalarBoundaryType::FIXED_VALUE)
                flux_mismatch_boundary_linf =
                    std::max(flux_mismatch_boundary_linf, mismatch);
        }
        // Normalized internal mismatch with an explicit physical floor:
        // E_phi = ||phi_RC - phi_U||_inf / max(||phi_RC||_inf,
        // rho * max|U| * max|Sf|). Never divides by a vanishing scale.
        const double flux_mismatch_denominator = std::max(
            flux_mismatch_scale,
            controls.density * max_velocity_abs * max_face_area);
        const double flux_mismatch_normalized =
            flux_mismatch_denominator > 0.0
                ? flux_mismatch_linf / flux_mismatch_denominator
                : 0.0;
        h.reconstructed_velocity_continuity_linf = reconstructed_linf;
        h.flux_velocity_mismatch_linf = flux_mismatch_linf;
        h.flux_velocity_mismatch_boundary_linf = flux_mismatch_boundary_linf;
        h.flux_velocity_mismatch_normalized = flux_mismatch_normalized;
        if (v20_relaxed_pressure_flux) {
            const auto decomposed=audit_v20_flux_velocity_decomposition(
                mesh,mass_flux,*v20_relaxed_pressure_flux,reconstructed_flux);
            std::cerr << std::scientific << std::setprecision(13)
                      << "V20_FLUX_VELOCITY_DECOMPOSITION iteration=" << iter
                      << " total_linf=" << decomposed.total_linf
                      << " pressure_split_linf=" << decomposed.pressure_split_linf
                      << " rc_vs_velocity_linf=" << decomposed.rhie_chow_linf
                      << " total_l1=" << decomposed.total_l1
                      << " pressure_split_l1=" << decomposed.pressure_split_l1
                      << " rc_vs_velocity_l1=" << decomposed.rhie_chow_l1
                      << " triangle_error_linf=" << decomposed.triangle_error_linf
                      << " internal_linf=" << decomposed.internal_total_linf
                      << " boundary_linf=" << decomposed.boundary_total_linf
                      << " total_worst_face=" << decomposed.total_worst_face
                      << " split_worst_face=" << decomposed.split_worst_face
                      << " rc_worst_face=" << decomposed.rc_worst_face
                      << " total_worst_boundary=" << decomposed.total_worst_boundary
                      << " qualified=0\n";
        }
        h.rhie_chow_identity_linf =
            rhie_chow_identity_evaluated ? rhie_chow_identity_audit_linf : 0.0;
        h.simplec_plain_rau_flux_gap_linf = simplec_plain_rau_gap;
        h.simplec_plain_rau_worst_face = simplec_plain_rau_worst;
        h.momentum_equation_residual_components =
            std::array<double, 3>{rx_diag.global, ry_diag.global, rz_diag.global};
        h.momentum_equation_residual_internal =
            std::max({rx_diag.interior, ry_diag.interior, rz_diag.interior});
        h.momentum_equation_residual_boundary =
            std::max({rx_diag.boundary, ry_diag.boundary, rz_diag.boundary});
        h.pressure_gradient_linf = pressure_gradient_linf;
        h.pressure_gradient_l2 = pressure_gradient_l2;
        const ResidualDiagnostic* worst_diag = &rx_diag;
        if (ry_diag.global > worst_diag->global) worst_diag = &ry_diag;
        if (rz_diag.global > worst_diag->global) worst_diag = &rz_diag;
        h.momentum_residual_cell = worst_diag->max_cell;
        h.momentum_residual_no_pressure = worst_diag->no_pressure;
        h.momentum_pressure_contribution = worst_diag->pressure_contribution;
        h.momentum_residual_patch = worst_diag->patch;

        // Optional forensic microscope for one cell. This is deliberately
        // disabled unless CFDX_DEBUG_CELL is set, so production/CI output is
        // unchanged. Use CFDX_DEBUG_CELL=<index> for a fixed cell or
        // CFDX_DEBUG_CELL=auto to inspect the independently reassembled
        // worst-residual cell. The dump is emitted once per algorithm at
        // convergence or at the configured iteration limit.
        const bool debug_cell_auto = controls.diagnostics.debug_cell_auto;
        std::size_t debug_cell = 0;
        const bool debug_cell_enabled =
            debug_cell_auto ||
            (controls.diagnostics.debug_cell >= 0 &&
             static_cast<std::size_t>(controls.diagnostics.debug_cell) < mesh.n_cells());
        if (debug_cell_auto)
            debug_cell = h.momentum_residual_cell;
        else if (debug_cell_enabled)
            debug_cell = static_cast<std::size_t>(controls.diagnostics.debug_cell);

        // The first-iteration-relative equation residual is diagnostic only:
        // if the initial residual is already at round-off, its denominator is
        // ill-conditioned and the ratio can grow while the actual equation
        // residual remains negligible. Use the RHS-normalized residual for the
        // nonlinear convergence monitor and keep the reference-relative value
        // available in history for diagnostics/adaptive-relaxation analysis.
        // The first-iteration normalization is useful only when its reference
        // is resolved relative to the momentum RHS. If the reference is at
        // round-off, keep it for diagnostics but do not let its ratio dominate
        // the nonlinear monitor.
        // The reference-relative metric is meaningful only when its
        // reference is resolved above the declared physical RHS-relative
        // tolerance (bounded below by the FP64 round-off floor). A tiny but
        // nonzero reference can otherwise pin the ratio near one despite an
        // already-satisfied physical residual.
        const double momentum_reference_resolution_threshold =
            std::max(
                100.0 * std::numeric_limits<double>::epsilon(),
                controls.convergence.relative_tolerance) *
            std::max(h.momentum_rhs_scale, tiny);
        const bool momentum_reference_ill_conditioned =
            h.momentum_equation_residual_reference <=
            momentum_reference_resolution_threshold;
        // COUPLED measures this residual after the monolithic solve, so its
        // first-iteration reference is not an independent pre-solve baseline.
        // Exclude that ratio for COUPLED regardless of its magnitude, while
        // keeping the independently RHS-normalized physical gate mandatory.
        const bool use_reference_relative_metric =
            !controls.convection_continuation.enabled &&
            !momentum_reference_ill_conditioned &&
            controls.algorithm != PressureVelocityAlgorithm::COUPLED;
        double nonlinear_metric = !use_reference_relative_metric
            ? std::max({
                  h.momentum_residual,
                  h.pressure_residual,
                  h.velocity_change_inf,
                  h.pressure_change_inf})
            : std::max({
                  h.momentum_residual,
                  h.momentum_equation_residual_relative,
                  h.pressure_residual,
                  h.velocity_change_inf,
                  h.pressure_change_inf});
        h.nonlinear_convergence_metric = nonlinear_metric;

        if (controls.diagnostics.iteration_trace &&
            (iter == 1 ||
             iter % std::max<std::size_t>(1, controls.diagnostics.iteration_trace_frequency) == 0)) {
            // Diagnostic-only trace: expose the authoritative nonlinear history
            // before the convergence gate. No criterion or tolerance is changed.
            std::cerr << "INCOMPRESSIBLE_ITER iter=" << iter
                      << " momentum_residual_rhs_relative=" << h.momentum_equation_residual_rhs_relative
                      << " momentum_residual_reference_relative=" << h.momentum_equation_residual_relative
                      << " momentum_residual_raw=" << h.momentum_equation_residual
                      << " momentum_rhs_scale=" << h.momentum_rhs_scale
                      << " momentum_residual_reference=" << h.momentum_equation_residual_reference
                      << " momentum_active_component_residual_max="
                      << h.momentum_equation_residual_rhs_active_component_max
                      << " momentum_rhs_scale_components=("
                      << h.momentum_rhs_scale_components[0] << ","
                      << h.momentum_rhs_scale_components[1] << ","
                      << h.momentum_rhs_scale_components[2] << ")"
                      << " momentum_rhs_relative_components=("
                      << h.momentum_equation_residual_rhs_relative_components[0] << ","
                      << h.momentum_equation_residual_rhs_relative_components[1] << ","
                      << h.momentum_equation_residual_rhs_relative_components[2] << ")"
                      << " continuity_normalized=" << h.continuity_normalized
                      << " velocity_change_inf=" << h.velocity_change_inf
                      << " pressure_change_inf=" << h.pressure_change_inf
                      << " alpha_u=" << h.effective_alpha_u
                      << " alpha_p=" << h.effective_alpha_p
                      << " pressure_linear_iterations=" << pressure_iterations
                      << " pressure_relative_residual=" << pressure_residual
                      << " momentum_linear_iterations="
                      << ((controls.algorithm == PressureVelocityAlgorithm::COUPLED)
                              ? 0 : std::max({rx.iterations, ry.iterations, rz.iterations}))
                      << " flux_velocity_mismatch_linf=" << h.flux_velocity_mismatch_linf
                      << " worst_momentum_cell=" << h.momentum_residual_cell
                      << " worst_momentum_patch=" << h.momentum_residual_patch
                      << " high_order_fraction=" << h.high_order_fraction
                      << " full_order_eligible=" << full_order_eligible
                      << "\n";
        }
        // Acceptance stays per-criterion. The monitor divides its synthetic
        // maximum by the maximum seen at the first iteration, so a metric whose
        // absolute scale sits far below that first maximum would be exempted
        // from the gate. nonlinear_metric is therefore a diagnostic summary and
        // never a replacement for the individual physical criteria.
        const bool momentum_residual_gate = momentum_component_certificate.converged;
        // COUPLED solves the monolithic momentum-continuity system before
        // this post-update residual is measured. Using that already-solved
        // residual as its own first-iteration reference makes the ratio
        // approximately one even when the independently recomputed true
        // residual is at the linear-solver floor. The RHS-normalized physical
        // residual gate above remains mandatory for COUPLED; the
        // post-update/reference ratio is not an independent convergence test
        // for this direct block solve.
        const bool momentum_equation_residual_gate =
            controls.algorithm == PressureVelocityAlgorithm::COUPLED ||
            // A pre-ramp reference is not a residual of the final operator.
            // During continuation the independently normalized per-component
            // physical momentum certificate is the mandatory convergence gate.
            controls.convection_continuation.enabled ||
            momentum_reference_ill_conditioned ||
            h.momentum_equation_residual_relative <= controls.convergence.relative_tolerance;
        const bool pressure_residual_gate =
            h.pressure_residual <= controls.convergence.relative_tolerance;
        const bool velocity_change_gate =
            h.velocity_change_inf <= controls.convergence.relative_tolerance;
        const bool pressure_change_gate =
            h.pressure_change_inf <= controls.convergence.relative_tolerance;
        // The convergence residual is always the RAW TVD flux defect; with the
        // defect memory enabled, its lag-filtered source must also have
        // vanished before the iteration may be declared converged.
        const bool tvd_memory_consistency_gate =
            !controls.tvd_deferred_memory.enabled ||
            !controls.tvd_deferred_memory.require_source_consistency ||
            h.tvd_memory_source_relative <= controls.convergence.relative_tolerance;
        const bool per_criterion_gates_ok =
            momentum_residual_gate && momentum_equation_residual_gate &&
            pressure_residual_gate && velocity_change_gate && pressure_change_gate &&
            mixed_modes_stable && tvd_memory_consistency_gate;

        // Reject nonfinite accepted-state physical metrics, even during a
        // continuation ramp where previous-operator ratios cannot be used.
        // The final per-component/face/continuity convergence gates remain
        // mandatory, regardless of this globalization policy.
        if (!std::isfinite(h.momentum_equation_residual_rhs_relative) ||
            !std::isfinite(h.continuity_normalized))
            throw NonlinearRetryableFailure(
                "solve_steady_incompressible: nonfinite nonlinear physical residual");

        // Multisecant candidates are only accepted after a complete SIMPLE
        // momentum/pressure/flux update, the full TVD residual reassembly and
        // independent continuity audit above. A poor predictor retries this
        // same outer iteration without mixing; rejected samples NEVER enter
        // the history or the final force certificate.
        if (multisecant_trial_used && !result.history.empty()) {
            const auto& previous=result.history.back();
            const double cap=controls.nonlinear_multisecant.maximum_physical_growth;
            const double mom_floor=10.0*controls.convergence.relative_tolerance;
            const double cont_floor=10.0*controls.convergence.continuity_tolerance;
            if (h.momentum_equation_residual_rhs_active_component_max >
                    cap*std::max(previous.momentum_equation_residual_rhs_active_component_max,
                                 mom_floor) ||
                h.continuity_normalized >
                    cap*std::max(previous.continuity_normalized,cont_floor)) {
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "MULTISECANT_REJECT iteration=" << iter
                              << " momentum=" << h.momentum_equation_residual_rhs_active_component_max
                              << " prior=" << previous.momentum_equation_residual_rhs_active_component_max
                              << " continuity=" << h.continuity_normalized
                              << " prior_continuity=" << previous.continuity_normalized << '\n';
                // A trial failure (including one before the monitor below)
                // has a single recovery path in the catch handler. The monitor
                // never observes a rejected candidate.
                throw NonlinearRetryableFailure(
                    "multisecant physical-residual growth; retry ordinary Picard");
            }
        }
        h.multisecant_proposed=multisecant_proposal.available;
        h.multisecant_accepted=multisecant_trial_used;
        if (multisecant_trial_used) {
            h.multisecant_relative_step=multisecant_proposal.relative_step;
            h.multisecant_estimated_objective=
                multisecant_proposal.estimated_objective;
        }

        // Segregated SIMPLE only: an opt-in physical-residual globalization
        // guard. The catch below restores U, p, mass_flux, monitor and SST
        // and shrinks the equation relaxation before replaying this iteration.
        // A residual from a different convection fraction is not a valid
        // reference for a growth ratio: perform this check only when the
        // previous ACCEPTED state used the same discrete operator. This
        // protects every UPWIND/full-TVD step, including after the ramp.
        // Do NOT apply independent post-solve U/p relaxation to COUPLED.
        if (controls.nonlinear_retry.reject_severe_physical_growth &&
            controls.algorithm == PressureVelocityAlgorithm::SIMPLE &&
            !result.history.empty() && same_convection_operator_as_previous) {
            const auto& previous = result.history.back();
            const bool momentum_growth = severe_physical_residual_growth(
                previous.momentum_equation_residual_rhs_relative,
                h.momentum_equation_residual_rhs_relative,
                controls.convergence.relative_tolerance,
                controls.nonlinear_retry);
            const bool continuity_growth = severe_physical_residual_growth(
                previous.continuity_normalized,
                h.continuity_normalized,
                controls.convergence.continuity_tolerance,
                controls.nonlinear_retry);
            if (momentum_growth || continuity_growth) {
                std::ostringstream reason;
                reason << "SIMPLE severe physical residual growth: momentum="
                       << previous.momentum_equation_residual_rhs_relative
                       << "->" << h.momentum_equation_residual_rhs_relative
                       << ", continuity=" << previous.continuity_normalized
                       << "->" << h.continuity_normalized;
                throw NonlinearRetryableFailure(reason.str());
            }
        }
        // Fail-closed globalization of the real *unrelaxed* physical operator.
        // The candidate above is fully reassembled with its current TVD face
        // reconstruction and its authoritative conservative Rhie-Chow flux.
        // Rejecting retries the SAME outer step with smaller momentum and
        // pressure equation relaxations, recomputing H/A and phi together.
        // Never interpolate phi independently of H/A and never accept a ramped
        // equation as a final-order convergence certificate.
        if (controls.picard_physical_globalization.enabled &&
            !result.history.empty() && same_convection_operator_as_previous &&
            iter >= controls.picard_physical_globalization.start_iteration &&
            full_order_eligible) {
            const auto& previous = result.history.back();
            const auto sample_of = [](const IncompressibleIteration& item) {
                PicardPhysicalSample sample;
                for (std::size_t d = 0; d < 3; ++d) {
                    sample.momentum_rhs_scale[d] = item.momentum_rhs_scale_components[d];
                    sample.momentum_absolute[d] =
                        item.momentum_equation_residual_rhs_relative_components[d] *
                        item.momentum_rhs_scale_components[d];
                }
                sample.continuity_l1 = item.continuity_l1;
                sample.conservative_face_flux_l1 = item.conservative_face_flux_l1;
                return sample;
            };
            const auto baseline = sample_of(previous);
            const auto current = sample_of(h);
            const auto after = picard_physical_merit(
                current, baseline, controls.convergence.relative_tolerance,
                controls.convergence.continuity_tolerance);
            PicardPhysicalMerit before;
            const std::size_t accepted_count = result.history.size();
            const std::size_t begin = accepted_count >
                controls.picard_physical_globalization.nonmonotone_window
                ? accepted_count - controls.picard_physical_globalization.nonmonotone_window : 0;
            for (std::size_t j = begin; j < accepted_count; ++j) {
                // Never combine residuals from two different discrete PDEs.
                if (result.history[j].high_order_fraction !=
                    previous.high_order_fraction) continue;
                const auto historic = picard_physical_merit(
                    sample_of(result.history[j]), baseline,
                    controls.convergence.relative_tolerance,
                    controls.convergence.continuity_tolerance);
                before.momentum = std::max(before.momentum, historic.momentum);
                before.continuity = std::max(before.continuity, historic.continuity);
                before.combined = std::max(before.combined, historic.combined);
            }
            if (!(before.combined > 0.0))
                throw std::runtime_error(
                    "Picard globalization: no comparable historical physical merit");
            const double step_scale = std::min(
                retry_controller.alpha_u(effective_alpha_u)/effective_alpha_u,
                retry_controller.alpha_p(effective_alpha_p)/effective_alpha_p);
            h.globalization_merit = after.combined;
            h.globalization_backtracks = retry_controller.retries();
            const bool accepted = picard_physical_accept(
                before, after, step_scale, controls.picard_physical_globalization);
            if (controls.diagnostics.iteration_trace || !accepted)
                std::cerr << "PICARD_PHYSICAL_GLOBALIZATION iteration=" << iter
                          << " trial=" << retry_controller.retries()
                          << " scale=" << step_scale
                          << " before=" << before.combined
                          << " after=" << after.combined
                          << " momentum_before=" << before.momentum
                          << " momentum_after=" << after.momentum
                          << " continuity_before=" << before.continuity
                          << " continuity_after=" << after.continuity
                          << " decision=" << (accepted ? "ACCEPT" : "REJECT") << '\n';
            if (!accepted)
                throw NonlinearRetryableFailure(
                    "Picard physical globalization: candidate increases or fails to "
                    "decrease the reassembled full-order residual");
        }
        // Reject severe adaptive overshoots before they reach the convergence
        // monitor. A retry is the same nonlinear iteration, so convergence
        // history must only observe the accepted state.
        // The same-operator contract applies to adaptive degradation as
        // well: do not punish the solver for changing its discretization.
        if (controls.adaptive_relaxation.enabled && !result.history.empty() &&
            same_convection_operator_as_previous) {
            const auto& previous = result.history.back();
            const bool severe_u =
                std::isfinite(previous.momentum_equation_residual_relative) &&
                std::isfinite(h.momentum_equation_residual_relative) &&
                previous.momentum_equation_residual_relative > 0.0 &&
                adaptive_relaxation_severe_degradation(
                    previous.momentum_equation_residual_relative,
                    h.momentum_equation_residual_relative,
                    previous.effective_alpha_u,
                    h.effective_alpha_u,
                    controls.adaptive_relaxation);
            const double previous_pressure_merit = v21_pressure_relaxation_merit(
                controls.v19_full_pressure_flux_step,
                previous.continuity_normalized,
                previous.momentum_equation_residual_rhs_relative,
                previous.flux_velocity_mismatch_normalized);
            const double current_pressure_merit = v21_pressure_relaxation_merit(
                controls.v19_full_pressure_flux_step,
                h.continuity_normalized,
                h.momentum_equation_residual_rhs_relative,
                h.flux_velocity_mismatch_normalized);
            const bool severe_p = previous_pressure_merit > 0.0 &&
                adaptive_relaxation_severe_degradation(
                    previous_pressure_merit,
                    current_pressure_merit,
                    previous.effective_alpha_p,
                    h.effective_alpha_p,
                    controls.adaptive_relaxation);

            if (severe_u || severe_p) {
                if (!retry_controller.can_retry())
                    throw NonlinearRetryableFailure(
                        "solve_steady_incompressible: adaptive relaxation severe overshoot "
                        "retry limit exhausted");
                if (severe_u) {
                    effective_alpha_u = std::max(
                        controls.adaptive_relaxation.min_alpha_u,
                        effective_alpha_u - controls.adaptive_relaxation.decrease_step);
                    adaptive_u_improvement_streak = 0;
                    adaptive_u_degradation_streak = 0;
                    adaptive_u_recovery_cooldown =
                        controls.adaptive_relaxation.recovery_cooldown_windows;
                }
                if (severe_p) {
                    effective_alpha_p = std::max(
                        controls.adaptive_relaxation.min_alpha_p,
                        effective_alpha_p - controls.adaptive_relaxation.decrease_step);
                    adaptive_p_improvement_streak = 0;
                    adaptive_p_degradation_streak = 0;
                    adaptive_p_recovery_cooldown =
                        controls.adaptive_relaxation.recovery_cooldown_windows;
                }
                transaction.reject();
                mass_flux = mass_flux_snapshot;
                pressure_low_state=pressure_low_snapshot;
                if (controls.sst_runtime && sst_attempt_snapshot)
                    restore_sst_production_state(*controls.sst_runtime->state, *sst_attempt_snapshot);
                if (controls.sst_runtime && v30_ptc_attempt_snapshot)
                    controls.sst_runtime->coupled_ptc_v30 = *v30_ptc_attempt_snapshot;
                frozen_state_valid = false;
                // The severe guard has already reduced the responsible effective
                // relaxation factor. Count the retry without applying the generic
                // retry-controller shrink a second time.
                retry_controller.retry_without_relaxation();
                continue;
            }
        }
        std::map<std::string, double> qoi_values;
        if (!controls.qoi_gates.empty()) {
            for (const auto& gate : controls.qoi_gates)
                qoi_values[gate.name] =
                    sample_incompressible_probe(
                        *std::find_if(
                            controls.probes.begin(), controls.probes.end(),
                            [&](const IncompressiblePointProbe& probe) {
                                return probe.name == gate.name;
                            }),
                        mesh, geometry, U, p);
        }
        SSTProductionStepResult sst_step;
        std::optional<SSTWallOperatorV17> post_wall_operator_v17;
        const std::vector<double>* post_wall_face_mu_v17=nullptr;
        std::optional<SstOpenFoamWallTripletV15> pending_wall_triplet_v15;
        MomentumResidualCertificate v8_refreshed_momentum =
            momentum_component_certificate;
        double sst_nu_t_relative_change = 0.0;
        if (controls.sst_runtime) {
            std::map<std::string, std::array<double, 3>> fixed_velocity_values;
            for (const auto& [patch, bc] : velocity_bcs) {
                if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    fixed_velocity_values.emplace(
                        patch, std::array<double, 3>{bc.value.x, bc.value.y, bc.value.z});
            }
            if (std::getenv("CFDX_SST_WALL_AUDIT") != nullptr &&
                !controls.sst_runtime->wall_omega_sst_patches.empty() &&
                (iter==1 || iter%5==0 || iter==controls.convergence.max_iterations)) {
                const auto wall_audit=audit_sst_wall_omega_discretization(
                    mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                    controls.sst_runtime->omega_face_values,U,
                    controls.sst_runtime->state->omega,fixed_velocity_values,
                    controls.sst_runtime->transport);
                std::cerr << std::scientific << std::setprecision(9)
                          << "SST_WALL_DISCRETE_AUDIT iteration=" << iter
                          << " faces=" << wall_audit.faces
                          << " d_min=" << wall_audit.min_distance
                          << " d_max=" << wall_audit.max_distance
                          << " lateral_normal_max=" << wall_audit.max_tangential_geometry_ratio
                          << " omega_face_min=" << wall_audit.min_face_value
                          << " omega_face_max=" << wall_audit.max_face_value
                          << " nasa_firstcell_min=" << wall_audit.min_nasa_first_cell_value
                          << " nasa_firstcell_max=" << wall_audit.max_nasa_first_cell_value
                          << " owner_to_nasa60_min=" << wall_audit.min_owner_to_first_cell_ratio
                          << " owner_to_nasa60_max=" << wall_audit.max_owner_to_first_cell_ratio
                          << " yplus_viscous_proxy_faces=" << wall_audit.velocity_proxy_faces
                          << " yplus_viscous_proxy_min=" << (wall_audit.velocity_proxy_faces
                              ? wall_audit.min_viscous_yplus_proxy : -1.0)
                          << " yplus_viscous_proxy_max=" << (wall_audit.velocity_proxy_faces
                              ? wall_audit.max_viscous_yplus_proxy : -1.0)
                          << " qualified=0\n";
                try {
                    const auto openfoam=audit_openfoam_omega_wall_reference(
                        mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                        controls.sst_runtime->omega_face_values,U,
                        controls.sst_runtime->state->k,fixed_velocity_values,
                        controls.sst_runtime->transport);
                    std::cerr << "OPENFOAM_OMEGA_WALL_REFERENCE iteration=" << iter
                              << " faces=" << openfoam.face_count
                              << " yplus_k_min=" << openfoam.minimum_k_based_y_plus
                              << " yplus_k_max=" << openfoam.maximum_k_based_y_plus
                              << " viscous_faces=" << openfoam.viscous_faces
                              << " buffer_faces=" << openfoam.buffer_faces
                              << " log_faces=" << openfoam.logarithmic_faces
                              << " max_face_omega_relative_gap=" << openfoam.worst_relative_omega_gap
                              << " worst_face=" << openfoam.worst_face
                              << " qualified=0\n";
                } catch (const std::invalid_argument& e) {
                    // This comparison is for qualification only. A moving
                    // wall or unknown velocity is not an OpenFOAM fixed-wall
                    // reference, but must not alter the production equation.
                    std::cerr << "OPENFOAM_OMEGA_WALL_REFERENCE iteration=" << iter
                              << " comparable=0 reason=" << e.what()
                              << " qualified=0\n";
                }
                // Reference-only triple audit: wall nut / omega / G are
                // linked. OpenFOAM's corner weighting and patch constraints
                // are deliberately NOT assumed equal to CFDX SST.
                try {
                    const auto triple = audit_openfoam_wall_coupling_v8(
                        mesh, geometry,
                        controls.sst_runtime->wall_omega_sst_patches,
                        U, *controls.sst_runtime->state,
                        fixed_velocity_values,
                        controls.sst_runtime->transport);
                    std::cerr << "SST_WALL_COUPLING_V8 iteration=" << iter
                              << " faces=" << triple.faces
                              << " omega_owner_gap=" << triple.max_omega_owner_gap
                              << " nut_spalding_gap=" << triple.max_nut_wall_gap
                              << " G_proxy_gap=" << triple.max_production_gap
                              << " yplus_k_max=" << triple.max_k_yplus
                              << " yplus_spalding_max=" << triple.max_spalding_yplus
                              << " worst_face=" << triple.worst_face
                              << " wall_owner_cells=" << triple.wall_owner_cells
                              << " corner_owner_cells=" << triple.multi_face_wall_owner_cells
                              << " omega_weighted_gap=" << triple.weighted_omega_owner_gap
                              << " nut_owner_proxy_weighted_gap=" << triple.weighted_nut_owner_proxy_gap
                              << " G_owner_proxy_weighted_gap=" << triple.weighted_G_owner_proxy_gap
                              << " openfoam_equivalent=" << triple.openfoam_equivalent
                              << " qualified=0\n";
                } catch (const std::exception& e) {
                    std::cerr << "SST_WALL_COUPLING_V8 iteration=" << iter
                              << " comparable=0 reason=" << e.what()
                              << " qualified=0\n";
                }
            }
            if (std::getenv("CFDX_BFS_STAGE_TRACE"))
                std::cerr << "BFS_PHASE_V10 iteration=" << iter
                          << " phase=SST_TRANSPORT_BEGIN budget="
                          << controls.sst_runtime->max_transport_iterations
                          << " target=" << controls.sst_runtime->transport_tolerance
                          << "\n";
            SSTWallRefreshV17 wall_refresh_v17;
            if (controls.sst_runtime->active_openfoam_wall_operator_v17) {
                wall_refresh_v17=[&](const auto& candidate_k,const auto& candidate_omega,
                                     ScalarBoundaryFaceValues& omega_faces,
                                     SSTWallProductionV17& production) {
                    SSTProductionState trial=*controls.sst_runtime->state;
                    trial.k=candidate_k;
                    trial.omega=candidate_omega;
                    cfdx::core::Field<double,cfdx::core::Location::CELL> F2(
                        mesh.n_cells(),"F2_wall_v17","1",1);
                    for (std::size_t cell=0;cell<mesh.n_cells();++cell) {
                        const auto blend=compute_sst_blending(
                            trial.k(cell),trial.omega(cell),trial.wall_distance(cell),
                            controls.kinematic_viscosity,
                            controls.sst_runtime->transport.beta_star,0.0,
                            controls.density,controls.sst_runtime->transport.sst_sigma_w2);
                        F2(cell)=blend.second;
                    }
                    refresh_sst_production_eddy_viscosity(
                        trial.k,trial.omega,trial.strain_rate,F2,
                        controls.sst_runtime->transport,trial.nu_t);
                    const auto ref=prepare_openfoam_wall_triplet_v15(
                        mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                        U,trial,fixed_velocity_values,controls.sst_runtime->transport,
                        controls.sst_runtime->omega_wall_corner_face_count_v19);
                    auto wall=make_sst_wall_operator_v17(
                        mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                        ref,controls.sst_runtime->transport,
                controls.sst_runtime->omega_wall_owner_constraints_v18);
                    omega_faces=std::move(wall.omega_faces);
                    production=std::move(wall.production);
                };
            }
            sst_step = advance_sst_production_iterate_transactionally(
                mesh, U, p, mass_flux, *controls.sst_runtime->state,
                controls.sst_runtime->transport, fixed_velocity_values,
                controls.sst_runtime->wall_patch_names,
                controls.sst_runtime->k_boundary_conditions,
                controls.sst_runtime->omega_boundary_conditions,
                controls.sst_runtime->max_transport_iterations,
                controls.sst_runtime->transport_tolerance,
                &controls.sst_runtime->omega_face_values,
                controls.sst_runtime->allow_bounded_provisional_transport,
                controls.sst_runtime->active_openfoam_wall_operator_v17
                    ? &wall_refresh_v17 : nullptr);
            if (std::getenv("CFDX_BFS_STAGE_TRACE"))
                std::cerr << "BFS_PHASE_V10 iteration=" << iter
                          << " phase=SST_TRANSPORT_END committed=" << sst_step.committed
                          << " converged=" << sst_step.transport.converged
                          << " iterations=" << sst_step.transport.iterations
                          << " k_residual=" << sst_step.transport.k_residual
                          << " omega_residual=" << sst_step.transport.second_residual
                          << " k_floor_events=" << sst_step.transport.k_floor_events
                          << " omega_floor_events=" << sst_step.transport.second_floor_events
                          << "\n";
            if (!sst_step.committed)
                throw NonlinearRetryableFailure(
                    "solve_steady_incompressible: SST update rejected: " + sst_step.reason);
            if (controls.sst_runtime->active_openfoam_wall_operator_v17) {
                const auto ref=prepare_openfoam_wall_triplet_v15(
                    mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                    U,*controls.sst_runtime->state,fixed_velocity_values,
                    controls.sst_runtime->transport,
                    controls.sst_runtime->omega_wall_corner_face_count_v19);
                post_wall_operator_v17=make_sst_wall_operator_v17(
                    mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                    ref,controls.sst_runtime->transport,
                controls.sst_runtime->omega_wall_owner_constraints_v18);
                post_wall_face_mu_v17=&post_wall_operator_v17->momentum_face_mu;
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "SST_WALL_OPERATOR_ACTIVE_V17 iteration=" << iter
                              << " faces=" << post_wall_operator_v17->faces
                              << " owners=" << post_wall_operator_v17->owners
                              << " momentum=1 omega=1 G=1"
                              << " openfoam_equivalent=0 qualified=0\n";
            }
            // One physical candidate -> one complete (nut,G,omega) reference.
            // This is stored only on outer transaction commit, never on a
            // rejected SST or rejected U/p/flux candidate. Reference-only:
            // none of the three FV boundary operators is changed in V15.
            if (controls.sst_runtime->stage_openfoam_wall_triplet_v15) {
                pending_wall_triplet_v15 = prepare_openfoam_wall_triplet_v15(
                    mesh,geometry,controls.sst_runtime->wall_omega_sst_patches,
                    U,*controls.sst_runtime->state,fixed_velocity_values,
                    controls.sst_runtime->transport,
                    controls.sst_runtime->omega_wall_corner_face_count_v19);
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "SST_WALL_TRIPLET_STAGED_V15 iteration=" << iter
                              << " faces=" << pending_wall_triplet_v15->wall_faces
                              << " owners=" << pending_wall_triplet_v15->wall_cells
                              << " corner_owners=" << pending_wall_triplet_v15->corner_cells
                              << " momentum_applied=0 sst_applied=0 qualified=0\n";
            }
            if (sst_step.provisional)
                std::cerr << "SST_PROVISIONAL_OUTER iteration=" << iter

                          << " inner_iterations=" << sst_step.transport.iterations
                          << " k_candidate_residual=" << sst_step.transport.k_residual
                          << " omega_candidate_residual=" << sst_step.transport.second_residual
                          << " qualified=0\n";
            // Turbulence convergence is an additional gate; it never relaxes
            // the existing carrier residual, continuity, or QoI criteria.
            if (!std::isfinite(sst_step.transport.k_residual) ||
                !std::isfinite(sst_step.transport.second_residual) ||
                sst_step.transport.nonfinite_bound_events != 0)
                throw NonlinearRetryableFailure(
                    "solve_steady_incompressible: SST residual/boundedness gate failed");
            // The transport solver's converged flag is based on its linear
            // solver status, bound audit, and field-update criterion. Its raw
            // equation residuals have dimensional scales and must not be mixed
            // directly with the normalized carrier convergence metric.
            if (!sst_attempt_snapshot)
                throw std::logic_error("SST production attempt is missing its state snapshot");
            sst_nu_t_relative_change = sst_relative_eddy_viscosity_change(
                *sst_attempt_snapshot, *controls.sst_runtime->state,
                controls.kinematic_viscosity);
            // The previous residual was assembled BEFORE SST changed nu_t.
            // Rebuild all three physical carrier equations with the NEW eddy
            // viscosity and the SAME accepted trial U, p, Rhie-Chow flux.
            // This is an independent nonlinear coupling gate, not a change
            // to the transport/SST equations or their convergence tolerance.
            if (controls.require_post_sst_momentum_reassembly ||
                controls.require_global_nonlinear_certificate_v8) {
                std::vector<double> trial_mu(mesh.n_cells());
                for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
                    trial_mu[cell]=controls.density*(
                        controls.kinematic_viscosity+
                        controls.sst_runtime->state->nu_t(cell));
                    if (!std::isfinite(trial_mu[cell]) || !(trial_mu[cell]>0.0))
                        throw NonlinearRetryableFailure(
                            "post-SST momentum reassembly: invalid updated viscosity");
                }
                // The finite-volume transpose of the stress uses nu_t too.
                // Replace its OLD accepted-U contribution, not the frozen
                // predictor contribution, before reassembling the equations.
                auto post_body_x=body_x,post_body_y=body_y,post_body_z=body_z;
                if (include_transpose_stress) {
                    const std::array<ScalarBoundaryConditions,3> bcs{
                        ubc_x,ubc_y,ubc_z};
                    const auto before=viscous_face_flux_divergence(
                        mesh,geometry,transpose_viscous_stress_flux(
                            mesh,geometry,U,mu_eff,bcs,
                            controls.density*controls.kinematic_viscosity,
                            active_wall_face_mu_v17));
                    const auto after=viscous_face_flux_divergence(
                        mesh,geometry,transpose_viscous_stress_flux(
                            mesh,geometry,U,trial_mu,bcs,
                            controls.density*controls.kinematic_viscosity,
                            post_wall_face_mu_v17));
                    for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
                        post_body_x(cell)+=after(cell,0)-before(cell,0);
                        post_body_y(cell)+=after(cell,1)-before(cell,1);
                        post_body_z(cell)+=after(cell,2)-before(cell,2);
                    }
                }
                auto ex_post=assemble_momentum_component(
                    mesh,geometry,mass_flux,final_grad_p,post_body_x,trial_mu,
                    ubc_x,0,controls.use_bounded_convection,
                    controls.convection_scheme,&final_ux_field,nullptr,
                    sst_molecular_wall_mu,outer_high_order_fraction,
                    final_face_velocity_conditions
                        ? &(*final_face_velocity_conditions)[0] : nullptr,
                    post_wall_face_mu_v17);
                auto ey_post=assemble_momentum_component(
                    mesh,geometry,mass_flux,final_grad_p,post_body_y,trial_mu,
                    ubc_y,1,controls.use_bounded_convection,
                    controls.convection_scheme,&final_uy_field,nullptr,
                    sst_molecular_wall_mu,outer_high_order_fraction,
                    final_face_velocity_conditions
                        ? &(*final_face_velocity_conditions)[1] : nullptr,
                    post_wall_face_mu_v17);
                auto ez_post=assemble_momentum_component(
                    mesh,geometry,mass_flux,final_grad_p,post_body_z,trial_mu,
                    ubc_z,2,controls.use_bounded_convection,
                    controls.convection_scheme,&final_uz_field,nullptr,
                    sst_molecular_wall_mu,outer_high_order_fraction,
                    final_face_velocity_conditions
                        ? &(*final_face_velocity_conditions)[2] : nullptr,
                    post_wall_face_mu_v17);
                if (controls.transient.enabled) {
                    const auto& prior=*controls.transient.previous;
                    const auto& older=controls.transient.history_valid
                        ? *controls.transient.previous_previous : prior;
                    const bool bdf2=controls.transient.scheme==
                        cfdx::core::TimeScheme::BDF2 && controls.transient.history_valid;
                    add_dual_time_mass_term(ex_post,prior,older,geometry,
                        controls.density,controls.transient.dt,
                        controls.transient.dt_previous,bdf2,0);
                    add_dual_time_mass_term(ey_post,prior,older,geometry,
                        controls.density,controls.transient.dt,
                        controls.transient.dt_previous,bdf2,1);
                    add_dual_time_mass_term(ez_post,prior,older,geometry,
                        controls.density,controls.transient.dt,
                        controls.transient.dt_previous,bdf2,2);
                }
                const std::array<const ScalarEquation*,3> equations{
                    &ex_post,&ey_post,&ez_post};
                const std::array<const Vector*,3> velocity{
                    &final_ux,&final_uy,&final_uz};
                std::array<double,3> absolute_residuals{},rhs_scales{};
                for (std::size_t component=0; component<3; ++component) {
                    const auto& eq=*equations[component];
                    rhs_scales[component]=1.e-300;
                    for (std::size_t cell=0; cell<mesh.n_cells(); ++cell)
                        rhs_scales[component]=std::max(
                            rhs_scales[component],std::abs(eq.rhs(cell)));
                    absolute_residuals[component]=scalar_equation_residual_inf(
                        eq,*velocity[component]);
                    if (!std::isfinite(absolute_residuals[component]))
                        throw NonlinearRetryableFailure(
                            "post-SST momentum reassembly: nonfinite residual");
                }
                // V45 opt-in microscope on the EXACT post-SST FV operator.
                // Never feed a compensated residual back into V7/V8, V13/V14,
                // the COUPLED solver or its convergence controls. Emit a
                // tentative marker; only GLOBAL_NONLINEAR_COMMITTED_V13
                // identifies states finally accepted by the full transaction.
                const char* post_sst_probe=std::getenv("CFDX_POST_SST_MOMENTUM_V45");
                if (post_sst_probe && std::string(post_sst_probe)=="1") {
                    for (std::size_t d=0;d<3;++d) {
                        const auto row=post_sst_momentum_hotspot_v45(
                            *equations[d],*velocity[d],d,rhs_scales[d]);
                        const auto cell=row.cell;
                        const double volume=geometry.cell_volumes[cell];
                        const double integrated_pressure_source=-volume*
                            final_grad_p.component_data(d)[cell];
                        const double integrated_body_source=volume*(d==0
                            ? post_body_x(cell) : d==1
                            ? post_body_y(cell) : post_body_z(cell));
                        std::cerr<<"POST_SST_MOMENTUM_ROW_V45 iteration="<<iter
                                 <<" component="<<d<<" cell="<<cell
                                 <<" original_signed="<<row.original_signed
                                 <<" original_abs="<<row.original_abs
                                 <<" rhs_scale="<<row.rhs_scale
                                 <<" rhs_relative="<<row.original_rhs_relative
                                 <<" diagonal_term="<<row.independent_row.diagonal_term
                                 <<" off_diagonal_terms="<<row.independent_row.off_diagonal_terms
                                 <<" rhs="<<row.independent_row.rhs
                                 <<" compensated_signed="<<row.independent_row.compensated_signed
                                 <<" pressure_source_integral="<<integrated_pressure_source
                                 <<" body_source_integral="<<integrated_body_source
                                 <<" pressure_body_rhs_remainder="
                                    <<(row.independent_row.rhs-
                                        integrated_pressure_source-integrated_body_source)
                                 <<" nu_t="<<controls.sst_runtime->state->nu_t(cell)
                                 <<" k="<<controls.sst_runtime->state->k(cell)
                                 <<" omega="<<controls.sst_runtime->state->omega(cell)
                                 <<" committed=0 qualified=0\n";
                    }
                }
                // V46: read-only face terms from the very same post-SST
                // operator. No gate, relaxation, RHS, flux or field is
                // changed; face balances are NOT used to accept a candidate.
                const char* face_probe_v46=std::getenv("CFDX_POST_SST_FACE_V46");
                if (face_probe_v46 && std::string(face_probe_v46)=="1") {
                    try {
                        const auto selection=post_sst_face_selection_v46(
                            std::getenv("CFDX_POST_SST_FACE_WINDOW_V46"),
                            std::getenv("CFDX_POST_SST_FACE_CELLS_V46"),mesh.n_cells());
                        if (selection.contains(iter)) {
                            if (controls.transient.enabled)
                                throw std::invalid_argument(
                                    "V46 steady face audit excludes physical time derivative");
                            const auto collocated=make_mass_flux(
                                mesh,geometry,U,controls.density,velocity_bcs,
                                pressure_bcs,controls.pressure_projected_face_interpolation,
                                final_face_velocity_conditions
                                    ? &*final_face_velocity_conditions : nullptr);
                            cfdx::core::Field<double,cfdx::core::Location::FACE>
                                new_transpose;
                            if (include_transpose_stress)
                                new_transpose=transpose_viscous_stress_flux(
                                    mesh,geometry,U,trial_mu,
                                    {ubc_x,ubc_y,ubc_z},sst_molecular_wall_mu,
                                    post_wall_face_mu_v17);
                            const std::array<const ScalarBoundaryConditions*,3> bcs{
                                &ubc_x,&ubc_y,&ubc_z};
                            const std::array<const cfdx::core::Field<double,
                                cfdx::core::Location::CELL>*,3> components{
                                &final_ux_field,&final_uy_field,&final_uz_field};
                            for (std::size_t d=0; d<2; ++d) {
                                cfdx::core::Field<double,cfdx::core::Location::CELL>
                                    pressure_part(mesh.n_cells(),"v46_grad_p","Pa/m",1);
                                cfdx::core::Field<double,cfdx::core::Location::CELL>
                                    physical_body(mesh.n_cells(),"v46_body","N/m3",1);
                                for (std::size_t c=0;c<mesh.n_cells();++c) {
                                    pressure_part(c)=final_grad_p.component_data(d)[c];
                                    const auto& acceleration=controls.cell_body_force.empty()
                                        ? controls.body_force : controls.cell_body_force[c];
                                    physical_body(c)=controls.density*(d==0
                                        ? acceleration.x : acceleration.y);
                                }
                                cfdx::core::Field<double,cfdx::core::Location::FACE>
                                    transpose_part;
                                const cfdx::core::Field<double,
                                    cfdx::core::Location::FACE>* trans_ptr=nullptr;
                                if (include_transpose_stress) {
                                    transpose_part=cfdx::core::Field<double,
                                        cfdx::core::Location::FACE>(
                                        mesh.n_faces(),"v46_transpose","N",1);
                                    for (std::size_t f=0;f<mesh.n_faces();++f)
                                        transpose_part(f)=new_transpose(f,d);
                                    trans_ptr=&transpose_part;
                                }
                                // Match the EXACT wall-face viscosity used
                                // by assemble_momentum_component. On resolved
                                // no-slip walls, this is molecular-only unless
                                // the V17 wall operator supplies an override.
                                const auto diagnostic_wall_face_mu=
                                    resolved_wall_face_viscosity(
                                        mesh,geometry,trial_mu,*bcs[d],
                                        sst_molecular_wall_mu,post_wall_face_mu_v17);
                                const auto audit=diagnose_post_sst_faces_v46(
                                    mesh,geometry,mass_flux,*components[d],
                                    pressure_part,physical_body,trial_mu,*bcs[d],
                                    controls.convection_scheme,outer_high_order_fraction,
                                    controls.use_bounded_convection,*equations[d],
                                    selection.cells,trans_ptr,&diagnostic_wall_face_mu,
                                    final_face_velocity_conditions
                                        ? &(*final_face_velocity_conditions)[d] : nullptr,
                                    &collocated);
                                for (const auto& record : audit.cells) {
                                    const auto& b=record.balance;
                                    std::ostringstream out;
                                    out<<std::scientific<<std::setprecision(17)
                                       <<"POST_SST_FACE_CELL_V46 iteration="<<iter
                                       <<" component="<<d<<" cell="<<record.cell
                                       <<" faces="<<record.faces.size()
                                       <<" u="<<(*components[d])(record.cell)
                                       <<" bounded_mode="<<(controls.use_bounded_convection?1:0)
                                       <<" conv="<<b.convection
                                       <<" visc="<<b.viscous
                                       <<" transpose="<<b.transpose_viscous
                                       <<" pressure="<<b.pressure
                                       <<" body="<<b.body_force
                                       <<" mass_div="<<b.mass_divergence
                                       <<" bounded="<<b.bounded_divergence
                                       <<" collocated_delta="<<b.collocated_difference
                                       <<" physical="<<b.physical_residual
                                       <<" face_operator="<<b.bounded_operator_residual
                                       <<" csr="<<b.algebraic_residual
                                       <<" closure="<<b.matrix_face_difference
                                       <<" closure_scale="<<b.matrix_face_scale
                                       <<" face_ok="<<(record.face_sum_consistent?1:0)
                                       <<" csr_ok="<<(record.matrix_face_consistent?1:0)
                                       <<" all_csr_ok="<<(audit.all_cells_matrix_face_consistent?1:0)
                                       <<" nu_t="<<controls.sst_runtime->state->nu_t(record.cell)
                                       <<" k="<<controls.sst_runtime->state->k(record.cell)
                                       <<" omega="<<controls.sst_runtime->state->omega(record.cell)
                                       <<" committed=0 qualified=0";
                                    std::cerr<<out.str()<<'\n';
                                    for (const auto& f : record.faces) {
                                        std::ostringstream flux;
                                        flux<<std::scientific<<std::setprecision(17)
                                            <<"POST_SST_FACE_V46 iteration="<<iter
                                            <<" component="<<d<<" cell="<<record.cell
                                            <<" face="<<f.face<<" owner="<<f.owner
                                            <<" neighbour="<<f.neighbour
                                            <<" sign="<<f.sign<<" patch="<<f.patch
                                            <<" mass="<<f.mass
                                            <<" conv="<<f.convection
                                            <<" visc="<<f.viscosity
                                            <<" transpose="<<f.transpose_stress
                                            <<" collocated_delta="<<f.collocated_convection_difference
                                            <<" committed=0 qualified=0";
                                        std::cerr<<flux.str()<<'\n';
                                    }
                                }
                            }
                        }
                    } catch (const std::exception& error) {
                        // Audit failures are explicit but never change the
                        // nonlinear transaction or silently certify anything.
                        std::cerr<<"POST_SST_FACE_ERROR_V46 iteration="<<iter
                                 <<" reason="<<error.what()<<" qualified=0\n";
                    }
                }
                // Use exactly the same active-component criterion as the
                // existing production momentum certificate. In particular,
                // inactive 2D/extrusion components do not acquire a different
                // artificial normalization in this second residual channel.
                v8_refreshed_momentum = certify_momentum_residual_components(
                    absolute_residuals,rhs_scales,
                    controls.convergence.relative_tolerance);
                h.post_sst_momentum_rhs_relative =
                    v8_refreshed_momentum.max_active_component_relative;
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "POST_SST_MOMENTUM iteration=" << iter
                              << " old_nut_operator="
                              << h.momentum_equation_residual_rhs_active_component_max
                              << " refreshed_nut_operator="
                              << h.post_sst_momentum_rhs_relative
                              << " relative_nut_change="
                              << sst_nu_t_relative_change << " qualified=0\n";
                // One nonlinear trial consists of both carrier AND SST. Do
                // not accept a severe cross-equation overshoot merely because
                // the old-viscosity momentum residual decreased. Keep the
                // existing retry policy and its unmodified physical tolerances.
                if (controls.nonlinear_retry.reject_severe_physical_growth &&
                    same_convection_operator_as_previous &&
                    !result.history.empty() &&
                    std::isfinite(result.history.back().post_sst_momentum_rhs_relative) &&
                    severe_physical_residual_growth(
                        result.history.back().post_sst_momentum_rhs_relative,
                        h.post_sst_momentum_rhs_relative,
                        controls.convergence.relative_tolerance,
                        controls.nonlinear_retry))
                    throw NonlinearRetryableFailure(
                        "post-SST physical momentum residual growth: retry full U/p/SST step");
            }
        }
        if (controls.require_global_nonlinear_certificate_v8) {
            double k_residual = std::numeric_limits<double>::quiet_NaN();
            double omega_residual = std::numeric_limits<double>::quiet_NaN();
            double viscosity_closure = std::numeric_limits<double>::quiet_NaN();
            if (controls.sst_runtime) {
                const auto& runtime = *controls.sst_runtime;
                const auto& state = *runtime.state;
                // Reassemble against the actual returned, coupled state:
                // this is independent of the frozen SST linear Krylov status.
                const auto sample = sst_reassembled_returned_state_residual(
                    mesh, geometry, mass_flux, state.k, state.omega,
                    state.strain_rate, state.wall_distance, runtime.transport,
                    runtime.k_boundary_conditions,
                    runtime.omega_boundary_conditions,
                    post_wall_operator_v17 ? &post_wall_operator_v17->omega_faces
                                           : &runtime.omega_face_values,
                    post_wall_operator_v17 ? &post_wall_operator_v17->production : nullptr);
                k_residual = sample.k_local_relative;
                omega_residual = sample.omega_local_relative;
                viscosity_closure = 0.;
                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                    const auto blend = compute_sst_blending(
                        state.k(cell), state.omega(cell), state.wall_distance(cell),
                        controls.kinematic_viscosity, runtime.transport.beta_star,
                        0., controls.density, runtime.transport.sst_sigma_w2);
                    const double expected = runtime.transport.a1 * state.k(cell) /
                        std::max(runtime.transport.a1 * state.omega(cell),
                                 state.strain_rate(cell) * blend.second);
                    const double scale = std::max({controls.kinematic_viscosity,
                        expected, state.nu_t(cell)});
                    viscosity_closure = std::max(viscosity_closure,
                        std::abs(state.nu_t(cell) - expected) / scale);
                }
            }
            h.global_nonlinear_v8 = certify_global_nonlinear_residual_v8(
                controls.algorithm, controls.coupling.schur_model,
                v8_refreshed_momentum, h.continuity_normalized,
                controls.convergence, controls.sst_runtime != nullptr,
                k_residual, omega_residual, viscosity_closure,
                !controls.sst_runtime ||
                    (sst_step.committed && !sst_step.provisional &&
                     sst_step.transport.converged),
                !controls.sst_runtime ||
                    (sst_step.transport.nonfinite_bound_events == 0 &&
                     sst_step.transport.k_floor_events == 0 &&
                     sst_step.transport.second_floor_events == 0));
            if (!h.global_nonlinear_v8.all_finite)
                throw NonlinearRetryableFailure(
                    "V8 global physical residual has nonfinite components");
            // The trial has already run full SST and refreshed U,p/phi.
            // Compare against the LAST COMMITTED history entry only, and
            // only with the same high-order convection operator. A rejection
            // propagates to the existing full U/p/phi/SST rollback.
            if (controls.joint_globalization_v13_experimental &&
                same_convection_operator_as_previous && !result.history.empty()) {
                GlobalNonlinearAcceptanceV13 joint;
                try {
                    joint = evaluate_global_nonlinear_trial_v13(
                        result.history.back().global_nonlinear_v8,
                        h.global_nonlinear_v8, controls.convergence);
                } catch (const std::invalid_argument& error) {
                    throw NonlinearRetryableFailure(
                        std::string("V13 joint state not comparable: ") + error.what());
                }
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "GLOBAL_JOINT_V13 iteration=" << iter
                              << " growth=" << joint.largest_growth
                              << " accepted=" << (joint.allowed ? 1 : 0)
                              << " failed_equation=" << joint.failed_equation
                              << " qualified=0\n";
                if (!joint.allowed)
                    throw NonlinearRetryableFailure(
                        "V13 joint post-SST physical residual growth: " +
                        joint.failed_equation);
            }
            if (controls.monotone_global_merit_v14_experimental) {
                const auto contraction=committed_global_merit_v14.evaluate(
                    h.global_nonlinear_v8,same_convection_operator_as_previous);
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "GLOBAL_CONTRACTION_TRIAL_V14 iteration=" << iter
                              << " baseline=" << contraction.baseline
                              << " old=" << contraction.old_merit
                              << " new=" << contraction.new_merit
                              << " ratio=" << contraction.ratio
                              << " allowed=" << contraction.allowed
                              << " qualified=0\n";
                if (!contraction.allowed)
                    throw NonlinearRetryableFailure(
                        "V14 independent full-state merit did not contract: retry U,p,SST");
            }
            nonlinear_metric = std::max(nonlinear_metric,
                h.global_nonlinear_v8.merit);
            if (controls.diagnostics.iteration_trace)
                std::cerr << "GLOBAL_NONLINEAR_V8 iteration=" << iter
                          << " schur=" << to_string(controls.coupling.schur_model)
                          << " method=" << static_cast<int>(controls.algorithm)
                          << " U=" << h.global_nonlinear_v8.momentum.max_active_component_relative
                          << " p_continuity=" << h.global_nonlinear_v8.continuity
                          << " k=" << h.global_nonlinear_v8.k
                          << " omega=" << h.global_nonlinear_v8.omega
                          << " nut_closure=" << h.global_nonlinear_v8.eddy_viscosity_closure
                          << " merit=" << h.global_nonlinear_v8.merit
                          << " accepted=" << h.global_nonlinear_v8.converged
                          << " qualified=0\n";
        }
        if (controls.convection_continuation.enabled &&
            outer_high_order_fraction == 1.0 && iter > 1 &&
            convection_high_order_fraction(controls.convection_continuation,
                                           iter - 1) < 1.0) {
            // Do not inherit a residual or stagnation history from the
            // incomplete UPWIND/deferred-correction operator.
            convergence_monitor = cfdx::core::ConvergenceMonitor(monitor_controls);
        }
        if (controls.require_post_sst_momentum_reassembly) {
            // Feed the cross-equation state into the nonlinear monitor too.
            // Otherwise the independent SST gate could reject convergence
            // while the reported global nonlinear merit silently tracked
            // only the old-viscosity momentum operator.
            nonlinear_metric=std::max(nonlinear_metric,
                h.post_sst_momentum_rhs_relative);
            h.nonlinear_convergence_metric=nonlinear_metric;
        }
        const auto convergence_report = convergence_monitor.update({
            iter,
            nonlinear_metric,
            h.continuity_normalized,
            std::move(qoi_values)});
        const bool sst_residual_gate = !controls.sst_runtime ||
            (sst_step.committed && !sst_step.provisional &&
             sst_step.transport.converged &&
             std::isfinite(sst_step.transport.k_residual) &&
             std::isfinite(sst_step.transport.second_residual) &&
             sst_step.transport.nonfinite_bound_events == 0 &&
             std::isfinite(sst_nu_t_relative_change) &&
             sst_nu_t_relative_change <= controls.convergence.relative_tolerance &&
             (!controls.require_post_sst_momentum_reassembly ||
              h.post_sst_momentum_rhs_relative <=
                  controls.convergence.relative_tolerance));
        const bool global_nonlinear_v8_gate =
            !controls.require_global_nonlinear_certificate_v8 ||
            h.global_nonlinear_v8.converged;
        const bool converged_now =
            iter >= minimum_outer_correctors &&
            convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED &&
            per_criterion_gates_ok && sst_residual_gate &&
            global_nonlinear_v8_gate && full_order_eligible;
        if (debug_cell_enabled && (converged_now || iter == controls.convergence.max_iterations)) {
            std::cerr << "\n=== CFDX MOMENTUM MICROSCOPE cell=" << debug_cell
                      << " cell_source=" << (debug_cell_auto ? "auto_worst" : "explicit")
                      << " algorithm=" << static_cast<int>(controls.algorithm)
                      << " iteration=" << iter << " ===\n";
            const auto dump_equation = [&](const char* name,
                                           const ScalarEquation& eq,
                                           const Vector& solution) {
                double ax = 0.0;
                const auto begin = eq.matrix.row_offsets_data()[debug_cell];
                const auto end = eq.matrix.row_offsets_data()[debug_cell + 1];
                std::cerr << name << " matrix row: rhs=" << eq.rhs(debug_cell)
                          << " diagonal=" << eq.diagonal[debug_cell]
                          << " nnz=" << (end - begin) << "\n";
                for (std::uint32_t k = begin; k < end; ++k) {
                    const auto col = eq.matrix.columns_data()[k];
                    const auto value = eq.matrix.values_data()[k];
                    ax += value * solution(col);
                    std::cerr << "  col=" << col << " a=" << value
                              << " x=" << solution(col)
                              << " ax=" << value * solution(col) << "\n";
                }
                const double residual = ax - eq.rhs(debug_cell);
                std::cerr << "  exact_row_ax=" << ax
                          << " exact_row_residual=" << residual
                          << " abs=" << std::abs(residual) << "\n";
            };
            dump_equation("Ux", final_ex, final_ux);
            dump_equation("Uy", final_ey, final_uy);
            dump_equation("Uz", final_ez, final_uz);

            // Reconstruct the same momentum predictor used by the pressure-
            // velocity coupling from the final, independently assembled
            // momentum rows. This is diagnostic only; it does not alter the
            // converged state. dAU = V/A_P is the exact pressure-response
            // coefficient used by the segregated Rhie-Chow operator.
            std::array<std::vector<double>, 3> diag_rAU;
            std::array<std::vector<double>, 3> hbyA;
            const ScalarEquation* final_eqs[3] = {&final_ex, &final_ey, &final_ez};
            const Vector* final_u[3] = {&final_ux, &final_uy, &final_uz};
            for (std::size_t d = 0; d < 3; ++d) {
                diag_rAU[d].resize(mesh.n_cells());
                hbyA[d].resize(mesh.n_cells());
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    const double ap = final_eqs[d]->diagonal[c];
                    if (!(ap > 0.0) || !std::isfinite(ap))
                        throw std::runtime_error(
                            "CFDX momentum microscope: invalid final momentum diagonal");
                    diag_rAU[d][c] = 1.0 / ap;
                    double h_explicit = final_eqs[d]->rhs(c) +
                        final_grad_p.component_data(d)[c] * geometry.cell_volumes[c];
                    const auto begin = final_eqs[d]->matrix.row_offsets_data()[c];
                    const auto end = final_eqs[d]->matrix.row_offsets_data()[c + 1];
                    for (std::uint32_t k = begin; k < end; ++k) {
                        const auto col = final_eqs[d]->matrix.columns_data()[k];
                        if (col != c)
                            h_explicit -=
                                final_eqs[d]->matrix.values_data()[k] * (*final_u[d])(col);
                    }
                    hbyA[d][c] = diag_rAU[d][c] * h_explicit;
                }
            }

            const auto reconstructed_flux =
                make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs,
                       {}, controls.pressure_projected_face_interpolation);
            double local_div_authoritative = 0.0;
            double local_div_reconstructed = 0.0;
            double max_delta_phi = 0.0;
            double max_delta_grad_p = 0.0;
            double max_delta_u_face = 0.0;
            double max_momentum_pressure_flux = 0.0;
            double max_correction_pressure_flux = 0.0;
            double max_rhie_chow_actual = 0.0;
            double max_rhie_chow_expected = 0.0;
            double max_rhie_chow_mismatch = 0.0;
            if (freeze_state && frozen_state_valid) {
                double phi_delta_linf = 0.0;
                double hbyA_delta_linf = 0.0;
                std::size_t hbyA_delta_cell = 0;
                for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                    phi_delta_linf = std::max(phi_delta_linf,
                        std::abs(mass_flux(f) - phi_used_in_momentum(f)));
                for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                    double delta = 0.0;
                    for (std::size_t d = 0; d < 3; ++d)
                        delta = std::max(delta, std::abs(hbyA[d][c] - frozen_hbyA[d][c]));
                    if (delta > hbyA_delta_linf) {
                        hbyA_delta_linf = delta;
                        hbyA_delta_cell = c;
                    }
                }
                const std::size_t comparison_cell =
                    debug_cell == 33 && mesh.n_cells() > 70 ? 70 : 33;
                auto dump_state_cell = [&](std::size_t cell, const char* label) {
                    if (cell >= mesh.n_cells()) return;
                    double phi_delta_l1 = 0.0;
                    double phi_used_l1 = 0.0;
                    double phi_reassembled_l1 = 0.0;
                    double hbyA_delta = 0.0;
                    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
                        const auto owner = mesh.ownership().owner(f);
                        const auto neighbour = mesh.ownership().neighbour(f);
                        if (owner != cell && (neighbour < 0 || static_cast<std::size_t>(neighbour) != cell)) continue;
                        const double sign = owner == cell ? 1.0 : -1.0;
                        const double used = sign * phi_used_in_momentum(f);
                        const double current = sign * mass_flux(f);
                        phi_used_l1 += std::abs(used);
                        phi_reassembled_l1 += std::abs(current);
                        phi_delta_l1 += std::abs(used-current);
                    }
                    for (std::size_t d = 0; d < 3; ++d)
                        hbyA_delta = std::max(hbyA_delta,
                            std::abs(hbyA[d][cell] - frozen_hbyA[d][cell]));
                    std::cerr << "FROZEN_STATE_CELL label=" << label
                              << " cell=" << cell
                              << " phi_used_in_momentum_assembly=" << phi_used_l1
                              << " phi_used_in_reassembly=" << phi_reassembled_l1
                              << " phi_delta_l1=" << phi_delta_l1
                              << " HbyA_delta_linf=" << hbyA_delta
                              << " frozen_rAU=(" << frozen_rAU[0][cell] << ","
                              << frozen_rAU[1][cell] << "," << frozen_rAU[2][cell] << ")"
                              << " final_rAU=(" << diag_rAU[0][cell] << ","
                              << diag_rAU[1][cell] << "," << diag_rAU[2][cell] << ")\n";
                };
                std::cerr << "FROZEN_STATE_COMPARE phi_global_linf=" << phi_delta_linf
                          << " HbyA_global_linf=" << hbyA_delta_linf
                          << " HbyA_max_cell=" << hbyA_delta_cell << "\n";
                dump_state_cell(debug_cell, "debug");
                if (comparison_cell != debug_cell)
                    dump_state_cell(comparison_cell, "comparison");
            }
            const Offset off = mesh.cells().offsets_data()[debug_cell];
            const Offset count = mesh.cells().offsets_data()[debug_cell + 1] - off;
            std::cerr << "faces=" << count << "\n";
            for (Offset k = 0; k < count; ++k) {
                const std::size_t face = mesh.cells().faces_data()[off + k];
                const auto owner = mesh.ownership().owner(face);
                const auto neighbour = mesh.ownership().neighbour(face);
                const bool owner_side = owner == debug_cell;
                const double sign = owner_side ? 1.0 : -1.0;
                const double phi_auth = sign * mass_flux(face);
                const double phi_reconstructed = sign * reconstructed_flux(face);
                local_div_authoritative += phi_auth;
                local_div_reconstructed += phi_reconstructed;
                const Vec3 Sf = geometry.face_area_vectors[face] * sign;
                std::size_t patch = geometry.face_patch[face];
                std::string patch_name = "INTERNAL";
                int patch_type = -1;
                if (neighbour < 0 && patch < mesh.boundary().n_patches()) {
                    patch_name = mesh.boundary().patch(patch).name;
                    patch_type = static_cast<int>(mesh.boundary().patch(patch).type);
                }

                Vec3 u_face_auth{0.0, 0.0, 0.0};
                Vec3 u_face_reconstructed{0.0, 0.0, 0.0};
                Vec3 grad_face_momentum{0.0, 0.0, 0.0};
                Vec3 grad_face_correction{0.0, 0.0, 0.0};
                Vec3 hbyA_face{0.0, 0.0, 0.0};
                double dAU_n = 0.0;
                double pressure_flux_correction = 0.0;
                double rhie_chow_expected = 0.0;
                double rhie_chow_actual = 0.0;

                if (neighbour >= 0) {
                    const std::size_t ncell = static_cast<std::size_t>(neighbour);
                    u_face_auth = Vec3{
                        0.5 * (final_ux(debug_cell) + final_ux(ncell)),
                        0.5 * (final_uy(debug_cell) + final_uy(ncell)),
                        0.5 * (final_uz(debug_cell) + final_uz(ncell))};
                    hbyA_face = Vec3{
                        0.5 * (hbyA[0][debug_cell] + hbyA[0][ncell]),
                        0.5 * (hbyA[1][debug_cell] + hbyA[1][ncell]),
                        0.5 * (hbyA[2][debug_cell] + hbyA[2][ncell])};
                    // Build the local face frame from the debug-cell side.
                    // For a neighbour-side debug cell, the global owner-oriented
                    // face vector and the centre-to-centre vector have opposite
                    // orientations. Mixing those frames was the source of the
                    // NaN/incorrect pressure-flux microscope records and could
                    // make a valid internal face look non-commuting.
                    const std::size_t other_cell = owner_side
                        ? ncell : owner;
                    const Vec3 dvec = geometry.cell_centres[other_cell] -
                                      geometry.cell_centres[debug_cell];
                    const double d = dvec.mag();
                    if (!(d > 0.0) || !std::isfinite(d))
                        throw std::runtime_error(
                            "CFDX momentum microscope: degenerate internal face distance");
                    const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                    const double cx = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[0][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[0][other_cell]);
                    const double cy = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[1][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[1][other_cell]);
                    const double cz = 0.5 * (
                        geometry.cell_volumes[debug_cell] * diag_rAU[2][debug_cell] +
                        geometry.cell_volumes[other_cell] * diag_rAU[2][other_cell]);
                    dAU_n = cx*e.x*e.x + cy*e.y*e.y + cz*e.z*e.z;
                    const Vec3 gp{
                        0.5 * (final_grad_p.component_data(0)[debug_cell] +
                               final_grad_p.component_data(0)[other_cell]),
                        0.5 * (final_grad_p.component_data(1)[debug_cell] +
                               final_grad_p.component_data(1)[other_cell]),
                        0.5 * (final_grad_p.component_data(2)[debug_cell] +
                               final_grad_p.component_data(2)[other_cell])};
                    // Compare the authoritative face-flux correction against
                    // the exact Rhie-Chow operator used by the pressure equation.
                    // This is deliberately evaluated from HbyA, not the final
                    // velocity, because phi = HbyA_f.Sf - RC(p) is the
                    // conservative pressure-velocity coupling operator.
                    // Use the exact inverse momentum diagonals reconstructed
                    // above from the final momentum rows.  The microscope is
                    // diagnostic only and must not introduce a separate
                    // SIMPLEC/SIMPLE coefficient path.
                    rhie_chow_expected = rhie_chow_pressure_flux_internal(
                        mesh, geometry, face, p, final_grad_p,
                        diag_rAU, controls.density);
                    const double hbyA_flux = controls.density * hbyA_face.dot(Sf);
                    rhie_chow_actual = hbyA_flux - phi_auth;
                    max_rhie_chow_actual = std::max(
                        max_rhie_chow_actual, std::abs(rhie_chow_actual));
                    max_rhie_chow_expected = std::max(
                        max_rhie_chow_expected, std::abs(rhie_chow_expected));
                    max_rhie_chow_mismatch = std::max(
                        max_rhie_chow_mismatch,
                        std::abs(rhie_chow_actual - rhie_chow_expected));
                    const Vec3 localSf = Sf;
                    const double normal_dp = (p(other_cell) - p(debug_cell)) / d;
                    const Vec3 dAU_face{cx, cy, cz};
                    grad_face_correction =
                        e * normal_dp + gp - e * gp.dot(e);
                    grad_face_momentum = gp;
                    pressure_flux_correction = controls.density * (
                        pressure_face_orthogonal_coefficient(
                            localSf, dvec, dAU_face) *
                            (p(other_cell) - p(debug_cell)) +
                        pressure_face_nonorthogonal_flux(
                            localSf, dvec, dAU_face, gp));
                    u_face_reconstructed = hbyA_face -
                        Vec3{cx * gp.x, cy * gp.y, cz * gp.z};
                } else {
                    if (patch < mesh.boundary().n_patches()) {
                        const auto& name = mesh.boundary().patch(patch).name;
                        const auto pit = pressure_bcs.find(name);
                        const auto uit = velocity_bcs.find(name);
                        if (uit != velocity_bcs.end() &&
                            uit->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                            u_face_auth = uit->second.value;
                        else
                            u_face_auth = Vec3{
                                final_ux(debug_cell), final_uy(debug_cell), final_uz(debug_cell)};
                        if (pit != pressure_bcs.end() &&
                            pit->second.type == ScalarBoundaryType::FIXED_VALUE) {
                            const Vec3 Sf_raw = geometry.face_area_vectors[face];
                            const Vec3 dvec = geometry.face_centres[face] -
                                              geometry.cell_centres[debug_cell];
                            const double d = dvec.mag();
                            if (!(d > 0.0) || !std::isfinite(d))
                                throw std::runtime_error(
                                    "CFDX momentum microscope: degenerate boundary face distance");
                            const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                            const Vec3 dAU_boundary{
                                geometry.cell_volumes[debug_cell] * diag_rAU[0][debug_cell],
                                geometry.cell_volumes[debug_cell] * diag_rAU[1][debug_cell],
                                geometry.cell_volumes[debug_cell] * diag_rAU[2][debug_cell]};
                            dAU_n = dAU_boundary.x*e.x*e.x +
                                    dAU_boundary.y*e.y*e.y +
                                    dAU_boundary.z*e.z*e.z;
                            const double normal_dp =
                                (pit->second.value - p(debug_cell)) / d;
                            const Vec3 gp{
                                final_grad_p.component_data(0)[debug_cell],
                                final_grad_p.component_data(1)[debug_cell],
                                final_grad_p.component_data(2)[debug_cell]};
                            pressure_flux_correction = controls.density * (
                                pressure_face_orthogonal_coefficient(
                                    Sf_raw, dvec, dAU_boundary) *
                                    (pit->second.value - p(debug_cell)) +
                                pressure_face_nonorthogonal_flux(
                                    Sf_raw, dvec, dAU_boundary, gp));
                            grad_face_correction = e * normal_dp + gp - e * gp.dot(e);
                        } else {
                            grad_face_correction = Vec3{
                                final_grad_p.component_data(0)[debug_cell],
                                final_grad_p.component_data(1)[debug_cell],
                                final_grad_p.component_data(2)[debug_cell]};
                        }
                        grad_face_momentum = Vec3{
                            final_grad_p.component_data(0)[debug_cell],
                            final_grad_p.component_data(1)[debug_cell],
                            final_grad_p.component_data(2)[debug_cell]};
                        hbyA_face = Vec3{
                            hbyA[0][debug_cell],
                            hbyA[1][debug_cell],
                            hbyA[2][debug_cell]};
                        const Vec3 gp = grad_face_momentum;
                        u_face_reconstructed = hbyA_face -
                            Vec3{dAU_n*gp.x, dAU_n*gp.y, dAU_n*gp.z};
                    }
                }

                const double momentum_pressure_flux = grad_face_momentum.dot(Sf);
                max_delta_phi = std::max(
                    max_delta_phi, std::abs(phi_auth - phi_reconstructed));
                max_delta_grad_p = std::max(
                    max_delta_grad_p,
                    (grad_face_momentum - grad_face_correction).mag());
                max_delta_u_face = std::max(
                    max_delta_u_face, (u_face_auth - u_face_reconstructed).mag());
                max_momentum_pressure_flux = std::max(
                    max_momentum_pressure_flux, std::abs(momentum_pressure_flux));
                max_correction_pressure_flux = std::max(
                    max_correction_pressure_flux, std::abs(pressure_flux_correction));

                std::cerr << "  face=" << face
                          << " owner=" << owner
                          << " neighbour=" << neighbour
                          << " patch=" << patch_name
                          << " patch_type=" << patch_type
                          << " Sf=(" << Sf.x << "," << Sf.y << "," << Sf.z << ")"
                          << " phi_authoritative=" << phi_auth
                          << " phi_reconstructed_from_Uface=" << phi_reconstructed
                          << " delta_phi=" << (phi_auth - phi_reconstructed)
                          << " U_face_authoritative=(" << u_face_auth.x << ","
                          << u_face_auth.y << "," << u_face_auth.z << ")"
                          << " U_face_reconstructed_HbyA_dAUgradp=("
                          << u_face_reconstructed.x << "," << u_face_reconstructed.y
                          << "," << u_face_reconstructed.z << ")"
                          << " delta_U_face=" << (u_face_auth-u_face_reconstructed).mag()
                          << " dAU=" << dAU_n
                          << " gradp_face_momentum=(" << grad_face_momentum.x << ","
                          << grad_face_momentum.y << "," << grad_face_momentum.z << ")"
                          << " gradp_face_correction=(" << grad_face_correction.x << ","
                          << grad_face_correction.y << "," << grad_face_correction.z << ")"
                          << " delta_gradp="
                          << (grad_face_momentum-grad_face_correction).mag()
                          << " momentum_pressure_flux=" << momentum_pressure_flux
                          << " correction_pressure_flux=" << pressure_flux_correction
                          << " rhie_chow_expected=" << rhie_chow_expected
                          << " rhie_chow_actual=" << rhie_chow_actual
                          << " delta_rhie_chow="
                          << (rhie_chow_actual-rhie_chow_expected)
                          << "\n";
            }

            std::cerr << "local_mass_flux_divergence=" << local_div_authoritative
                      << " local_reconstructed_flux_divergence=" << local_div_reconstructed
                      << " reconstructed_velocity_continuity_linf="
                      << h.reconstructed_velocity_continuity_linf
                      << " delta_phi_linf=" << max_delta_phi
                      << " delta_u_face_linf=" << max_delta_u_face
                      << " delta_gradp_linf=" << max_delta_grad_p
                      << " momentum_pressure_flux_linf=" << max_momentum_pressure_flux
                      << " correction_pressure_flux_linf=" << max_correction_pressure_flux
                      << " rhie_chow_expected_linf=" << max_rhie_chow_expected
                      << " rhie_chow_actual_linf=" << max_rhie_chow_actual
                      << " delta_rhie_chow_linf=" << max_rhie_chow_mismatch
                      << "\n";
            std::cerr << "diagnostic_global=(" << h.momentum_equation_residual_components[0]
                      << "," << h.momentum_equation_residual_components[1]
                      << "," << h.momentum_equation_residual_components[2] << ")"
                      << " internal=" << h.momentum_equation_residual_internal
                      << " boundary=" << h.momentum_equation_residual_boundary
                      << " worst_cell=" << h.momentum_residual_cell
                      << " worst_patch=" << h.momentum_residual_patch << "\n";
            std::cerr << "no_pressure=" << h.momentum_residual_no_pressure
                      << " pressure_contribution=" << h.momentum_pressure_contribution
                      << " gradp_linf=" << h.pressure_gradient_linf
                      << "\n";
            std::cerr << "=== END CFDX MOMENTUM MICROSCOPE ===\n";
        }
        {
            const auto conservation =
                cfdx::core::audit_face_flux_conservation(mesh, mass_flux);
            const double scale =
                std::max(conservation.global_abs_boundary_flux,
                         std::max(conservation.l1_local_imbalance, 1.0));
            h.mass_boundary_flux = conservation.global_boundary_flux;
            h.mass_global_cell_balance = conservation.global_cell_balance;
            h.mass_local_l1 = conservation.l1_local_imbalance;
            h.mass_local_linf = conservation.max_local_imbalance;
            h.mass_local_l2 = conservation.l2_local_imbalance;
            h.mass_normalized_imbalance =
                conservation.l1_local_imbalance / scale;
            h.mass_worst_cell = conservation.worst_cell;
            h.mass_nonfinite_faces = conservation.nonfinite_faces;
            if (controls.diagnostics.flux_patch_audit) {
                const auto mass_audit=audit_v17_authoritative_flux(
                    mesh,geometry,mass_flux);
                std::cerr << std::scientific << std::setprecision(13)
                          << "V17_FLUX_AUDIT iteration=" << iter
                          << " cell_l1=" << mass_audit.cell_l1
                          << " cell_linf=" << mass_audit.cell_linf
                          << " cell_worst=" << mass_audit.worst_cell
                          << " worst_x=" << geometry.cell_centres[mass_audit.worst_cell].x
                          << " worst_y=" << geometry.cell_centres[mass_audit.worst_cell].y
                          << " worst_div=" << mass_audit.cell_divergence[mass_audit.worst_cell]
                          << " boundary_signed=" << mass_audit.boundary_signed
                          << " boundary_abs=" << mass_audit.boundary_abs
                          << " cell_signed=" << mass_audit.cell_signed
                          << " topology_error=" << mass_audit.topology_error
                          << " rhie_chow_identity_linf=" << h.rhie_chow_identity_linf
                          << " qualified=0\n";
                for (const auto& [patch_name,patch] : mass_audit.patches)
                    std::cerr << "V17_FLUX_PATCH iteration=" << iter
                              << " name=" << patch_name
                              << " faces=" << patch.faces
                              << " signed=" << patch.signed_flux
                              << " absolute=" << patch.abs_flux << " qualified=0\n";
                const auto* cf=mesh.cells().faces_data();
                const auto* offs=mesh.cells().offsets_data();
                const auto cell=mass_audit.worst_cell;
                for (auto j=offs[cell];j<offs[cell+1];++j) {
                    const std::size_t f=cf[j];
                    const auto nb=mesh.ownership().neighbour(f);
                    const double outward=(mesh.ownership().owner(f)==cell
                        ? mass_flux(f) : -mass_flux(f));
                    std::string name="internal";
                    if(nb<0) {
                        const auto pi=geometry.face_patch[f];
                        name=mesh.boundary().patch(pi).name;
                    }
                    std::cerr << "V17_FLUX_FACE iteration=" << iter
                              << " cell=" << cell << " face=" << f
                              << " patch=" << name
                              << " outward=" << outward
                              << " owner_flux=" << mass_flux(f)
                              << " qualified=0\n";
                }
            }
            for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
                for (std::size_t d = 0; d < 3; ++d) {
                    if (!std::isfinite(U.component_data(d)[c]))
                        ++h.boundedness_nonfinite_velocity;
                }
            }

            // Independent post-solve momentum conservation reconstruction.
            // The flux is rebuilt from the accepted final field, independently
            // of the assembled matrix residual. Internal faces therefore have
            // one owner-oriented flux and are audited separately for balance.
            Field<double, Location::FACE> audit_transpose_flux;
            if (include_transpose_stress)
                audit_transpose_flux = transpose_viscous_stress_flux(
                    mesh, geometry, U, mu_eff, {ubc_x, ubc_y, ubc_z},
                controls.sst_runtime
                    ? controls.density * controls.kinematic_viscosity : -1.0,
                    post_wall_face_mu_v17);
            for (std::size_t d = 0; d < 3; ++d) {
                cfdx::core::Field<double, cfdx::core::Location::CELL> component(
                    mesh.n_cells(), "U_reconstructed", "m/s", 1);
                cfdx::core::Field<double, cfdx::core::Location::CELL> source(
                    mesh.n_cells(), "momentum_source_reconstructed", "N/m3", 1);
                ScalarBoundaryConditions scalar_bcs;
                for (const auto& [name, bc] : velocity_bcs) {
                    if (bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                        scalar_bcs[name] = {ScalarBoundaryType::FIXED_VALUE, bc.value.x, 0.0};
                    else
                        scalar_bcs[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
                    if (d == 1 && bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                        scalar_bcs[name].value = bc.value.y;
                    if (d == 2 && bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                        scalar_bcs[name].value = bc.value.z;
                }
                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                    component(cell) = U.component_data(d)[cell];
                    // Same channel selection as the momentum assembly
                    // (cellwise replaces the uniform vector when
                    // supplied), preserving this audit's own convention
                    // for the source term.
                    const bool cellwise_body_force =
                        !controls.cell_body_force.empty();
                    const double body_component =
                        cellwise_body_force
                            ? (d == 0
                                   ? controls.cell_body_force[cell].x
                                   : (d == 1
                                          ? controls.cell_body_force[cell].y
                                          : controls.cell_body_force[cell].z))
                            : (d == 0
                                   ? controls.body_force.x
                                   : (d == 1
                                          ? controls.body_force.y
                                          : controls.body_force.z));
                    // audit_integrated_balance sums outward flux + source.
                    // The physical equation is div(convective-viscous flux)
                    // + grad(p*) - rho*f = 0.
                    source(cell) = final_grad_p.component_data(d)[cell] -
                                   controls.density * body_component;
                }
                const auto wall_face_mu = controls.sst_runtime
                    ? resolved_wall_face_viscosity(
                        mesh, geometry, mu_eff, scalar_bcs,
                        sst_molecular_wall_mu,post_wall_face_mu_v17)
                    : std::vector<double>{};
                auto flux = reconstruct_scalar_transport_flux(
                    mesh, geometry, mass_flux, component, 0.0, scalar_bcs,
                    controls.use_bounded_convection, controls.convection_scheme,
                    final_face_velocity_conditions ? &(*final_face_velocity_conditions)[d] : nullptr,
                    &mu_eff, nullptr,
                    wall_face_mu.empty() ? nullptr : &wall_face_mu,
                    outer_high_order_fraction);
                // Transport flux is convection minus viscous stress. Rebuild
                // the transpose contribution from final U independently of
                // the predictor source and assembled matrix residual.
                if (include_transpose_stress)
                    for (std::size_t face=0; face<mesh.n_faces(); ++face)
                        flux(face) -= audit_transpose_flux(face,d);
                // Independent face balance must include the real physical
                // acceleration in transient mode. A quasi-steady or pseudo-
                // time diagonal is not part of this evidence channel.
                std::optional<cfdx::core::Field<double, cfdx::core::Location::CELL>>
                    negative_time_accumulation;
                if (controls.transient.enabled) {
                    const auto& prior = *controls.transient.previous;
                    const auto& older = controls.transient.history_valid
                        ? *controls.transient.previous_previous : prior;
                    const bool bdf2 =
                        controls.transient.scheme == cfdx::core::TimeScheme::BDF2 &&
                        controls.transient.history_valid;
                    negative_time_accumulation.emplace(
                        mesh.n_cells(), "negative_momentum_accumulation", "N", 1);
                    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                        (*negative_time_accumulation)(cell) =
                            momentum_audit_negative_accumulation(
                                component(cell), prior.component_data(d)[cell],
                                older.component_data(d)[cell], controls.density,
                                geometry.cell_volumes[cell], controls.transient.dt,
                                controls.transient.dt_previous, bdf2);
                    }
                }
                const auto* time_accumulation = negative_time_accumulation
                    ? &*negative_time_accumulation : nullptr;
                const auto balance = cfdx::core::audit_integrated_balance(
                    mesh, flux, source, geometry.cell_volumes, time_accumulation);
                const auto physical = certify_nonlinear_face_residual(
                    balance, flux, source, geometry.cell_volumes, time_accumulation);
                h.momentum_face_local_linf[d] = physical.local_linf;
                h.momentum_face_local_relative[d] = physical.local_relative;
                h.momentum_conservation_residual[d] = balance.residual;
                h.momentum_conservation_normalized[d] = balance.normalized_residual;
                h.momentum_conservation_worst_cell[d] = balance.worst_cell;
                if (controls.momentum_face_budget_callback &&
                    controls.momentum_face_budget_every > 0 &&
                    (iter == 1 || iter % controls.momentum_face_budget_every == 0)) {
                    // A/B diagnostic: exact accepted U, pressure gradient,
                    // actual Rhie-Chow mass flux and same transport operator.
                    // `collocated` is an observation, not a replacement for
                    // the conservative flux used in the momentum equation.
                    const auto collocated = make_mass_flux(
                        mesh, geometry, U, controls.density, velocity_bcs,
                        pressure_bcs,
                        controls.pressure_projected_face_interpolation,
                        final_face_velocity_conditions ? &*final_face_velocity_conditions : nullptr);
                    cfdx::core::Field<double, cfdx::core::Location::CELL>
                        pressure_component(mesh.n_cells(), "gradp_d", "Pa/m", 1);
                    cfdx::core::Field<double, cfdx::core::Location::CELL>
                        momentum_body(mesh.n_cells(), "body_d", "N/m3", 1);
                    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
                        pressure_component(cell) = final_grad_p.component_data(d)[cell];
                        const auto& force = controls.cell_body_force.empty()
                            ? controls.body_force : controls.cell_body_force[cell];
                        momentum_body(cell) = controls.density *
                            (d == 0 ? force.x : d == 1 ? force.y : force.z);
                    }
                    cfdx::core::Field<double, cfdx::core::Location::FACE>
                        transpose_component;
                    const auto* transpose_ptr =
                        static_cast<const cfdx::core::Field<double, cfdx::core::Location::FACE>*>(nullptr);
                    if (include_transpose_stress) {
                        transpose_component = cfdx::core::Field<double, cfdx::core::Location::FACE>(
                            mesh.n_faces(), "transpose_face_d", "N", 1);
                        for (std::size_t f = 0; f < mesh.n_faces(); ++f)
                            transpose_component(f) = audit_transpose_flux(f,d);
                        transpose_ptr = &transpose_component;
                    }
                    auto budget = diagnose_momentum_face_budget(
                        mesh, geometry, mass_flux, component,
                        pressure_component, momentum_body, mu_eff, scalar_bcs,
                        controls.convection_scheme, outer_high_order_fraction,
                        controls.use_bounded_convection, &collocated,
                        transpose_ptr,
                        wall_face_mu.empty() ? nullptr : &wall_face_mu,
                        time_accumulation,
                        final_face_velocity_conditions ? &(*final_face_velocity_conditions)[d] : nullptr);
                    // The independent face residual must match the freshly
                    // reassembled nonlinear Picard operator for EACH cell.
                    // Comparing maxima alone cannot detect displaced errors.
                    const ScalarEquation* final_equations[] = {
                        &final_ex, &final_ey, &final_ez};
                    reconcile_momentum_budget_with_matrix(
                        budget, *final_equations[d], component);
                    if (!momentum_matrix_face_closure_consistent(budget))
                        throw std::runtime_error(
                            "momentum face/CSR operator disagrees cellwise at cell " +
                            std::to_string(budget.matrix_closure_worst_cell));
                    // The budget and the independent face certificate have
                    // to agree: missing pressure/viscous terms must fail.
                    const double agreement = std::abs(
                        budget.worst_absolute_residual - physical.local_linf);
                    if (agreement > 1e-10 * std::max({1.0,
                            budget.worst_absolute_residual, physical.local_linf}))
                        throw std::runtime_error(
                            "momentum face budget disagrees with independent residual");
                    controls.momentum_face_budget_callback(iter, d, budget, geometry);
                }

            }
        }
        if (controls.pseudo_momentum.enabled)
            std::cerr << std::scientific << std::setprecision(9)
                      << "V17_PSEUDO_MOMENTUM iteration=" << iter
                      << " cfl=" << h.v17_pseudo_cfl
                      << " mass_min=" << h.v17_pseudo_mass_min
                      << " mass_max=" << h.v17_pseudo_mass_max
                      << " capped_cells=" << h.v17_pseudo_capped_cells
                      << " physical_momentum_rhs_relative="
                      << h.momentum_equation_residual_rhs_relative
                      << " continuity=" << h.continuity_normalized
                      << " qualified=0\n";
        const bool independent_face_residual_ok =
            !controls.require_independent_face_residual_for_convergence ||
            independent_momentum_face_residual_gate(
                h.momentum_face_local_relative,
                controls.convergence.relative_tolerance);
        result.history.push_back(h);

        // Samples are retained regardless of probe_callback. External observers
        // are notified only after the nonlinear transaction commits below, so a
        // retry cannot leak a sample from a candidate state that was rolled back.
        for (const auto& probe : controls.probes) {
            if (probe.name.empty())
                throw std::invalid_argument("solve_steady_incompressible: probe name must not be empty");
            const double value = sample_incompressible_probe(probe, mesh, geometry, U, p);
            if (!std::isfinite(value))
                throw std::runtime_error("solve_steady_incompressible: probe value is not finite");
            result.probe_samples.push_back(
                IncompressibleProbeSample{probe.name, iter, 0.0, value});
        }

        result.convergence_status = convergence_report.status;
        result.convergence_reason = convergence_report.reason;
        if ((convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED ||
             convergence_report.status == cfdx::core::ConvergenceStatus::MAX_ITERATIONS) &&
            !sst_residual_gate) {
            // The carrier monitor does not know about the additional SST
            // closure/transport gate. Do not publish a converged status while
            // the accepted momentum residual still corresponds to an earlier
            // effective-viscosity field.
            result.convergence_status =
                iter == controls.convergence.max_iterations
                    ? cfdx::core::ConvergenceStatus::MAX_ITERATIONS
                    : cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason =
                convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED
                    ? "carrier criteria satisfied but SST transport/closure gate is not satisfied"
                    : "max iterations reached; carrier and/or SST transport/closure gates unsatisfied";
        }
        if (controls.convection_continuation.enabled &&
            !full_order_eligible &&
            convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED) {
            result.convergence_status = cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason =
                "high-order convection continuation: full operator hold not completed";
        }
        if (controls.require_global_nonlinear_certificate_v8 &&
            !global_nonlinear_v8_gate &&
            convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED) {
            result.convergence_status =
                iter == controls.convergence.max_iterations
                    ? cfdx::core::ConvergenceStatus::MAX_ITERATIONS
                    : cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason =
                "V8 true U-p-k-omega residual certificate is not satisfied";
        }
        if (convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED &&
            (!per_criterion_gates_ok || iter < minimum_outer_correctors)) {
            result.convergence_status = cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason = !per_criterion_gates_ok
                ? "global nonlinear monitor met but one or more per-component physical gates remain open"
                : "minimum required outer pressure correctors not yet completed";
        }
        if (!independent_face_residual_ok &&
            (convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED ||
             convergence_report.status == cfdx::core::ConvergenceStatus::MAX_ITERATIONS)) {
            result.convergence_status =
                iter == controls.convergence.max_iterations
                    ? cfdx::core::ConvergenceStatus::MAX_ITERATIONS
                    : cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason =
                convergence_report.status == cfdx::core::ConvergenceStatus::CONVERGED
                    ? "carrier monitor satisfied, independent local face momentum residual gate not met"
                    : "maximum nonlinear iterations; independent face momentum and/or carrier gates not met";
        }
        if (converged_now && independent_face_residual_ok) {
            result.converged = true;
            result.iterations = iter;
            stop_after_iteration = true;
        }
        if (controls.convection_continuation.enabled &&
            !full_order_eligible &&
            convergence_report.status == cfdx::core::ConvergenceStatus::STAGNATED) {
            // The operator is changing; a plateau before completion of the
            // ramp is NOT stagnation of the final physical discretization.
            convergence_monitor = cfdx::core::ConvergenceMonitor(monitor_controls);
            result.convergence_status = cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason =
                "high-order convection continuation: intermediate-operator plateau";
        }
        if (convergence_report.status == cfdx::core::ConvergenceStatus::STAGNATED &&
            full_order_eligible && controls.keep_iterating_after_stagnation) {
            result.convergence_status = cfdx::core::ConvergenceStatus::CONTINUE;
            result.convergence_reason =
                "stagnation observed; continue for full physical residual gates";
        }
        if (convergence_report.status == cfdx::core::ConvergenceStatus::DIVERGED ||
            (convergence_report.status == cfdx::core::ConvergenceStatus::STAGNATED &&
             full_order_eligible && !controls.keep_iterating_after_stagnation)) {

            result.iterations = iter;
            stop_after_iteration = true;
        }
        result.iterations = iter;
            // The gate runs on a complete trial state *before* any mutation
            // of committed TVD histories, before transaction.commit(), and
            // before any irreversible output or checkpoint is published.
            if (controls.tentative_iteration_gate_callback)
                controls.tentative_iteration_gate_callback(iter,0.0,mesh,U,p);
            // A retry must never reuse a deferred TVD flux from a rejected
            // attempt: the filtered memory is only promoted on commit.
            if (controls.tvd_deferred_memory.enabled || controls.tvd_secant_implicit.enabled) {
                if (committed_tvd_deferred_valid) {
                    committed_tvd_two_steps_ago = std::move(committed_tvd_raw);
                    committed_tvd_two_steps_valid = true;
                }
                committed_tvd_raw = std::move(pending_raw_tvd_deferred);
                // Raw-memory mode stabilizes the two-cycle without creating
                // an extra recursive time scale. The recursive variant is
                // retained for A/B diagnostics, never silently selected.
                if (controls.tvd_deferred_memory.enabled &&
                    controls.tvd_deferred_memory.recursive)
                    committed_tvd_deferred = std::move(pending_effective_tvd_deferred);
                else
                    committed_tvd_deferred = committed_tvd_raw;
                committed_tvd_input_velocity = std::move(pending_tvd_input_velocity);
                committed_tvd_deferred_valid = true;
            }
            std::optional<SSTV31OuterPTC> v31_outer_proposal;
            std::optional<SSTV32CouplingBudget> v32_balance_proposal;
            if (controls.outer_ptc_v31_experimental) {
                const auto& global=h.global_nonlinear_v8;
                if (!global.all_finite || !global.sst_transport_committed ||
                    !global.bounded)
                    throw NonlinearRetryableFailure(
                        "V31 outer PTC requires a complete bounded SST certificate");
                const std::array<double,5> metrics{{
                    global.momentum.max_active_component_relative,
                    global.continuity,global.k,global.omega,
                    global.eddy_viscosity_closure}};
                v31_outer_proposal=controls.sst_runtime->outer_ptc_v31;
                // Preview before irreversible mutation; nothing is
                // published if full physical contraction fails.
                try {
                    v31_outer_proposal->committed(
                        metrics,same_convection_operator_as_previous,
                        global.converged);
                    if (controls.outer_ptc_balance_v32_experimental) {
                        v32_balance_proposal=
                            controls.sst_runtime->outer_ptc_balance_v32;
                        const double old_merit=controls.sst_runtime->outer_ptc_v31
                            .baseline_valid ? SSTV31OuterPTC::max_merit(
                                controls.sst_runtime->outer_ptc_v31.last)
                            : std::numeric_limits<double>::infinity();
                        // The first accepted outer state establishes a
                        // baseline; subsequent ones must contract. Do NOT
                        // adapt the budget on an uncomparable baseline.
                        if (std::isfinite(old_merit) &&
                            same_convection_operator_as_previous)
                            v32_balance_proposal->committed(
                                old_merit,SSTV31OuterPTC::max_merit(metrics),
                                global.converged);
                    }
                } catch (const std::invalid_argument& error) {
                    throw NonlinearRetryableFailure(error.what());
                }
            }
            // Validate/snapshot the candidate fields BEFORE irreversible
            // commit. If a field observer refuses a nonfinite candidate, the
            // U/p/phi/SST transaction remains eligible for full rollback.
            std::optional<OuterFieldPeriodReportV43> field_period_report_v44;
            if (v43_field_probe) {
                const double pressure_gauge=p(controls.pressure_reference_cell);
                const auto& state=*controls.sst_runtime->state;
                field_period_report_v44=v43_field_probe->record(
                    [&](std::size_t channel,std::size_t index) -> double {
                        if (channel<3) return U.component_data(channel)[index];
                        if (channel==3) return p(index)-pressure_gauge;
                        if (channel==4) return mass_flux(index);
                        if (channel==5) return state.k(index);
                        if (channel==6) return state.omega(index);
                        return state.nu_t(index);
                    });
            }
            transaction.commit();
            // Advance the V45 window ONLY after the complete U/p/phi/SST
            // transaction commits. Retries of one outer iteration must retain
            // the same trigger and vary only the safeguarded trial beta.
            if (controls.coupled_period_four_damping_v45_experimental > 0.0)
                coupled_phase_escape_v45.commit(phase_escape_trial_v45,iter,
                    controls.coupled_period_four_damping_steps_v45);
            if (field_period_report_v44 && field_period_report_v44->has_four_back) {
                const auto& report=*field_period_report_v44;
                std::cerr<<"OUTER_FIELD_PERIOD_V43 iteration="<<iter
                             <<" Ux_four="<<report.relative_four_back[0]
                             <<" Uy_four="<<report.relative_four_back[1]
                             <<" Uz_four="<<report.relative_four_back[2]
                             <<" p_gauge_four="<<report.relative_four_back[3]
                             <<" phi_four="<<report.relative_four_back[4]
                             <<" k_four="<<report.relative_four_back[5]
                             <<" omega_four="<<report.relative_four_back[6]
                             <<" nut_four="<<report.relative_four_back[7]
                             <<" max_four="<<report.maximum_four_back
                             <<" max_previous="<<report.maximum_previous
                             <<" worst_channel="<<report.worst_channel_four_back
                             <<" worst_index="<<report.worst_index_four_back
                             <<" same_phase_count="<<report.consecutive_four_back
                             <<" period_four_candidate="<<report.period_four_candidate
                             <<" physical_qualified=0\n";
            }
            if (v32_balance_proposal)
                controls.sst_runtime->outer_ptc_balance_v32=*v32_balance_proposal;
            if (v31_outer_proposal) {
                controls.sst_runtime->outer_ptc_v31=*v31_outer_proposal;
                if (controls.diagnostics.iteration_trace)
                    std::cerr<<"OUTER_PTC_V31 iteration="<<iter
                             <<" cfl_committed="<<v31_outer_proposal->cfl
                             <<" committed="<<v31_outer_proposal->committed_outer
                             <<" qualified=0\n";
            }
            if (pending_wall_triplet_v15 && controls.sst_runtime) {
                controls.sst_runtime->committed_wall_triplet_v15 =
                    std::move(pending_wall_triplet_v15);
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "SST_WALL_TRIPLET_COMMITTED_V15 iteration=" << iter
                              << " faces=" << controls.sst_runtime->committed_wall_triplet_v15->wall_faces
                              << " physical_operator_modified=0 qualified=0\n";
            }
            if (controls.monotone_global_merit_v14_experimental) {
                committed_global_merit_v14.commit(
                    h.global_nonlinear_v8,same_convection_operator_as_previous);
                if (controls.diagnostics.iteration_trace)
                    std::cerr << "GLOBAL_CONTRACTION_COMMITTED_V14 iteration=" << iter
                              << " committed=" << committed_global_merit_v14.committed_steps()
                              << " contractions=" << committed_global_merit_v14.contraction_steps()
                              << " merit=" << h.global_nonlinear_v8.merit
                              << " converged=" << h.global_nonlinear_v8.converged
                              << " qualified=0\n";
            }
            // Explicit COMMITTED evidence: GLOBAL_NONLINEAR_V8 is printed
            // while the candidate is still tentative and may subsequently
            // roll back. Only this marker proves the complete U,p,phi,SST
            // transaction has been accepted. No physical qualification is
            // implied by a successful nonlinear update.
            if (controls.require_global_nonlinear_certificate_v8 &&
                controls.diagnostics.iteration_trace) {
                const auto& global = h.global_nonlinear_v8;
                std::cerr << "GLOBAL_NONLINEAR_COMMITTED_V13 iteration=" << iter
                          << " U=" << global.momentum.max_active_component_relative
                          << " p_continuity=" << global.continuity
                          << " k=" << global.k
                          << " omega=" << global.omega
                          << " nut_closure=" << global.eddy_viscosity_closure
                          << " merit=" << global.merit
                          << " converged=" << (global.converged ? 1 : 0)
                          << " qualified=0\n";
            }
            if (controls.nonlinear_multisecant.method !=
                NonlinearMultisecantMethod::OFF) {
                const auto output=pack_nonlinear_state(U,p,mass_flux);
                std::vector<double> physical;
                physical.reserve(4*nc_multisecant);
                const std::array<const ScalarEquation*,3> equations{
                    &final_ex,&final_ey,&final_ez};
                const std::array<const Vector*,3> solutions{
                    &final_ux,&final_uy,&final_uz};
                if (!multisecant_scales_valid) {
                    multisecant_rhs_reference=h.momentum_rhs_scale_components;
                    multisecant_flux_reference=
                        std::max(h.conservative_face_flux_l1,1e-30);
                    multisecant_scales_valid=true;
                }
                for (std::size_t d=0;d<3;++d) {
                    const auto& matrix=equations[d]->matrix;
                    for (std::size_t cell=0;cell<nc_multisecant;++cell) {
                        double signed_res=-equations[d]->rhs(cell);
                        for (auto k=matrix.row_offsets_data()[cell];
                             k<matrix.row_offsets_data()[cell+1];++k)
                            signed_res+=matrix.values_data()[k]*
                                (*solutions[d])(matrix.columns_data()[k]);
                        physical.push_back(signed_res /
                            std::max(multisecant_rhs_reference[d],1e-30));
                    }
                }
                std::vector<double> div(nc_multisecant,0.0);
                for (std::size_t f=0;f<nf_multisecant;++f) {
                    const auto o=mesh.ownership().owner(f);
                    div[o]+=mass_flux(f);
                    const auto n=mesh.ownership().neighbour(f);
                    if (n>=0) div[static_cast<std::size_t>(n)]-=mass_flux(f);
                }
                for (const double r:div)
                    physical.push_back(r/multisecant_flux_reference);
                double u_scale=1.0, p_scale=1e-6;
                for (std::size_t j=0;j<3*nc_multisecant;++j)
                    u_scale=std::max(u_scale,std::abs(output[j]));
                for (std::size_t j=3*nc_multisecant;j<4*nc_multisecant;++j)
                    p_scale=std::max(p_scale,std::abs(output[j]));
                const double rho_scale=controls.density*u_scale*u_scale;
                p_scale=std::max(p_scale,rho_scale);
                double phi_scale=std::max(1e-12,
                    h.conservative_face_flux_l1 / static_cast<double>(std::max<std::size_t>(1,nf_multisecant)));
                std::vector<double> scales(output.size(),u_scale);
                for (std::size_t j=3*nc_multisecant;j<4*nc_multisecant;++j)
                    scales[j]=p_scale;
                for (std::size_t j=4*nc_multisecant;j<output.size();++j)
                    scales[j]=phi_scale;
                multisecant_history.append(
                    MultisecantSample{multisecant_input,output,std::move(physical)},
                    scales);
                if (controls.diagnostics.iteration_trace &&
                    (multisecant_trial_used || iter%25==0))
                    std::cerr << "MULTISECANT_PREDICTOR iteration=" << iter
                              << " method=" << (controls.nonlinear_multisecant.method==
                                  NonlinearMultisecantMethod::ANDERSON?"ANDERSON":"NGMRES")
                              << " applied=" << multisecant_trial_used
                              << " history=" << multisecant_history.size()
                              << " predictor_step=" << multisecant_proposal.relative_step
                              << " physical_momentum=" << h.momentum_equation_residual_rhs_active_component_max
                              << " continuity=" << h.continuity_normalized << '\n';
            }
            if (controls.sst_runtime) {
                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                    mu_eff[cell] = controls.density *
                        (controls.kinematic_viscosity +
                         controls.sst_runtime->state->nu_t(cell));
            }
            // Only publish a face flux after the nonlinear transaction has
            // accepted the corresponding U/p state. mass_flux is mutable during
            // a tentative iteration and is not part of NonlinearStateRollback;
            // publishing it after a rejected attempt would therefore expose a
            // flux belonging to a state that was rolled back.
            result.authoritative_mass_flux = mass_flux;
            if (controls.pressure_high_low_production)
                result.authoritative_pressure_low = pressure_low_state;
            if (controls.sst_runtime &&
                controls.sst_runtime->transport.sst_coupled_ptc_v30_experimental &&
                sst_step.committed && !sst_step.provisional &&
                sst_step.transport.converged)
                controls.sst_runtime->coupled_ptc_v30.outer_commit();
            iteration_completed = true;
            } catch (const NonlinearRetryableFailure& error) {
                transaction.reject();
                mass_flux = mass_flux_snapshot;
                pressure_low_state=pressure_low_snapshot;
                if (controls.sst_runtime && sst_attempt_snapshot)
                    restore_sst_production_state(*controls.sst_runtime->state, *sst_attempt_snapshot);
                if (controls.sst_runtime && v30_ptc_attempt_snapshot)
                    controls.sst_runtime->coupled_ptc_v30 = *v30_ptc_attempt_snapshot;
                if (controls.sst_runtime) {
                    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
                        mu_eff[cell] = controls.density *
                            (controls.kinematic_viscosity +
                             controls.sst_runtime->state->nu_t(cell));
                }
                // The monitor and exported iteration evidence are tentative
                // until the carrier/SST transaction commits. Replaying a retry
                // with the same iteration index must not see stale samples or
                // a previously published converged flag.
                convergence_monitor = convergence_monitor_snapshot;
                result.history.resize(result_history_size);
                result.probe_samples.resize(probe_samples_size);
                result.converged = previous_result_converged;
                result.convergence_status = previous_result_status;
                result.convergence_reason = previous_result_reason;
                result.iterations = previous_result_iterations;
                stop_after_iteration = false;
                frozen_state_valid = false;
                // V37: retrying identical A,b,PC with smaller update damping
                // cannot fix the failed Krylov solve. The full rollback above
                // is still mandatory. SIMPLEC can change its Schur PC when
                // alpha_u changes; experimental outer PTC can change A.
                const CoupledLinearRetryV37 linear_retry_policy{
                    controls.coupling.schur_model == CoupledSchurModel::SIMPLEC ||
                        controls.outer_ptc_v31_experimental,
                    coupled_linear_failed_v37,
                    multisecant_trial_used && !picard_fallback_after_multisecant};
                if (linear_retry_policy.fail_fast()) {
                    std::cerr << "COUPLED_LINEAR_UNCHANGED_RETRY_V37 iteration=" << iter
                              << " schur=" << to_string(controls.coupling.schur_model)
                              << " extra_linear_solves_skipped=" <<
                                  (controls.nonlinear_retry.max_retries - retry_controller.retries())
                              << " committed=0 qualified=0\n";
                    throw;
                }
                if (!retry_controller.can_retry())
                    throw;
                if (multisecant_trial_used &&
                    !picard_fallback_after_multisecant) {
                    // Do NOT shrink equation relaxation on an Anderson/NGMRES
                    // failure. The predictor changed the input state, not the
                    // convergence tolerances or the fixed-point map.
                    picard_fallback_after_multisecant = true;
                    multisecant_history.reset();
                    retry_controller.retry_without_relaxation();
                    std::cerr << "MULTISECANT_FALLBACK iteration=" << iter
                              << " reason=" << error.what()
                              << " alpha_u=" << retry_controller.alpha_u(effective_alpha_u)
                              << " alpha_p=" << retry_controller.alpha_p(effective_alpha_p)
                              << " history_reset=1\n";
                    continue;
                }
                retry_controller.reject();
                std::cerr << "CFDX nonlinear iteration " << iter
                          << " rejected; retry=" << retry_controller.retries()
                          << " alpha_u=" << retry_controller.alpha_u(effective_alpha_u)
                          << " alpha_p=" << retry_controller.alpha_p(effective_alpha_p)
                          << " reason=" << error.what() << "\n";
            } catch (const std::runtime_error& error) {
                // A post-commit diagnostics/I/O failure is NOT a nonlinear
                // trial rejection. Never roll back an already committed state.
                if (!transaction.active()) throw;
                transaction.reject();
                mass_flux = mass_flux_snapshot;
                pressure_low_state=pressure_low_snapshot;
                if (controls.sst_runtime && sst_attempt_snapshot)
                    restore_sst_production_state(*controls.sst_runtime->state, *sst_attempt_snapshot);
                if (controls.sst_runtime && v30_ptc_attempt_snapshot)
                    controls.sst_runtime->coupled_ptc_v30 = *v30_ptc_attempt_snapshot;
                if (controls.sst_runtime)
                    for (std::size_t cell=0;cell<mesh.n_cells();++cell)
                        mu_eff[cell]=controls.density*(
                            controls.kinematic_viscosity+
                            controls.sst_runtime->state->nu_t(cell));
                convergence_monitor = convergence_monitor_snapshot;
                result.history.resize(result_history_size);
                result.probe_samples.resize(probe_samples_size);
                result.converged = previous_result_converged;
                result.convergence_status = previous_result_status;
                result.convergence_reason = previous_result_reason;
                result.iterations = previous_result_iterations;
                frozen_state_valid = false;
                if (multisecant_trial_used &&
                    !picard_fallback_after_multisecant &&
                    retry_controller.can_retry()) {
                    picard_fallback_after_multisecant = true;
                    multisecant_history.reset();
                    retry_controller.retry_without_relaxation();
                    std::cerr << "MULTISECANT_FALLBACK iteration=" << iter
                              << " before_physical_gate=1"
                              << " reason=" << error.what()
                              << " alpha_u=" << retry_controller.alpha_u(effective_alpha_u)
                              << " alpha_p=" << retry_controller.alpha_p(effective_alpha_p)
                              << " history_reset=1\n";
                    continue;
                }
                throw;
            } catch (...) {
                if (transaction.active()) transaction.reject();
                mass_flux = mass_flux_snapshot;
                pressure_low_state=pressure_low_snapshot;
                if (controls.sst_runtime && sst_attempt_snapshot)
                    restore_sst_production_state(*controls.sst_runtime->state, *sst_attempt_snapshot);
                if (controls.sst_runtime && v30_ptc_attempt_snapshot)
                    controls.sst_runtime->coupled_ptc_v30 = *v30_ptc_attempt_snapshot;
                throw;
            }
            // Restart observers must never serialize U,p before commit, nor
            // reconstruct phi from U on a collocated grid. Keep I/O outside the
            // rollback try/catch: a failed write is not a nonlinear retry.
            if (iteration_completed && controls.accepted_iteration_callback)
                controls.accepted_iteration_callback(
                    iter + controls.initial_accepted_iteration,
                    mesh, U, p, result.authoritative_mass_flux,
                    result.reference_momentum_residual,
                    controls.initial_face_mass_flux_from_paired_checkpoint
                        ? controls.initial_nonlinear_metric_reference
                        : result.history.front().nonlinear_convergence_metric);
            if (iteration_completed && controls.pressure_high_low_production &&
                controls.accepted_compensated_iteration_callback)
                controls.accepted_compensated_iteration_callback(
                    iter + controls.initial_accepted_iteration,
                    mesh,U,p,result.authoritative_mass_flux,
                    result.authoritative_pressure_low);
            if (iteration_completed && controls.iteration_output_callback &&
                !controls.iteration_output_callback(iter, 0.0, mesh, U, p)) {
                result.iterations = iter;
                stop_after_iteration = true;
            }
            if (iteration_completed && controls.momentum_linear_relaxation_callback) {
                const double used_tolerance = result.history.empty()
                    ? pending_linear_tolerance : result.history.back().linear_tolerance_used;
                for (std::size_t d=0;d<3;++d)
                    if (pending_linear_relaxation[d])
                        controls.momentum_linear_relaxation_callback(
                            iter,d,*pending_linear_relaxation[d],
                            pending_alpha_u,used_tolerance);
            }
            if (iteration_completed && pending_picard_tvd_trace &&
                controls.momentum_picard_tvd_face_callback)
                controls.momentum_picard_tvd_face_callback(
                    iter, *pending_picard_tvd_trace);
            if (iteration_completed && pending_same_state_trace &&
                controls.momentum_same_state_face_callback)
                controls.momentum_same_state_face_callback(
                    iter, *pending_same_state_trace);
            if (iteration_completed && cell_tracker) {
                // These rows and authoritative face fluxes were reassembled for
                // the same physical state whose nonlinear transaction committed.
                cell_tracker->observe(
                    iter + controls.initial_accepted_iteration,
                    cell_trace_momentum, cell_trace_continuity,
                    pending_cell_rhs_scales, cell_trace_centers,
                    same_convection_operator_as_previous,
                    stop_after_iteration || iter == controls.convergence.max_iterations);
                last_committed_cell_iteration = iter + controls.initial_accepted_iteration;
            }
            if (iteration_completed && controls.probe_callback) {
                for (std::size_t sample_index = probe_samples_size;
                     sample_index < result.probe_samples.size(); ++sample_index)
                    controls.probe_callback(result.probe_samples[sample_index]);
            }
        }
        if (stop_after_iteration)
            break;
    }

    if (cell_tracker)
        cell_tracker->finish(last_committed_cell_iteration);

    // authoritative_mass_flux is updated only after an accepted nonlinear
    // transaction above. In particular, a rejected retry cannot overwrite the
    // evidence for the last accepted U/p state. Do not reconstruct it from U.
    if (pressure_context)
        result.pressure_linear_context = pressure_context->stats();
    return result;
}

} // namespace cfdx::physics

