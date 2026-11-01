#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/linear_solver_context.h"
#include "cfdx/core/linalg/mgr_preconditioner.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/linear_solver_dispatch.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/physics/advanced_convergence.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/momentum_viscous_stress.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/pressure_velocity_system.h"
#include "cfdx/physics/pressure_continuity_line_search.h"
#include "cfdx/physics/pressure_affine_flux_trials.h"
#include "cfdx/physics/solver_control.h"
#include "cfdx/physics/sst_production_coupling.h"
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
};

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
};

struct IncompressibleSolverControls {
    InitializationControls initialization;
    // Opt-in only. The production RANS selector remains fail-closed; supplying
    // this runtime explicitly opts the caller into the SST-coupled path.
    SSTProductionRuntime* sst_runtime = nullptr;
    PressureVelocityAlgorithm algorithm = PressureVelocityAlgorithm::SIMPLE;
    IncompressibleTransientControls transient;
    CouplingControls coupling;
    ConvergenceCriteria convergence;
    std::size_t linear_max_iterations = 1000;
    double linear_tolerance = 1e-10;
    ConvergenceAccelerationControls acceleration;
    double density = 1.0;
    double kinematic_viscosity = 1.0e-3;
    double turbulent_viscosity = 0.0;
    // Optional cellwise kinematic eddy viscosity (m^2/s). Replaces the scalar
    // value when supplied; the production RANS guard remains active.
    std::vector<double> cell_turbulent_viscosity;
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
    cfdx::core::LinearSolverRequest coupled_linear_solver;
    // PCD pressure-side solves are nested inside the coupled preconditioner.
    // They deliberately have independent solver controls so the outer coupled
    // iteration budget can never become an accidental inner Krylov budget.
    cfdx::core::LinearSolverRequest pcd_laplacian_linear_solver{
        cfdx::core::KrylovModel::GMRES,
        cfdx::core::PreconditionerModel::NativeAMG,
        40,
        false,
        cfdx::core::NullSpaceModel::None};
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
        cfdx::core::NullSpaceModel::None};
    std::size_t pcd_laplacian_max_iterations = 100;
    std::size_t pcd_convection_diffusion_max_iterations = 100;
    double pcd_laplacian_tolerance = 1e-8;
    double pcd_convection_diffusion_tolerance = 1e-8;
    bool use_bounded_convection = true;
    ConvectionScheme convection_scheme = ConvectionScheme::UPWIND;
    DiagnosticsControls diagnostics;
    NonlinearRetryControls nonlinear_retry;
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
    // Evidence-only hook: invoked at the exact Krylov call point of the
    // coupled solve - after assembly, gauge treatment, finite checks and
    // preconditioner setup - so the snapshot matrix, RHS and initial guess
    // are the very objects handed to the solver. Null by default: no
    // capture, no overhead, no behavior change. It exists so an isolated
    // benchmark can replay the REAL solved system instead of a divergent
    // replica - see tests/validation/test_pcd_schur_action_benchmark.cpp.
    std::function<void(const CoupledAssemblySnapshot&)> coupled_assembly_callback;
};

struct IncompressibleIteration {
    std::size_t iteration = 0;
    double momentum_residual = std::numeric_limits<double>::infinity();
    double pressure_residual = std::numeric_limits<double>::infinity();
    double continuity_l1 = std::numeric_limits<double>::infinity();
    double continuity_linf = std::numeric_limits<double>::infinity();
    double continuity_normalized = std::numeric_limits<double>::infinity();
    // Explicit normalization evidence: raw residual / max(abs(momentum RHS))
    // and raw residual / first-iteration raw residual are different metrics.
    double momentum_rhs_scale = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_reference = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_rhs_relative = std::numeric_limits<double>::infinity();
    std::array<double, 3> momentum_rhs_scale_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    std::array<double, 3> momentum_equation_residual_rhs_relative_components{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity()};
    double momentum_equation_residual = std::numeric_limits<double>::infinity();
    double momentum_equation_residual_relative = std::numeric_limits<double>::infinity();
    double velocity_change_inf = std::numeric_limits<double>::infinity();
    double pressure_change_inf = std::numeric_limits<double>::infinity();
    double effective_alpha_u = 0.0;
    double effective_alpha_p = 0.0;
    double nonlinear_convergence_metric = std::numeric_limits<double>::infinity();
    std::size_t momentum_linear_iterations = 0;
    std::size_t pressure_linear_iterations = 0;
    std::size_t pressure_correctors_used = 0;
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
    validate_coupling_controls(c.coupling);
    validate_convergence_criteria(c.convergence);
    validate_convergence_acceleration_controls(c.acceleration);
    validate_adaptive_relaxation_controls(c.adaptive_relaxation);
    validate_nonlinear_retry_controls(c.nonlinear_retry);
    validate_pressure_continuity_line_search(c.pressure_continuity_line_search);
    // The monolithic U-p block still assembles its analytical G operator
    // with the historical face weights. Until that block is augmented with
    // the skew-correction Jacobian, activating this option would make the
    // coupled matrix inconsistent with the reconstructed Rhie-Chow flux.
    if (c.pressure_skew_corrected_gradient &&
        c.algorithm == PressureVelocityAlgorithm::COUPLED)
        throw std::invalid_argument(
            "skew-corrected pressure gradient is not qualified for COUPLED");
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
    if (c.linear_max_iterations == 0 || !(c.linear_tolerance > 0.0))
        throw std::invalid_argument("invalid linear solver controls");
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

inline cfdx::core::Field<double, cfdx::core::Location::CELL>
gauss_gradient_with_boundary(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& bcs,
    bool skew_correct = false)
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
                const double dc = (geometry.face_centres[f] -
                                   geometry.cell_centres[c]).mag();
                const double dn = (geometry.cell_centres[other] -
                                   geometry.face_centres[f]).mag();
                if (!(dc + dn > 0.0))
                    throw std::runtime_error(
                        "gauss_gradient_with_boundary: degenerate face distance");
                const double w = dn / (dc + dn);
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
                    else if (skew_correct && it != bcs.end() &&
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

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {})
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
            Uf=(Uf+Un)*0.5;
        } else {
            const std::size_t p=geometry.face_patch[f];
            if (p < mesh.boundary().n_patches() &&
                mesh.boundary().patch(p).type == PatchType::EMPTY) {
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
    double rho)
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
    const double cx = 0.5*(geometry.cell_volumes[o] * rAU[0][o] +
                           geometry.cell_volumes[n] * rAU[0][n]);
    const double cy = 0.5*(geometry.cell_volumes[o] * rAU[1][o] +
                           geometry.cell_volumes[n] * rAU[1][n]);
    const double cz = 0.5*(geometry.cell_volumes[o] * rAU[2][o] +
                           geometry.cell_volumes[n] * rAU[2][n]);
    const Vec3 gp{
        0.5*(grad_p.component_data(0)[o] + grad_p.component_data(0)[n]),
        0.5*(grad_p.component_data(1)[o] + grad_p.component_data(1)[n]),
        0.5*(grad_p.component_data(2)[o] + grad_p.component_data(2)[n])};
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
    bool skew_correct_gradient = false)
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
    bool skew_correct_gradient)
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
        p, mesh, geometry, pressure_bcs, skew_correct_gradient);
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
            Hf = (Hf + Hn) * 0.5;
            correction = rhie_chow_pressure_flux_internal(
                mesh, geometry, f, p, grad_p, rAU, rho);
            if (coupled_physical_velocity) {
                // In the monolithic block U is the physical velocity, not H/A.
                // Add the interpolated dAU*grad(p) compensation so an affine
                // pressure field is not subtracted twice from the face flux.
                const double cx = 0.5 * (
                    geometry.cell_volumes[o] * rAU[0][o] +
                    geometry.cell_volumes[n] * rAU[0][n]);
                const double cy = 0.5 * (
                    geometry.cell_volumes[o] * rAU[1][o] +
                    geometry.cell_volumes[n] * rAU[1][n]);
                const double cz = 0.5 * (
                    geometry.cell_volumes[o] * rAU[2][o] +
                    geometry.cell_volumes[n] * rAU[2][n]);
                const Vec3 gp{
                    0.5 * (grad_p.component_data(0)[o] + grad_p.component_data(0)[n]),
                    0.5 * (grad_p.component_data(1)[o] + grad_p.component_data(1)[n]),
                    0.5 * (grad_p.component_data(2)[o] + grad_p.component_data(2)[n])};
                const Vec3 Sf = geometry.face_area_vectors[f];
                correction -= rho * (cx * gp.x * Sf.x +
                                     cy * gp.y * Sf.y +
                                     cz * gp.z * Sf.z);
            }
        } else {
            const std::size_t patch = geometry.face_patch[f];
            if (patch < mesh.boundary().n_patches() &&
                mesh.boundary().patch(patch).type == PatchType::EMPTY) {
                flux(f) = 0.0;
                continue;
            }
            if (patch < mesh.boundary().n_patches()) {
                const auto& name = mesh.boundary().patch(patch).name;
                const auto ubc = bcs.find(name);
                if (ubc != bcs.end() &&
                    ubc->second.type == VelocityBoundaryCondition::Type::FIXED_VALUE)
                    Hf = ubc->second.value;

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
    const cfdx::core::Field<double, cfdx::core::Location::CELL>* convected_field = nullptr)
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
        convection_scheme, convected_field);
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
