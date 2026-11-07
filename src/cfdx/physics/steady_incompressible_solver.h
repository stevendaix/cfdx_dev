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
#include "cfdx/core/linalg/vector.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/io/restart/dat_restart.h"
#include "cfdx/physics/initial_face_flux_contract.h"
#include "cfdx/physics/advanced_convergence.h"
#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/momentum_viscous_stress.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/pressure_velocity_system.h"
#include "cfdx/physics/pressure_continuity_line_search.h"
#include "cfdx/physics/convection_continuation.h"
#include "cfdx/physics/pressure_affine_flux_trials.h"
#include "cfdx/physics/solver_control.h"
#include "cfdx/physics/sst_production_coupling.h"
#include "cfdx/physics/momentum_residual_certification.h"
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
};

inline void validate_v17_pseudo_momentum_controls(const V17PseudoMomentumControls& c) {
    if (!std::isfinite(c.initial_cfl) || !std::isfinite(c.minimum_cfl) ||
        !std::isfinite(c.maximum_cfl) || !std::isfinite(c.growth) ||
        !std::isfinite(c.shrink) || !(c.minimum_cfl > 0.0) ||
        !(c.minimum_cfl <= c.initial_cfl) ||
        !(c.initial_cfl <= c.maximum_cfl) || !(c.growth >= 1.0) ||
        !(c.shrink > 0.0 && c.shrink <= 1.0))
        throw std::invalid_argument("V17 pseudo-momentum controls: invalid CFL schedule");
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
};

struct IncompressibleSolverControls {
    InitializationControls initialization;
    // Opt-in only. The production RANS selector remains fail-closed; supplying
    // this runtime explicitly opts the caller into the SST-coupled path.
    SSTProductionRuntime* sst_runtime = nullptr;
    PressureVelocityAlgorithm algorithm = PressureVelocityAlgorithm::SIMPLE;
    IncompressibleTransientControls transient;
    V17PseudoMomentumControls pseudo_momentum;
    CouplingControls coupling;
    ConvergenceCriteria convergence;
    std::size_t linear_max_iterations = 1000;
    double linear_tolerance = 1e-10;
    // Optional independent pressure linear budget; zero inherits the global one.
    // The nonlinear mass and momentum conservation gates are unaffected.
    double pressure_linear_tolerance = 0.0;
    std::size_t pressure_linear_max_iterations = 0;
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
    if (c.coupling.coupled_schwarz_experimental &&
        (c.algorithm != PressureVelocityAlgorithm::COUPLED ||
         c.coupling.schur_model != CoupledSchurModel::BlockLocal))
        throw std::invalid_argument(
            "experimental Schwarz requires COUPLED with BlockLocal Schur");
    validate_convergence_criteria(c.convergence);
    validate_convergence_acceleration_controls(c.acceleration);
    validate_adaptive_relaxation_controls(c.adaptive_relaxation);
    validate_nonlinear_retry_controls(c.nonlinear_retry);
    validate_convection_continuation(c.convection_continuation,
                                     c.convergence.max_iterations);
    if (c.convection_continuation.enabled &&
        ((c.convection_scheme != ConvectionScheme::TVD &&
          c.convection_scheme != ConvectionScheme::TVD_VANLEER &&
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

inline cfdx::core::Field<double, cfdx::core::Location::FACE>
make_mass_flux(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
    double rho,
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {},
    bool projected_face_interpolation = false)
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
