#pragma once

#include <array>
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <string>

namespace cfdx::core {

enum class LinearProblemKind {
    General,
    PressurePoisson,
    Diffusion,
    Momentum,
    ScalarTransport,
    CoupledPressureVelocity
};

enum class KrylovModel {
    Auto,
    CG,
    BiCGStab,
    GMRES,
    FGMRES,
    LGMRES,
    PipeCG,
    PipeFGMRES
};

enum class PreconditionerModel {
    Auto,
    None,
    Jacobi,
    GaussSeidel,
    ILU0,
    ILUT,
    NativeAMG,
    SmoothedAggregationAMG,
    AnisotropicLineAMG,
    FSAI,
    RAS,
    NativeFieldSplit,
    CoupledBlockSchur,
    PCD,
    LSC,
    MGR,
    AnisotropicLine,
    AnisotropicBoomerAMG,
    AnisotropicSmoothedAggregationAMG,
    DirectionalSemiCoarseningAMG,
    EnergyMinimizingAMG,
    VankaPressurePatch,
    BraessSarazin
};

enum class ModelAvailability { Available, Planned };

enum class NullSpaceModel { None, Constant };

enum class LinearScalingModel { None, Row, Column, RowColumn, SymmetricDiagonal };

inline const char* to_string(LinearScalingModel value) {
    switch (value) {
        case LinearScalingModel::None: return "none";
        case LinearScalingModel::Row: return "row";
        case LinearScalingModel::Column: return "column";
        case LinearScalingModel::RowColumn: return "row_column";
        case LinearScalingModel::SymmetricDiagonal: return "symmetric_diagonal";
    }
    return "unknown";
}

struct SolverModelDescriptor {
    const char* name;
    ModelAvailability availability;
    bool supports_gpu;
    bool supports_mpi;
};

inline const std::array<SolverModelDescriptor, 8>& krylov_model_catalog() {
    static const std::array<SolverModelDescriptor, 8> models{{
        {"auto", ModelAvailability::Available, false, true},
        {"cg", ModelAvailability::Available, false, true},
        {"bicgstab", ModelAvailability::Available, false, true},
        {"gmres", ModelAvailability::Available, false, true},
        {"fgmres", ModelAvailability::Available, false, true},
        {"lgmres", ModelAvailability::Planned, false, false},
        {"pipecg", ModelAvailability::Planned, false, true},
        {"pipefgmres", ModelAvailability::Planned, false, true}
    }};
    return models;
}

inline const std::array<SolverModelDescriptor, 23>& preconditioner_model_catalog() {
    static const std::array<SolverModelDescriptor, 23> models{{
        {"auto", ModelAvailability::Available, false, true},
        {"none", ModelAvailability::Available, false, true},
        {"jacobi", ModelAvailability::Available, false, true},
        {"gauss_seidel", ModelAvailability::Available, false, false},
        {"ilu0", ModelAvailability::Available, false, false},
        {"ilut", ModelAvailability::Planned, false, false},
        {"native_amg", ModelAvailability::Available, false, false},
        {"smoothed_aggregation_amg", ModelAvailability::Available, false, false},
        {"anisotropic_line_amg", ModelAvailability::Available, false, false},
        {"fsai", ModelAvailability::Planned, false, false},
        {"ras", ModelAvailability::Planned, false, true},
        {"native_fieldsplit", ModelAvailability::Available, false, false},
        {"coupled_block_schur", ModelAvailability::Available, false, false},
        // PCD is requestable so the production coupled branch it is wired into
        // is reachable and can be exercised. This flag is the dispatch gate read
        // by select_linear_solver; it is not a qualification claim. PCD maturity
        // is tracked separately in numerical_method_registry.h
        // (preconditioner.pcd = Planned) and NUMERICAL_METHOD_CAPABILITY_MATRIX.json
        // (status: partial), and neither is advanced by this entry.
        {"pcd", ModelAvailability::Available, false, false},
        {"lsc", ModelAvailability::Planned, false, false},
        {"mgr", ModelAvailability::Available, false, false},
        {"anisotropic_line", ModelAvailability::Available, false, false},
        {"anisotropic_boomer_amg", ModelAvailability::Available, false, false},
        {"anisotropic_smoothed_aggregation_amg", ModelAvailability::Available, false, false},
        {"directional_semi_coarsening_amg", ModelAvailability::Available, false, false},
        {"energy_minimizing_amg", ModelAvailability::Available, false, false},
        {"vanka_pressure_patch", ModelAvailability::Available, false, false},
        {"braess_sarazin", ModelAvailability::Available, false, false}
    }};
    return models;
}

inline const std::array<SolverModelDescriptor, 2>& null_space_model_catalog() {
    static const std::array<SolverModelDescriptor, 2> models{{
        {"none", ModelAvailability::Available, false, false},
        {"constant", ModelAvailability::Available, false, false}
    }};
    return models;
}

inline const char* to_string(LinearProblemKind value) {
    switch (value) {
        case LinearProblemKind::General: return "general";
        case LinearProblemKind::PressurePoisson: return "pressure_poisson";
        case LinearProblemKind::Diffusion: return "diffusion";
        case LinearProblemKind::Momentum: return "momentum";
        case LinearProblemKind::ScalarTransport: return "scalar_transport";
        case LinearProblemKind::CoupledPressureVelocity: return "coupled_pressure_velocity";
    }
    return "unknown";
}

inline const char* to_string(KrylovModel value) {
    return krylov_model_catalog()[static_cast<std::size_t>(value)].name;
}

inline const char* to_string(PreconditionerModel value) {
    return preconditioner_model_catalog()[static_cast<std::size_t>(value)].name;
}

inline const char* to_string(NullSpaceModel value) {
    return null_space_model_catalog()[static_cast<std::size_t>(value)].name;
}

struct AnisotropyGeometry;

struct LinearSolverRequest {
    KrylovModel krylov = KrylovModel::Auto;
    PreconditionerModel preconditioner = PreconditionerModel::Auto;
    int gmres_restart = 40;
    bool allow_fallback = false;
    NullSpaceModel null_space = NullSpaceModel::None;
    // Optional algebraic scaling. None preserves the historical production path.
    LinearScalingModel scaling = LinearScalingModel::None;
    // V26 opt-in: correction solves on the ORIGINAL matrix after a stalled
    // Krylov solve. The final status still requires the original FP64
    // true-relative-residual threshold. Disabled by default.
    std::size_t residual_refinement_steps = 0;
    // Opt-in observability of the EXACT original assembled CSR matrix, before
    // scaling and preconditioning. Produces JSON in LinearSolveReport.
    bool diagnose_anisotropy = false;
    std::string anisotropy_equation{};
    // Explicit opt-in: write exact, unscaled CSR/RHS/initial guess for replay.
    // Caller owns naming/rotation; repeated use of a path overwrites it.
    std::string linear_capture_path{};
    // Borrowed for the duration of solve_linear_system(); optional, no
    // fabricated mesh metrics when absent.
    const AnisotropyGeometry* anisotropy_geometry = nullptr;
    // Only used with BraessSarazin: M=alpha*diag(Auu). No auto tuning until
    // a qualified pilot contraction is available for the assembled operator.
    double braess_sarazin_alpha = 1.0;
};

struct MatrixCharacteristics {
    std::size_t equations = 0;
    double average_nnz_per_row = 0.0;
    double coefficient_range = 1.0;
    bool square = true;
    bool numerically_symmetric = false;
    bool diagonally_dominant = false;
    bool strongly_scaled = false;
    // Large diagonal coefficient spread; this is not a geometric anisotropy metric.
    bool large_diagonal_dynamic_range = false;
    // Algebraic proxy: in at least one row, negative off-diagonal couplings
    // differ by >= 10^3. This is not a mesh aspect-ratio measurement.
    bool strongly_directional_coupling = false;
    bool positive_diagonal = true;
    bool saddle_point = false;
};

struct LinearSolverPlan {
    LinearProblemKind problem = LinearProblemKind::General;
    KrylovModel krylov = KrylovModel::GMRES;
    PreconditionerModel preconditioner = PreconditionerModel::Jacobi;
    NullSpaceModel null_space = NullSpaceModel::None;
    bool automatic_krylov = true;
    bool automatic_preconditioner = true;
    std::string reason;
};

inline bool is_available(KrylovModel model) {
    return krylov_model_catalog()[static_cast<std::size_t>(model)].availability ==
           ModelAvailability::Available;
}

inline bool is_available(PreconditionerModel model) {
    return preconditioner_model_catalog()[static_cast<std::size_t>(model)].availability ==
           ModelAvailability::Available;
}

inline LinearSolverPlan select_linear_solver(LinearProblemKind problem,
                                             const MatrixCharacteristics& characteristics,
                                             const LinearSolverRequest& request = {}) {
    if (request.gmres_restart <= 0)
        throw std::invalid_argument("GMRES restart must be positive");
    if (!(request.braess_sarazin_alpha >= 0.25 &&
          request.braess_sarazin_alpha <= 4.0) ||
        !std::isfinite(request.braess_sarazin_alpha))
        throw std::invalid_argument("Braess-Sarazin alpha must be finite in [0.25,4]");
    if (!characteristics.square)
        throw std::invalid_argument("automatic linear-solver selection requires a square matrix");
    if (!is_available(request.krylov))
        throw std::invalid_argument(
            std::string("requested Krylov model is not implemented: ") +
            to_string(request.krylov));
    if (!is_available(request.preconditioner))
        throw std::invalid_argument(
            std::string("requested preconditioner model is not implemented: ") +
            to_string(request.preconditioner));

    const std::size_t equations = characteristics.equations;
    LinearSolverPlan plan;
    plan.problem = problem;
    plan.automatic_krylov = request.krylov == KrylovModel::Auto;
    plan.automatic_preconditioner = request.preconditioner == PreconditionerModel::Auto;
    plan.null_space = request.null_space;

    switch (problem) {
        case LinearProblemKind::PressurePoisson:
        case LinearProblemKind::Diffusion:
            if (characteristics.numerically_symmetric &&
                characteristics.diagonally_dominant &&
                characteristics.positive_diagonal &&
                characteristics.strongly_directional_coupling &&
                request.null_space != NullSpaceModel::Constant) {
                // The existing certificate validates damping, not the
                // complete SPD character of the multilevel V-cycle.
                // Use flexible right-preconditioned Krylov rather than CG
                // for the difficult directionally-coupled pressure path.
                plan.krylov = KrylovModel::FGMRES;
                plan.preconditioner = equations < 24 ? PreconditionerModel::Jacobi
                                                     : PreconditionerModel::NativeAMG;
                plan.reason = "strong directional off-diagonal coupling: flexible GMRES and strength-aware algebraic multigrid; true residual certified";
            } else if (characteristics.numerically_symmetric &&
                characteristics.diagonally_dominant &&
                characteristics.positive_diagonal &&
                !characteristics.large_diagonal_dynamic_range &&
                !characteristics.strongly_scaled) {
                plan.krylov = KrylovModel::CG;
                plan.preconditioner = equations < 24 ? PreconditionerModel::Jacobi
                                                     : PreconditionerModel::NativeAMG;
                plan.reason = "symmetric diagonally-dominant, moderately scaled elliptic matrix: CG with an SPD-qualified preconditioner";
            } else {
                // Symmetry alone is not enough to justify CG with NativeAMG:
                // high diagonal dynamic range / coefficient scaling can make
                // the current AMG application non-SPD in practice. Keep the
                // physical matrix unchanged and choose a nonsymmetric-safe
                // Krylov/preconditioner pair for this unqualified regime.
                plan.krylov = KrylovModel::GMRES;
                plan.preconditioner = equations < 24 ? PreconditionerModel::Jacobi
                                                     : PreconditionerModel::ILU0;
                plan.reason = (characteristics.large_diagonal_dynamic_range || characteristics.strongly_scaled)
                    ? "elliptic matrix has a large diagonal dynamic range or coefficient spread: GMRES with a nonsymmetric-safe preconditioner"
                    : "elliptic matrix lacks a verified SPD signature: GMRES with a nonsymmetric-safe preconditioner";
            }
            break;
        case LinearProblemKind::Momentum:
        case LinearProblemKind::ScalarTransport:
            plan.krylov = characteristics.strongly_scaled ? KrylovModel::FGMRES
                                                           : KrylovModel::BiCGStab;
            plan.preconditioner = equations < 24 ? PreconditionerModel::Jacobi
                                                 : PreconditionerModel::ILU0;
            plan.reason = characteristics.strongly_scaled
                ? "strong coefficient scaling: flexible GMRES preserves a conservative nonsymmetric policy"
                : "nonsymmetric transport matrix: BiCGStab with local factorization";
            break;
        case LinearProblemKind::CoupledPressureVelocity:
            if (!characteristics.saddle_point)
                throw std::invalid_argument(
                    "coupled selection requires an explicit saddle-point matrix classification");
            plan.krylov = KrylovModel::FGMRES;
            plan.preconditioner = PreconditionerModel::CoupledBlockSchur;
            plan.reason = "saddle-point block system: flexible GMRES with a CFD Schur approximation";
            break;
        case LinearProblemKind::General:
            plan.krylov = KrylovModel::GMRES;
            plan.preconditioner = PreconditionerModel::Jacobi;
            plan.reason = "general operator: conservative GMRES/Jacobi policy";
            break;
    }

    if (characteristics.large_diagonal_dynamic_range)
        plan.reason += "; large diagonal dynamic range detected — AMG qualification remains required";
    if (characteristics.strongly_scaled)
        plan.reason += "; coefficient scaling detected";

    if (!plan.automatic_krylov) plan.krylov = request.krylov;
    if (!plan.automatic_preconditioner) plan.preconditioner = request.preconditioner;
    // Line-smoothed AMG is not a certified SPD preconditioner. With an
    // explicit line-AMG request and Krylov=auto, choose flexible GMRES,
    // rather than selecting CG first and rejecting a valid scalar path.
    if (plan.automatic_krylov &&
        plan.preconditioner == PreconditionerModel::AnisotropicLineAMG) {
        plan.krylov = KrylovModel::FGMRES;
        plan.reason += "; experimental line-AMG requires flexible nonsymmetric Krylov";
    }

    // The optional algebraic-line AMG smoothers are qualified with a flexible
    // outer Krylov method. They must not inherit CG from a benign SPD-looking
    // matrix solely because the user specified a preconditioner and left
    // Krylov on Auto. This is an explicit algorithmic pairing contract.
    const bool anisotropic_amg =
        plan.preconditioner == PreconditionerModel::AnisotropicBoomerAMG ||
        plan.preconditioner == PreconditionerModel::AnisotropicSmoothedAggregationAMG ||
        plan.preconditioner == PreconditionerModel::DirectionalSemiCoarseningAMG ||
        plan.preconditioner == PreconditionerModel::EnergyMinimizingAMG;
    if (plan.automatic_krylov &&
        (anisotropic_amg ||
         plan.preconditioner == PreconditionerModel::AnisotropicLine)) {
        plan.krylov = KrylovModel::FGMRES;
        plan.reason += "; anisotropic line preconditioning uses flexible GMRES";
    }

    // Row/column equilibration is generally NOT a symmetry-preserving
    // congruence transformation. Reject explicit CG and use nonsymmetric-safe
    // Krylov for auto selection; symmetric-diagonal scaling is admissible.
    const bool nonsymmetric_scaling =
        request.scaling == LinearScalingModel::Row ||
        request.scaling == LinearScalingModel::Column ||
        request.scaling == LinearScalingModel::RowColumn;
    if (plan.krylov == KrylovModel::CG && nonsymmetric_scaling) {
        if (!plan.automatic_krylov)
            throw std::invalid_argument("CG cannot use nonsymmetric algebraic scaling");
        plan.krylov = KrylovModel::FGMRES;
        plan.reason += "; nonsymmetric scaling requires flexible GMRES";
    }

    if (plan.null_space == NullSpaceModel::Constant && plan.automatic_preconditioner)
        plan.reason += "; projected constant null space is enforced across Krylov and compatible preconditioning";

    const bool cg_problem = problem == LinearProblemKind::PressurePoisson ||
                            problem == LinearProblemKind::Diffusion;
    if (plan.krylov == KrylovModel::CG &&
        (!cg_problem || !characteristics.numerically_symmetric ||
         !characteristics.positive_diagonal))
        throw std::invalid_argument(
            "CG requires a symmetric, positive-diagonal pressure/diffusion operator");
    if (plan.krylov == KrylovModel::CG &&
        plan.preconditioner != PreconditionerModel::None &&
        plan.preconditioner != PreconditionerModel::Jacobi &&
        plan.preconditioner != PreconditionerModel::NativeAMG &&
        plan.preconditioner != PreconditionerModel::SmoothedAggregationAMG)
        throw std::invalid_argument(
            "CG requires an SPD-qualified preconditioner (none, Jacobi, or a native AMG variant)");
    if (plan.preconditioner == PreconditionerModel::AnisotropicLineAMG &&
        (plan.krylov == KrylovModel::CG ||
         problem == LinearProblemKind::CoupledPressureVelocity))
        throw std::invalid_argument(
            "anisotropic line AMG requires nonsymmetric scalar Krylov");
    if ((plan.preconditioner == PreconditionerModel::NativeAMG ||
         plan.preconditioner == PreconditionerModel::SmoothedAggregationAMG ||
         anisotropic_amg) && !cg_problem)
        throw std::invalid_argument(
            "native AMG is currently qualified only for SPD elliptic problems");
    if (anisotropic_amg &&
        (!characteristics.numerically_symmetric || !characteristics.positive_diagonal))
        throw std::invalid_argument(
            "anisotropic line-AMG currently requires a symmetric positive-diagonal operator");
    if (plan.preconditioner == PreconditionerModel::BraessSarazin &&
        (problem != LinearProblemKind::CoupledPressureVelocity ||
         plan.krylov != KrylovModel::FGMRES))
        throw std::invalid_argument(
            "Braess-Sarazin is opt-in for coupled FGMRES only");
    if (plan.preconditioner == PreconditionerModel::VankaPressurePatch &&
        (problem != LinearProblemKind::CoupledPressureVelocity ||
         plan.krylov != KrylovModel::FGMRES))
        throw std::invalid_argument(
            "pressure-patch Vanka is opt-in for coupled FGMRES only");
    if (plan.preconditioner == PreconditionerModel::NativeFieldSplit &&
        problem != LinearProblemKind::CoupledPressureVelocity)
        throw std::invalid_argument("native FieldSplit requires a coupled pressure-velocity profile");
    if ((plan.preconditioner == PreconditionerModel::CoupledBlockSchur ||
         plan.preconditioner == PreconditionerModel::PCD ||
         plan.preconditioner == PreconditionerModel::MGR) &&
        problem != LinearProblemKind::CoupledPressureVelocity)
        throw std::invalid_argument(
            "coupled block Schur and MGR require a pressure-velocity profile");
    if (plan.null_space == NullSpaceModel::Constant && !cg_problem)
        throw std::invalid_argument(
            "constant null space is currently qualified only for pressure/diffusion problems");
    if (plan.null_space == NullSpaceModel::Constant &&
        plan.krylov != KrylovModel::CG &&
        plan.krylov != KrylovModel::GMRES &&
        plan.krylov != KrylovModel::FGMRES)
        throw std::invalid_argument(
            "constant null space requires CG, GMRES or FGMRES with explicit projection");
    return plan;
}

inline LinearSolverPlan select_linear_solver(LinearProblemKind problem,
                                             std::size_t equations,
                                             const LinearSolverRequest& request = {}) {
    MatrixCharacteristics characteristics;
    characteristics.equations = equations;
    characteristics.square = true;
    characteristics.numerically_symmetric =
        problem == LinearProblemKind::PressurePoisson ||
        problem == LinearProblemKind::Diffusion;
    characteristics.diagonally_dominant =
        problem == LinearProblemKind::PressurePoisson ||
        problem == LinearProblemKind::Diffusion;
    characteristics.saddle_point =
        problem == LinearProblemKind::CoupledPressureVelocity;
    return select_linear_solver(problem, characteristics, request);
}

} // namespace cfdx::core
