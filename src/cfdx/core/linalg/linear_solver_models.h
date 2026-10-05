#pragma once

#include <array>
#include <cstddef>
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
    FSAI,
    RAS,
    NativeFieldSplit,
    CoupledBlockSchur,
    PCD,
    LSC,
    MGR
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

inline const std::array<SolverModelDescriptor, 15>& preconditioner_model_catalog() {
    static const std::array<SolverModelDescriptor, 15> models{{
        {"auto", ModelAvailability::Available, false, true},
        {"none", ModelAvailability::Available, false, true},
        {"jacobi", ModelAvailability::Available, false, true},
        {"gauss_seidel", ModelAvailability::Available, false, false},
        {"ilu0", ModelAvailability::Available, false, false},
        {"ilut", ModelAvailability::Planned, false, false},
        {"native_amg", ModelAvailability::Available, false, false},
        {"smoothed_aggregation_amg", ModelAvailability::Available, false, false},
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
        {"mgr", ModelAvailability::Available, false, false}
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

struct LinearSolverRequest {
    KrylovModel krylov = KrylovModel::Auto;
    PreconditionerModel preconditioner = PreconditionerModel::Auto;
    int gmres_restart = 40;
    bool allow_fallback = false;
    NullSpaceModel null_space = NullSpaceModel::None;
    // Optional algebraic scaling. None preserves the historical production path.
    LinearScalingModel scaling = LinearScalingModel::None;
};

struct MatrixCharacteristics {
    std::size_t equations = 0;
    double average_nnz_per_row = 0.0;
    double coefficient_range = 1.0;
    bool square = true;
    bool numerically_symmetric = false;
    bool diagonally_dominant = false;
    bool strongly_scaled = false;
    bool anisotropic = false;
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
            if (characteristics.numerically_symmetric && characteristics.diagonally_dominant) {
                plan.krylov = KrylovModel::CG;
                plan.preconditioner = equations < 24 ? PreconditionerModel::Jacobi
                                                     : PreconditionerModel::NativeAMG;
                plan.reason = "symmetric diagonally-dominant elliptic matrix: CG with an SPD-qualified preconditioner";
            } else {
                plan.krylov = KrylovModel::GMRES;
                plan.preconditioner = equations < 24 ? PreconditionerModel::Jacobi
                                                     : PreconditionerModel::ILU0;
                plan.reason = "elliptic matrix lacks a verified SPD signature: GMRES with a nonsymmetric-safe preconditioner";
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

    if (characteristics.anisotropic)
        plan.reason += "; anisotropy detected — dedicated AMG anisotropy qualification remains required";
    if (characteristics.strongly_scaled)
        plan.reason += "; coefficient scaling detected";

    if (!plan.automatic_krylov) plan.krylov = request.krylov;
    if (!plan.automatic_preconditioner) plan.preconditioner = request.preconditioner;

    if (plan.null_space == NullSpaceModel::Constant && plan.automatic_preconditioner)
        plan.reason += "; projected constant null space is enforced across Krylov and compatible preconditioning";

    const bool cg_problem = problem == LinearProblemKind::PressurePoisson ||
                            problem == LinearProblemKind::Diffusion;
    if (plan.krylov == KrylovModel::CG && !cg_problem)
        throw std::invalid_argument("CG requires an SPD pressure/diffusion problem profile");
    if (plan.krylov == KrylovModel::CG &&
        plan.preconditioner != PreconditionerModel::None &&
        plan.preconditioner != PreconditionerModel::Jacobi &&
        plan.preconditioner != PreconditionerModel::NativeAMG &&
        plan.preconditioner != PreconditionerModel::SmoothedAggregationAMG)
        throw std::invalid_argument(
            "CG requires an SPD-qualified preconditioner (none, Jacobi, or a native AMG variant)");
    if ((plan.preconditioner == PreconditionerModel::NativeAMG ||
         plan.preconditioner == PreconditionerModel::SmoothedAggregationAMG) && !cg_problem)
        throw std::invalid_argument(
            "native AMG is currently qualified only for SPD elliptic problems");
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
    if (plan.null_space == NullSpaceModel::Constant && plan.krylov != KrylovModel::CG)
        throw std::invalid_argument("constant null space currently requires projected CG");
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
