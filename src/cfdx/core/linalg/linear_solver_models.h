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
    LSC,
    MGR
};

enum class ModelAvailability { Available, Planned };

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

inline const std::array<SolverModelDescriptor, 14>& preconditioner_model_catalog() {
    static const std::array<SolverModelDescriptor, 14> models{{
        {"auto", ModelAvailability::Available, false, true},
        {"none", ModelAvailability::Available, false, true},
        {"jacobi", ModelAvailability::Available, false, true},
        {"gauss_seidel", ModelAvailability::Available, false, false},
        {"ilu0", ModelAvailability::Available, false, false},
        {"ilut", ModelAvailability::Planned, false, false},
        {"native_amg", ModelAvailability::Available, false, false},
        {"smoothed_aggregation_amg", ModelAvailability::Planned, false, false},
        {"fsai", ModelAvailability::Planned, false, false},
        {"ras", ModelAvailability::Planned, false, true},
        {"native_fieldsplit", ModelAvailability::Available, false, false},
        {"coupled_block_schur", ModelAvailability::Available, false, false},
        {"lsc", ModelAvailability::Planned, false, false},
        {"mgr", ModelAvailability::Planned, false, false}
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

struct LinearSolverRequest {
    KrylovModel krylov = KrylovModel::Auto;
    PreconditionerModel preconditioner = PreconditionerModel::Auto;
    int gmres_restart = 40;
    bool allow_fallback = false;
};

struct LinearSolverPlan {
    LinearProblemKind problem = LinearProblemKind::General;
    KrylovModel krylov = KrylovModel::GMRES;
    PreconditionerModel preconditioner = PreconditionerModel::Jacobi;
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
                                             std::size_t equations,
                                             const LinearSolverRequest& request = {}) {
    if (request.gmres_restart <= 0)
        throw std::invalid_argument("GMRES restart must be positive");
    if (!is_available(request.krylov))
        throw std::invalid_argument(
            std::string("requested Krylov model is not implemented: ") +
            to_string(request.krylov));
    if (!is_available(request.preconditioner))
        throw std::invalid_argument(
            std::string("requested preconditioner model is not implemented: ") +
            to_string(request.preconditioner));

    LinearSolverPlan plan;
    plan.problem = problem;
    plan.automatic_krylov = request.krylov == KrylovModel::Auto;
    plan.automatic_preconditioner =
        request.preconditioner == PreconditionerModel::Auto;

    switch (problem) {
        case LinearProblemKind::PressurePoisson:
        case LinearProblemKind::Diffusion:
            plan.krylov = KrylovModel::CG;
            plan.preconditioner = equations < 24
                ? PreconditionerModel::Jacobi
                : PreconditionerModel::NativeAMG;
            plan.reason = "SPD elliptic operator: CG with an SPD preconditioner";
            break;
        case LinearProblemKind::Momentum:
        case LinearProblemKind::ScalarTransport:
            plan.krylov = KrylovModel::BiCGStab;
            plan.preconditioner = equations < 24
                ? PreconditionerModel::Jacobi
                : PreconditionerModel::ILU0;
            plan.reason = "nonsymmetric transport operator: BiCGStab with local factorization";
            break;
        case LinearProblemKind::CoupledPressureVelocity:
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

    if (!plan.automatic_krylov) plan.krylov = request.krylov;
    if (!plan.automatic_preconditioner)
        plan.preconditioner = request.preconditioner;

    const bool cg_problem = problem == LinearProblemKind::PressurePoisson ||
                            problem == LinearProblemKind::Diffusion;
    if (plan.krylov == KrylovModel::CG && !cg_problem)
        throw std::invalid_argument("CG requires an SPD pressure/diffusion problem profile");
    if (plan.krylov == KrylovModel::CG &&
        plan.preconditioner != PreconditionerModel::None &&
        plan.preconditioner != PreconditionerModel::Jacobi &&
        plan.preconditioner != PreconditionerModel::NativeAMG)
        throw std::invalid_argument(
            "CG requires an SPD-qualified preconditioner (none, Jacobi, or native AMG)");
    if (plan.preconditioner == PreconditionerModel::NativeAMG && !cg_problem)
        throw std::invalid_argument("native AMG is currently qualified only for SPD elliptic problems");
    if (plan.preconditioner == PreconditionerModel::NativeFieldSplit &&
        problem != LinearProblemKind::CoupledPressureVelocity)
        throw std::invalid_argument("native FieldSplit requires a coupled pressure-velocity profile");
    if (plan.preconditioner == PreconditionerModel::CoupledBlockSchur &&
        problem != LinearProblemKind::CoupledPressureVelocity)
        throw std::invalid_argument("coupled block Schur requires a pressure-velocity profile");
    return plan;
}

} // namespace cfdx::core
