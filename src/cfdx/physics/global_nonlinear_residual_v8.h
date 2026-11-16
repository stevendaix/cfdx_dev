#pragma once

#include "cfdx/physics/momentum_residual_certification.h"
#include "cfdx/physics/pressure_velocity_algorithms.h"
#include "cfdx/physics/solver_control.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

// Certificate of the returned *physical* finite-volume state. No Krylov,
// AMG, Schur or PCD estimate enters this object. It can be used unchanged
// by segregated and all six production COUPLED Schur approximations.
//
// The momentum certificate is from the three reassembled equations, using
// the accepted face flux and, with SST, the FINAL eddy viscosity. The k and
// omega components must come from independently reassembled returned-state
// SST equations, not the linear scalar solver's reported residuals.
struct GlobalNonlinearResidualV8 {
    PressureVelocityAlgorithm algorithm = PressureVelocityAlgorithm::SIMPLE;
    CoupledSchurModel schur = CoupledSchurModel::BlockLocal;
    MomentumResidualCertificate momentum{};
    double continuity = std::numeric_limits<double>::infinity();
    double k = std::numeric_limits<double>::quiet_NaN();
    double omega = std::numeric_limits<double>::quiet_NaN();
    double eddy_viscosity_closure = std::numeric_limits<double>::quiet_NaN();
    double merit = std::numeric_limits<double>::infinity();
    bool sst = false;
    bool sst_transport_committed = false;
    bool bounded = false;
    bool all_finite = false;
    bool converged = false;
};

inline GlobalNonlinearResidualV8 certify_global_nonlinear_residual_v8(
    PressureVelocityAlgorithm algorithm, CoupledSchurModel schur,
    const MomentumResidualCertificate& momentum,
    double continuity,
    const ConvergenceCriteria& gates,
    bool has_sst = false,
    double reassembled_k = std::numeric_limits<double>::quiet_NaN(),
    double reassembled_omega = std::numeric_limits<double>::quiet_NaN(),
    double closure = std::numeric_limits<double>::quiet_NaN(),
    bool sst_transport_committed = false,
    bool sst_bounded = false)
{
    if (!(std::isfinite(gates.relative_tolerance) && gates.relative_tolerance > 0.) ||
        !(std::isfinite(gates.continuity_tolerance) && gates.continuity_tolerance > 0.))
        throw std::invalid_argument("global nonlinear V8: invalid original physical gates");
    GlobalNonlinearResidualV8 result;
    result.algorithm = algorithm;
    result.schur = schur;
    result.momentum = momentum;
    result.continuity = continuity;
    result.sst = has_sst;
    result.sst_transport_committed = sst_transport_committed;
    result.bounded = has_sst ? sst_bounded : true;
    result.k = reassembled_k;
    result.omega = reassembled_omega;
    result.eddy_viscosity_closure = closure;
    result.all_finite = std::isfinite(continuity) && continuity >= 0. &&
        std::isfinite(momentum.global_relative) && momentum.global_relative >= 0. &&
        std::isfinite(momentum.max_active_component_relative) &&
        momentum.max_active_component_relative >= 0.;
    for (std::size_t d = 0; d < 3; ++d)
        result.all_finite = result.all_finite &&
            std::isfinite(momentum.component_relative[d]) &&
            momentum.component_relative[d] >= 0.;
    if (has_sst)
        result.all_finite = result.all_finite &&
            std::isfinite(reassembled_k) && reassembled_k >= 0. &&
            std::isfinite(reassembled_omega) && reassembled_omega >= 0. &&
            std::isfinite(closure) && closure >= 0.;
    if (!result.all_finite) return result;
    result.merit = std::max({momentum.global_relative,
        momentum.max_active_component_relative, continuity});
    if (has_sst) result.merit = std::max({result.merit, reassembled_k,
        reassembled_omega, closure});
    const bool sst_ok = !has_sst ||
        (sst_transport_committed && sst_bounded &&
         reassembled_k <= gates.relative_tolerance &&
         reassembled_omega <= gates.relative_tolerance &&
         closure <= gates.relative_tolerance);
    result.converged = momentum.converged &&
        continuity <= gates.continuity_tolerance && sst_ok;
    return result;
}

} // namespace cfdx::physics
