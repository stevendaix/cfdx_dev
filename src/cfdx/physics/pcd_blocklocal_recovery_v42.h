#pragma once
// Certify an explicitly selected BlockLocal retry against the ORIGINAL A,b.
// A true CSR residual checks both momentum and pressure/gauge equations;
// this is a linear-system certificate, never a nonlinear CFD qualification.
#include "cfdx/physics/coupled_block_residual_v16.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {
struct PcdBlockLocalRecoveryDecisionV42 {
    bool adopt = false;
    double true_relative = 0.0;
    CoupledBlockResidualV16 blocks;
};
inline PcdBlockLocalRecoveryDecisionV42 certify_pcd_blocklocal_recovery_v42(
    const cfdx::core::SparseMatrix& A, const cfdx::core::Vector& b,
    const cfdx::core::Vector& candidate, std::size_t velocity_rows,
    double unchanged_tolerance, bool solver_reported_converged)
{
    if (!(unchanged_tolerance > 0.0) || !std::isfinite(unchanged_tolerance))
        throw std::invalid_argument("V42 PCD retry: invalid original tolerance");
    PcdBlockLocalRecoveryDecisionV42 d;
    d.blocks=independent_coupled_block_residual_v16(A,b,candidate,velocity_rows);
    d.true_relative=static_cast<double>(d.blocks.relative_to_global_rhs());
    d.adopt=solver_reported_converged && std::isfinite(d.true_relative) &&
        d.true_relative<=unchanged_tolerance;
    return d;
}
} // namespace cfdx::physics
