#pragma once
// V80 experimental NATIVE high/low algebraic refinement used before the
// pre-existing real-face Rhie-Chow refinement. The CSR is only a preconditioner;
// no numerical CSR certificate may authorize a CFD state transaction.
#include "cfdx/core/linalg/pressure_compensated_refinement.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {
struct NativeAlgebraicDefectV80 {
    double initial_relative=std::numeric_limits<double>::infinity();
    double final_relative=std::numeric_limits<double>::infinity();
    double final_double_only_relative=std::numeric_limits<double>::infinity();
    std::size_t attempted=0,accepted=0;
    bool converged=false;
    bool production_commit_authorized=false;
};
// InnerSolve: bool(const Vector& defect, Vector& delta) called with zero delta.
// An inexact FGMRES step is permitted only when its TRUE residual improves.
// All updates are high/low. The independent physical-face gate is in the caller.
template<class InnerSolve>
NativeAlgebraicDefectV80 refine_native_pressure_algebraic_v80(
    const cfdx::core::SparseMatrix& a,const cfdx::core::Vector& b,
    cfdx::core::Vector& high,cfdx::core::Vector& low,
    double original_tolerance,std::size_t max_passes,InnerSolve&& inner_solve)
{
    if (max_passes>8 || !(original_tolerance>0.) ||
        !std::isfinite(original_tolerance))
        throw std::invalid_argument("V80: invalid algebraic refinement controls");
    const auto r=cfdx::core::refine_pressure_high_low(
        a,b,high,low,original_tolerance,max_passes,
        [&](const cfdx::core::SparseMatrix&,const cfdx::core::Vector& defect,
            cfdx::core::Vector& delta) {
            if (!inner_solve(defect,delta)) return false;
            for (std::size_t i=0;i<delta.size();++i)
                if (!std::isfinite(delta(i))) return false;
            return true;
        });
    NativeAlgebraicDefectV80 result;
    result.initial_relative=r.initial_relative_residual;
    result.final_relative=r.final_relative_residual;
    result.final_double_only_relative=r.final_double_only_residual;
    result.attempted=r.passes_attempted;
    result.accepted=r.passes_accepted;
    result.converged=r.algebraic_converged;
    return result;
}
} // namespace cfdx::physics