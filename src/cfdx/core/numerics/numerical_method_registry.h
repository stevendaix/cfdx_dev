#pragma once

#include "cfdx/core/numerics/numerical_method_contract.h"

#include <vector>

namespace cfdx::core {

inline std::vector<NumericalMethodContract> numerical_method_registry()
{
    using F = NumericalMethodFamily;
    using S = VerificationStatus;
    using C = ConservationContract;
    return {
        {"gradient.gauss_cell","Green-Gauss cell gradient",F::Gradient,S::Implemented,C::LocalFaceConservative,false,false,true,1,0,
         "grad(phi)=1/V sum_f phi_f Sf","numerics.gradient.gauss",{"constant","linear","quadratic","skewed","non_orthogonal"}},
        {"gradient.least_squares","Weighted least-squares gradient",F::Gradient,S::Implemented,C::NotApplicable,false,false,true,1,0,
         "argmin_g sum_N w_N(phi_N-phi_P-g.dx)^2","numerics.gradient.least_squares",{"constant","linear","rank_detection","skewed"}},
        {"interpolation.linear","Cell-to-face linear interpolation",F::Interpolation,S::Implemented,C::LocalFaceConservative,false,false,true,2,0,
         "phi_f=(1-w)phi_P+w phi_N","numerics.interpolation.linear",{"linear_exactness","conservation"}},
        {"interpolation.upwind","First-order upwind interpolation",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "phi_f=phi_upwind","numerics.interpolation.upwind",{"boundedness","monotonicity","conservation"}},
        {"convection.tvd.minmod","TVD MinMod",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "MUSCL with psi=minmod","numerics.convection.tvd.minmod",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.vanleer","TVD Van Leer",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "MUSCL with Van Leer limiter","numerics.convection.tvd.vanleer",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.superbee","TVD Superbee",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "MUSCL with Superbee limiter","numerics.convection.tvd.superbee",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.vanalbada","TVD Van Albada",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "MUSCL with Van Albada limiter","numerics.convection.tvd.vanalbada",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.mc","TVD monotonized central",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "MUSCL with MC limiter","numerics.convection.tvd.mc",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"diffusion.orthogonal","Orthogonal two-point diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,true,2,0,
         "Gamma Sf.d/|d|^2 delta_phi","numerics.diffusion.orthogonal",{"MMS","refinement"}},
        {"diffusion.corrected","Corrected non-orthogonal diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,false,1,0,
         "orthogonal + tangential gradient correction","numerics.diffusion.corrected",{"MMS","non_orthogonal_refinement"}},
        {"diffusion.limited","Limited non-orthogonal diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,false,1,0,
         "orthogonal + bounded correction","numerics.diffusion.limited",{"MMS","non_orthogonal_robustness"}},
        {"temporal.euler_explicit","Explicit Euler",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,1,
         "phi(n+1)=phi(n)+dt RHS(n)","numerics.temporal.euler_explicit",{"analytical_decay","temporal_refinement"}},
        {"temporal.euler_implicit","Implicit Euler",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,1,
         "backward Euler","numerics.temporal.euler_implicit",{"analytical_decay","temporal_refinement"}},
        {"temporal.crank_nicolson","Crank-Nicolson",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,2,
         "trapezoidal time integration","numerics.temporal.crank_nicolson",{"analytical_decay","temporal_refinement"}},
        {"temporal.bdf2","BDF2",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,2,
         "3phi(n+1)-4phi(n)+phi(n-1) over 2dt","numerics.temporal.bdf2",{"transient_MMS","temporal_refinement","restart_history"}},
        {"linear.cg","Conjugate Gradient",F::LinearSolver,S::Verified,C::NotApplicable,false,false,false,0,0,
         "Krylov CG for SPD systems","linear.cg",{"SPD_exact_solve","true_residual"}},
        {"linear.bicgstab","BiCGStab",F::LinearSolver,S::Verified,C::NotApplicable,false,false,false,0,0,
         "BiCGStab Krylov iteration","linear.bicgstab",{"nonsymmetric_exact_solve","true_residual"}},
        {"linear.gmres","GMRES",F::LinearSolver,S::Verified,C::NotApplicable,false,false,false,0,0,
         "Restarted/unrestarted GMRES operator iteration","linear.gmres",{"nonsymmetric_exact_solve","true_residual"}},
        {"pressure_velocity.simple","SIMPLE",F::PressureVelocity,S::Implemented,C::GloballyConservative,false,false,false,0,0,
         "segregated pressure correction","pressure_velocity.simple",{"Couette","Poiseuille","cavity","skew_mesh"}},
        {"pressure_velocity.simplec","SIMPLEC",F::PressureVelocity,S::Implemented,C::GloballyConservative,false,false,false,0,0,
         "consistent SIMPLE pressure correction","pressure_velocity.simplec",{"Couette","Poiseuille","cavity","skew_mesh"}},
        {"pressure_velocity.piso","PISO",F::PressureVelocity,S::Implemented,C::GloballyConservative,false,false,false,0,0,
         "multiple pressure corrections","pressure_velocity.piso",{"Couette","Poiseuille","cavity","skew_mesh"}},
        {"pressure_velocity.pimple","PIMPLE",F::PressureVelocity,S::Implemented,C::GloballyConservative,false,false,false,0,0,
         "outer-loop PISO/SIMPLE hybrid","pressure_velocity.pimple",{"Couette","Poiseuille","cavity","skew_mesh"}},
        {"pressure_velocity.fractional_step","Fractional step",F::PressureVelocity,S::Implemented,C::GloballyConservative,false,false,false,0,0,
         "projection/fractional-step coupling","pressure_velocity.fractional_step",{"Couette","Poiseuille","cavity","skew_mesh"}},
        {"pressure_velocity.coupled","Fully coupled U-p",F::PressureVelocity,S::Implemented,C::GloballyConservative,false,false,false,0,0,
         "block coupled momentum-pressure system","pressure_velocity.coupled",{"Couette","Poiseuille","cavity","skew_mesh"}}
    };
}

inline void validate_numerical_method_registry()
{
    for (const auto& method : numerical_method_registry())
        validate_numerical_method_contract(method);
}

} // namespace cfdx::core
