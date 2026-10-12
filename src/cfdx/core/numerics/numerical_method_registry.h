#pragma once

#include "cfdx/core/numerics/numerical_method_contract.h"

#include <vector>
#include <unordered_set>

namespace cfdx::core {

inline std::vector<NumericalMethodContract> numerical_method_registry()
{
    using F = NumericalMethodFamily;
    using S = VerificationStatus;
    using C = ConservationContract;
    return {
        {"gradient.gauss_cell","Green-Gauss cell gradient",F::Gradient,S::Implemented,C::LocalFaceConservative,false,false,true,1,0,
         "grad(phi)=1/V sum_f phi_f Sf","numerics.gradient.gauss",{"constant","linear","quadratic","skewed","non_orthogonal"}},
        {"gradient.gauss_vertex","Vertex-based Green-Gauss gradient",F::Gradient,S::Implemented,C::LocalFaceConservative,false,false,false,1,0,
         "inverse-distance vertex average; Green-Gauss over vertex-derived face values","numerics.gradient.gauss_vertex",{"constant","linear","refinement"}},
        {"gradient.gauss_point","Point-linear (skew-corrected) Green-Gauss gradient",F::Gradient,S::Implemented,C::LocalFaceConservative,false,false,false,2,0,
         "face value = signed-area average of inverse-distance vertex values and the centre-line face-point value","numerics.gradient.gauss_point",{"linear_exactness","skewed","polyhedral","refinement"}},
        {"gradient.least_squares","Weighted least-squares gradient",F::Gradient,S::Implemented,C::NotApplicable,false,false,true,1,0,
         "argmin_g sum_N w_N(phi_N-phi_P-g.dx)^2","numerics.gradient.least_squares",{"constant","linear","rank_detection","skewed"}},
        {"gradient.least_squares_weighted","Explicit weighted least-squares gradient",F::Gradient,S::Implemented,C::NotApplicable,false,false,false,2,0,
         "argmin_g sum_N w_N(phi_N-phi_P-g.dx)^2 with selectable distance weighting","numerics.gradient.weighted_least_squares",{"constant","linear","rank_detection","conditioning","skewed","polyhedral"}},
        {"gradient.least_squares_weighted_extended","Extended weighted least-squares gradient",F::Gradient,S::Implemented,C::NotApplicable,false,false,false,3,0,
         "same weighted LS fit over a two-ring cell stencil","numerics.gradient.weighted_least_squares_extended",{"constant","linear","smooth","tetrahedral","polyhedral","conditioning"}},
        {"gradient.least_squares_quadratic","Quadratic-basis least-squares gradient",F::Gradient,S::Implemented,C::NotApplicable,false,false,false,2,0,
         "linear + Hessian fit over two rings of face neighbours; best polyhedral accuracy","numerics.gradient.least_squares_quadratic",{"constant","linear","skewed","polyhedral","refinement"}},
        {"interpolation.linear","Cell-to-face linear interpolation",F::Interpolation,S::Implemented,C::LocalFaceConservative,false,false,true,2,0,
         "phi_f=(1-w)phi_P+w phi_N","numerics.interpolation.linear",{"linear_exactness","conservation"}},
        {"interpolation.upwind","First-order upwind interpolation",F::Interpolation,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "phi_f=phi_upwind","numerics.interpolation.upwind",{"boundedness","monotonicity","conservation"}},
        {"convection.upwind","First-order upwind convection",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "F phi_f with phi_f=phi_upwind","numerics.convection.upwind",{"boundedness","monotonicity","conservation"}},
        {"convection.second_order_upwind","Bounded second-order upwind reconstruction",F::Convection,S::Implemented,C::LocalFaceConservative,true,false,false,2,0,
         "phi_U + grad(phi)_U dot (Cf-C_U), bounded to [phi_U,phi_D]","numerics.convection.second_order_upwind",{"smooth_refinement","boundedness","conservation"}},
        {"convection.tvd.minmod","TVD MinMod",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,1,0,
         "MUSCL with psi=minmod","numerics.convection.tvd.minmod",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.vanleer","TVD Van Leer",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,2,0,
         "MUSCL with psi=vanleer","numerics.convection.tvd.vanleer",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.superbee","TVD Superbee",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,2,0,
         "MUSCL with psi=superbee","numerics.convection.tvd.superbee",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.vanalbada","TVD Van Albada",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,2,0,
         "MUSCL with psi=vanalbada","numerics.convection.tvd.vanalbada",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.tvd.mc","TVD MC",F::Convection,S::Implemented,C::LocalFaceConservative,true,true,false,2,0,
         "MUSCL with psi=mc","numerics.convection.tvd.mc",{"boundedness","TVD","smooth_refinement","steep_gradient"}},
        {"convection.central","Central (average) convection",F::Convection,S::Implemented,C::LocalFaceConservative,false,false,false,2,0,
         "phi_f=0.5(phi_P+phi_N)","numerics.convection.central",{"smooth_refinement","conservation"}},
        {"convection.blended","Blended linear/upwind convection",F::Convection,S::Implemented,C::LocalFaceConservative,false,false,false,1,0,
         "phi_f=beta*0.5(phi_P+phi_N)+(1-beta)*upwind","numerics.convection.blended",{"refinement","skewed","bounded"}},
        {"convection.quick","Multidimensional quadratic upwind (QUICK-equivalent)",F::Convection,S::Implemented,C::LocalFaceConservative,false,false,false,3,0,
         "quadratic least-squares fit evaluated at Cf from the upwind cell","numerics.convection.quick",{"smooth_refinement","conservation","polyhedral","steep_gradient"}},
        {"convection.quick_bounded","Bounded multidimensional quadratic upwind",F::Convection,S::Implemented,C::LocalFaceConservative,true,false,false,2,0,
         "QUICK-equivalent reconstruction clipped to the adjacent-cell envelope","numerics.convection.quick_bounded",{"boundedness","conservation","polyhedral","steep_gradient"}},
        {"interpolation.limiter.minmod","MinMod limiter coefficient",F::Interpolation,S::Implemented,C::NotApplicable,true,true,false,0,0,
         "psi=max(0,min(1,r))","numerics.limiter.minmod",{"coefficient_exactness","boundedness"}},
        {"interpolation.limiter.vanleer","Van Leer limiter coefficient",F::Interpolation,S::Implemented,C::NotApplicable,true,true,false,0,0,
         "psi=(r+|r|)/(1+|r|)","numerics.limiter.vanleer",{"coefficient_exactness","boundedness"}},
        {"interpolation.limiter.superbee","Superbee limiter coefficient",F::Interpolation,S::Implemented,C::NotApplicable,true,true,false,0,0,
         "psi=max(0,max(min(1,2r),min(2,r)))","numerics.limiter.superbee",{"coefficient_exactness","boundedness"}},
        {"interpolation.limiter.vanalbada","Van Albada limiter coefficient",F::Interpolation,S::Implemented,C::NotApplicable,true,true,false,0,0,
         "psi=(r^2+r)/(r^2+1)","numerics.limiter.vanalbada",{"coefficient_exactness","boundedness"}},
        {"interpolation.limiter.mc","MC limiter coefficient",F::Interpolation,S::Implemented,C::NotApplicable,true,true,false,0,0,
         "psi=max(0,min(2r,(1+r)/2,2))","numerics.limiter.mc",{"coefficient_exactness","boundedness"}},
        {"diffusion.uncorrected","Uncorrected non-orthogonal diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,true,1,0,
         "orthogonal projection only","numerics.diffusion.uncorrected",{"MMS","non_orthogonal_refinement"}},
        {"diffusion.orthogonal","Orthogonal two-point diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,true,2,0,
         "Gamma Sf.d/|d|^2 delta_phi","numerics.diffusion.orthogonal",{"MMS","refinement"}},
        {"diffusion.corrected","Corrected non-orthogonal diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,false,1,0,
         "orthogonal + tangential gradient correction","numerics.diffusion.corrected",{"MMS","non_orthogonal_refinement"}},
        {"diffusion.limited","Limited non-orthogonal diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,false,1,0,
         "orthogonal + bounded correction","numerics.diffusion.limited",{"MMS","non_orthogonal_robustness"}},
        {"diffusion.over_relaxed","Over-relaxed non-orthogonal diffusion",F::Diffusion,S::Implemented,C::LocalFaceConservative,false,false,false,1,0,
         "alpha=|Sf|^2/(Sf.d), correction (Sf-alpha d).grad_f with Sf-alpha d orthogonal to Sf","numerics.diffusion.over_relaxed",{"linear_exactness","non_orthogonal_refinement","conservation"}},
        {"temporal.euler_explicit","Explicit Euler",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,1,
         "phi(n+1)=phi(n)+dt RHS(n)","numerics.temporal.euler_explicit",{"analytical_decay","temporal_refinement"}},
        {"temporal.euler_implicit","Implicit Euler",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,1,
         "backward Euler","numerics.temporal.euler_implicit",{"analytical_decay","temporal_refinement"}},
        {"temporal.crank_nicolson","Crank-Nicolson",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,2,
         "trapezoidal time integration","numerics.temporal.crank_nicolson",{"analytical_decay","temporal_refinement"}},
        {"temporal.bdf2","BDF2",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,2,
         "3phi(n+1)-4phi(n)+phi(n-1) over 2dt","numerics.temporal.bdf2",{"transient_MMS","temporal_refinement","restart_history"}},
        {"temporal.rk2","Explicit RK2 (midpoint)",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,2,
         "two-stage explicit Runge-Kutta","numerics.temporal.rk2",{"analytical_decay","temporal_refinement"}},
        {"temporal.rk3","Explicit RK3 (Williamson low-storage)",F::Temporal,S::Implemented,C::NotApplicable,false,false,false,0,3,
         "3-register low-storage Runge-Kutta","numerics.temporal.rk3",{"analytical_decay","temporal_refinement"}},
        {"linear.cg","Conjugate Gradient",F::LinearSolver,S::Verified,C::NotApplicable,false,false,false,0,0,
         "Krylov CG for SPD systems","linear.cg",{"SPD_exact_solve","true_residual"}},
        {"linear.bicgstab","BiCGStab",F::LinearSolver,S::Verified,C::NotApplicable,false,false,false,0,0,
         "BiCGStab Krylov iteration","linear.bicgstab",{"nonsymmetric_exact_solve","true_residual"}},
        {"linear.gmres","GMRES",F::LinearSolver,S::Verified,C::NotApplicable,false,false,false,0,0,
         "Restarted/unrestarted GMRES operator iteration","linear.gmres",{"nonsymmetric_exact_solve","true_residual"}},
        {"linear.fgmres","Flexible GMRES",F::LinearSolver,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Flexible Krylov iteration for varying preconditioners","linear.fgmres",{"nonsymmetric_exact_solve","true_residual","variable_preconditioner"}},
        {"preconditioner.native_amg","Native AMG",F::Preconditioner,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Smoothed aggregation / native AMG hierarchy for SPD elliptic systems","preconditioner.native_amg",{"SPD","anisotropy","true_residual","reuse"}},
        {"preconditioner.smoothed_aggregation_amg","Smoothed aggregation AMG",F::Preconditioner,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Smoothed aggregation AMG for SPD elliptic systems","preconditioner.smoothed_aggregation_amg",{"SPD","true_residual","reuse"}},
        {"preconditioner.native_fieldsplit","Native FieldSplit",F::Preconditioner,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Block field split for coupled pressure-velocity systems","preconditioner.native_fieldsplit",{"block_exact_small_system","Schur_variants"}},
        {"preconditioner.coupled_block_schur","Coupled block Schur",F::Preconditioner,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Velocity block plus pressure Schur approximation","preconditioner.coupled_block_schur",{"block_exact_small_system","true_residual","Schur_variants"}},
        {"preconditioner.pcd","PCD Schur",F::Preconditioner,S::Planned,C::NotApplicable,false,false,false,0,0,
         "Pressure-convection-diffusion Schur approximation","preconditioner.pcd",{"algebraic_reference","coupled_block_integration","pressure_gauge"}},
        {"schur.block_local","Block-local Schur",F::Schur,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Cell-local velocity-block Schur approximation","schur.block_local",{"coupled_block_integration","true_residual"}},
        {"schur.pcd","PCD Schur",F::Schur,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Pressure-convection-diffusion Schur inverse approximation","schur.pcd",{"algebraic_reference","coupled_block_integration","pressure_gauge"}},
        {"schur.lsc","LSC Schur",F::Schur,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Least-squares commutator Schur inverse approximation","schur.lsc",{"algebraic_reference","coupled_block_integration","true_residual"}},
        {"schur.bfbt","BFBT Schur",F::Schur,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "BFBT Schur inverse approximation","schur.bfbt",{"algebraic_reference","coupled_block_integration","true_residual"}},
        {"schur.simple","SIMPLE Schur",F::Schur,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Diagonal momentum SIMPLE Schur operator approximation","schur.simple",{"algebraic_reference","coupled_block_integration","action_contract"}},
        {"schur.simplec","SIMPLEC Schur",F::Schur,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "Consistent momentum SIMPLEC Schur operator approximation","schur.simplec",{"algebraic_reference","coupled_block_integration","action_contract"}},
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
         "block coupled momentum-pressure system","pressure_velocity.coupled",{"Couette","Poiseuille","cavity","skew_mesh"}},
        {"time_step.adaptive_cfl","Adaptive CFL controller",F::TimeStep,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "dt_new=dt*clip(CFL_target/CFL_current,min_factor,max_factor)","adaptive_cfl.adaptive_time_step",{"analytical_controller","growth","shrink","finite_inputs"}},
        {"time_step.pseudo_transient_cfl","Pseudo-transient CFL controller",F::TimeStep,S::Implemented,C::NotApplicable,false,false,false,0,0,
         "bounded CFL growth/saturation schedule","adaptive_cfl.pseudo_transient_cfl",{"growth","saturation","boundedness"}}
    };
}

inline void validate_numerical_method_registry()
{
    std::unordered_set<std::string> ids;
    for (const auto& method : numerical_method_registry()) {
        validate_numerical_method_contract(method);
        if (!ids.insert(method.id).second)
            throw std::invalid_argument("duplicate numerical method id: " + method.id);
    }
}

} // namespace cfdx::core
