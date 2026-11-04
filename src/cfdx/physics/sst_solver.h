#pragma once

#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/scalar_transport_gradient.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cfdx::physics {

inline std::vector<double> sst_boundary_gradient_dot(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarBoundaryConditions& k_bcs,
    const ScalarBoundaryConditions& omega_bcs,
    const ScalarBoundaryFaceValues* omega_face_values = nullptr)
{
    const auto grad_k = compute_boundary_aware_wls_gradient(
        k, mesh, geometry, k_bcs);
    const auto grad_omega = compute_boundary_aware_wls_gradient(
        omega, mesh, geometry, omega_bcs, omega_face_values);
    std::vector<double> result(mesh.n_cells());
    for (std::size_t cell = 0; cell < result.size(); ++cell) {
        result[cell] = 0.0;
        for (std::size_t component = 0; component < 3; ++component)
            result[cell] += grad_k(cell, component) * grad_omega(cell, component);
        if (!std::isfinite(result[cell]))
            throw std::invalid_argument("SST boundary gradient contraction is non-finite");
    }
    return result;
}

// Diagnostic-only scale-aware update. k [m2/s2] and omega [1/s] cannot be
// compared to a single absolute Picard tolerance as dimensionless quantities.
// Use one reference scale per FIELD; global-relative change is informative but
// is NOT a substitute for a cellwise equation defect or conservation gate.
struct SSTFieldUpdateDiagnostic {
    double absolute_inf = 0.0;
    double relative_global_inf = 0.0;
    double reference_scale = 0.0;
};

inline SSTFieldUpdateDiagnostic sst_field_update_diagnostic(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& previous,
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& current,
    double minimum_physical_scale)
{
    if (previous.size() == 0 || previous.size() != current.size() ||
        !std::isfinite(minimum_physical_scale) || minimum_physical_scale <= 0.0)
        throw std::invalid_argument("SST update diagnostic: invalid field size or scale");
    SSTFieldUpdateDiagnostic diagnostic;
    diagnostic.reference_scale = minimum_physical_scale;
    for (std::size_t cell = 0; cell < current.size(); ++cell) {
        if (!std::isfinite(previous(cell)) || !std::isfinite(current(cell)))
            throw std::invalid_argument("SST update diagnostic: non-finite field");
        diagnostic.reference_scale = std::max(
            diagnostic.reference_scale,
            std::max(std::abs(previous(cell)), std::abs(current(cell))));
        diagnostic.absolute_inf = std::max(
            diagnostic.absolute_inf, std::abs(current(cell) - previous(cell)));
    }
    diagnostic.relative_global_inf =
        diagnostic.absolute_inf / diagnostic.reference_scale;
    return diagnostic;
}

// Residual in solution units. Scaling each row by sum_j |a_ij| makes this
// comparable with the existing absolute iterate-change tolerance. Unlike a
// RHS-relative residual, it remains meaningful for a homogeneous equation.
// The unscaled physical residual is still reported in the transport result.
inline double sst_equation_solution_defect(
    const ScalarEquation& equation, const cfdx::core::Vector& solution)
{
    if (solution.size() != equation.rhs.size())
        throw std::invalid_argument("SST residual: field size mismatch");
    double maximum = 0.0;
    for (std::size_t cell = 0; cell < solution.size(); ++cell) {
        double residual = -equation.rhs(cell), scale = 0.0;
        for (auto entry = equation.matrix.row_offsets_data()[cell];
             entry < equation.matrix.row_offsets_data()[cell + 1]; ++entry) {
            const double coefficient = equation.matrix.values_data()[entry];
            residual += coefficient * solution(equation.matrix.columns_data()[entry]);
            scale += std::abs(coefficient);
        }
        if (!std::isfinite(residual) || !std::isfinite(scale) || !(scale > 0.0))
            return std::numeric_limits<double>::infinity();
        maximum = std::max(maximum, std::abs(residual) / scale);
    }
    return maximum;
}

// Independent CELL-LOCAL residual certificate. A global max(k) or max(omega)
// can hide an O(1) defect in near-wall cells whose field is much smaller.
// For each row, measure the local solution-unit defect |A*x-b| / ||A_i||_1,
// then compare it with the largest *local stencil* field magnitude, bounded
// below by the explicitly specified field floor. This is NOT a substitute for
// mass/flux conservation or for the full nonlinear operator reassembly.
struct SSTLocalEquationDefect {
    double maximum_relative = 0.0;
    std::size_t worst_cell = 0;
};

inline SSTLocalEquationDefect sst_local_equation_defect(
    const ScalarEquation& equation,
    const cfdx::core::Vector& solution,
    double minimum_field_scale)
{
    if (solution.empty() || solution.size() != equation.rhs.size() ||
        !std::isfinite(minimum_field_scale) || minimum_field_scale <= 0.0)
        throw std::invalid_argument("SST local defect: invalid field or scale");
    SSTLocalEquationDefect result;
    const auto* offsets = equation.matrix.row_offsets_data();
    const auto* cols = equation.matrix.columns_data();
    const auto* vals = equation.matrix.values_data();
    for (std::size_t row = 0; row < solution.size(); ++row) {
        double residual = -equation.rhs(row);
        double matrix_row_norm = 0.0;
        double local_field_scale = std::max(minimum_field_scale,
                                             std::abs(solution(row)));
        for (std::uint32_t k=offsets[row]; k<offsets[row+1]; ++k) {
            if (cols[k] >= solution.size() || !std::isfinite(vals[k]) ||
                !std::isfinite(solution(cols[k])))
                throw std::invalid_argument("SST local defect: invalid row coefficient");
            residual += vals[k]*solution(cols[k]);
            matrix_row_norm += std::abs(vals[k]);
            local_field_scale = std::max(local_field_scale,
                                         std::abs(solution(cols[k])));
        }
        if (!(matrix_row_norm > 0.0) || !std::isfinite(matrix_row_norm) ||
            !std::isfinite(residual) || !std::isfinite(local_field_scale))
            throw std::invalid_argument("SST local defect: invalid assembled equation");
        const double relative = (std::abs(residual) / matrix_row_norm) /
                                local_field_scale;
        if (!std::isfinite(relative))
            throw std::invalid_argument("SST local defect: normalization overflow");
        if (relative > result.maximum_relative) {
            result.maximum_relative = relative;
            result.worst_cell = row;
        }
    }
    return result;
}

// Dimensionless SST convergence certificate: the same numeric threshold
// applies to the k and omega *relative* updates and equation defects. The
// equation defect remains an independent true residual of A*x-b, not a
// difference of iterates. Report the actual raw residuals separately.
struct SSTRelativeConvergenceCertificate {
    double k_update_relative = 0.0;
    double omega_update_relative = 0.0;
    double k_equation_relative = 0.0;
    double omega_equation_relative = 0.0;
    bool converged = false;
};

inline SSTRelativeConvergenceCertificate sst_relative_convergence_certificate(
    const SSTFieldUpdateDiagnostic& k_change,
    const SSTFieldUpdateDiagnostic& omega_change,
    double k_row_scaled_equation_defect,
    double omega_row_scaled_equation_defect,
    double tolerance)
{
    if (!(tolerance > 0.0) || !std::isfinite(tolerance) ||
        !(k_change.reference_scale > 0.) ||
        !(omega_change.reference_scale > 0.) ||
        !std::isfinite(k_row_scaled_equation_defect) ||
        !std::isfinite(omega_row_scaled_equation_defect) ||
        k_row_scaled_equation_defect < 0. ||
        omega_row_scaled_equation_defect < 0.)
        throw std::invalid_argument("SST relative convergence certificate: invalid input");
    SSTRelativeConvergenceCertificate result;
    result.k_update_relative = k_change.relative_global_inf;
    result.omega_update_relative = omega_change.relative_global_inf;
    result.k_equation_relative =
        k_row_scaled_equation_defect / k_change.reference_scale;
    result.omega_equation_relative =
        omega_row_scaled_equation_defect / omega_change.reference_scale;
    result.converged = std::isfinite(result.k_update_relative) &&
        std::isfinite(result.omega_update_relative) &&
        result.k_update_relative <= tolerance &&
        result.omega_update_relative <= tolerance &&
        result.k_equation_relative <= tolerance &&
        result.omega_equation_relative <= tolerance;
    return result;
}

// Translate the absolute solution-unit gate into a sufficient relative
// Krylov tolerance. Otherwise a warm start can satisfy ||r||2/||b||2 while
// failing the row-scaled SST gate, so every solve returns the same rejected
// iterate forever. This only tightens the requested linear tolerance; the
// external tolerance, iteration budgets and relaxation are unchanged.
inline double sst_linear_tolerance_for_solution_defect(
    const ScalarEquation& equation, double tolerance)
{
    if (!std::isfinite(tolerance) || tolerance<=0.0 || equation.rhs.size()==0)
        throw std::invalid_argument("SST linear tolerance: invalid requested tolerance or empty equation");
    double rhs_norm=0.0;
    double minimum_row_norm=std::numeric_limits<double>::infinity();
    for (std::size_t cell=0; cell<equation.rhs.size(); ++cell) {
        if (!std::isfinite(equation.rhs(cell)))
            throw std::invalid_argument("SST linear tolerance: non-finite right-hand side");
        rhs_norm=std::hypot(rhs_norm,equation.rhs(cell));
        double row_norm=0.0;
        for (auto entry=equation.matrix.row_offsets_data()[cell];
             entry<equation.matrix.row_offsets_data()[cell+1]; ++entry)
            row_norm+=std::abs(equation.matrix.values_data()[entry]);
        if (!std::isfinite(row_norm) || !(row_norm>0.0))
            throw std::invalid_argument("SST linear tolerance: zero or non-finite matrix row norm");
        minimum_row_norm=std::min(minimum_row_norm,row_norm);
    }
    if (!std::isfinite(rhs_norm))
        throw std::invalid_argument("SST linear tolerance: right-hand side norm overflow");
    // Native CG/BiCGStab use max(||b||2,1e-15). GMRES/dense fallbacks use
    // smaller RHS scales, so this is sufficient for those paths as well.
    const double ratio=std::min(1.0,minimum_row_norm/std::max(rhs_norm,1e-15));
    const double translated=tolerance*ratio;
    if (!std::isfinite(translated) || !(translated>0.0))
        throw std::domain_error("SST linear tolerance: requested absolute defect cannot be represented");
    // A global min-row-norm bound can demand a relative Krylov target far
    // below FP64 precision on stretched grids. Prevent impossible linear
    // oversolving. The independent physical SST solution-defect and iterate
    // convergence tests below remain REQUIRED before accepting a solution.
    const double fp64_floor = 64.0 * std::numeric_limits<double>::epsilon();
    return std::max(translated, fp64_floor);
}

inline std::pair<double,double> compute_sst_blending(
    double k, double omega, double wall_distance,
    double molecular_viscosity, double beta_star,
    double grad_k_dot_grad_omega, double density, double sigma_w2)
{
    if(!std::isfinite(k) || !std::isfinite(omega) || !std::isfinite(wall_distance) ||
       !std::isfinite(molecular_viscosity) || !std::isfinite(beta_star) ||
       !std::isfinite(grad_k_dot_grad_omega) ||
       !(k >= 0.0) || !(omega > 0.0) || !(wall_distance > 0.0) ||
       !(molecular_viscosity >= 0.0) || !(beta_star > 0.0) ||
       !(density > 0.0) || !(sigma_w2 > 0.0))
        throw std::invalid_argument("compute_sst_blending: invalid inputs");
    const double ki=k, wi=omega, y=wall_distance;
    const double cd_kw=std::max(2.0*density*sigma_w2/wi*grad_k_dot_grad_omega,1e-10);
    const double arg1=std::min(
        std::max(std::sqrt(ki)/(beta_star*wi*y),500.0*molecular_viscosity/(y*y*wi)),
        4.0*density*sigma_w2*ki/(cd_kw*y*y));
    const double arg2=std::max(2.0*std::sqrt(ki)/(beta_star*wi*y),
                                500.0*molecular_viscosity/(y*y*wi));
    return {std::clamp(std::tanh(std::pow(arg1,4.0)),0.0,1.0),
            std::clamp(std::tanh(arg2*arg2),0.0,1.0)};
}

inline std::pair<double,double> compute_sst_blending(
    double k, double omega, double wall_distance,
    double molecular_viscosity, double beta_star)
{
    return compute_sst_blending(k,omega,wall_distance,molecular_viscosity,beta_star,
                                0.0,1.0,0.856);
}

struct SSTCellSources {
    double sk=0.0, spk=0.0, sw=0.0, spw=0.0;
    double gamma_k=0.0, gamma_w=0.0, nut=0.0;
};

inline SSTCellSources sst_cell_sources(
    double k, double omega, double strain, double F1, double F2,
    double grad_k_dot_grad_omega, const TurbulenceTransportControls& c)
{
    // std::clamp(NaN, 0, 1) returns NaN, while std::max(finite, NaN)
    // returns the finite argument on common C++ implementations. Thus a
    // corrupted F2 or strain could silently select the omega branch and
    // manufacture a plausible eddy viscosity. Validate before any clamp.
    if (!std::isfinite(k) || !std::isfinite(omega) ||
        !std::isfinite(strain) || strain < 0.0 || !std::isfinite(F1) ||
        !std::isfinite(F2) || !std::isfinite(grad_k_dot_grad_omega) ||
        F1 < 0.0 || F1 > 1.0 || F2 < 0.0 || F2 > 1.0)
        throw std::invalid_argument("SST cell sources: invalid state or blending factors");
    const double ki=std::max(k,c.k_min), wi=std::max(omega,c.omega_min);
    const double S=std::max(strain,0.0), f1=F1, f2=F2;
    const auto blend=[f1](double a,double b){return f1*a+(1.0-f1)*b;};
    SSTCellSources s;
    s.nut=c.a1*ki/std::max(c.a1*wi,S*f2);
    const double Pk=std::min(c.density*s.nut*S*S,
                             c.sst_production_limiter*c.beta_star*c.density*ki*wi);
    s.sk=Pk;
    s.spk=-c.density*c.beta_star*wi;
    const double cross=2.0*(1.0-f1)*c.density*c.sst_sigma_w2/wi*grad_k_dot_grad_omega;
    // SST-2003 uses the limited k-production in the omega production term.
    // If nut is zero, the limited production is also zero; avoid a 0/0 path.
    const double omega_production = s.nut > 0.0 ? Pk / s.nut : 0.0;
    const double beta = blend(c.beta1,c.beta2);
    const double sink_jacobian = c.sst_omega_newton_sink ? 2.0 : 1.0;
    // Newton linearization of -rho*beta*omega^2 around omega_old:
    // source = +rho*beta*omega_old^2, Sp = -2*rho*beta*omega_old.
    // At omega=omega_old the net sink remains EXACTLY -rho*beta*omega_old^2.
    // The positive explicit compensation preserves the M-matrix sink sign.
    s.sw=blend(c.gamma1,c.gamma2)*omega_production+std::max(cross,0.0);
    if (c.sst_omega_newton_sink)
        s.sw += c.density*beta*wi*wi;
    s.spw=-c.density*beta*wi*sink_jacobian+std::min(cross,0.0)/wi;
    s.gamma_k=c.density*(c.molecular_viscosity+
        blend(c.sst_sigma_k1,c.sst_sigma_k2)*s.nut);
    s.gamma_w=c.density*(c.molecular_viscosity+
        blend(c.sst_sigma_w1,c.sst_sigma_w2)*s.nut);
    if (!std::isfinite(s.nut) || s.nut < 0.0 ||
        !std::isfinite(s.sk) || !std::isfinite(s.spk) ||
        !std::isfinite(s.sw) || !std::isfinite(s.spw) ||
        !std::isfinite(s.gamma_k) || !std::isfinite(s.gamma_w) ||
        s.gamma_k < 0.0 || s.gamma_w < 0.0)
        throw std::domain_error("SST cell sources: non-finite transport coefficient");
    return s;
}

inline TurbulenceTransportResult solve_sst_transport(
    const cfdx::core::Mesh& mesh,
    const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& F1,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& F2,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs = {},
    const ScalarBoundaryConditions& omega_bcs = {},
    std::size_t max_iterations = 100,
    double tolerance = 1e-8,
    const ScalarBoundaryFaceValues* omega_face_values = nullptr)
{
    validate_turbulence_controls(controls);
    if(controls.model!=TurbulenceModel::SST)
        throw std::invalid_argument("solve_sst_transport requires SST model");
    const std::size_t n=mesh.n_cells();
    if(F1.size()!=n || F2.size()!=n || strain_rate.size()!=n ||
       k.size()!=n || omega.size()!=n || mass_flux.size()!=mesh.n_faces() ||
       k.dimension()!=1 || omega.dimension()!=1 || strain_rate.dimension()!=1 ||
       F1.dimension()!=1 || F2.dimension()!=1 || mass_flux.dimension()!=1)
        throw std::invalid_argument("SST blending field size mismatch");
    if (max_iterations == 0 || !std::isfinite(tolerance) || tolerance <= 0.0)
        throw std::invalid_argument("SST transport requires valid iteration controls");
    if (omega_face_values)
        validate_transport_boundary_face_values(mesh, geometry, *omega_face_values, true);

    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter) {
        const auto cross = sst_boundary_gradient_dot(
            k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
        auto oldk=k;
        auto oldw=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),sw(n,"Sw","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(n,"Spk","kg/m3/s",1),spw(n,"Spw","kg/m3/s",1);
        std::vector<double> gamma_k(n),gamma_w(n);
        for(std::size_t i=0;i<n;++i) {
            const auto q=sst_cell_sources(k(i),omega(i),strain_rate(i),F1(i),F2(i),cross[i],controls);
            sk(i)=q.sk; spk(i)=q.spk; sw(i)=q.sw; spw(i)=q.spw;
            gamma_k[i]=q.gamma_k; gamma_w[i]=q.gamma_w;
        }

        auto eqk=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
            nullptr,nullptr,nullptr,&gamma_k, ConvectionScheme::UPWIND, &k);
        auto eqw=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
            omega_face_values,nullptr,nullptr,&gamma_w, ConvectionScheme::UPWIND, &omega);

        // Relative mode converts the requested dimensionless tolerance into
        // each field's physical solution units before the Krylov solve.
        const double k_scale = sst_field_update_diagnostic(
            oldk, oldk, controls.k_min).reference_scale;
        const double omega_scale = sst_field_update_diagnostic(
            oldw, oldw, controls.omega_min).reference_scale;
        const double k_linear_target = controls.sst_scale_aware_convergence
            ? tolerance * k_scale : tolerance;
        const double omega_linear_target = controls.sst_scale_aware_convergence
            ? tolerance * omega_scale : tolerance;
        const ScalarSolveControls k_sc{
            2000,sst_linear_tolerance_for_solution_defect(eqk,k_linear_target),0.7};
        const ScalarSolveControls omega_sc{
            2000,sst_linear_tolerance_for_solution_defect(eqw,omega_linear_target),0.7};
        cfdx::core::Vector k_solution(n, 0.0);
        cfdx::core::Vector omega_solution(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            k_solution(i) = k(i);
            omega_solution(i) = omega(i);
        }
        cfdx::core::SolverResult rk, rw;
        if (controls.sst_sequential_omega_first) {
            // Reuse the existing scalar transport assembly: no second
            // discretization of gradients, wall values or mass fluxes.
            rw = solve_scalar_equation(eqw, omega_solution, omega_sc);
            for (std::size_t cell=0; cell<n; ++cell)
                omega(cell)=omega_solution(cell);
            record_turbulence_bounds(result,
                enforce_turbulence_bounds(k,omega,controls));
            for (std::size_t cell=0; cell<n; ++cell)
                omega_solution(cell)=omega(cell);
            const auto updated_cross = sst_boundary_gradient_dot(
                k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
            for (std::size_t i=0;i<n;++i) {
                const auto q=sst_cell_sources(k(i),omega(i),strain_rate(i),
                    F1(i),F2(i),updated_cross[i],controls);
                sk(i)=q.sk;
                spk(i)=q.spk;
                gamma_k[i]=q.gamma_k;
            }
            eqk=assemble_scalar_equation(
                mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
                nullptr,nullptr,nullptr,&gamma_k,ConvectionScheme::UPWIND,&k);
            const ScalarSolveControls updated_k_sc{
                2000,sst_linear_tolerance_for_solution_defect(eqk,k_linear_target),0.7};
            rk=solve_scalar_equation(eqk,k_solution,updated_k_sc);
        } else {
            rk=solve_scalar_equation(eqk,k_solution,k_sc);
            rw=solve_scalar_equation(eqw,omega_solution,omega_sc);
        }
        for (std::size_t i = 0; i < n; ++i) {
            k(i) = k_solution(i);
            omega(i) = omega_solution(i);
        }
        record_turbulence_bounds(result, enforce_turbulence_bounds(k,omega,controls));

        double dk=0.0,dw=0.0;
        for(std::size_t i=0;i<n;++i) {
            dk=std::max(dk,std::abs(k(i)-oldk(i)));
            dw=std::max(dw,std::abs(omega(i)-oldw(i)));
            // Residuals must describe the state returned to the caller, after
            // both linear relaxation and turbulence positivity enforcement.
            k_solution(i)=k(i);
            omega_solution(i)=omega(i);
        }
        result.k_residual=scalar_equation_residual_inf(eqk,k_solution);
        result.second_residual=scalar_equation_residual_inf(eqw,omega_solution);
        result.iterations=iter;
        const double k_defect = sst_equation_solution_defect(eqk,k_solution);
        const double omega_defect = sst_equation_solution_defect(eqw,omega_solution);
        const auto relative_certificate = sst_relative_convergence_certificate(
            sst_field_update_diagnostic(oldk,k,controls.k_min),
            sst_field_update_diagnostic(oldw,omega,controls.omega_min),
            k_defect, omega_defect, tolerance);
        const bool local_equations_ok = !controls.sst_scale_aware_convergence ||
            (sst_local_equation_defect(eqk,k_solution,controls.k_min).maximum_relative <= tolerance &&
             sst_local_equation_defect(eqw,omega_solution,controls.omega_min).maximum_relative <= tolerance);
        const bool iterate_and_equations_ok = controls.sst_scale_aware_convergence
            ? (relative_certificate.converged && local_equations_ok)
            : std::max(dk,dw)<=tolerance &&
              k_defect<=tolerance && omega_defect<=tolerance;
        if (result.nonfinite_bound_events==0 && rk.status==cfdx::core::SolverStatus::CONVERGED &&
           rw.status==cfdx::core::SolverStatus::CONVERGED &&
           std::isfinite(result.k_residual) && std::isfinite(result.second_residual) &&
           iterate_and_equations_ok) {
            result.converged=true;
            break;
        }
    }
    return result;
}

// Reassemble the *nonlinear* SST operator for the state actually returned
// after Picard damping. The one-step linear residual is computed against the
// previous k/omega/F1/F2 and cannot certify a damped candidate. This is a
// measurement path using the SAME production assembly, not a new discretizer.
struct SSTNonlinearResidualSample {
    double k_raw_inf = std::numeric_limits<double>::infinity();
    double omega_raw_inf = std::numeric_limits<double>::infinity();
    double k_row_scaled = std::numeric_limits<double>::infinity();
    double omega_row_scaled = std::numeric_limits<double>::infinity();
    double k_local_relative = std::numeric_limits<double>::infinity();
    double omega_local_relative = std::numeric_limits<double>::infinity();
    std::size_t k_worst_cell = 0;
    std::size_t omega_worst_cell = 0;
};

inline SSTNonlinearResidualSample sst_reassembled_returned_state_residual(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs,
    const ScalarBoundaryConditions& omega_bcs,
    const ScalarBoundaryFaceValues* omega_face_values)
{
    const std::size_t n = mesh.n_cells();
    const auto cross = sst_boundary_gradient_dot(
        k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
    cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),
        sw(n,"Sw","W/m3",1),spk(n,"Spk","kg/m3/s",1),
        spw(n,"Spw","kg/m3/s",1);
    std::vector<double> gamma_k(n),gamma_w(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto [f1,f2] = compute_sst_blending(
            k(i),omega(i),wall_distance(i),controls.molecular_viscosity,
            controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
        const auto q=sst_cell_sources(
            k(i),omega(i),strain_rate(i),f1,f2,cross[i],controls);
        sk(i)=q.sk; sw(i)=q.sw; spk(i)=q.spk; spw(i)=q.spw;
        gamma_k[i]=q.gamma_k; gamma_w[i]=q.gamma_w;
    }
    const auto eqk=assemble_scalar_equation(
        mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
        nullptr,nullptr,nullptr,&gamma_k,ConvectionScheme::UPWIND,&k);
    const auto eqw=assemble_scalar_equation(
        mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
        omega_face_values,nullptr,nullptr,&gamma_w,ConvectionScheme::UPWIND,&omega);
    cfdx::core::Vector ks(n), ws(n);
    for(std::size_t i=0;i<n;++i) { ks(i)=k(i); ws(i)=omega(i); }
    const auto k_local=sst_local_equation_defect(eqk,ks,controls.k_min);
    const auto omega_local=sst_local_equation_defect(eqw,ws,controls.omega_min);
    return {
        scalar_equation_residual_inf(eqk,ks),
        scalar_equation_residual_inf(eqw,ws),
        sst_equation_solution_defect(eqk,ks),
        sst_equation_solution_defect(eqw,ws),
        k_local.maximum_relative,omega_local.maximum_relative,
        k_local.worst_cell,omega_local.worst_cell};
}

inline TurbulenceTransportResult solve_sst_transport_dynamic_blending(
    const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,const ScalarBoundaryConditions& k_bcs={},
    const ScalarBoundaryConditions& omega_bcs={},std::size_t max_iterations=100,double tolerance=1e-8,
    const ScalarBoundaryFaceValues* omega_face_values=nullptr)
{
    validate_turbulence_controls(controls);
    const std::size_t n=mesh.n_cells();
    if(n==0||wall_distance.size()!=n||k.size()!=n||omega.size()!=n||strain_rate.size()!=n||mass_flux.size()!=mesh.n_faces()||
       wall_distance.dimension()!=1||k.dimension()!=1||omega.dimension()!=1||strain_rate.dimension()!=1||mass_flux.dimension()!=1)
        throw std::invalid_argument("SST dynamic-blending field size mismatch");
    if (controls.model!=TurbulenceModel::SST || max_iterations==0 ||
        !std::isfinite(tolerance) || tolerance<=0.0)
        throw std::invalid_argument("SST dynamic blending requires SST and valid iteration controls");
    TurbulenceTransportResult result;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old_k=k,old_w=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> F1(n,"F1","1",1),F2(n,"F2","1",1);
        // Compute the full-cell gradient contraction once per outer iteration.
        // Recomputing this n-cell vector inside the n-cell loop made this O(n^2).
        const auto cross=sst_boundary_gradient_dot(
            k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
        for(std::size_t i=0;i<n;++i){
            const auto b=compute_sst_blending(k(i),omega(i),wall_distance(i),controls.molecular_viscosity,controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
            F1(i)=b.first; F2(i)=b.second;
        }
        const auto inner=solve_sst_transport(
            mesh,geometry,mass_flux,k,omega,strain_rate,F1,F2,controls,
            k_bcs,omega_bcs,1,tolerance,omega_face_values);
        // When the dynamic-blending Picard map oscillates, damp the update
        // without changing the physical equation, positivity bounds or the
        // exact final convergence gate. If the candidate is fully certified,
        // apply it unmodified; a damped intermediate can NEVER be certified
        // based on residuals of a different (undamped) field.
        if (!inner.converged && controls.sst_picard_damping < 1.0) {
            const double alpha=controls.sst_picard_damping;
            for(std::size_t i=0;i<n;++i) {
                k(i)=old_k(i)+alpha*(k(i)-old_k(i));
                omega(i)=old_w(i)+alpha*(omega(i)-old_w(i));
            }
        }
        const auto k_change = sst_field_update_diagnostic(
            old_k, k, controls.k_min);
        const auto omega_change = sst_field_update_diagnostic(
            old_w, omega, controls.omega_min);
        const double dk = k_change.absolute_inf;
        const double dw = omega_change.absolute_inf;
        accumulate_turbulence_diagnostics(result, inner);
        result.k_residual = inner.k_residual;
        result.second_residual = inner.second_residual;
        result.iterations = iter;
        // Compare the residual of the CURRENT (possibly damped) state to the
        // residual of the undamped candidate that the inner solver returned.
        // Sampling is enabled when needed for physical acceptance or forensic
        // diagnostics; otherwise legacy non-diagnostic performance is intact.
        const bool must_audit_returned = inner.converged ||
            controls.sst_scale_aware_convergence ||
            controls.sst_picard_damping < 1.0 ||
            std::getenv("CFDX_SST_DIAGNOSTICS") != nullptr;
        SSTNonlinearResidualSample returned;
        if (must_audit_returned) {
            returned = sst_reassembled_returned_state_residual(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,omega_face_values);
            result.k_residual = returned.k_raw_inf;
            result.second_residual = returned.omega_raw_inf;
        }
        // Opt-in diagnostics: dk/domega describe the DAMPED returned state,
        // whereas k_raw/omega_raw are the PRE-DAMPING candidate equation
        // residuals reported by the one-step scalar solver. They must not be
        // presented as a certificate of the returned damped iterate.
        if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
            (iter == 1 || iter % 100 == 0 || iter == max_iterations)) {
            std::cerr << std::scientific << std::setprecision(9)
                      << "SST_TRANSPORT_DIAGNOSTIC iteration=" << iter
                      << " picard_damping=" << controls.sst_picard_damping
                      << " dk_inf=" << dk << " domega_inf=" << dw
                      << " dk_relative_global=" << k_change.relative_global_inf
                      << " domega_relative_global=" << omega_change.relative_global_inf
                      << " k_reference_scale=" << k_change.reference_scale
                      << " omega_reference_scale=" << omega_change.reference_scale
                      << " k_raw=" << inner.k_residual
                      << " omega_raw=" << inner.second_residual
                      << " relative_criterion=" << controls.sst_scale_aware_convergence
                      << " k_returned_nonlinear_raw=" << result.k_residual
                      << " omega_returned_nonlinear_raw=" << result.second_residual
                      << " k_returned_row_scaled=" << returned.k_row_scaled
                      << " omega_returned_row_scaled=" << returned.omega_row_scaled
                      << " k_local_relative=" << returned.k_local_relative
                      << " omega_local_relative=" << returned.omega_local_relative
                      << " k_local_worst_cell=" << returned.k_worst_cell
                      << " omega_local_worst_cell=" << returned.omega_worst_cell
                      << " k_inner_converged=" << inner.converged
                      << " finite_events=" << result.nonfinite_bound_events
                      << " floor_k=" << result.k_floor_events
                      << " floor_omega=" << result.second_floor_events
                      << '\n';
        }
        const bool outer_update_ok = controls.sst_scale_aware_convergence
            ? (k_change.relative_global_inf <= tolerance &&
               omega_change.relative_global_inf <= tolerance)
            : std::max(dk, dw) <= tolerance;
        const bool returned_equations_ok = !must_audit_returned ||
            (controls.sst_scale_aware_convergence
                ? returned.k_row_scaled / k_change.reference_scale <= tolerance &&
                  returned.omega_row_scaled / omega_change.reference_scale <= tolerance &&
                  returned.k_local_relative <= tolerance &&
                  returned.omega_local_relative <= tolerance
                : returned.k_row_scaled <= tolerance &&
                  returned.omega_row_scaled <= tolerance);
        if (inner.converged && result.nonfinite_bound_events == 0 &&
            outer_update_ok && returned_equations_ok) {
            result.converged = true;
            break;
        }
    }
    return result;
}

} // namespace cfdx::physics
