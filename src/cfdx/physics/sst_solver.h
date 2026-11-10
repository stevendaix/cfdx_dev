#pragma once

#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/sst_jfnk_v22.h"
#include "cfdx/physics/sst_anderson_v20.h"
#include "cfdx/physics/sst_fp64_stall_v14.h"
#include "cfdx/physics/sst_fp64_stall_v14.h"
#include "cfdx/physics/scalar_transport_gradient.h"
#include <algorithm>
#include <functional>
#include <array>
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
    const ScalarBoundaryFaceValues* omega_face_values = nullptr,
    bool use_wls = false,
    int intrinsic_dimension = 3)
{
    // Merged contract: the default keeps the boundary-aware certified WLS
    // contraction; the wall-distance/gradient audit opt-in selects the
    // dimension-aware transport-boundaries WLS variant instead.
    const auto grad_k = use_wls
        ? scalar_gradient_with_transport_boundaries_wls(
              k, mesh, geometry, k_bcs, nullptr, intrinsic_dimension)
        : compute_boundary_aware_wls_gradient(k, mesh, geometry, k_bcs);
    const auto grad_omega = use_wls
        ? scalar_gradient_with_transport_boundaries_wls(
              omega, mesh, geometry, omega_bcs, omega_face_values,
              intrinsic_dimension)
        : compute_boundary_aware_wls_gradient(
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

// V16: the turbulent eddy viscosity vanishes at a smooth resolved no-slip
// wall, whereas its owner-cell value generally does not. The generic scalar
// assembler defaults to Gamma_owner at Dirichlet faces; keep that historical
// choice by default, but enable a physically distinct wall-face Gamma=rho*nu
// without changing any interior face, inlet, outlet or symmetry condition.
inline std::vector<double> sst_face_diffusion_coefficients(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const std::vector<double>& gamma_cell,
    const ScalarBoundaryConditions& boundary_conditions,
    const ScalarBoundaryFaceValues* face_values,
    const TurbulenceTransportControls& controls)
{
    if (gamma_cell.size() != mesh.n_cells() ||
        geometry.face_patch.size() != mesh.n_faces())
        throw std::invalid_argument("SST wall diffusion: field/mesh mismatch");
    const double wall_gamma = controls.density * controls.molecular_viscosity;
    if (!(wall_gamma > 0.0) || !std::isfinite(wall_gamma))
        throw std::invalid_argument("SST wall diffusion: invalid molecular coefficient");
    std::vector<double> gamma_face(mesh.n_faces());
    const auto& ownership = mesh.ownership();
    for (std::size_t face = 0; face < mesh.n_faces(); ++face) {
        const auto owner = ownership.owner(face);
        if (owner >= gamma_cell.size() || !std::isfinite(gamma_cell[owner]) ||
            gamma_cell[owner] < wall_gamma)
            throw std::invalid_argument("SST wall diffusion: invalid owner coefficient");
        gamma_face[face] = gamma_cell[owner];
        if (ownership.neighbour(face) >= 0) continue;
        const auto patch = geometry.face_patch[face];
        if (patch >= mesh.boundary().n_patches())
            throw std::invalid_argument("SST wall diffusion: unassigned boundary face");
        if (mesh.boundary().patch(patch).type != cfdx::core::PatchType::WALL)
            continue;
        const auto bc = resolve_scalar_boundary_condition(
            mesh,geometry,face,boundary_conditions,face_values);
        if (bc.type == ScalarBoundaryType::FIXED_VALUE)
            gamma_face[face] = wall_gamma;
    }
    return gamma_face;
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

// Diagnostic on EXACTLY the same sparse A, b and returned state. The old
// FP64 summation is retained for comparison. Long-double Neumaier summation
// reduces cancellation on platforms with extended precision; on MSVC where
// long double == double, compensation still reduces summation error. This is
// not an a-posteriori physical error bound, and must not be used to relax
// conservation or the configured tolerance.
struct SSTResidualPrecisionAudit {
    SSTLocalEquationDefect native;
    SSTLocalEquationDefect compensated;
    double maximum_row_disagreement = 0.0;
    std::size_t disagreement_cell = 0;
};

inline SSTResidualPrecisionAudit sst_local_residual_precision_audit(
    const ScalarEquation& equation,
    const cfdx::core::Vector& solution,
    double minimum_field_scale)
{
    SSTResidualPrecisionAudit audit;
    audit.native = sst_local_equation_defect(equation, solution, minimum_field_scale);
    const auto* offsets = equation.matrix.row_offsets_data();
    const auto* cols = equation.matrix.columns_data();
    const auto* vals = equation.matrix.values_data();
    if (offsets == nullptr || cols == nullptr || vals == nullptr)
        throw std::invalid_argument("SST precision audit: missing sparse matrix storage");
    for (std::size_t row = 0; row < solution.size(); ++row) {
        const double rhs = equation.rhs(row);
        if (!std::isfinite(rhs))
            throw std::invalid_argument("SST precision audit: nonfinite rhs");
        long double sum = -static_cast<long double>(rhs);
        long double correction = 0.0L;
        long double row_norm = 0.0L;
        long double local_scale = std::max(
            static_cast<long double>(minimum_field_scale),
            std::abs(static_cast<long double>(solution(row))));
        double native_residual = -rhs;
        double native_norm = 0.0;
        double native_scale = std::max(minimum_field_scale, std::abs(solution(row)));
        for (std::uint32_t entry = offsets[row]; entry < offsets[row+1]; ++entry) {
            const std::size_t col = cols[entry];
            if (col >= solution.size() || !std::isfinite(vals[entry]) ||
                !std::isfinite(solution(col)))
                throw std::invalid_argument("SST precision audit: invalid matrix row");
            const long double term = static_cast<long double>(vals[entry]) *
                                     static_cast<long double>(solution(col));
            const long double next = sum + term;
            if (std::abs(sum) >= std::abs(term))
                correction += (sum-next)+term;
            else
                correction += (term-next)+sum;
            sum = next;
            row_norm += std::abs(static_cast<long double>(vals[entry]));
            local_scale = std::max(local_scale,
                std::abs(static_cast<long double>(solution(col))));
            native_residual += vals[entry] * solution(col);
            native_norm += std::abs(vals[entry]);
            native_scale = std::max(native_scale, std::abs(solution(col)));
        }
        if (!(row_norm > 0.0L) || !(local_scale > 0.0L) ||
            !std::isfinite(sum+correction))
            throw std::invalid_argument("SST precision audit: invalid row norm or sum");
        const long double precise = std::abs(sum+correction) / (row_norm*local_scale);
        if (!std::isfinite(precise) || precise > std::numeric_limits<double>::max())
            throw std::domain_error("SST precision audit: residual overflow");
        const double compensated = static_cast<double>(precise);
        if (compensated > audit.compensated.maximum_relative) {
            audit.compensated.maximum_relative = compensated;
            audit.compensated.worst_cell = row;
        }
        const double native = std::abs(native_residual) / native_norm / native_scale;
        const double gap = std::abs(native-compensated);
        if (gap > audit.maximum_row_disagreement) {
            audit.maximum_row_disagreement = gap;
            audit.disagreement_cell = row;
        }
    }
    return audit;
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

// V19 numerical improvement: iterative refinement of the *unchanged*
// frozen SST A*x=b in solution units. A large Dirichlet omega wall RHS can
// dominate ||b||_2 and force a relative Krylov stop while a small interior
// row is still wrong. Solve D^{-1} A delta = D^{-1}(b-A*x) with D_ii =
// sum_j |A_ij|, then verify A*(x+delta)-b on ORIGINAL coefficients.
// This is an opt-in local linear correction, not a nonlinear or physical
// convergence certificate. Never change the accepted iterate on failure.
struct SSTEquilibratedDefectCorrectionV19 {
    bool attempted = false;
    bool accepted = false;
    double before = 0.0;
    double after = 0.0;
    std::size_t iterations = 0;
};

inline SSTEquilibratedDefectCorrectionV19 sst_equilibrated_defect_correction_v19(
    const ScalarEquation& original, cfdx::core::Vector& solution,
    double defect_tolerance)
{
    if (!std::isfinite(defect_tolerance) || !(defect_tolerance > 0.0) ||
        solution.size() != original.rhs.size() || solution.empty() ||
        original.diagonal.size() != solution.size() ||
        original.matrix.n_rows() != solution.size() ||
        original.matrix.n_cols() != solution.size())
        throw std::invalid_argument("V19 SST equilibrated correction: invalid equation");
    SSTEquilibratedDefectCorrectionV19 audit;
    audit.before = sst_equation_solution_defect(original,solution);
    audit.after = audit.before;
    if (audit.before <= defect_tolerance) return audit;
    audit.attempted=true;
    ScalarEquation scaled=original;
    const auto* offsets=original.matrix.row_offsets_data();
    const auto* columns=original.matrix.columns_data();
    const auto* values=original.matrix.values_data();
    auto* scaled_values=scaled.matrix.values_data();
    for (std::size_t row=0;row<solution.size();++row) {
        long double rhs=static_cast<long double>(original.rhs(row));
        long double row_norm=0.0L;
        for (auto j=offsets[row];j<offsets[row+1];++j) {
            const auto col=columns[j];
            if (col>=solution.size() || !std::isfinite(values[j]) ||
                !std::isfinite(solution(col)))
                throw std::invalid_argument("V19 SST correction: non-finite original row");
            rhs-=static_cast<long double>(values[j]) *
                 static_cast<long double>(solution(col));
            row_norm+=std::abs(static_cast<long double>(values[j]));
        }
        if (!(row_norm>0.0L) || !std::isfinite(rhs) || !std::isfinite(row_norm))
            throw std::invalid_argument("V19 SST correction: singular or non-finite row");
        const long double scaled_rhs=rhs/row_norm;
        if (!std::isfinite(scaled_rhs) ||
            std::abs(scaled_rhs)>std::numeric_limits<double>::max())
            throw std::domain_error("V19 SST correction: unrepresentable scaled RHS");
        scaled.rhs(row)=static_cast<double>(scaled_rhs);
        for (auto j=offsets[row];j<offsets[row+1];++j)
            scaled_values[j]=static_cast<double>(
                static_cast<long double>(values[j])/row_norm);
        scaled.diagonal[row]=original.diagonal[row]/static_cast<double>(row_norm);
    }
    cfdx::core::Vector delta(solution.size(),0.0);
    const ScalarSolveControls stricter{2000,
        std::min(1e-12,defect_tolerance*0.01),1.0};
    const auto linear=solve_scalar_equation(scaled,delta,stricter);
    audit.iterations=linear.iterations;
    if (linear.status!=cfdx::core::SolverStatus::CONVERGED) return audit;
    cfdx::core::Vector candidate=solution;
    bool changed=false;
    for (std::size_t i=0;i<solution.size();++i) {
        if (!std::isfinite(delta(i)) || !std::isfinite(candidate(i)+delta(i)))
            return audit;
        const double next=candidate(i)+delta(i);
        changed |= next!=candidate(i);
        candidate(i)=next;
    }
    if (!changed) return audit;
    const double after=sst_equation_solution_defect(original,candidate);
    // Decline even a converged auxiliary correction unless it proves a
    // substantial decrease on the original matrix and RHS.
    if (std::isfinite(after) && after<audit.before*0.9) {
        solution=std::move(candidate);
        audit.after=after;
        audit.accepted=true;
    }
    return audit;
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


}
    return report;
    production = std::move(staged_production);
    omega_faces = std::move(staged_faces);
    omega = std::move(staged_omega);
    }
        staged_production = std::move(settled_production);
        staged_faces = std::move(settled_faces);
            throw std::domain_error("V21 SST owner synchronization: wall target depends on projected omega");
            64.0*std::numeric_limits<double>::epsilon())
        if (report.relative_after >
        // an implicit wall Jacobian; fail closed instead of manufacturing one.
        // resynchronizations. A target that depends materially on omega needs
        // Do not hide a coupled wall fixed point with an arbitrary number of
        }
                gap/std::max({omega_floor,std::abs(staged_omega(row)),std::abs(target)}));
            report.relative_after = std::max(report.relative_after,
            report.maximum_after = std::max(report.maximum_after,gap);
            const double gap = std::abs(staged_omega(row)-target);
                throw std::domain_error("V21 SST owner synchronization: invalid refreshed target");
            if (!std::isfinite(target) || target < omega_floor)
            const double target = settled_production.omega_owner_by_cell[row];
            if (!owner_mask[row]) continue;
        for (std::size_t row = 0; row < omega.size(); ++row) {
            throw std::domain_error("V21 SST owner synchronization: unstable wall ownership");
            settled_production.has_omega_owner_cell != owner_mask)
        if (!settled_production.has_owner_omega_constraints() ||
        settled_production.validate(k.size());
        callback(k, staged_omega, settled_faces, settled_production);
        SSTWallProductionV17 settled_production;
        ScalarBoundaryFaceValues settled_faces;
        // must correspond to the *projected* omega, not the original trial.
        // A second evaluation is essential: the production source G and nu_t
        }
            ++report.owners;
            staged_omega(row) = value;
            }
                report.worst_cell = row;
                report.maximum_before = gap;
            if (gap >= report.maximum_before) {
            const double gap = std::abs(previous-value);
                throw std::domain_error("V21 SST owner synchronization: invalid wall value");
                !std::isfinite(value) || value < omega_floor)
            if (!std::isfinite(previous) || previous < omega_floor ||
            const double previous = staged_omega(row);
            const double value = staged_production.omega_owner_by_cell[row];
            if (!owner_mask[row]) continue;
        for (std::size_t row = 0; row < omega.size(); ++row) {
        const auto owner_mask = staged_production.has_omega_owner_cell;
    if (staged_production.has_owner_omega_constraints()) {
    SSTWallOwnerSyncV21 report;
    staged_production.validate(k.size());
    callback(k, staged_omega, staged_faces, staged_production);
    SSTWallProductionV17 staged_production;
    ScalarBoundaryFaceValues staged_faces;
    auto staged_omega = omega;
        throw std::invalid_argument("V21 SST owner synchronization: invalid inputs");
        !(omega_floor > 0.0) || !std::isfinite(omega_floor))
    if (!callback || k.size() == 0 || k.size() != omega.size() ||
{
    double omega_floor)
    SSTWallProductionV17& production,
    ScalarBoundaryFaceValues& omega_faces,
    const SSTWallRefreshV17& callback,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
inline SSTWallOwnerSyncV21 sst_synchronize_wall_owner_v21(

};
    double relative_after = 0.0;
    double maximum_after = 0.0;
    double maximum_before = 0.0;
    std::size_t worst_cell = 0;
    std::size_t owners = 0;
struct SSTWallOwnerSyncV21 {
// Publish no field or wall source until the complete consistency check passes.
// wall callback whose omega target changes when the owner value changes.
// project all owner cells on a STAGED field, then refresh again to detect a
// of any convergence criterion. Refresh k/omega/nu_t/G/face omega together,
// This optional projection is a constrained nonlinear trial, NOT a relaxation
//
// against the original refreshed constraint, never a frozen approximation.
// after k is updated. The returned-state residual MUST still be evaluated
// k, hence its omega solution is not generally on the constraint manifold
// omega_i = omega_wall_i(k,U,...). A frozen scalar solve uses the previous
// V21: on omegaWallFunction owner cells the PDE row is *replaced* by

}
    return audit;
    }
        ++audit.owners;
        }
            audit.worst_cell=i;
            audit.max_absolute_mismatch=delta;
        if (delta>audit.max_absolute_mismatch) {
            throw std::invalid_argument("V18 omega owner audit: non-finite state");
        if (!std::isfinite(delta))
        const double delta=std::abs(omega(i)-walls->omega_owner_by_cell[i]);
        if (!walls->has_omega_owner_cell[i]) continue;
    for (std::size_t i=0;i<omega.size();++i) {
    walls->validate(omega.size());
    if (!walls || !walls->has_owner_omega_constraints()) return audit;
    SSTWallOwnerConstraintAuditV18 audit;
{
    const SSTWallProductionV17* walls)
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
inline SSTWallOwnerConstraintAuditV18 sst_wall_owner_audit_v18(
};
    std::size_t worst_cell=0;
    double max_absolute_mismatch=0.0;
    std::size_t owners=0;
struct SSTWallOwnerConstraintAuditV18 {

}
    return constrained;
    equation=std::move(updated);
    }
        updated.diagonal[row]=1.0;
        updated.rhs(row)=walls->omega_owner_by_cell[row];
            values[j]=cols[j]==row ? 1.0 : 0.0;
        for (auto j=offsets[row];j<offsets[row+1];++j)
        if (!walls->has_omega_owner_cell[row]) continue;
    for (std::size_t row=0;row<n;++row) {
    auto* values=updated.matrix.values_data();
    auto updated=equation;
    }
        ++constrained;
            throw std::invalid_argument("V18 omega owner constraint: missing diagonal");
        if (!diagonal_found)
            diagonal_found |= cols[j]==row;
        for (auto j=offsets[row];j<offsets[row+1];++j)
        bool diagonal_found=false;
        if (!walls->has_omega_owner_cell[row]) continue;
    for (std::size_t row=0;row<n;++row) {
    std::size_t constrained=0;
    const auto* cols=equation.matrix.columns_data();
    const auto* offsets=equation.matrix.row_offsets_data();
    // not permission to manufacture a new graph without informing Krylov.
    // BEFORE modifying any coefficient. A missing diagonal is an error,
    // Transactional: check diagonal presence for every constrained row
        throw std::invalid_argument("V18 omega owner constraint: malformed equation");
        equation.diagonal.size()!=n)
    if (equation.matrix.n_rows()!=n || equation.matrix.n_cols()!=n ||
    walls->validate(n);
    const auto n=equation.rhs.size();
    if (!walls || !walls->has_owner_omega_constraints()) return 0;
{
    ScalarEquation& equation, const SSTWallProductionV17* walls)
inline std::size_t sst_apply_omega_owner_constraints_v18(
// different operator could incorrectly certify the returned SST state.
// Reassembly must call this same routine: otherwise a small residual of a
// CURRENT omega equation, retaining every other row and the original RHS.
// only Dirichlet values at wall faces. Apply an exact row constraint to the
// The OpenFOAM omegaWallFunction manipulates the owner-cell equations, not
    }
        return !has_omega_owner_cell.empty();
    bool has_owner_omega_constraints() const noexcept {
                throw std::invalid_argument("V18 omega owner constraint: invalid value");
                 !std::isfinite(omega_owner_by_cell[i])))
                (!has_wall_cell[i] || !(omega_owner_by_cell[i]>0.) ||
            if (has_omega_owner_cell[i] &&
        for (std::size_t i=0;owner_enabled && i<n;++i)
            throw std::invalid_argument("V18 omega owner constraint: invalid size");
            (omega_owner_by_cell.size()!=n || has_omega_owner_cell.size()!=n))
        if (owner_enabled &&
            !has_omega_owner_cell.empty();
        const bool owner_enabled = !omega_owner_by_cell.empty() ||
    std::vector<unsigned char> has_omega_owner_cell;
    std::vector<double> omega_owner_by_cell;
    // Empty vectors mean the V17 behaviour is preserved exactly.
    // Dirichlet value, this replaces the owner ROW of the omega equation.
    // V18 experimental owner-cell omega constraint. Unlike a boundary face
// V17: a wall-production constraint is expressed per unit carrier mass
// (m2/s3), not as the volume-integrated source rho*G. It is intentionally
// optional; the legacy SST operator stays unchanged when this is absent.
struct SSTWallProductionV17 {
    std::vector<double> g_by_cell;
    std::vector<unsigned char> has_wall_cell;

    void validate(std::size_t n) const {
        if (g_by_cell.size() != n || has_wall_cell.size() != n)
            throw std::invalid_argument("V17 SST wall G field size mismatch");
        for (std::size_t i=0; i<n; ++i)
            if (has_wall_cell[i] &&
                (!(g_by_cell[i]>=0.) || !std::isfinite(g_by_cell[i])))
                throw std::invalid_argument("V17 SST wall G is nonfinite or negative");
    }
    double for_cell(std::size_t i) const noexcept {
        return has_wall_cell[i] ? g_by_cell[i] : -1.0;
    }
};
using SSTWallRefreshV17 = std::function<void(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>&,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>&,
    ScalarBoundaryFaceValues&, SSTWallProductionV17&)>;

struct SSTCellSources {
    double sk=0.0, spk=0.0, sw=0.0, spw=0.0;
    double gamma_k=0.0, gamma_w=0.0, nut=0.0;
};

inline SSTCellSources sst_cell_sources(
    double k, double omega, double strain, double F1, double F2,
    double grad_k_dot_grad_omega, const TurbulenceTransportControls& c,
    double wall_g_kinematic = -1.0)
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
    if (!std::isfinite(wall_g_kinematic))
        throw std::invalid_argument("SST wall G is nonfinite");
    const double production_per_mass = wall_g_kinematic >= 0.0
        ? wall_g_kinematic : s.nut*S*S;
    const double Pk=std::min(c.density*production_per_mass,
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
    const ScalarBoundaryFaceValues* omega_face_values = nullptr,
    const SSTWallProductionV17* wall_production = nullptr)
{
    validate_turbulence_controls(controls);
    if (wall_production) wall_production->validate(mesh.n_cells());
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
        // Keep the cross-diffusion source and outer blending on the same
        // gradient scheme. Previously opt-in WLS only affected F1/F2 in the
        // outer loop while the inner omega equation still used Gauss.
        const auto cross = sst_boundary_gradient_dot(
            k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values,
            controls.sst_use_wls_boundary_gradients, controls.sst_gradient_dimension);
        auto oldk=k;
        auto oldw=omega;
        cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),sw(n,"Sw","W/m3",1);
        cfdx::core::Field<double,cfdx::core::Location::CELL> spk(n,"Spk","kg/m3/s",1),spw(n,"Spw","kg/m3/s",1);
        std::vector<double> gamma_k(n),gamma_w(n);
        for(std::size_t i=0;i<n;++i) {
            const auto q=sst_cell_sources(k(i),omega(i),strain_rate(i),F1(i),F2(i),cross[i],controls,
                wall_production ? wall_production->for_cell(i) : -1.0);
            sk(i)=q.sk; spk(i)=q.spk; sw(i)=q.sw; spw(i)=q.spw;
            gamma_k[i]=q.gamma_k; gamma_w[i]=q.gamma_w;
        }

        const auto k_wall_gamma = controls.sst_molecular_wall_diffusion
            ? sst_face_diffusion_coefficients(mesh,geometry,gamma_k,k_bcs,nullptr,controls)
            : std::vector<double>{};
        const auto omega_wall_gamma = controls.sst_molecular_wall_diffusion
            ? sst_face_diffusion_coefficients(mesh,geometry,gamma_w,omega_bcs,
                                              omega_face_values,controls)
            : std::vector<double>{};
        auto eqk=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
            nullptr,nullptr,nullptr,&gamma_k, ConvectionScheme::UPWIND, &k,
            nullptr,true,controls.sst_molecular_wall_diffusion ? &k_wall_gamma : nullptr);
        auto eqw=assemble_scalar_equation(
            mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
            omega_face_values,nullptr,nullptr,&gamma_w, ConvectionScheme::UPWIND, &omega,
            nullptr,true,controls.sst_molecular_wall_diffusion ? &omega_wall_gamma : nullptr);
        sst_apply_omega_owner_constraints_v18(eqw,wall_production);

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
            2000,sst_linear_tolerance_for_solution_defect(eqk,k_linear_target),controls.sst_linear_relaxation};
        const ScalarSolveControls omega_sc{
            2000,sst_linear_tolerance_for_solution_defect(eqw,omega_linear_target),controls.sst_linear_relaxation};
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
                k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values,
                controls.sst_use_wls_boundary_gradients, controls.sst_gradient_dimension);
            for (std::size_t i=0;i<n;++i) {
                const auto q=sst_cell_sources(k(i),omega(i),strain_rate(i),
                    F1(i),F2(i),updated_cross[i],controls,
                    wall_production ? wall_production->for_cell(i) : -1.0);
                sk(i)=q.sk;
                spk(i)=q.spk;
                gamma_k[i]=q.gamma_k;
            }
            const auto updated_k_wall_gamma = controls.sst_molecular_wall_diffusion
                ? sst_face_diffusion_coefficients(mesh,geometry,gamma_k,k_bcs,nullptr,controls)
                : std::vector<double>{};
            eqk=assemble_scalar_equation(
                mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
                nullptr,nullptr,nullptr,&gamma_k,ConvectionScheme::UPWIND,&k,
                nullptr,true,controls.sst_molecular_wall_diffusion ? &updated_k_wall_gamma : nullptr);
            const ScalarSolveControls updated_k_sc{
                2000,sst_linear_tolerance_for_solution_defect(eqk,k_linear_target),controls.sst_linear_relaxation};
            rk=solve_scalar_equation(eqk,k_solution,updated_k_sc);
        } else {
            rk=solve_scalar_equation(eqk,k_solution,k_sc);
            rw=solve_scalar_equation(eqw,omega_solution,omega_sc);
        }
        if (controls.sst_equilibrated_defect_correction_v19) {
            double omega_step=0.0;
            for (std::size_t i=0;i<n;++i)
                omega_step=std::max(omega_step,std::abs(omega_solution(i)-oldw(i)));
            if (omega_step<=1e-5) {
                const auto before=omega_solution;
                const auto correction=sst_equilibrated_defect_correction_v19(
                    eqw,omega_solution,tolerance);
                if (correction.accepted) {
                    for (std::size_t i=0;i<n;++i) {
                        if (!(omega_solution(i)>=controls.omega_min) ||
                            !std::isfinite(omega_solution(i))) {
                            omega_solution=before;
                            break;
                        }
                    }
                }
                if (std::getenv("CFDX_SST_DIAGNOSTICS") && correction.attempted)
                    std::cerr << "SST_EQUILIBRATED_DEFECT_V19"
                              << " before=" << correction.before
                              << " after=" << sst_equation_solution_defect(eqw,omega_solution)
                              << " accepted=" << (correction.accepted &&
                                   sst_equation_solution_defect(eqw,omega_solution)<correction.before)
                              << " iterations=" << correction.iterations << '\n';
            }
        }
        result.last_k_linear_status = rk.status;
        result.last_second_linear_status = rw.status;
        result.last_k_linear_before = rk.residual_relative_before_relaxation;
        result.last_k_linear_after = rk.residual_relative_after_relaxation;
        result.last_second_linear_before = rw.residual_relative_before_relaxation;
        result.last_second_linear_after = rw.residual_relative_after_relaxation;
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
        const auto k_local = controls.sst_compensated_residual_gate
            ? sst_local_residual_precision_audit(eqk,k_solution,controls.k_min).compensated
            : sst_local_equation_defect(eqk,k_solution,controls.k_min);
        const auto omega_local = controls.sst_compensated_residual_gate
            ? sst_local_residual_precision_audit(eqw,omega_solution,controls.omega_min).compensated
            : sst_local_equation_defect(eqw,omega_solution,controls.omega_min);
        const bool local_equations_ok = !controls.sst_scale_aware_convergence ||
            (k_local.maximum_relative <= tolerance &&
             omega_local.maximum_relative <= tolerance);
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
    double k_local_compensated = std::numeric_limits<double>::quiet_NaN();
    double omega_local_compensated = std::numeric_limits<double>::quiet_NaN();
    double k_local_roundoff_gap = std::numeric_limits<double>::quiet_NaN();
    double omega_local_roundoff_gap = std::numeric_limits<double>::quiet_NaN();
};

// V22: a single FV assembly is the source of BOTH the scalar convergence
// certificate and the matrix-free coupled SST nonlinear residual. No alternate
// source formula, wall closure or gradient operator is introduced here.
inline std::pair<ScalarEquation,ScalarEquation> sst_reassembled_equations_v22(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs,
    const ScalarBoundaryConditions& omega_bcs,
    const ScalarBoundaryFaceValues* omega_face_values,
    const SSTWallProductionV17* wall_production = nullptr)
{
    const std::size_t n = mesh.n_cells();
    if (wall_production) wall_production->validate(n);
    const auto cross = sst_boundary_gradient_dot(
        k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values,
        controls.sst_use_wls_boundary_gradients, controls.sst_gradient_dimension);
    cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),
        sw(n,"Sw","W/m3",1),spk(n,"Spk","kg/m3/s",1),
        spw(n,"Spw","kg/m3/s",1);
    std::vector<double> gamma_k(n),gamma_w(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto [f1,f2] = compute_sst_blending(
            k(i),omega(i),wall_distance(i),controls.molecular_viscosity,
            controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
        const auto q=sst_cell_sources(
            k(i),omega(i),strain_rate(i),f1,f2,cross[i],controls,
            wall_production ? wall_production->for_cell(i) : -1.0);
        sk(i)=q.sk; sw(i)=q.sw; spk(i)=q.spk; spw(i)=q.spw;
        gamma_k[i]=q.gamma_k; gamma_w[i]=q.gamma_w;
    }
    const auto k_wall_gamma = controls.sst_molecular_wall_diffusion
        ? sst_face_diffusion_coefficients(mesh,geometry,gamma_k,k_bcs,nullptr,controls)
        : std::vector<double>{};
    const auto omega_wall_gamma = controls.sst_molecular_wall_diffusion
        ? sst_face_diffusion_coefficients(mesh,geometry,gamma_w,omega_bcs,
                                          omega_face_values,controls)
        : std::vector<double>{};
    auto eqk=assemble_scalar_equation(
        mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
        nullptr,nullptr,nullptr,&gamma_k,ConvectionScheme::UPWIND,&k,
        nullptr,true,controls.sst_molecular_wall_diffusion ? &k_wall_gamma : nullptr);
    auto eqw=assemble_scalar_equation(
        mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
        omega_face_values,nullptr,nullptr,&gamma_w,ConvectionScheme::UPWIND,&omega,
        nullptr,true,controls.sst_molecular_wall_diffusion ? &omega_wall_gamma : nullptr);
    sst_apply_omega_owner_constraints_v18(eqw,wall_production);
    return {std::move(eqk),std::move(eqw)};
}

// V22 signed, per-ROW nonlinear residual. Unlike a scalar max certificate this
// retains the full vector [Rk,Romega] needed by J*v. The fixed nonlinear
// state is assembled afresh at each trial, using the production FV equations.
struct SSTCoupledResidualV22 {
    std::vector<double> values;
    std::vector<double> row_norms;
    std::vector<double> local_scales;
    std::vector<double> diagonal;
};

inline SSTCoupledResidualV22 sst_coupled_residual_v22(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs,
    const ScalarBoundaryConditions& omega_bcs,
    const ScalarBoundaryFaceValues* omega_face_values,
    const SSTWallProductionV17* wall_production = nullptr)
{
    const auto eq = sst_reassembled_equations_v22(
        mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,controls,
        k_bcs,omega_bcs,omega_face_values,wall_production);
    const std::size_t n=mesh.n_cells();
    SSTCoupledResidualV22 out;
    out.values.resize(2*n);
    out.row_norms.resize(2*n);
    out.local_scales.resize(2*n);
    out.diagonal.resize(2*n);
    for (std::size_t part=0;part<2;++part) {
        const auto& equation=part==0?eq.first:eq.second;
        const auto& field=part==0?k:omega;
        const double floor=part==0?controls.k_min:controls.omega_min;
        const auto* offsets=equation.matrix.row_offsets_data();
        const auto* cols=equation.matrix.columns_data();
        const auto* vals=equation.matrix.values_data();
        for (std::size_t i=0;i<n;++i) {
            long double signed_sum=-static_cast<long double>(equation.rhs(i));
            double row_norm=0.0,diag=0.0,scale=std::max(floor,std::abs(field(i)));
            for (auto entry=offsets[i];entry<offsets[i+1];++entry) {
                const auto col=cols[entry];
                if (col>=n || !std::isfinite(vals[entry]) ||
                    !std::isfinite(field(col)))
                    throw std::domain_error("V22 SST vector residual: invalid FV row");
                signed_sum+=static_cast<long double>(vals[entry])*field(col);
                row_norm+=std::abs(vals[entry]);
                if (col==i) diag+=vals[entry];
                scale=std::max(scale,std::abs(field(col)));
            }
            const auto residual=static_cast<double>(signed_sum);
            if (!(row_norm>0.) || !std::isfinite(row_norm) ||
                !(scale>0.) || !std::isfinite(scale) || !std::isfinite(residual))
                throw std::domain_error("V22 SST vector residual: nonfinite physical defect");
            out.values[part*n+i]=residual;
            out.row_norms[part*n+i]=row_norm;
            out.local_scales[part*n+i]=scale;
            out.diagonal[part*n+i]=diag;
        }
    }
    return out;
}

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
    const ScalarBoundaryFaceValues* omega_face_values,
    const SSTWallProductionV17* wall_production = nullptr)
{
    const std::size_t n = mesh.n_cells();
    const auto equations=sst_reassembled_equations_v22(
        mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,controls,
        k_bcs,omega_bcs,omega_face_values,wall_production);
    const auto& eqk=equations.first;
    const auto& eqw=equations.second;
    cfdx::core::Vector ks(n), ws(n);
    for(std::size_t i=0;i<n;++i) { ks(i)=k(i); ws(i)=omega(i); }
    const auto k_audit=sst_local_residual_precision_audit(eqk,ks,controls.k_min);
    const auto omega_audit=sst_local_residual_precision_audit(eqw,ws,controls.omega_min);
    const auto& k_local = controls.sst_compensated_residual_gate
        ? k_audit.compensated : k_audit.native;
    const auto& omega_local = controls.sst_compensated_residual_gate
        ? omega_audit.compensated : omega_audit.native;
    return {
        scalar_equation_residual_inf(eqk,ks),
        scalar_equation_residual_inf(eqw,ws),
        sst_equation_solution_defect(eqk,ks),
        sst_equation_solution_defect(eqw,ws),
        k_local.maximum_relative,omega_local.maximum_relative,
        k_local.worst_cell,omega_local.worst_cell,
        k_audit.compensated.maximum_relative,omega_audit.compensated.maximum_relative,
        k_audit.maximum_row_disagreement,omega_audit.maximum_row_disagreement};
}

// V16 opt-in operator budget. Split the exact assembled finite-volume row
// into source-free advection and diffusion plus the evaluated nonlinear source
// and sink, using the production assembler for all three checks. The sum is
// compared with the full equation A*x-b, preventing unverified hand-derived
// flux signs or hidden deferred nonorthogonal corrections.
struct SSTCellDiscreteBudget {
    std::size_t cell = 0;
    double convection = 0.0;
    double diffusion = 0.0;
    double explicit_source = 0.0;
    double implicit_source = 0.0;
    double recomposed = 0.0;
    double full_residual = 0.0;
    double closure_error = 0.0;
    double relative_local = 0.0;
    // Physical source split evaluated on the CURRENT state, independent of
    // Newton / sign-split linearization (row-residual units).
    double physical_production = 0.0;
    double physical_destruction = 0.0;
    double physical_cross_diffusion = 0.0;
    double nonlinear_source_identity_error = 0.0;
    double continuity_divergence = 0.0;
    double wall_orthogonal_flux = 0.0;
    double wall_owner_diffusion_flux = 0.0;
    double wall_molecular_diffusion_flux = 0.0;
    std::size_t dirichlet_wall_faces = 0;
    bool owner_constraint_row = false;
    double owner_constraint_target = 0.0;
};

inline double sst_signed_equation_row_residual(
    const ScalarEquation& eq, const cfdx::core::Vector& x, std::size_t row)
{
    if (row >= eq.rhs.size() || x.size() != eq.rhs.size())
        throw std::invalid_argument("SST row balance: invalid row or solution size");
    long double sum = -static_cast<long double>(eq.rhs(row));
    const auto* offsets = eq.matrix.row_offsets_data();
    for (auto entry=offsets[row]; entry<offsets[row+1]; ++entry) {
        const auto col=eq.matrix.columns_data()[entry];
        if (col >= x.size())
            throw std::invalid_argument("SST row balance: invalid column");
        sum += static_cast<long double>(eq.matrix.values_data()[entry]) *
               static_cast<long double>(x(col));
    }
    const double value=static_cast<double>(sum);
    if (!std::isfinite(value))
        throw std::invalid_argument("SST row balance: nonfinite residual");
    return value;
}

inline SSTCellDiscreteBudget sst_cell_discrete_budget(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& field,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& explicit_source,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& implicit_source,
    const std::vector<double>& gamma,
    const ScalarBoundaryConditions& bcs,
    const ScalarBoundaryFaceValues* face_values,
    const TurbulenceTransportControls& controls,
    std::size_t row,
    double minimum_scale)
{
    const auto nc=mesh.n_cells();
    if (row>=nc || field.size()!=nc || explicit_source.size()!=nc ||
        implicit_source.size()!=nc || gamma.size()!=nc ||
        mass_flux.size()!=mesh.n_faces() ||
        !(minimum_scale>0.0) || !std::isfinite(minimum_scale))
        throw std::invalid_argument("SST row budget: invalid arguments");
    cfdx::core::Field<double,cfdx::core::Location::CELL> zero_source(
        nc,"zero_source","",1);
    zero_source.fill(0.0);
    cfdx::core::Field<double,cfdx::core::Location::FACE> zero_flux(
        mesh.n_faces(),"zero_flux","",1);
    zero_flux.fill(0.0);
    const auto gamma_faces=controls.sst_molecular_wall_diffusion
        ? sst_face_diffusion_coefficients(mesh,geometry,gamma,bcs,face_values,controls)
        : std::vector<double>{};
    const auto* boundary_gamma=controls.sst_molecular_wall_diffusion
        ? &gamma_faces : nullptr;
    const auto full=assemble_scalar_equation(
        mesh,geometry,mass_flux,0.0,explicit_source,implicit_source,bcs,true,
        face_values,nullptr,nullptr,&gamma,ConvectionScheme::UPWIND,&field,
        nullptr,true,boundary_gamma);
    const auto diff=assemble_scalar_equation(
        mesh,geometry,zero_flux,0.0,zero_source,zero_source,bcs,true,
        face_values,nullptr,nullptr,&gamma,ConvectionScheme::UPWIND,&field,
        nullptr,true,boundary_gamma);
    cfdx::core::Vector x(nc);
    for (std::size_t i=0;i<nc;++i) x(i)=field(i);
    SSTCellDiscreteBudget budget;
    budget.cell=row;
    // Reconstruct bounded-upwind advection directly face by face. A
    // convection-only matrix can have a zero diagonal at stagnant cells and
    // therefore cannot be passed through the production assembler's positive
    // diagonal contract.
    budget.convection=0.0;
    budget.diffusion=sst_signed_equation_row_residual(diff,x,row);
    budget.explicit_source=-explicit_source(row)*geometry.cell_volumes[row];
    budget.implicit_source=-implicit_source(row)*field(row)*geometry.cell_volumes[row];
    budget.full_residual=sst_signed_equation_row_residual(full,x,row);
    double row_norm=0.0;
    double stencil_scale=std::max(minimum_scale,std::abs(field(row)));
    const auto* offsets=full.matrix.row_offsets_data();
    for (auto entry=offsets[row];entry<offsets[row+1];++entry) {
        row_norm+=std::abs(full.matrix.values_data()[entry]);
        stencil_scale=std::max(stencil_scale,std::abs(x(full.matrix.columns_data()[entry])));
    }
    if (!(row_norm>0.0) || !std::isfinite(row_norm) || !(stencil_scale>0.0))
        throw std::invalid_argument("SST row budget: invalid row normalization");
    budget.relative_local=std::abs(budget.full_residual)/(row_norm*stencil_scale);
    const auto& ownership=mesh.ownership();
    for (std::size_t face=0;face<mesh.n_faces();++face) {
        const auto owner=ownership.owner(face);
        const auto neighbour=ownership.neighbour(face);
        if (owner!=row && (neighbour < 0 || static_cast<std::size_t>(neighbour)!=row))
            continue;
        if (neighbour>=0) {
            const auto other=static_cast<std::size_t>(neighbour);
            const double F=mass_flux(face);
            const double transported = F >= 0.0 ? field(owner) : field(other);
            const double sign=owner==row ? 1.0 : -1.0;
            budget.continuity_divergence+=sign*F;
            budget.convection+=sign*F*transported;
            continue;
        }
        const auto patch=geometry.face_patch[face];
        if (patch>=mesh.boundary().n_patches())
            throw std::invalid_argument("SST row budget: unassigned patch");
        if (mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
            continue;
        budget.continuity_divergence+=mass_flux(face);
        const auto bc=resolve_scalar_boundary_condition(mesh,geometry,face,bcs,face_values);
        // Fixed-gradient discretization uses its prescribed boundary state
        // for either face-flux sign; other types use upwind or zero-gradient.
        double psi_face=field(row);
        if (bc.type==ScalarBoundaryType::FIXED_VALUE && mass_flux(face)<0.0)
            psi_face=bc.value;
        else if (bc.type==ScalarBoundaryType::FIXED_GRADIENT)
            psi_face=field(row)+bc.gradient*boundary_normal_distance(geometry,face,row);
        budget.convection+=mass_flux(face)*psi_face;
        if (mesh.boundary().patch(patch).type!=cfdx::core::PatchType::WALL ||
            bc.type!=ScalarBoundaryType::FIXED_VALUE) continue;
        ++budget.dirichlet_wall_faces;
        const double area=geometry.face_area_vectors[face].mag();
        const double dn=boundary_normal_distance(geometry,face,row);
        const double owner_flux=gamma[row]*area/dn*(field(row)-bc.value);
        const double molecular_flux=controls.density*controls.molecular_viscosity*
            area/dn*(field(row)-bc.value);
        budget.wall_owner_diffusion_flux+=owner_flux;
        budget.wall_molecular_diffusion_flux+=molecular_flux;
        budget.wall_orthogonal_flux+=controls.sst_molecular_wall_diffusion
            ? molecular_flux : owner_flux;
    }
    budget.convection-=budget.continuity_divergence*field(row);
    budget.recomposed=budget.convection+budget.diffusion+
        budget.explicit_source+budget.implicit_source;
    budget.closure_error=budget.full_residual-budget.recomposed;
    if (!std::isfinite(budget.recomposed) || !std::isfinite(budget.full_residual))
        throw std::invalid_argument("SST row budget: non-finite result");
    return budget;
}

inline std::pair<SSTCellDiscreteBudget,SSTCellDiscreteBudget>
sst_reassembled_discrete_cell_budgets(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs,
    const ScalarBoundaryConditions& omega_bcs,
    const ScalarBoundaryFaceValues* omega_face_values,
    std::size_t row, const SSTWallProductionV17* wall_production = nullptr)
{
    if (wall_production) wall_production->validate(mesh.n_cells());
    const auto n=mesh.n_cells();
    if (row>=n) throw std::invalid_argument("SST row budget: invalid row");
    const auto cross=sst_boundary_gradient_dot(
        k,omega,mesh,geometry,k_bcs,omega_bcs,omega_face_values,
        controls.sst_use_wls_boundary_gradients, controls.sst_gradient_dimension);
    cfdx::core::Field<double,cfdx::core::Location::CELL>
        sk(n,"Sk","",1),spk(n,"Spk","",1),sw(n,"Sw","",1),spw(n,"Spw","",1);
    std::vector<double> gamma_k(n),gamma_w(n);
    for (std::size_t i=0;i<n;++i) {
        const auto [f1,f2]=compute_sst_blending(
            k(i),omega(i),wall_distance(i),controls.molecular_viscosity,
            controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
        const auto s=sst_cell_sources(
            k(i),omega(i),strain_rate(i),f1,f2,cross[i],controls,
            wall_production ? wall_production->for_cell(i) : -1.0);
        sk(i)=s.sk;spk(i)=s.spk;sw(i)=s.sw;spw(i)=s.spw;
        gamma_k[i]=s.gamma_k;gamma_w[i]=s.gamma_w;
    }
    auto k_budget=sst_cell_discrete_budget(mesh,geometry,mass_flux,k,sk,spk,gamma_k,
                                 k_bcs,nullptr,controls,row,controls.k_min);
    auto omega_budget=sst_cell_discrete_budget(mesh,geometry,mass_flux,omega,sw,spw,gamma_w,
                                 omega_bcs,omega_face_values,controls,row,controls.omega_min);
    }
        omega_budget.wall_molecular_diffusion_flux=0.;
        omega_budget.wall_owner_diffusion_flux=0.;
        omega_budget.wall_orthogonal_flux=0.;
        omega_budget.dirichlet_wall_faces=0;
                      std::abs(omega_budget.owner_constraint_target)});
            std::max({controls.omega_min,std::abs(omega(row)),
        omega_budget.relative_local=std::abs(omega_budget.full_residual)/
        omega_budget.closure_error=0.;
        omega_budget.recomposed=omega_budget.full_residual;
        omega_budget.implicit_source=omega(row);
        omega_budget.explicit_source=-omega_budget.owner_constraint_target;
        omega_budget.diffusion=0.;
        omega_budget.convection=0.;
        omega_budget.full_residual=omega(row)-omega_budget.owner_constraint_target;
        omega_budget.owner_constraint_target=wall_production->omega_owner_by_cell[row];
        omega_budget.owner_constraint_row=true;
        // enforced row; they remain separately diagnosable by the raw budget.
        // Never report the unmodified PDE sources as if they decomposed this
        // The transport equation is replaced by an owner value constraint.
        wall_production->has_omega_owner_cell[row]) {
    if (wall_production && wall_production->has_owner_omega_constraints() &&
    const auto [f1,f2]=compute_sst_blending(k(row),omega(row),wall_distance(row),
        controls.molecular_viscosity,controls.beta_star,cross[row],controls.density,
        controls.sst_sigma_w2);
    const auto source=sst_cell_sources(
        k(row),omega(row),strain_rate(row),f1,f2,cross[row],controls,
        wall_production ? wall_production->for_cell(row) : -1.0);
    const double volume=geometry.cell_volumes[row];
    const double blended_gamma=f1*controls.gamma1+(1.0-f1)*controls.gamma2;
    const double blended_beta=f1*controls.beta1+(1.0-f1)*controls.beta2;
    const double omega_production=source.nut>0.0 ? source.sk/source.nut : 0.0;
    const double cross_source=2.0*(1.0-f1)*controls.density*controls.sst_sigma_w2/
                              omega(row)*cross[row];
    k_budget.physical_production=-source.sk*volume;
    k_budget.physical_destruction=controls.density*controls.beta_star*
                                  omega(row)*k(row)*volume;
    k_budget.nonlinear_source_identity_error=
        k_budget.explicit_source+k_budget.implicit_source-
        (k_budget.physical_production+k_budget.physical_destruction);
    omega_budget.physical_production=-blended_gamma*omega_production*volume;
    omega_budget.physical_destruction=controls.density*blended_beta*
                                     omega(row)*omega(row)*volume;
    omega_budget.physical_cross_diffusion=-cross_source*volume;
    omega_budget.nonlinear_source_identity_error=
        omega_budget.explicit_source+omega_budget.implicit_source-
        (omega_budget.physical_production+omega_budget.physical_destruction+
         omega_budget.physical_cross_diffusion);
    if (omega_budget.owner_constraint_row) {
        // Replaced equation: no PDE source decomposition is claimed here.
        // `owner_constraint_row` tells consumers to read this as an
        // enforced boundary constraint, not a volume-source balance.
        omega_budget.physical_production=0.;
        omega_budget.physical_destruction=0.;
        omega_budget.physical_cross_diffusion=0.;
        omega_budget.nonlinear_source_identity_error=0.;
    }
    return {k_budget,omega_budget};
}

// Conservative merit for an OPTIONAL nonlinear Picard backtracking step.
// This is the fully reassembled dimensionless cell-local SST equation defect,
// not the tiny Krylov residual of the frozen linear system. No line-search
// acceptance can mark an inner SST or outer carrier step as converged.
inline double sst_nonlinear_merit(const SSTNonlinearResidualSample& sample)
{
    const double merit = std::max(sample.k_local_relative, sample.omega_local_relative);
    if (!(merit >= 0.0) || !std::isfinite(merit))
        throw std::invalid_argument("SST line search: non-finite nonlinear merit");
    return merit;
}

inline bool sst_nonlinear_step_reduces_merit(double old_merit,
                                             double next_merit, double alpha)
{
    if (!std::isfinite(old_merit) || !std::isfinite(next_merit) ||
        old_merit < 0.0 || next_merit < 0.0 ||
        !(alpha > 0.0) || alpha > 1.0)
        return false;
    // A one-ulp fluctuation cannot be used to claim progress.
    return next_merit < old_merit * (1.0 - 1e-9*alpha);
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
    const ScalarBoundaryFaceValues* omega_face_values=nullptr,
    const SSTWallRefreshV17* wall_refresh = nullptr)
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
    std::size_t consecutive_line_search_stalls=0;
    };
        return sync;
        }
            refresh_wall(kk,ww);
        } else {
                controls.omega_min);
                kk,ww,*wall_refresh,active_omega_values,active_wall_g,
            sync=sst_synchronize_wall_owner_v21(
        if (controls.sst_owner_consistent_projection_v21) {
        SSTWallOwnerSyncV21 sync;
    const auto refresh_projected = [&](const auto& kk, auto& ww) {
            "V21 SST owner projection requires a dynamic wall refresh callback");
        throw std::invalid_argument(
    if (controls.sst_owner_consistent_projection_v21 && !wall_refresh)
    SSTAndersonV20 anderson_v20;
    SSTFp64StallV14 fp64_stall;
    // Refresh exactly one (wall omega, wall G) closure for every candidate
    // k/omega pair. No pre-SST frozen wall source may be passed off as the
    // residual of the returned state.
    ScalarBoundaryFaceValues active_omega_values;
    SSTWallProductionV17 active_wall_g;
    const auto refresh_wall = [&](const auto& kk, const auto& ww) {
        if (wall_refresh) {
            active_omega_values = ScalarBoundaryFaceValues{};
            active_wall_g = SSTWallProductionV17{};
            (*wall_refresh)(kk, ww, active_omega_values, active_wall_g);
            active_wall_g.validate(n);
        }
    };
    const auto face_values_now = [&]() -> const ScalarBoundaryFaceValues* {
        return wall_refresh ? &active_omega_values : omega_face_values;
    };
    const auto sources_now = [&]() -> const SSTWallProductionV17* {
        return wall_refresh ? &active_wall_g : nullptr;
    };
    SSTFp64StallV14 fp64_stall;
    for(std::size_t iter=1;iter<=max_iterations;++iter){
        auto old_k=k,old_w=omega;
        refresh_wall(k,omega);
        cfdx::core::Field<double,cfdx::core::Location::CELL> F1(n,"F1","1",1),F2(n,"F2","1",1);
        // Compute the full-cell gradient contraction once per outer iteration.
        // Recomputing this n-cell vector inside the n-cell loop made this O(n^2).
        const auto cross=sst_boundary_gradient_dot(
            k, omega, mesh, geometry, k_bcs,omega_bcs,face_values_now(),
            controls.sst_use_wls_boundary_gradients, controls.sst_gradient_dimension);
        for(std::size_t i=0;i<n;++i){
            const auto b=compute_sst_blending(k(i),omega(i),wall_distance(i),controls.molecular_viscosity,controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
            F1(i)=b.first; F2(i)=b.second;
        }
        const auto inner=solve_sst_transport(
            mesh,geometry,mass_flux,k,omega,strain_rate,F1,F2,controls,
            raw_wall_sync=refresh_projected(k,omega);
        if (controls.sst_owner_consistent_projection_v21)
        SSTWallOwnerSyncV21 raw_wall_sync;
        // Re-evaluation later uses the same refreshed omega/G operator.
        // line search. The old state remains untouched as the merit baseline.
        // Assemble a wall-consistent *raw Picard proposal* for Anderson and
            k_bcs,omega_bcs,1,tolerance,face_values_now(),sources_now());
        // Optional Armijo-like backtracking on the *returned nonlinear SST
        // equations*, with the pressure/velocity/face flux held fixed.
        // It never changes the global conservation or convergence criteria.
        double line_alpha = 1.0;
        bool line_accepted = true;
        double line_old_merit = std::numeric_limits<double>::quiet_NaN();
        double line_new_merit = std::numeric_limits<double>::quiet_NaN();
        } else if (controls.sst_nonlinear_line_search && !inner.converged) {
                          << " qualified=0\n";
                          << " accepted=" << line_accepted
                          << " new_merit=" << line_new_merit
                          << " old_merit=" << line_old_merit
                          << " alpha=" << line_alpha
                          << " theta=" << best_theta
                          << "SST_ANDERSON_V20 iteration=" << iter
                std::cerr << std::scientific << std::setprecision(10)
                (iter==1||iter%10==0||!line_accepted))
            if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
            refresh_wall(k,omega);
            }
                ++consecutive_line_search_stalls;
                line_alpha=0.;line_new_merit=line_old_merit;
                k=old_k;omega=old_w;
            } else {
                consecutive_line_search_stalls=0;
                line_alpha=best_alpha;line_new_merit=best_merit;
                k=std::move(best_k);omega=std::move(best_w);
            if (line_accepted) {
            line_accepted=std::isfinite(best_merit);
            }
                }
                    }
                        best_alpha=alpha;best_theta=choice==0?0.:accelerated.theta;
                        best_merit=merit;best_k=k;best_w=omega;
                        merit<best_merit) {
                    if (sst_nonlinear_step_reduces_merit(line_old_merit,merit,alpha)&&
                        controls,k_bcs,omega_bcs,face_values_now(),sources_now()));
                        mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                    const double merit=sst_nonlinear_merit(sst_reassembled_returned_state_residual(
                    refresh_projected(k,omega);
                    if (!admissible) continue;
                    }
                        k(i)=ck;omega(i)=cw;
                        }
                            admissible=false;break;
                            ck<controls.k_min||cw<controls.omega_min) {
                        if (!std::isfinite(ck)||!std::isfinite(cw)||
                        const double cw=old_w(i)+alpha*(tw[i]-old_w(i));
                        const double ck=old_k(i)+alpha*(tk[i]-old_k(i));
                    for (std::size_t i=0;i<n;++i) {
                    bool admissible=true;
                for (double alpha:{1.,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125}) {
                const auto& tw=choice==0?gw:accelerated.omega;
                const auto& tk=choice==0?gk:accelerated.k;
            for (int choice=0;choice<(accelerated.available?2:1);++choice) {
            double best_alpha=0.,best_theta=0.;
            auto best_k=old_k,best_w=old_w;
            double best_merit=std::numeric_limits<double>::infinity();
                controls,k_bcs,omega_bcs,face_values_now(),sources_now()));
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
            line_old_merit=sst_nonlinear_merit(sst_reassembled_returned_state_residual(
            refresh_wall(k,omega);
            k=old_k;omega=old_w;
            anderson_v20.record(xk,xw,gk,gw);
                xk,xw,gk,gw,controls.k_min,controls.omega_min);
            const auto accelerated=anderson_v20.propose(
            }
                gk[i]=raw_k(i);gw[i]=raw_w(i);
                xk[i]=old_k(i);xw[i]=old_w(i);
            for (std::size_t i=0;i<n;++i) {
            std::vector<double> xk(n),xw(n),gk(n),gw(n);
            const auto raw_k=k,raw_w=omega;
            // never controls acceptance.
            // reassembled SST physical residual. The frozen Krylov residual
            // Both trial families use the *current* wall omega/G and a fully
            // Safeguard: the two-depth extrapolation is only a proposal.

        } else if (controls.sst_anderson_m2_v20 && !inner.converged) {
                          << " qualified=0\n";
                          << " reason=" << newton.reason
                          << " accepted=" << line_accepted
                              (best_choice==0 ? "picard" : "jfnk") : "none")
                          << " selected=" << (line_accepted ?
                          << " new_merit=" << line_new_merit
                          << " old_merit=" << line_old_merit
                          << " alpha=" << line_alpha
                          << " evaluations=" << newton.function_evaluations
                          << " linear_true_relative=" << newton.linear_true_relative
                          << " proposed=" << newton.available
                          << "SST_JFNK_V22 iteration=" << iter
                std::cerr << std::scientific << std::setprecision(10)
                (iter==1 || iter%10==0 || !line_accepted))
            if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
            }
                ++consecutive_line_search_stalls;
                line_alpha=0.;line_new_merit=line_old_merit;
                k=old_k;omega=old_w;
            } else {
                consecutive_line_search_stalls=0;
                line_alpha=best_alpha;line_new_merit=best_merit;
                k=std::move(best_k);omega=std::move(best_w);
            if (line_accepted) {
            line_accepted=std::isfinite(best_merit);
            }
                }
                    }
                        // nonlinear field is committed before acceptance.
                        // Reject this trial only. Neither wall state nor
                    } catch (const std::exception&) {
                        }
                            best_choice=choice;
                            best_alpha=alpha;
                            best_k=k;best_w=omega;
                            best_merit=merit;
                                line_old_merit,merit,alpha) && merit<best_merit) {
                        if (sst_nonlinear_step_reduces_merit(
                                face_values_now(),sources_now()));
                                wall_distance,controls,k_bcs,omega_bcs,
                                mesh,geometry,mass_flux,k,omega,strain_rate,
                            sst_reassembled_returned_state_residual(
                        const double merit=sst_nonlinear_merit(
                        refresh_projected(k,omega);
                    try {
                    if (!feasible) continue;
                    }
                        k(i)=ck;omega(i)=cw;
                        }
                            feasible=false;break;
                            ck<controls.k_min||cw<controls.omega_min) {
                        if (!std::isfinite(ck)||!std::isfinite(cw)||
                        const double cw=start_w(i)+alpha*(target_w-start_w(i));
                        const double ck=start_k(i)+alpha*(target_k-start_k(i));
                        const double target_w=choice==0?raw_w(i):newton.fields[n+i];
                        const double target_k=choice==0?raw_k(i):newton.fields[i];
                    for (std::size_t i=0;i<n;++i) {
                    bool feasible=true;
                                  0.015625,0.0078125}) {
                for (double alpha:{1.,0.5,0.25,0.125,0.0625,0.03125,
            for (std::size_t choice=0;choice<choices;++choice) {
            const std::size_t choices=newton.available?2:1;
            std::size_t best_choice=0;
            double best_alpha=0.;
            auto best_k=old_k,best_w=old_w;
            double best_merit=std::numeric_limits<double>::infinity();
                full_residual,jctl,diagonal_preconditioner);
                physical,fixed_scales,n,controls.k_min,controls.omega_min,
            const auto newton=sst_jfnk_propose_v22(
            jctl.maximum_backtracks=7;
            jctl.finite_difference_step=1e-7;
            jctl.linear_relative_tolerance=0.25; // linear only, NOT SST convergence
            jctl.krylov_dimension=32;
            CoupledNewtonKrylovControls jctl;
            };
                    controls,k_bcs,omega_bcs,f,g).values;
                    mesh,geometry,mass_flux,kk,ww,strain_rate,wall_distance,
                return sst_coupled_residual_v22(
                }
                    f=&staging_faces;g=&staging_wall;
                    staging_wall.validate(n);
                    (*wall_refresh)(kk,ww,staging_faces,staging_wall);
                if (wall_refresh) {
                const SSTWallProductionV17* g=nullptr;
                const ScalarBoundaryFaceValues* f=omega_face_values;
                SSTWallProductionV17 staging_wall;
                ScalarBoundaryFaceValues staging_faces;
                }
                    kk(i)=state[i];ww(i)=state[n+i];
                for (std::size_t i=0;i<n;++i) {
                auto kk=start_k,ww=start_w;
            const SSTPhysicalResidualV22 full_residual=[&](const auto& state) {
            }
                diagonal_preconditioner[i]=std::clamp(frozen_log_diagonal,0.02,50.0);
                    throw std::domain_error("V22 SST JFNK: invalid block diagonal");
                if (!std::isfinite(frozen_log_diagonal))
                    base_full.diagonal[i]*physical[i]/fixed_scales[i]);
                const double frozen_log_diagonal=std::abs(
                // It does not approximate away k/omega derivatives in J*v.
                // preconditioner for the log-variable Newton system.
                // Frozen scalar-transport Jacobi is only a LEFT
                    throw std::domain_error("V22 SST JFNK: unusable fixed FV scale");
                if (!(fixed_scales[i]>0.) || !std::isfinite(fixed_scales[i]))
                fixed_scales[i]=base_full.row_norms[i]*base_full.local_scales[i];
            for (std::size_t i=0;i<2*n;++i) {
            }
                physical[n+i]=start_w(i);
                physical[i]=start_k(i);
            for (std::size_t i=0;i<n;++i) {
                                diagonal_preconditioner(2*n);
            std::vector<double> physical(2*n),fixed_scales(2*n),
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
            const auto base_full=sst_coupled_residual_v22(
            line_old_merit=sst_nonlinear_merit(base_sample);
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
            const auto base_sample=sst_reassembled_returned_state_residual(
            const auto start_k=k,start_w=omega;
            refresh_projected(k,omega);
            k=old_k;omega=old_w;
            const auto raw_k=k,raw_w=omega;
            // candidate with all wall face values and G refreshed.
            // system: the final nonlinear merit is always recomputed on the
            // SST certificate. Fixed row scales only condition the linear
            // finite differences of the SAME FV assembly as the returned
            // V22: joint [Rk,Romega] JFNK. The Krylov action is built by
        if (controls.sst_jfnk_v22_experimental && !inner.converged) {
            auto trial_k=k, trial_w=omega;
            k=old_k; omega=old_w;
            refresh_wall(k,omega);
            refresh_wall(k,omega);
            const auto starting=sst_reassembled_returned_state_residual(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
            line_old_merit=sst_nonlinear_merit(starting);
            line_accepted=false;
            // Largest reducing step wins; decline must be present on the
            // newly assembled operator, not only on the linear predictor.
            for (const double alpha : {1.0,0.5,0.25,0.125,0.0625,0.03125}) {
                for (std::size_t i=0; i<n; ++i) {
                    k(i)=old_k(i)+alpha*(trial_k(i)-old_k(i));
                    omega(i)=old_w(i)+alpha*(trial_w(i)-old_w(i));
                }
                refresh_projected(k,omega);
                refresh_wall(k,omega);
                const auto check=sst_reassembled_returned_state_residual(
                    mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                    controls,k_bcs,omega_bcs,face_values_now(),sources_now());
                const double merit=sst_nonlinear_merit(check);
                if (sst_nonlinear_step_reduces_merit(line_old_merit,merit,alpha)) {
                    line_alpha=alpha;
                    line_new_merit=merit;
                    line_accepted=true;
                    break;
                }
            }
            if (!line_accepted) {
                k=old_k; omega=old_w;
                line_alpha=0.0;
                line_new_merit=line_old_merit;
                ++consecutive_line_search_stalls;
            } else {
                consecutive_line_search_stalls=0;
            }
            if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
                (iter==1 || iter%100==0 || iter==max_iterations || !line_accepted))
                std::cerr << std::scientific << std::setprecision(9)
                          << "SST_NONLINEAR_LINESEARCH iteration=" << iter
                          << " alpha=" << line_alpha
                          << " old_merit=" << line_old_merit
                          << " new_merit=" << line_new_merit
                          << " accepted=" << line_accepted
                          << " qualified=0\n";
        } else if (!inner.converged && controls.sst_picard_damping < 1.0) {
            // Historical fixed damping remains unchanged when not using
            // line search, and never certifies an undamped candidate.
            const double alpha=controls.sst_picard_damping;
            for(std::size_t i=0;i<n;++i) {
                k(i)=old_k(i)+alpha*(k(i)-old_k(i));
                omega(i)=old_w(i)+alpha*(omega(i)-old_w(i));
            }
        }
                      << " qualified=0\n";
                      << " line_accepted=" << line_accepted
                      << " worst_cell=" << raw_wall_sync.worst_cell
                      << " final_relative=" << final_wall_sync.relative_after
                      << " final_gap=" << final_wall_sync.maximum_after
                      << " trial_gap=" << final_wall_sync.maximum_before
                      << " initial_gap=" << raw_wall_sync.maximum_before
                      << " owners=" << final_wall_sync.owners
                      << "SST_OWNER_SYNC_V21 iteration=" << iter
            std::cerr << std::scientific << std::setprecision(12)
            (iter==1 || iter%100==0 || iter==max_iterations))
            controls.sst_owner_consistent_projection_v21 &&
        if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
            refresh_wall(k,omega);
        else
            final_wall_sync=refresh_projected(k,omega);
        if (controls.sst_owner_consistent_projection_v21 && line_accepted)
        SSTWallOwnerSyncV21 final_wall_sync;
        // again to document the exact state used by the residual certificate.
        // lies on the owner constraint manifold; refresh both G and omega
        // Never modify a rejected old state. Any accepted candidate already
        refresh_wall(k,omega);
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
        result.last_k_linear_status = inner.last_k_linear_status;
        result.last_second_linear_status = inner.last_second_linear_status;
        result.last_k_linear_before = inner.last_k_linear_before;
        result.last_k_linear_after = inner.last_k_linear_after;
        result.last_second_linear_before = inner.last_second_linear_before;
        result.last_second_linear_after = inner.last_second_linear_after;
        // Compare the residual of the CURRENT (possibly damped) state to the
        // residual of the undamped candidate that the inner solver returned.
        // Sampling is enabled when needed for physical acceptance or forensic
        // diagnostics; otherwise legacy non-diagnostic performance is intact.
        const bool must_audit_returned = inner.converged ||
            controls.sst_scale_aware_convergence ||
            controls.sst_picard_damping < 1.0 ||
            controls.sst_nonlinear_line_search ||
            controls.sst_jfnk_v22_experimental ||
            controls.sst_anderson_m2_v20 ||
            std::getenv("CFDX_SST_DIAGNOSTICS") != nullptr;
        SSTNonlinearResidualSample returned;
        if (must_audit_returned) {
            returned = sst_reassembled_returned_state_residual(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
            result.k_residual = returned.k_raw_inf;
            result.second_residual = returned.omega_raw_inf;
        }
        }
                      << " qualified=0\n";
                      << " worst_cell=" << audit.worst_cell
                      << " max_absolute_mismatch=" << audit.max_absolute_mismatch
                      << " owners=" << audit.owners
                      << "SST_WALL_OWNER_CONSTRAINT_V18 iteration=" << iter
            std::cerr << std::scientific << std::setprecision(12)
            const auto audit=sst_wall_owner_audit_v18(omega,sources_now());
            (iter==1 || iter%100==0 || iter==max_iterations)) {
            std::getenv("CFDX_SST_DIAGNOSTICS") &&
        if (sources_now() && sources_now()->has_owner_omega_constraints() &&
        // OpenFOAM or experimental BFS qualification certificate.
        // pair, including corner area weights. Never count this as an
        // unconstrained PDE. The target is refreshed on the SAME trial k/w
        // V18: distinguish the exact owner equation from residuals of the
        // V16: expensive on-demand physical row budget, with no effect on the
        // solver path or any convergence gate. Compare exact full A*x-b with
        // separate advection, diffusion, production and destruction terms.
        const char* row_budget_flag=std::getenv("CFDX_SST_ROW_BUDGET");
        if (row_budget_flag != nullptr && must_audit_returned &&
            (iter==1 || iter%100==0 || iter==max_iterations)) {
            char* end=nullptr;
            const unsigned long requested=std::strtoul(row_budget_flag,&end,10);
            if (end==row_budget_flag || *end!='\0' || requested>=n)
                throw std::invalid_argument("CFDX_SST_ROW_BUDGET requires a valid zero-based cell index");
            const auto budget=sst_reassembled_discrete_cell_budgets(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                // The budget must use exactly the face values used by this
                // nonlinear candidate; V17 dynamically refreshes omega wall
                // faces, so the static pre-SST pointer would silently audit a
                // DIFFERENT equation.
                controls,k_bcs,omega_bcs,face_values_now(),
                static_cast<std::size_t>(requested),sources_now());
            for (const auto& item : {std::pair<const char*,SSTCellDiscreteBudget>{"k",budget.first},
                                     {"omega",budget.second}}) {
                const auto& b=item.second;
                std::cerr << std::scientific << std::setprecision(12)
                    << "SST_CELL_BUDGET iteration=" << iter
                    << " field=" << item.first
                    << " cell=" << b.cell
                    << " convection=" << b.convection
                    << " diffusion=" << b.diffusion
                    << " production=" << b.explicit_source
                    << " destruction=" << b.implicit_source
                    << " physical_production=" << b.physical_production
                    << " physical_destruction=" << b.physical_destruction
                    << " physical_cross=" << b.physical_cross_diffusion
                    << " source_identity_error=" << b.nonlinear_source_identity_error
                    << " local_relative=" << b.relative_local
                    << " sum=" << b.recomposed
                    << " residual=" << b.full_residual
                    << " decomposition_error=" << b.closure_error
                    << " div_phi=" << b.continuity_divergence
                    << " wall_faces=" << b.dirichlet_wall_faces
                    << " wall_orth_flux=" << b.wall_orthogonal_flux
                    << " wall_owner_flux=" << b.wall_owner_diffusion_flux
                    << " wall_molecular_flux=" << b.wall_molecular_diffusion_flux
                    << " qualified=0\n";
            }
        }
        // V20 read-only row-forensic comparison: the same discrete production
        // operator and dynamically refreshed omega/G wall state are used for
        // every sample. Guarded for small meshes to prevent accidental O(n^2)
        // diagnostics on the 20k-cell BFS qualification mesh.
        if (std::getenv("CFDX_SST_ASSEMBLY_TRACE_V20") && n<=4096 &&
            must_audit_returned &&
            (iter==1 || iter==2 || iter%100==0 || iter==max_iterations)) {
            // V21: print the currently worst reassembled k/omega cells as
            // well as the historical 77/357/369/426 probes. This matters
            // after owner projection: the dominant defect can move into an
            // adjacent *unconstrained* transport row (e.g. BFS cell 363).
            // Diagnostic only, with no effect on assembly or acceptance.
            const std::array<std::size_t, 7> audit_rows{{
                77,357,363,369,426,returned.k_worst_cell,returned.omega_worst_cell}};
            for (std::size_t j=0; j<audit_rows.size(); ++j) {
                const std::size_t row=audit_rows[j];
                if (row>=n || std::find(audit_rows.begin(),
                   audit_rows.begin()+static_cast<std::ptrdiff_t>(j),row)
                   !=audit_rows.begin()+static_cast<std::ptrdiff_t>(j)) continue;
                const auto budgets=sst_reassembled_discrete_cell_budgets(
                    mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                    controls,k_bcs,omega_bcs,face_values_now(),row,sources_now());
                for (const auto& item : {
                        std::pair<const char*,SSTCellDiscreteBudget>{"k",budgets.first},
                        std::pair<const char*,SSTCellDiscreteBudget>{"omega",budgets.second}}) {
                    const auto& b=item.second;
                    std::cerr << std::scientific << std::setprecision(12)
                              << "SST_ASSEMBLY_V20 iteration=" << iter
                              << " field=" << item.first << " cell=" << row
                              << " convection=" << b.convection
                              << " diffusion=" << b.diffusion
                              << " production=" << b.physical_production
                              << " destruction=" << b.physical_destruction
                              << " cross=" << b.physical_cross_diffusion
                              << " balance=" << b.full_residual
                              << " closure_error=" << b.closure_error
                              << " source_identity=" << b.nonlinear_source_identity_error
                              << " constrained=" << b.owner_constraint_row
                              << " qualified=0\n";
                }
            }
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
                      << " linear_relaxation=" << controls.sst_linear_relaxation
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
                      << " compensated_gate=" << controls.sst_compensated_residual_gate
                      << " k_local_compensated=" << returned.k_local_compensated
                      << " omega_local_compensated=" << returned.omega_local_compensated
                      << " k_local_roundoff_gap=" << returned.k_local_roundoff_gap
                      << " omega_local_roundoff_gap=" << returned.omega_local_roundoff_gap
                      << " k_local_worst_cell=" << returned.k_worst_cell
                      << " omega_local_worst_cell=" << returned.omega_worst_cell
                      << " k_inner_converged=" << inner.converged
                      << " k_linear_status=" << cfdx::core::to_string(inner.last_k_linear_status)
                      << " omega_linear_status=" << cfdx::core::to_string(inner.last_second_linear_status)
                      << " k_linear_relative_before_relaxation=" << inner.last_k_linear_before
                      << " k_linear_relative_after_relaxation=" << inner.last_k_linear_after
                      << " omega_linear_relative_before_relaxation=" << inner.last_second_linear_before
                      << " omega_linear_relative_after_relaxation=" << inner.last_second_linear_after
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
        // With zero accepted step the old state is returned unchanged. A
        // continued repeat will reassemble the same operator indefinitely.
        // Exit early as UNCONVERGED; the carrier transaction will roll back.
        if ((controls.sst_nonlinear_line_search || controls.sst_anderson_m2_v20 ||
             controls.sst_jfnk_v22_experimental) &&
            consecutive_line_search_stalls >= 8 && !inner.converged) {
            if (std::getenv("CFDX_SST_DIAGNOSTICS"))
                std::cerr << "SST_LINESEARCH_STAGNATED iteration=" << iter
                          << " consecutive_rejections=" << consecutive_line_search_stalls
                          << " qualified=0\n";
            break;
        }
        // The raw k/omega residuals alone do not certify convergence. In a
        // stiff BFS, they can be < tolerance while the *returned* equations
        // or field change remain above the same tolerance. Expose each gate
        // without ever changing the acceptance predicate below.
        if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
            (iter == 1 || iter % 100 == 0 || iter == max_iterations ||
             (inner.converged && outer_update_ok && returned_equations_ok))) {
            std::cerr << std::scientific << std::setprecision(12)
                      << "SST_ACCEPTANCE_GATES iteration=" << iter
                      << " inner_converged=" << inner.converged
                      << " line_accepted=" << line_accepted
                      << " finite=" << (result.nonfinite_bound_events == 0)
                      << " outer_update_ok=" << outer_update_ok
                      << " returned_equations_ok=" << returned_equations_ok
                      << " dk_absolute=" << dk
                      << " domega_absolute=" << dw
                      << " k_returned_row_scaled=" << returned.k_row_scaled
                      << " omega_returned_row_scaled=" << returned.omega_row_scaled
                      << " tolerance=" << tolerance
                      << " qualified_bfs=0\n";
        }
        const bool all_sst_gates_ok = inner.converged && line_accepted &&
            result.nonfinite_bound_events == 0 && outer_update_ok &&
            returned_equations_ok;
        if (fp64_stall.observe(dk,dw,all_sst_gates_ok)) {
            // Fail CLOSED. At a stall the inner linearized solve can be
            // unconverged, so must_audit_returned may be false. In that case
            // the old log emitted "inf" placeholders instead of the actual
            // defects and worst cells. Always sample the RETURNED nonlinear
            // equations here; this is read-only and cannot promote a state.
            if (!must_audit_returned)
                returned = sst_reassembled_returned_state_residual(
                    mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                    controls,k_bcs,omega_bcs,face_values_now(),sources_now());
            std::cerr << "SST_FP64_PRECISION_STALL_V14 iteration=" << iter
                      << " identical_steps=" << fp64_stall.identical_unconverged
                      << " k_returned_row_scaled=" << returned.k_row_scaled
                      << " omega_returned_row_scaled=" << returned.omega_row_scaled
                      << " k_local_relative=" << returned.k_local_relative
                      << " omega_local_relative=" << returned.omega_local_relative
                      << " k_worst_cell=" << returned.k_worst_cell
                      << " omega_worst_cell=" << returned.omega_worst_cell
                      << " k_compensated=" << returned.k_local_compensated
                      << " omega_compensated=" << returned.omega_local_compensated
                      << " target=" << tolerance
                      << " accepted=0 qualified=0\n";
            break;
        }
        if (all_sst_gates_ok) {
        // The raw k/omega residuals alone do not certify convergence. In a
        // stiff BFS, they can be < tolerance while the *returned* equations
        // or field change remain above the same tolerance. Expose each gate
        // without ever changing the acceptance predicate below.
        if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
            (iter == 1 || iter % 100 == 0 || iter == max_iterations ||
             (inner.converged && outer_update_ok && returned_equations_ok))) {
            std::cerr << std::scientific << std::setprecision(12)
                      << "SST_ACCEPTANCE_GATES iteration=" << iter
                      << " inner_converged=" << inner.converged
                      << " line_accepted=" << line_accepted
                      << " finite=" << (result.nonfinite_bound_events == 0)
                      << " outer_update_ok=" << outer_update_ok
                      << " returned_equations_ok=" << returned_equations_ok
                      << " dk_absolute=" << dk
                      << " domega_absolute=" << dw
                      << " k_returned_row_scaled=" << returned.k_row_scaled
                      << " omega_returned_row_scaled=" << returned.omega_row_scaled
                      << " tolerance=" << tolerance
                      << " qualified_bfs=0\n";
        }
        const bool all_sst_gates_ok = inner.converged && line_accepted &&
            result.nonfinite_bound_events == 0 && outer_update_ok &&
            returned_equations_ok;
        if (fp64_stall.observe(dk,dw,all_sst_gates_ok)) {
            // Fail CLOSED: a repeated FP64 state with a nonzero discrete
            // equation defect cannot converge on subsequent identical steps.
            // This only saves work; the outer U,p,phi,SST transaction rolls back.
            std::cerr << "SST_FP64_PRECISION_STALL_V14 iteration=" << iter
                      << " identical_steps=" << fp64_stall.identical_unconverged
                      << " k_returned_row_scaled=" << returned.k_row_scaled
                      << " omega_returned_row_scaled=" << returned.omega_row_scaled
                      << " target=" << tolerance
                      << " accepted=0 qualified=0\n";
            break;
        }
        if (all_sst_gates_ok) {
            result.converged = true;
            if (std::getenv("CFDX_SST_DIAGNOSTICS"))
                std::cerr << "SST_CERTIFIED_INNER"
                          << " iteration=" << iter
                          << " k_nonlinear_raw=" << result.k_residual
                          << " omega_nonlinear_raw=" << result.second_residual
                          << " k_local_relative=" << returned.k_local_relative
                          << " omega_local_relative=" << returned.omega_local_relative
                          << " k_linear_status=" << cfdx::core::to_string(result.last_k_linear_status)
                          << " omega_linear_status=" << cfdx::core::to_string(result.last_second_linear_status)
                          << " physical_qualification=NOT_QUALIFIED"
                          << '\n';
            break;
        }
    }
    return result;
}

} // namespace cfdx::physics
