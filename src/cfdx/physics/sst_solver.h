#pragma once

#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/sst_fp64_stall_v14.h"
#include "cfdx/physics/sst_anderson_v20.h"
#include "cfdx/physics/sst_jfnk_v22.h"
#include "cfdx/physics/sst_block_preconditioner_v23.h"
#include "cfdx/physics/sst_spatial_block_preconditioner_v24.h"
#include "cfdx/physics/sst_globalization_v25.h"
#include "cfdx/physics/sst_globalization_v26.h"
#include "cfdx/physics/sst_trust_subspace_v27.h"
#include "cfdx/physics/sst_adaptive_trust_v28.h"
#include "cfdx/physics/sst_pseudo_transient_v29.h"
#include "cfdx/physics/sst_hotspot_derivatives_v31.h"
#include "cfdx/physics/sst_row_trace_v33.h"
#include "cfdx/physics/sst_stall_forensics_v42.h"
#include "cfdx/physics/sst_two_cycle_v41.h"
#include "cfdx/physics/sst_branch_preconditioner_v34.h"
#include "cfdx/physics/sst_openfoam_weighted_owner_v29.h"
#include "cfdx/physics/scalar_transport_gradient.h"
#include <algorithm>
#include <functional>
#include <array>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <sstream>
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
    double rms_relative = 0.; // V26 diagnostic, does not replace maximum_relative
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
    double rss=0.;
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
        rss=std::hypot(rss,relative);
        if (relative > result.maximum_relative) {
            result.maximum_relative = relative;
            result.worst_cell = row;
        }
    }
    result.rms_relative=rss/std::sqrt(static_cast<double>(solution.size()));
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
    double compensated_rss=0.;
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
        compensated_rss=std::hypot(compensated_rss,compensated);
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
    audit.compensated.rms_relative=compensated_rss/
        std::sqrt(static_cast<double>(solution.size()));
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

// V17: a wall-production constraint is expressed per unit carrier mass
// (m2/s3), not as the volume-integrated source rho*G. It is intentionally
// optional; the legacy SST operator stays unchanged when this is absent.
struct SSTWallProductionV17 {
    std::vector<double> g_by_cell;
    std::vector<unsigned char> has_wall_cell;
    // V18 experimental owner-cell omega constraint. Unlike a boundary face
    // Dirichlet value, this replaces the owner ROW of the omega equation.
    // Empty vectors mean the V17 behaviour is preserved exactly.
    std::vector<double> omega_owner_by_cell;
    std::vector<unsigned char> has_omega_owner_cell;

    void validate(std::size_t n) const {
        if (g_by_cell.size() != n || has_wall_cell.size() != n)
            throw std::invalid_argument("V17 SST wall G field size mismatch");
        for (std::size_t i=0; i<n; ++i)
            if (has_wall_cell[i] &&
                (!(g_by_cell[i]>=0.) || !std::isfinite(g_by_cell[i])))
                throw std::invalid_argument("V17 SST wall G is nonfinite or negative");
        const bool owner_enabled = !omega_owner_by_cell.empty() ||
            !has_omega_owner_cell.empty();
        if (owner_enabled &&
            (omega_owner_by_cell.size()!=n || has_omega_owner_cell.size()!=n))
            throw std::invalid_argument("V18 omega owner constraint: invalid size");
        for (std::size_t i=0;owner_enabled && i<n;++i)
            if (has_omega_owner_cell[i] &&
                (!has_wall_cell[i] || !(omega_owner_by_cell[i]>0.) ||
                 !std::isfinite(omega_owner_by_cell[i])))
                throw std::invalid_argument("V18 omega owner constraint: invalid value");
    }
    double for_cell(std::size_t i) const noexcept {
        return has_wall_cell[i] ? g_by_cell[i] : -1.0;
    }
    bool has_owner_omega_constraints() const noexcept {
        return !has_omega_owner_cell.empty();
    }
};

// The OpenFOAM omegaWallFunction manipulates the owner-cell equations, not
// only Dirichlet values at wall faces. Apply an exact row constraint to the
// CURRENT omega equation, retaining every other row and the original RHS.
// Reassembly must call this same routine: otherwise a small residual of a
// different operator could incorrectly certify the returned SST state.
inline std::size_t sst_apply_omega_owner_constraints_v18(
    ScalarEquation& equation, const SSTWallProductionV17* walls)
{
    if (!walls || !walls->has_owner_omega_constraints()) return 0;
    const auto n=equation.rhs.size();
    walls->validate(n);
    if (equation.matrix.n_rows()!=n || equation.matrix.n_cols()!=n ||
        equation.diagonal.size()!=n)
        throw std::invalid_argument("V18 omega owner constraint: malformed equation");
    // Transactional: check diagonal presence for every constrained row
    // BEFORE modifying any coefficient. A missing diagonal is an error,
    // not permission to manufacture a new graph without informing Krylov.
    const auto* offsets=equation.matrix.row_offsets_data();
    const auto* cols=equation.matrix.columns_data();
    std::size_t constrained=0;
    for (std::size_t row=0;row<n;++row) {
        if (!walls->has_omega_owner_cell[row]) continue;
        bool diagonal_found=false;
        for (auto j=offsets[row];j<offsets[row+1];++j)
            diagonal_found |= cols[j]==row;
        if (!diagonal_found)
            throw std::invalid_argument("V18 omega owner constraint: missing diagonal");
        ++constrained;
    }
    auto updated=equation;
    auto* values=updated.matrix.values_data();
    for (std::size_t row=0;row<n;++row) {
        if (!walls->has_omega_owner_cell[row]) continue;
        for (auto j=offsets[row];j<offsets[row+1];++j)
            values[j]=cols[j]==row ? 1.0 : 0.0;
        updated.rhs(row)=walls->omega_owner_by_cell[row];
        updated.diagonal[row]=1.0;
    }
    equation=std::move(updated);
    return constrained;
}

// V29 opt-in: on conventional fully-covered wall patches, mimic OpenFOAM
// fvMatrix::setValues(weights=1): fix the owner AND eliminate neighboring
// column contributions. Incomplete/partial face coverage is intentionally
// rejected by this production adapter pending patch-level geometry wiring.
inline std::size_t sst_apply_omega_owner_constraints_v29(
    ScalarEquation& equation, const SSTWallProductionV17* walls)
{
    if (!walls || !walls->has_owner_omega_constraints()) return 0;
    const auto n=equation.rhs.size();
    walls->validate(n);
    const auto audit=sst_apply_weighted_owner_matrix_v29(
        equation,walls->has_omega_owner_cell,walls->omega_owner_by_cell,
        std::vector<double>(n,1.),true);
    return audit.full_rows;
}

struct SSTWallOwnerConstraintAuditV18 {
    std::size_t owners=0;
    double max_absolute_mismatch=0.0;
    std::size_t worst_cell=0;
};
inline SSTWallOwnerConstraintAuditV18 sst_wall_owner_audit_v18(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const SSTWallProductionV17* walls)
{
    SSTWallOwnerConstraintAuditV18 audit;
    if (!walls || !walls->has_owner_omega_constraints()) return audit;
    walls->validate(omega.size());
    for (std::size_t i=0;i<omega.size();++i) {
        if (!walls->has_omega_owner_cell[i]) continue;
        const double delta=std::abs(omega(i)-walls->omega_owner_by_cell[i]);
        if (!std::isfinite(delta))
            throw std::invalid_argument("V18 omega owner audit: non-finite state");
        if (delta>audit.max_absolute_mismatch) {
            audit.max_absolute_mismatch=delta;
            audit.worst_cell=i;
        }
        ++audit.owners;
    }
    return audit;
}
using SSTWallRefreshV17 = std::function<void(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>&,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>&,
    ScalarBoundaryFaceValues&, SSTWallProductionV17&)>;

// V21: on omegaWallFunction owner cells the PDE row is *replaced* by
// omega_i = omega_wall_i(k,U,...). A frozen scalar solve uses the previous
// k, hence its omega solution is not generally on the constraint manifold
// after k is updated. The returned-state residual MUST still be evaluated
// against the original refreshed constraint, never a frozen approximation.
//
// This optional projection is a constrained nonlinear trial, NOT a relaxation
// of any convergence criterion. Refresh k/omega/nu_t/G/face omega together,
// project all owner cells on a STAGED field, then refresh again to detect a
// wall callback whose omega target changes when the owner value changes.
// Publish no field or wall source until the complete consistency check passes.
struct SSTWallOwnerSyncV21 {
    std::size_t owners = 0;
    std::size_t worst_cell = 0;
    double maximum_before = 0.0;
    double maximum_after = 0.0;
    double relative_after = 0.0;
};

inline SSTWallOwnerSyncV21 sst_synchronize_wall_owner_v21(
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const SSTWallRefreshV17& callback,
    ScalarBoundaryFaceValues& omega_faces,
    SSTWallProductionV17& production,
    double omega_floor)
{
    if (!callback || k.size() == 0 || k.size() != omega.size() ||
        !(omega_floor > 0.0) || !std::isfinite(omega_floor))
        throw std::invalid_argument("V21 SST owner synchronization: invalid inputs");
    auto staged_omega = omega;
    ScalarBoundaryFaceValues staged_faces;
    SSTWallProductionV17 staged_production;
    callback(k, staged_omega, staged_faces, staged_production);
    staged_production.validate(k.size());
    SSTWallOwnerSyncV21 report;
    if (staged_production.has_owner_omega_constraints()) {
        const auto owner_mask = staged_production.has_omega_owner_cell;
        for (std::size_t row = 0; row < omega.size(); ++row) {
            if (!owner_mask[row]) continue;
            const double value = staged_production.omega_owner_by_cell[row];
            const double previous = staged_omega(row);
            if (!std::isfinite(previous) || previous < omega_floor ||
                !std::isfinite(value) || value < omega_floor)
                throw std::domain_error("V21 SST owner synchronization: invalid wall value");
            const double gap = std::abs(previous-value);
            if (gap >= report.maximum_before) {
                report.maximum_before = gap;
                report.worst_cell = row;
            }
            staged_omega(row) = value;
            ++report.owners;
        }
        // A second evaluation is essential: the production source G and nu_t
        // must correspond to the *projected* omega, not the original trial.
        ScalarBoundaryFaceValues settled_faces;
        SSTWallProductionV17 settled_production;
        callback(k, staged_omega, settled_faces, settled_production);
        settled_production.validate(k.size());
        if (!settled_production.has_owner_omega_constraints() ||
            settled_production.has_omega_owner_cell != owner_mask)
            throw std::domain_error("V21 SST owner synchronization: unstable wall ownership");
        for (std::size_t row = 0; row < omega.size(); ++row) {
            if (!owner_mask[row]) continue;
            const double target = settled_production.omega_owner_by_cell[row];
            if (!std::isfinite(target) || target < omega_floor)
                throw std::domain_error("V21 SST owner synchronization: invalid refreshed target");
            const double gap = std::abs(staged_omega(row)-target);
            report.maximum_after = std::max(report.maximum_after,gap);
            report.relative_after = std::max(report.relative_after,
                gap/std::max({omega_floor,std::abs(staged_omega(row)),std::abs(target)}));
        }
        // Do not hide a coupled wall fixed point with an arbitrary number of
        // resynchronizations. A target that depends materially on omega needs
        // an implicit wall Jacobian; fail closed instead of manufacturing one.
        if (report.relative_after >
            64.0*std::numeric_limits<double>::epsilon())
            throw std::domain_error("V21 SST owner synchronization: wall target depends on projected omega");
        staged_faces = std::move(settled_faces);
        staged_production = std::move(settled_production);
    }
    omega = std::move(staged_omega);
    omega_faces = std::move(staged_faces);
    production = std::move(staged_production);
    return report;
}

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
        if (controls.sst_weighted_wall_matrix_v29_experimental)
            sst_apply_omega_owner_constraints_v29(eqw,wall_production);
        else
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
    // Equal-weight, per-equation RMS of the same dynamic local FV rows.
    // These are SECONDARY diagnostics; the existing Linf gate is unchanged.
    double k_local_l2 = std::numeric_limits<double>::infinity();
    double omega_local_l2 = std::numeric_limits<double>::infinity();
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
    if (controls.sst_weighted_wall_matrix_v29_experimental)
        sst_apply_omega_owner_constraints_v29(eqw,wall_production);
    else
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
    const SSTWallProductionV17* wall_production = nullptr,
    std::pair<ScalarEquation,ScalarEquation>* frozen_v53_transport = nullptr)
{
    auto eq = sst_reassembled_equations_v22(
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
    if(frozen_v53_transport)
        *frozen_v53_transport=std::move(eq);
    return out;
}

// V23: construct a CHEAP frozen local 2x2 SST source Jacobian. The base
// FV matrix provides advection, diffusion and the frozen sink diagonal;
// centered LOCAL source differences add the missing k/omega cross terms.
// The matrix-free Jv still contains ALL gradients, F1/F2, wall refresh and
// cross-cell diffusion dependencies: these blocks are never the true J.
// G, grad(k).grad(omega), and wall targets are fixed for this approximate P.
inline SSTBlockPreconditionerV23 sst_local_source_blocks_v23(
    const cfdx::core::Mesh& mesh, const FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& strain_rate,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>& wall_distance,
    const TurbulenceTransportControls& controls,
    const ScalarBoundaryConditions& k_bcs,
    const ScalarBoundaryConditions& omega_bcs,
    const ScalarBoundaryFaceValues* omega_face_values,
    const SSTWallProductionV17* wall_production,
    const SSTCoupledResidualV22& base,
    const std::vector<double>& frozen_scales,
    const std::vector<double>& frozen_left)
{
    const std::size_t n=mesh.n_cells();
    if (!n || geometry.cell_volumes.size()!=n || k.size()!=n || omega.size()!=n ||
        strain_rate.size()!=n || wall_distance.size()!=n ||
        frozen_scales.size()!=2*n || frozen_left.size()!=2*n ||
        base.diagonal.size()!=2*n)
        throw std::invalid_argument("V23 SST blocks: inconsistent FV fields");
    if (wall_production) wall_production->validate(n);
    const auto cross=sst_boundary_gradient_dot(
        k,omega,mesh,geometry,k_bcs,omega_bcs,omega_face_values,
        controls.sst_use_wls_boundary_gradients,controls.sst_gradient_dimension);
    std::vector<SSTBlockPreconditionerV23::Block> blocks;
    blocks.reserve(n);
    for (std::size_t i=0;i<n;++i) {
        const double ki=k(i),wi=omega(i),V=geometry.cell_volumes[i];
        if (!(V>0.) || !(ki>controls.k_min) || !(wi>controls.omega_min) ||
            !std::isfinite(V))
            throw std::invalid_argument("V23 SST blocks: invalid local physical state");
        const auto q=[&](double kval,double wval) {
            const auto blend=compute_sst_blending(
                kval,wval,wall_distance(i),controls.molecular_viscosity,
                controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
            return sst_cell_sources(
                kval,wval,strain_rate(i),blend.first,blend.second,cross[i],
                controls,wall_production ? wall_production->for_cell(i) : -1.);
        };
        // The caller has already established that Newton is away from floors.
        // A scale-relative perturbation avoids differentiating dimensional
        // k and omega with the same absolute step.
        const double dk=std::min(1e-5*ki,0.49*(ki-controls.k_min));
        const double dw=std::min(1e-5*wi,0.49*(wi-controls.omega_min));
        if (!(dk>0.) || !(dw>0.) || !std::isfinite(dk) || !std::isfinite(dw))
            throw std::invalid_argument("V23 SST blocks: invalid source perturbation");
        const auto kp=q(ki+dk,wi),km=q(ki-dk,wi);
        const auto wp=q(ki,wi+dw),wm=q(ki,wi-dw);
        const auto partial=[&](double p,double m,double step) {
            return (p-m)/(2.0*step);
        };
        // R = A(state)*x-b(state), with A_ii already including -V*Sp.
        // Only derivatives of Sk and Sp beyond this frozen A are added.
        const double rkk=base.diagonal[i]-V*(
            partial(kp.sk,km.sk,dk)+ki*partial(kp.spk,km.spk,dk));
        const double rkw=-V*(partial(wp.sk,wm.sk,dw)+
                             ki*partial(wp.spk,wm.spk,dw));
        double rwk=-V*(partial(kp.sw,km.sw,dk)+
                       wi*partial(kp.spw,km.spw,dk));
        double rww=base.diagonal[n+i]-V*(
            partial(wp.sw,wm.sw,dw)+wi*partial(wp.spw,wm.spw,dw));
        // The replaced owner omega row is NOT a transport source equation.
        // Its unknown wall-target derivative is retained by the true Jv,
        // not fabricated from transport source derivatives in this P.
        if (wall_production && wall_production->has_owner_omega_constraints() &&
            wall_production->has_omega_owner_cell[i]) {
            rwk=0.;rww=base.diagonal[n+i];
        }
        const double rowk=frozen_scales[i]*frozen_left[i];
        const double roww=frozen_scales[n+i]*frozen_left[n+i];
        if (!(rowk>0.) || !(roww>0.) || !std::isfinite(rowk) ||
            !std::isfinite(roww))
            throw std::invalid_argument("V23 SST blocks: invalid fixed row scaling");
        blocks.push_back({rkk*ki/rowk,rkw*wi/rowk,
                          rwk*ki/roww,rww*wi/roww});
    }
    return SSTBlockPreconditionerV23(blocks);
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
        k_audit.maximum_row_disagreement,omega_audit.maximum_row_disagreement,
        k_local.rms_relative,omega_local.rms_relative};
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
    if (wall_production && wall_production->has_owner_omega_constraints() &&
        wall_production->has_omega_owner_cell[row]) {
        // The transport equation is replaced by an owner value constraint.
        // Never report the unmodified PDE sources as if they decomposed this
        // enforced row; they remain separately diagnosable by the raw budget.
        omega_budget.owner_constraint_row=true;
        omega_budget.owner_constraint_target=wall_production->omega_owner_by_cell[row];
        omega_budget.full_residual=omega(row)-omega_budget.owner_constraint_target;
        omega_budget.convection=0.;
        omega_budget.diffusion=0.;
        omega_budget.explicit_source=-omega_budget.owner_constraint_target;
        omega_budget.implicit_source=omega(row);
        omega_budget.recomposed=omega_budget.full_residual;
        omega_budget.closure_error=0.;
        omega_budget.relative_local=std::abs(omega_budget.full_residual)/
            std::max({controls.omega_min,std::abs(omega(row)),
                      std::abs(omega_budget.owner_constraint_target)});
        omega_budget.dirichlet_wall_faces=0;
        omega_budget.wall_orthogonal_flux=0.;
        omega_budget.wall_owner_diffusion_flux=0.;
        omega_budget.wall_molecular_diffusion_flux=0.;
    }
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

inline double sst_nonlinear_smooth_merit_v26(const SSTNonlinearResidualSample& sample)
{
    if (!(sample.k_local_l2>=0.) || !(sample.omega_local_l2>=0.) ||
        !std::isfinite(sample.k_local_l2) ||
        !std::isfinite(sample.omega_local_l2))
        throw std::invalid_argument("V26 SST smooth merit: invalid RMS certificate");
    return std::hypot(sample.k_local_l2,sample.omega_local_l2)/std::sqrt(2.);
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
    SSTFp64StallV14 fp64_stall;
    // Defaults to the exact V40 trajectory. Allocations occur only when the
    // user explicitly enables the V41 diagnostic/rescue policy.
    double picard_damping_v41=controls.sst_picard_damping;
    std::unique_ptr<SSTTwoCycleAuditV41> two_cycle_v41;
    std::vector<double> two_cycle_snapshot_k_v41, two_cycle_snapshot_omega_v41;
    if (controls.sst_two_cycle_rescue_v41_experimental) {
        two_cycle_v41=std::make_unique<SSTTwoCycleAuditV41>(n);
        // Reuse sampling buffers for all SST iterations. The experimental
        // observation adds no per-iteration O(N) allocations or default-on
        // O(N) storage cost.
        two_cycle_snapshot_k_v41.resize(n);
        two_cycle_snapshot_omega_v41.resize(n);
    }
    SSTAndersonV20 anderson_v20;
    SSTV28AdaptiveTrustState adaptive_trust_v28;
    SSTV29PseudoTimeState pseudo_time_v29;
    SSTV30ContinuationState ptc_v30 = controls.sst_continuation_state_v30
        ? *controls.sst_continuation_state_v30 : SSTV30ContinuationState{};
    ptc_v30.validate();
    // Strict, once-per-call validation before mutating k, omega or wall data.
    // V33 forensic mode is opt-in and cannot certify an outer iteration.
    const auto v33_trace_cells=sst_v33_trace_cells(
        std::getenv("CFDX_SST_ROW_TRACE_V33"),n);
    std::vector<SSTV33SourceBranches> v33_previous_branches(v33_trace_cells.size());
    if (controls.sst_branch_preconditioner_v34_experimental && n>4096)
        throw std::invalid_argument(
            "V34 full-FV hotspot correction is explicitly limited to <=4096 cells; disable for large BFS");
    if (controls.sst_owner_consistent_projection_v21 && !wall_refresh)
        throw std::invalid_argument(
            "V21 SST owner projection requires a dynamic wall refresh callback");
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
    const auto refresh_projected = [&](const auto& kk, auto& ww) {
        SSTWallOwnerSyncV21 sync;
        if (controls.sst_owner_consistent_projection_v21) {
            sync=sst_synchronize_wall_owner_v21(
                kk,ww,*wall_refresh,active_omega_values,active_wall_g,
                controls.omega_min);
        } else {
            refresh_wall(kk,ww);
        }
        return sync;
    };
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
            k_bcs,omega_bcs,1,tolerance,face_values_now(),sources_now());
        // Assemble a wall-consistent *raw Picard proposal* for Anderson and
        // line search. The old state remains untouched as the merit baseline.
        // Re-evaluation later uses the same refreshed omega/G operator.
        SSTWallOwnerSyncV21 raw_wall_sync;
        if (controls.sst_owner_consistent_projection_v21)
            raw_wall_sync=refresh_projected(k,omega);
        // Optional Armijo-like backtracking on the *returned nonlinear SST
        // equations*, with the pressure/velocity/face flux held fixed.
        // It never changes the global conservation or convergence criteria.
        double line_alpha = 1.0;
        bool line_accepted = true;
        double line_old_merit = std::numeric_limits<double>::quiet_NaN();
        double line_new_merit = std::numeric_limits<double>::quiet_NaN();
        if (controls.sst_jfnk_v22_experimental && !inner.converged) {
            // V22: joint [Rk,Romega] JFNK. The Krylov action is built by
            // finite differences of the SAME FV assembly as the returned
            // SST certificate. Fixed row scales only condition the linear
            // system: the final nonlinear merit is always recomputed on the
            // candidate with all wall face values and G refreshed.
            const auto raw_k=k,raw_w=omega;
            k=old_k;omega=old_w;
            refresh_projected(k,omega);
            const auto start_k=k,start_w=omega;
            const auto base_sample=sst_reassembled_returned_state_residual(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
            line_old_merit=sst_nonlinear_merit(base_sample);
            const auto base_full=sst_coupled_residual_v22(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
            std::vector<double> physical(2*n),fixed_scales(2*n),
                                diagonal_preconditioner(2*n);
            for (std::size_t i=0;i<n;++i) {
                physical[i]=start_k(i);
                physical[n+i]=start_w(i);
            }
            for (std::size_t i=0;i<2*n;++i) {
                fixed_scales[i]=base_full.row_norms[i]*base_full.local_scales[i];
                if (!(fixed_scales[i]>0.) || !std::isfinite(fixed_scales[i]))
                    throw std::domain_error("V22 SST JFNK: unusable fixed FV scale");
                // Frozen scalar-transport Jacobi is only a LEFT
                // preconditioner for the log-variable Newton system.
                // It does not approximate away k/omega derivatives in J*v.
                const double frozen_log_diagonal=std::abs(
                    base_full.diagonal[i]*physical[i]/fixed_scales[i]);
                if (!std::isfinite(frozen_log_diagonal))
                    throw std::domain_error("V22 SST JFNK: invalid block diagonal");
                diagonal_preconditioner[i]=std::clamp(frozen_log_diagonal,0.02,50.0);
            }
            const SSTPhysicalResidualV22 full_residual=[&](const auto& state) {
                auto kk=start_k,ww=start_w;
                for (std::size_t i=0;i<n;++i) {
                    kk(i)=state[i];ww(i)=state[n+i];
                }
                ScalarBoundaryFaceValues staging_faces;
                SSTWallProductionV17 staging_wall;
                const ScalarBoundaryFaceValues* f=omega_face_values;
                const SSTWallProductionV17* g=nullptr;
                if (wall_refresh) {
                    (*wall_refresh)(kk,ww,staging_faces,staging_wall);
                    staging_wall.validate(n);
                    f=&staging_faces;g=&staging_wall;
                }
                return sst_coupled_residual_v22(
                    mesh,geometry,mass_flux,kk,ww,strain_rate,wall_distance,
                    controls,k_bcs,omega_bcs,f,g).values;
            };
            CoupledNewtonKrylovControls jctl;
            jctl.krylov_dimension=32;
            jctl.linear_relative_tolerance=0.25; // linear only, NOT SST convergence
            jctl.finite_difference_step=1e-7;
            jctl.maximum_backtracks=7;
            // Always retain the original V22 Jacobi proposal. The new block
            // is a separate candidate, sampled at a bounded cadence to avoid
            // doubling costly FV Jv evaluations at every SST iteration.
            // Only the reassembled physical merit can choose a direction.
            const auto newton=sst_jfnk_propose_v22(
                physical,fixed_scales,n,controls.k_min,controls.omega_min,
                full_residual,jctl,diagonal_preconditioner);
            std::unique_ptr<SSTBlockPreconditionerV23> local_blocks;
            std::unique_ptr<SSTSpatialBlockPreconditionerV24> spatial_blocks;
            SSTJFNKV22Candidate block_newton,spatial_newton,v25_newton,v26_newton,v27_newton,v28_newton;
            SSTJFNKV22Candidate pseudo_v29_newton, pseudo_v30_newton;
            bool pseudo_v30_probe=false;
            double pseudo_v30_cfl_before=ptc_v30.cfl;
            std::size_t pseudo_v30_shifted=0;
            std::string pseudo_v30_reason="disabled";
            bool pseudo_v29_probe=false;
            std::size_t pseudo_v29_shifted=0, pseudo_v29_skipped_wall=0;
            double pseudo_v29_cfl_before=pseudo_time_v29.cfl;
            double pseudo_v29_proposal_linf=std::numeric_limits<double>::infinity();
            double pseudo_v29_proposal_l2=std::numeric_limits<double>::infinity();
            std::string pseudo_v29_reason="disabled";
            SSTV25GlobalizationResult v25_search;
            SSTV26SmoothSearch v26_search;
            SSTV27TrustSearch v27_search;
            SSTV28BlendSearch v28_search;
            double v28_radius_before=adaptive_trust_v28.radius;
            bool v28_committed=false;
            double v28_commit_log_step=0.;
            double v28_commit_rho=0.;
            SSTV26DirectionalAudit v26_jvp;
            std::size_t v25_wall_corrected=0;
            double v25_wall_error=0.;
            std::string v25_wall_sampled_ids="none";
            std::string preconditioner_reason;
            const bool block_probe=controls.sst_block_preconditioner_v23_experimental &&
                                   (iter==1 || iter%10==0);
            const bool v28_rescue_probe=sst_v28_stagnation_probe(
                controls.sst_adaptive_trust_v28_experimental,iter,
                consecutive_line_search_stalls);
            const bool spatial_probe=controls.sst_spatial_block_preconditioner_v24_experimental &&
                                     (iter==1 || iter%10==0 || v28_rescue_probe);
            if (block_probe || spatial_probe) {
                try {
                    local_blocks=std::make_unique<SSTBlockPreconditionerV23>(
                        sst_local_source_blocks_v23(
                            mesh,geometry,k,omega,strain_rate,wall_distance,
                            controls,k_bcs,omega_bcs,face_values_now(),
                            sources_now(),base_full,fixed_scales,
                            diagonal_preconditioner));
                    if (std::getenv("CFDX_SST_JACOBIAN_AUDIT_V24") &&
                        (iter==1 || iter==10) && n<=5000) {
                        for (const std::size_t cell:{std::size_t{357},
                             std::size_t{363},std::size_t{426}}) {
                            if (cell>=n) continue;
                            try {
                                const auto audit=sst_local_jacobian_audit_v24(
                                    cell,physical,fixed_scales,
                                    diagonal_preconditioner,*local_blocks,full_residual);
                                std::cerr << std::scientific << std::setprecision(10)
                                          << "SST_JACOBIAN_AUDIT_V24 iteration=" << iter
                                          << " cell=" << cell
                                          << " approx_kk=" << audit.approximate[0]
                                          << " approx_kw=" << audit.approximate[1]
                                          << " approx_wk=" << audit.approximate[2]
                                          << " approx_ww=" << audit.approximate[3]
                                          << " true_kk=" << audit.actual[0]
                                          << " true_kw=" << audit.actual[1]
                                          << " true_wk=" << audit.actual[2]
                                          << " true_ww=" << audit.actual[3]
                                          << " max_relative_difference=" << audit.maximum_relative_difference
                                          << " qualified=0\n";
                            } catch (const std::exception& e) {
                                std::cerr << "SST_JACOBIAN_AUDIT_V24 iteration=" << iter
                                          << " cell=" << cell
                                          << " status=unavailable reason=" << e.what()
                                          << " qualified=0\n";
                            }
                        }
                    }
                    if (controls.sst_log_globalization_v25_experimental &&
                        spatial_probe && wall_refresh &&
                        active_wall_g.has_owner_omega_constraints()) {
                        // The owner cells with the largest *k*-row defect are
                        // sampled first. Do not special-case mesh indices 357
                        // or 363: the selection is independent of numbering.
                        std::vector<std::size_t> selected_wall;
                        for (std::size_t i=0;i<n;++i)
                            if (active_wall_g.has_omega_owner_cell[i])
                                selected_wall.push_back(i);
                        std::stable_sort(selected_wall.begin(),selected_wall.end(),
                            [&](std::size_t a,std::size_t b) {
                                return std::abs(base_full.values[a]/fixed_scales[a]) >
                                       std::abs(base_full.values[b]/fixed_scales[b]);
                            });
                        // A reproducible audit may prioritize a specifically
                        // identified wall owner (e.g. the V24 hotspot 357).
                        // It is never enabled from defaults, never edits F,
                        // and rejects invalid IDs rather than silently using
                        // a different cell. The normal path is index-agnostic.
                        if (const char* requested=std::getenv("CFDX_SST_V25_WALL_CELL")) {
                            std::size_t consumed=0;
                            const auto chosen=std::stoull(requested,&consumed,10);
                            if (consumed!=std::string(requested).size() ||
                                chosen>=n ||
                                !active_wall_g.has_omega_owner_cell[chosen])
                                throw std::invalid_argument(
                                    "V25 wall diagnostic cell is not a valid owner");
                            const auto pos=std::find(selected_wall.begin(),
                                                     selected_wall.end(),chosen);
                            if (pos==selected_wall.end())
                                throw std::invalid_argument("V25 wall owner not in sorted set");
                            selected_wall.erase(pos);
                            selected_wall.insert(selected_wall.begin(),chosen);
                        }
                        {
                            std::ostringstream ids;
                            for (std::size_t slot=0;slot<std::min<std::size_t>(4,selected_wall.size());++slot) {
                                if (slot) ids<<',';
                                ids<<selected_wall[slot];
                            }
                            v25_wall_sampled_ids=ids.str();
                        }
                        const auto correction=sst_correct_wall_blocks_v25(
                            *local_blocks,selected_wall,physical,fixed_scales,
                            diagonal_preconditioner,full_residual,4);
                        v25_wall_corrected=correction.corrected;
                        v25_wall_error=correction.max_relative_error;
                        local_blocks=std::make_unique<SSTBlockPreconditionerV23>(
                            correction.preconditioner);
                    }
                    if (controls.sst_branch_preconditioner_v34_experimental &&
                        spatial_probe && n<=4096) {
                        // The same full FV residual used by Newton is the
                        // authority for this bounded local P correction.
                        // Select by the signed, scale-normalized *physical*
                        // row defect. This does not hard-code BFS numbering.
                        auto hotspots=sst_v34_hotspot_cells(
                            base_full.values,fixed_scales,n,1);
                        // Controlled forensic override for comparing the
                        // same physical cells across solver versions. It is
                        // not a numerically privileged mesh cell index.
                        if (const auto* explicit_ids=std::getenv("CFDX_SST_V34_CELL_IDS")) {
                            hotspots=sst_v33_trace_cells(explicit_ids,n);
                            if (hotspots.size()>2)
                                throw std::invalid_argument(
                                    "V34 diagnostic override permits at most two cells");
                        }
                        const auto cross_now=sst_boundary_gradient_dot(
                            k,omega,mesh,geometry,k_bcs,omega_bcs,
                            face_values_now(),controls.sst_use_wls_boundary_gradients,
                            controls.sst_gradient_dimension);
                        std::vector<SSTV33SourceBranches> branches;
                        branches.reserve(hotspots.size());
                        for (const auto cell:hotspots) {
                            const auto blending=compute_sst_blending(
                                k(cell),omega(cell),wall_distance(cell),
                                controls.molecular_viscosity,controls.beta_star,
                                cross_now[cell],controls.density,controls.sst_sigma_w2);
                            branches.push_back(sst_v33_source_branches(
                                k(cell),omega(cell),strain_rate(cell),blending.second,
                                sources_now() ? sources_now()->for_cell(cell) : -1.,
                                controls.density,controls.a1,controls.beta_star,
                                controls.sst_production_limiter,controls.k_min,
                                controls.omega_min));
                        }
                        const auto correction=sst_v34_correct_hotspot_blocks(
                            *local_blocks,hotspots,branches,physical,fixed_scales,
                            diagonal_preconditioner,full_residual);
                        local_blocks=std::make_unique<SSTBlockPreconditionerV23>(
                            correction.preconditioner);
                        std::cerr << std::scientific << std::setprecision(10)
                            << "SST_BRANCH_PRECOND_V34 iteration=" << iter
                            << " evaluated=" << correction.residual_evaluations
                            << " corrected=" << correction.corrected_cells.size()
                            << " skipped_switch=" << correction.skipped_near_switch
                            << " skipped_nonsmooth=" << correction.skipped_nonsmooth
                            << " skipped_singular=" << correction.skipped_singular
                            << " max_difference=" << correction.maximum_relative_difference
                            << " max_unilateral_gap=" << correction.maximum_one_sided_disagreement
                            << " cells=";
                        for (std::size_t id=0;id<correction.sampled_cells.size();++id) {
                            if (id) std::cerr << ',';
                            std::cerr << correction.sampled_cells[id];
                        }
                        std::cerr << " qualified=0\n";
                    }
                    if (block_probe) {
                        const CoupledRightPreconditioner block_right_action=[&local_blocks](const auto& rhs) {
                            return local_blocks->apply(rhs);
                        };
                        block_newton=sst_jfnk_propose_v22(
                            physical,fixed_scales,n,controls.k_min,controls.omega_min,
                            full_residual,jctl,diagonal_preconditioner,
                            block_right_action);
                    }
                    if (spatial_probe) {
                        // Freeze the EXACT same production FV transport CSR at
                        // the current state. Approximate only the Jacobian used
                        // by P; true Jv and the nonlinear merit remain unchanged.
                        auto equations=sst_reassembled_equations_v22(
                            mesh,geometry,mass_flux,k,omega,strain_rate,
                            wall_distance,controls,k_bcs,omega_bcs,
                            face_values_now(),sources_now());
                        spatial_blocks=std::make_unique<SSTSpatialBlockPreconditionerV24>(
                            equations.first.matrix,equations.second.matrix,
                            physical,fixed_scales,diagonal_preconditioner,*local_blocks);
                        if (const char* pair_text=std::getenv("CFDX_SST_V34_CROSS_AUDIT")) {
                            if (!controls.sst_branch_preconditioner_v34_experimental)
                                throw std::invalid_argument("V34 cross-cell audit requires V34 opt-in");
                            const auto pair=sst_v34_parse_cross_pair(pair_text,n);
                            const double kk=sst_v34_frozen_offcell_coefficient(
                                equations.first.matrix,pair[0],pair[1]);
                            const double ww=sst_v34_frozen_offcell_coefficient(
                                equations.second.matrix,pair[0],pair[1]);
                            const auto cross_audit=sst_v34_cross_cell_jacobian_audit(
                                pair[0],pair[1],physical,fixed_scales,
                                diagonal_preconditioner,full_residual,kk,ww);
                            std::cerr << std::scientific << std::setprecision(10)
                                << "SST_CROSS_JACOBIAN_V34 iteration=" << iter
                                << " row=" << pair[0] << " column=" << pair[1]
                                << " approximate_kk=" << cross_audit.approximate[0]
                                << " approximate_kw=" << cross_audit.approximate[1]
                                << " approximate_wk=" << cross_audit.approximate[2]
                                << " approximate_ww=" << cross_audit.approximate[3]
                                << " true_kk=" << cross_audit.actual[0]
                                << " true_kw=" << cross_audit.actual[1]
                                << " true_wk=" << cross_audit.actual[2]
                                << " true_ww=" << cross_audit.actual[3]
                                << " relative_difference=" << cross_audit.maximum_relative_difference
                                << " evaluations=" << cross_audit.residual_evaluations
                                << " qualified=0\n";
                        }
                        const CoupledRightPreconditioner spatial_action=[&spatial_blocks](const auto& rhs) {
                            return spatial_blocks->apply(rhs);
                        };
                        spatial_newton=sst_jfnk_propose_v22(
                            physical,fixed_scales,n,controls.k_min,controls.omega_min,
                            full_residual,jctl,diagonal_preconditioner,
                            spatial_action);
                    }
                } catch (const std::exception& error) {
                    // Preconditioning is never a reason to accept a bad step.
                    preconditioner_reason=error.what();
                }
            }
            // V25: search in the actual LOG Newton direction after the true
            // GMRES gate, even when the frozen-scale inner line search failed.
            // Each trial uses a LOCAL staged wall triplet, not caller state.
            if (controls.sst_log_globalization_v25_experimental &&
                spatial_probe && !spatial_newton.log_direction.empty()) {
                const SSTV25PhysicalEvaluator physical_trial=[&](const auto& candidate) {
                    auto kk=start_k,ww=start_w;
                    for (std::size_t i=0;i<n;++i) {
                        kk(i)=candidate[i];ww(i)=candidate[n+i];
                    }
                    ScalarBoundaryFaceValues trial_faces;
                    SSTWallProductionV17 trial_production;
                    const ScalarBoundaryFaceValues* f=omega_face_values;
                    const SSTWallProductionV17* g=nullptr;
                    if (wall_refresh) {
                        if (controls.sst_owner_consistent_projection_v21)
                            sst_synchronize_wall_owner_v21(kk,ww,*wall_refresh,
                                trial_faces,trial_production,controls.omega_min);
                        else {
                            (*wall_refresh)(kk,ww,trial_faces,trial_production);
                            trial_production.validate(n);
                        }
                        f=&trial_faces;g=&trial_production;
                    }
                    const auto sample=sst_reassembled_returned_state_residual(
                        mesh,geometry,mass_flux,kk,ww,strain_rate,wall_distance,
                        controls,k_bcs,omega_bcs,f,g);
                    SSTV25PhysicalTrial out;
                    out.merit=sst_nonlinear_merit(sample);
                    out.smooth_l2=sst_nonlinear_smooth_merit_v26(sample);
                    out.k_local=sample.k_local_relative;
                    out.omega_local=sample.omega_local_relative;
                    out.k_worst=sample.k_worst_cell;
                    out.omega_worst=sample.omega_worst_cell;
                    out.fields.resize(2*n);
                    for (std::size_t i=0;i<n;++i) {
                        out.fields[i]=kk(i);out.fields[n+i]=ww(i);
                    }
                    return out;
                };
                if (controls.sst_smooth_globalization_v26_experimental) {
                    // Three log perturbation amplitudes, measured on the SAME
                    // unprojected FV residual as Jv; informational only.
                    v26_jvp=sst_audit_log_jvp_v26(
                        physical,spatial_newton.log_direction,n,
                        controls.k_min,controls.omega_min,full_residual);
                    // ONE set of reassembled projected FV trials yields both
                    // V25 (minimum physical Linf) and V26 (balanced RMS) picks.
                    v26_search=sst_globalize_log_direction_v26(
                        physical,spatial_newton.log_direction,n,controls.k_min,
                        controls.omega_min,line_old_merit,
                        sst_nonlinear_smooth_merit_v26(base_sample),physical_trial);
                    v25_search=v26_search.linf;
                    if (v26_search.smooth_accepted) {
                        v26_newton.available=true;
                        v26_newton.fields=v26_search.smooth_fields;
                        v26_newton.reason="V26 balanced L2 candidate passed physical Linf gate";
                    }
                } else {
                    v25_search=sst_globalize_log_direction_v25(
                        physical,spatial_newton.log_direction,n,controls.k_min,
                        controls.omega_min,line_old_merit,physical_trial);
                }
                if (v25_search.accepted) {
                    v25_newton.available=true;
                    v25_newton.fields=v25_search.fields;
                    v25_newton.reason="V25 physical log line search accepted";
                }
                if (controls.sst_subspace_trust_v27_experimental &&
                    controls.sst_smooth_globalization_v26_experimental) {
                    std::vector<double> picard_physical(2*n);
                    for (std::size_t i=0;i<n;++i) {
                        picard_physical[i]=raw_k(i);
                        picard_physical[n+i]=raw_w(i);
                    }
                    v27_search=sst_cauchy_subspace_v27(
                        physical,spatial_newton.log_direction,picard_physical,
                        n,controls.k_min,controls.omega_min,
                        line_old_merit,sst_nonlinear_smooth_merit_v26(base_sample),
                        base_sample.k_worst_cell,base_sample.omega_worst_cell,
                        physical_trial,
                        controls.sst_adaptive_trust_v28_experimental ?
                           adaptive_trust_v28.radius : 0.5);
                    if (v27_search.accepted) {
                        v27_newton.available=true;
                        v27_newton.fields=v27_search.fields;
                        v27_newton.reason="V27 true projected SST subspace Cauchy candidate";
                    }
                    if (controls.sst_adaptive_trust_v28_experimental &&
                        v27_search.accepted && newton.available) {
                        v28_search=sst_blend_adaptive_trust_v28(
                            physical,newton.fields,v27_search.fields,n,
                            controls.k_min,controls.omega_min,
                            line_old_merit,
                            sst_nonlinear_smooth_merit_v26(base_sample),
                            adaptive_trust_v28.radius,physical_trial);
                        if (v28_search.accepted) {
                            v28_newton.available=true;
                            v28_newton.fields=v28_search.fields;
                            v28_newton.reason=
                                "V28 projected logarithmic Newton-Cauchy blend";
                        }
                        // A losing V28 blend must not shrink the radius before
                        // the final FV comparison: V27/Jacobi may still win.
                        // Only the final rejected state shrinks it below.
                    }
                }
            }
            // V30: use the SAME joint FV residual as Newton, but shift its
            // full k-omega Jacobian, not two segregated frozen equations.
            // Only F(x)+D*(x-x_old) is passed to JFNK; the accepted state
            // must still improve F(x) with its original steady wall closure.
            pseudo_v30_probe=controls.sst_coupled_ptc_v30_experimental &&
                (iter==1 || iter%10==0 || consecutive_line_search_stalls>0);
            if (pseudo_v30_probe) {
                try {
                    std::vector<double> face_abs(n,0.);
                    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
                        const double flux=std::abs(mass_flux(f));
                        if (!std::isfinite(flux))
                            throw std::invalid_argument("V30 nonfinite face flux");
                        face_abs[mesh.ownership().owner(f)]+=flux;
                        const auto neighbour=mesh.ownership().neighbour(f);
                        if (neighbour>=0)
                            face_abs[static_cast<std::size_t>(neighbour)]+=flux;
                    }
                    std::vector<double> flux2(2*n);
                    for (std::size_t i=0;i<n;++i)
                        flux2[i]=flux2[n+i]=face_abs[i];
                    std::vector<unsigned char> fixed(2*n,0);
                    if (sources_now() && sources_now()->has_owner_omega_constraints())
                        for (std::size_t i=0;i<n;++i)
                            fixed[n+i]=sources_now()->has_omega_owner_cell[i];
                    const auto mass=sst_v30_coupled_pseudo_mass(
                        base_full.diagonal,base_full.row_norms,flux2,fixed,ptc_v30.cfl,
                        ptc_v30.v32_max_diagonal_fraction);
                    for (double value:mass)
                        pseudo_v30_shifted+=(value>0.);
                    SSTPhysicalResidualV22 augmented=[&](const auto& state) {
                        return sst_v30_shift_residual(
                            full_residual(state),state,physical,mass);
                    };
                    // Left-scale the shifted Newton rows with A_ii + D_tau,
                    // rather than the unshifted A_ii (V22 row scaling).
                    // Full J*v and every accepted FV merit remain untouched.
                    const auto augmented_scaling=sst_v30_augmented_row_scaling(
                        base_full.diagonal,mass,physical,fixed_scales);
                    pseudo_v30_newton=sst_jfnk_propose_v22(
                        physical,fixed_scales,n,controls.k_min,controls.omega_min,
                        augmented,jctl,augmented_scaling);
                    pseudo_v30_reason=pseudo_v30_newton.available
                        ? "jfnk_shifted_candidate" : pseudo_v30_newton.reason;
                } catch (const std::exception& e) {
                    pseudo_v30_reason=e.what();
                }
                refresh_projected(k,omega);
            }
            // V29: independent local pseudo-time continuation proposal. The
            // pseudo-mass is used ONLY in a frozen candidate solve. The
            // unshifted FV operator with refreshed wall G and omega is used
            // for every acceptance test, exactly as for Newton/Picard.
            pseudo_v29_probe=controls.sst_pseudo_transient_v29_experimental &&
                (iter==1 || iter%10==0 || consecutive_line_search_stalls>0);
            if (pseudo_v29_probe) {
                try {
                    auto frozen=sst_reassembled_equations_v22(
                        mesh,geometry,mass_flux,start_k,start_w,strain_rate,
                        wall_distance,controls,k_bcs,omega_bcs,
                        face_values_now(),sources_now());
                    cfdx::core::Vector xk(n,0.),xw(n,0.);
                    for (std::size_t i=0;i<n;++i) {
                        xk(i)=start_k(i);xw(i)=start_w(i);
                    }
                    const auto ak=sst_add_local_pseudo_time_v29(
                        frozen.first,xk,pseudo_time_v29.cfl);
                    const auto aw=sst_add_local_pseudo_time_v29(
                        frozen.second,xw,pseudo_time_v29.cfl,
                        sources_now() && sources_now()->has_owner_omega_constraints()
                            ? sources_now()->has_omega_owner_cell
                            : std::vector<unsigned char>{});
                    pseudo_v29_shifted=ak.shifted_rows+aw.shifted_rows;
                    pseudo_v29_skipped_wall=aw.skipped_wall_rows;
                    // Only these provisional linear systems have pseudo-time
                    // stabilization. Never call a shifted residual 'physical'.
                    const ScalarSolveControls sc_k{
                        1000,sst_linear_tolerance_for_solution_defect(
                            frozen.first,tolerance),1.0};
                    const ScalarSolveControls sc_w{
                        1000,sst_linear_tolerance_for_solution_defect(
                            frozen.second,tolerance),1.0};
                    const auto rk=solve_scalar_equation(frozen.first,xk,sc_k);
                    const auto rw=solve_scalar_equation(frozen.second,xw,sc_w);
                    if (rk.status!=cfdx::core::SolverStatus::CONVERGED ||
                        rw.status!=cfdx::core::SolverStatus::CONVERGED)
                        throw std::runtime_error("V29 pseudo-time Krylov incomplete");
                    auto ck=start_k,cw=start_w;
                    for (std::size_t i=0;i<n;++i) {
                        if (!(xk(i)>=controls.k_min) ||
                            !(xw(i)>=controls.omega_min) ||
                            !std::isfinite(xk(i)) || !std::isfinite(xw(i)))
                            throw std::runtime_error("V29 pseudo-time positivity");
                        ck(i)=xk(i);cw(i)=xw(i);
                    }
                    // This mutates a staging wall cache, not the accepted k/w.
                    refresh_projected(ck,cw);
                    const auto sample=sst_reassembled_returned_state_residual(
                        mesh,geometry,mass_flux,ck,cw,strain_rate,wall_distance,
                        controls,k_bcs,omega_bcs,face_values_now(),sources_now());
                    pseudo_v29_proposal_linf=sst_nonlinear_merit(sample);
                    pseudo_v29_proposal_l2=sst_nonlinear_smooth_merit_v26(sample);
                    // Strict dual objective for the pseudo-time candidate;
                    // neither criterion is allowed to override final SST gates.
                    if (!sst_nonlinear_step_reduces_merit(
                            line_old_merit,pseudo_v29_proposal_linf,1.) ||
                        !(pseudo_v29_proposal_l2 <
                          sst_nonlinear_smooth_merit_v26(base_sample)))
                        throw std::runtime_error("V29 pseudo-time not dual-descending");
                    pseudo_v29_newton.available=true;
                    pseudo_v29_newton.fields.resize(2*n);
                    for (std::size_t i=0;i<n;++i) {
                        pseudo_v29_newton.fields[i]=ck(i);
                        pseudo_v29_newton.fields[n+i]=cw(i);
                    }
                    pseudo_v29_newton.reason="V29 local pseudo-time dual FV decrease";
                    pseudo_v29_reason="dual_descending";
                } catch (const std::exception& e) {
                    pseudo_v29_reason=e.what();
                }
                // Restore wall closure for the committed baseline; the final
                // selection reprojects each trial independently.
                refresh_projected(k,omega);
            }
            // Safeguard: Picard plus the original V22 Jacobi always compete.
            // Neither a V23 nor a V24/V25 failed candidate can mask V22.
            std::vector<const SSTJFNKV22Candidate*> proposals={nullptr};
            if (newton.available) proposals.push_back(&newton);
            if (block_newton.available) proposals.push_back(&block_newton);
            if (spatial_newton.available) proposals.push_back(&spatial_newton);
            if (v25_newton.available) proposals.push_back(&v25_newton);
            if (v26_newton.available) proposals.push_back(&v26_newton);
            if (v27_newton.available) proposals.push_back(&v27_newton);
            if (v28_newton.available) proposals.push_back(&v28_newton);
            if (pseudo_v29_newton.available) proposals.push_back(&pseudo_v29_newton);
            if (pseudo_v30_newton.available) proposals.push_back(&pseudo_v30_newton);
            double best_merit=std::numeric_limits<double>::infinity();
            auto best_k=old_k,best_w=old_w;
            double best_alpha=0.;
            std::size_t best_choice=0;
            for (std::size_t choice=0;choice<proposals.size();++choice) {
                for (double alpha:{1.,0.5,0.25,0.125,0.0625,0.03125,
                                  0.015625,0.0078125}) {
                    bool feasible=true;
                    for (std::size_t i=0;i<n;++i) {
                        const double target_k=proposals[choice]
                            ? proposals[choice]->fields[i] : raw_k(i);
                        const double target_w=proposals[choice]
                            ? proposals[choice]->fields[n+i] : raw_w(i);
                        const double ck=start_k(i)+alpha*(target_k-start_k(i));
                        const double cw=start_w(i)+alpha*(target_w-start_w(i));
                        if (!std::isfinite(ck)||!std::isfinite(cw)||
                            ck<controls.k_min||cw<controls.omega_min) {
                            feasible=false;break;
                        }
                        k(i)=ck;omega(i)=cw;
                    }
                    if (!feasible) continue;
                    try {
                        refresh_projected(k,omega);
                        const double merit=sst_nonlinear_merit(
                            sst_reassembled_returned_state_residual(
                                mesh,geometry,mass_flux,k,omega,strain_rate,
                                wall_distance,controls,k_bcs,omega_bcs,
                                face_values_now(),sources_now()));
                        if (sst_nonlinear_step_reduces_merit(
                                line_old_merit,merit,alpha) && merit<best_merit) {
                            best_merit=merit;
                            best_k=k;best_w=omega;
                            best_alpha=alpha;
                            best_choice=choice;
                        }
                    } catch (const std::exception&) {
                        // Reject this trial only. Neither wall state nor
                        // nonlinear field is committed before acceptance.
                    }
                }
            }
            line_accepted=std::isfinite(best_merit);
            if (line_accepted) {
                k=std::move(best_k);omega=std::move(best_w);
                line_alpha=best_alpha;line_new_merit=best_merit;
                consecutive_line_search_stalls=0;
                if (controls.sst_coupled_ptc_v30_experimental) {
                    ptc_v30.accept_inner(line_old_merit,line_new_merit);
                    ptc_v30.record_selected_source(
                        pseudo_v30_newton.available &&
                        proposals[best_choice]==&pseudo_v30_newton);
                }
                if (controls.sst_pseudo_transient_v29_experimental &&
                    pseudo_v29_newton.available &&
                    proposals[best_choice]==&pseudo_v29_newton)
                    pseudo_time_v29.selected(line_old_merit,line_new_merit);
                // Adapt the radius only for a candidate ACTUALLY selected
                // by the strict physical FV line search, not one that merely
                // passed an internal V28 trial but lost to V22/V27.
                if (controls.sst_adaptive_trust_v28_experimental &&
                    best_alpha==1.) {
                    // Both V27 and V28 are measured on the same CURRENT,
                    // wall-projected FV operator. The shared controller can
                    // update from either one, but ONLY when it wins the
                    // unchanged final global Linf comparison.
                    const bool chose_v28=v28_newton.available &&
                        proposals[best_choice]==&v28_newton &&
                        std::abs(line_new_merit-v28_search.merit)<=
                           1e-8*std::max(1.,line_new_merit);
                    const bool chose_v27=v27_newton.available &&
                        proposals[best_choice]==&v27_newton &&
                        std::abs(line_new_merit-v27_search.merit)<=
                           1e-8*std::max(1.,line_new_merit);
                    if (chose_v28) {
                        adaptive_trust_v28.commit(
                            v28_search.actual_over_predicted,
                            v28_search.used_log_radius);
                        v28_committed=true;
                        v28_commit_log_step=v28_search.used_log_radius;
                        v28_commit_rho=v28_search.actual_over_predicted;
                    } else if (chose_v27) {
                        double log_size=0.;
                        for (std::size_t i=0;i<physical.size();++i)
                            log_size=std::max(log_size,std::abs(
                                std::log(v27_search.fields[i]/physical[i])));
                        adaptive_trust_v28.commit(v27_search.ratio,log_size);
                        v28_committed=true;
                        v28_commit_log_step=log_size;
                        v28_commit_rho=v27_search.ratio;
                    }
                }
            } else {
                k=old_k;omega=old_w;
                line_alpha=0.;line_new_merit=line_old_merit;
                ++consecutive_line_search_stalls;
                if (controls.sst_coupled_ptc_v30_experimental)
                    ptc_v30.reject_inner();
                if (controls.sst_pseudo_transient_v29_experimental)
                    pseudo_time_v29.rejected();
                // Only the final physical FV decision can shrink the radius.
                if (controls.sst_adaptive_trust_v28_experimental)
                    adaptive_trust_v28.reject();
            }
            refresh_wall(k,omega);
            if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
                (iter==1 || iter%10==0 || !line_accepted ||
                 v28_rescue_probe || v28_committed))
                std::cerr << std::scientific << std::setprecision(10)
                          << "SST_JFNK_V22 iteration=" << iter
                          << " proposed=" << newton.available
                          << " linear_true_relative=" << newton.linear_true_relative
                          << " evaluations=" << newton.function_evaluations
                          << " preconditioner=row_jacobi"
                          << " block_probe=" << block_probe
                          << " block_proposed=" << block_newton.available
                          << " block_linear_true_relative=" << block_newton.linear_true_relative
                          << " block_evaluations=" << block_newton.function_evaluations
                          << " spatial_probe=" << spatial_probe
                          << " spatial_proposed=" << spatial_newton.available
                          << " spatial_linear_true_relative=" << spatial_newton.linear_true_relative
                          << " spatial_evaluations=" << spatial_newton.function_evaluations
                          << " spatial_reason=" << (spatial_newton.reason.empty() ? "none" :
                              (spatial_newton.available ? "accepted" : "rejected"))
                          << " spatial_offcell=" << (spatial_blocks ? spatial_blocks->off_cell_coefficients() : 0)
                          << " v25_active=" << controls.sst_log_globalization_v25_experimental
                          << " v25_direction=" << !spatial_newton.log_direction.empty()
                          << " v25_attempted=" << v25_search.attempted
                          << " v25_valid=" << v25_search.valid
                          << " v25_accepted=" << v25_search.accepted
                          << " v25_alpha=" << v25_search.alpha
                          << " v25_first_merit=" << v25_search.full_step_merit
                          << " v25_best_trial_merit=" << v25_search.best_trial_merit
                          << " v25_best_merit=" << v25_search.merit
                          << " v25_max_log_direction=" << v25_search.maximum_absolute_log_direction
                          << " v25_kworst=" << v25_search.k_worst
                          << " v25_wworst=" << v25_search.omega_worst
                          << " v26_active=" << controls.sst_smooth_globalization_v26_experimental
                          << " v26_jvp_available=" << v26_jvp.available
                          << " v26_jvp_evaluations=" << v26_jvp.evaluations
                          << " v26_jvp_h1=" << v26_jvp.steps[0]
                          << " v26_jvp_h2=" << v26_jvp.steps[1]
                          << " v26_jvp_h3=" << v26_jvp.steps[2]
                          << " v26_jvp_norm_h1=" << v26_jvp.derivative_norms[0]
                          << " v26_jvp_norm_h2=" << v26_jvp.derivative_norms[1]
                          << " v26_jvp_norm_h3=" << v26_jvp.derivative_norms[2]
                          << " v26_jvp_one_h1=" << v26_jvp.one_sided_disagreement[0]
                          << " v26_jvp_one_h2=" << v26_jvp.one_sided_disagreement[1]
                          << " v26_jvp_one_h3=" << v26_jvp.one_sided_disagreement[2]
                          << " v26_jvp_adj_12=" << v26_jvp.adjacent_scale_disagreement[0]
                          << " v26_jvp_adj_23=" << v26_jvp.adjacent_scale_disagreement[1]
                          << " v26_jvp_max_disagreement=" << v26_jvp.maximum_disagreement
                          << " v26_l2_before=" << v26_search.baseline_l2
                          << " v26_l2_best_trial=" << v26_search.best_trial_l2
                          << " v26_l2_accepted=" << v26_search.smooth_accepted
                          << " v26_l2_selected=" << v26_search.selected_l2
                          << " v26_linf_selected=" << v26_search.selected_linf
                          << " v26_alpha=" << v26_search.smooth_alpha
                          << " v27_active=" << controls.sst_subspace_trust_v27_experimental
                          << " v27_bases=" << v27_search.bases
                          << " v27_evaluations=" << v27_search.evaluations
                          << " v27_valid=" << v27_search.valid_trials
                          << " v27_accepted=" << v27_search.accepted
                          << " v27_linf=" << v27_search.merit
                          << " v27_l2=" << v27_search.smooth_l2
                          << " v27_radius=" << v27_search.chosen_radius
                          << " v27_ratio=" << v27_search.ratio
                          << " v27_slope=" << v27_search.directional_l2_slope
                          << " v27_kworst=" << v27_search.k_active
                          << " v27_wworst=" << v27_search.omega_active
                          << " v27_choice=" << v27_search.direction
                          << " v30_active=" << controls.sst_coupled_ptc_v30_experimental
                          << " v30_probe=" << pseudo_v30_probe
                          << " v30_cfl_before=" << pseudo_v30_cfl_before
                          << " v30_cfl_after=" << ptc_v30.cfl
                          << " v30_selected_own=" << ptc_v30.selected_own
                          << " v30_selected_other=" << ptc_v30.selected_other
                          << " v30_shifted=" << pseudo_v30_shifted
                          << " v30_available=" << pseudo_v30_newton.available
                          << " v30_linear_true_relative=" << pseudo_v30_newton.linear_true_relative
                          << " v30_reason=" << pseudo_v30_reason
                          << " v29_active=" << controls.sst_pseudo_transient_v29_experimental
                          << " v29_probe=" << pseudo_v29_probe
                          << " v29_cfl_before=" << pseudo_v29_cfl_before
                          << " v29_cfl_after=" << pseudo_time_v29.cfl
                          << " v29_shifted=" << pseudo_v29_shifted
                          << " v29_wall_skip=" << pseudo_v29_skipped_wall
                          << " v29_linf=" << pseudo_v29_proposal_linf
                          << " v29_l2=" << pseudo_v29_proposal_l2
                          << " v29_available=" << pseudo_v29_newton.available
                          << " v29_reason=" << pseudo_v29_reason
                          << " v28_active=" << controls.sst_adaptive_trust_v28_experimental
                          << " v28_rescue_probe=" << v28_rescue_probe
                          << " v28_accepted=" << v28_search.accepted
                          << " v28_candidates=" << v28_search.admissible_trials
                          << " v28_evaluations=" << v28_search.evaluations
                          << " v28_linf=" << v28_search.merit
                          << " v28_l2=" << v28_search.smooth_l2
                          << " v28_rho=" << v28_search.actual_over_predicted
                          << " v28_jacobi_weight=" << v28_search.weight_jacobi
                          << " v28_step=" << v28_search.used_log_radius
                          << " v28_radius_before=" << v28_radius_before
                          << " v28_radius_after=" << adaptive_trust_v28.radius
                          << " v28_committed=" << v28_committed
                          << " v28_commit_log_step=" << v28_commit_log_step
                          << " v28_commit_rho=" << v28_commit_rho
                          << " v28_commit_source=" << (v28_committed ?
                              (proposals[best_choice]==&v27_newton ? "v27" : "v28"):
                              "none")
                          << " v25_wall_sampled_ids=" << v25_wall_sampled_ids
                          << " v25_wall_corrected=" << v25_wall_corrected
                          << " v25_wall_relative_error=" << v25_wall_error
                          << " coupled_blocks=" << (local_blocks ? local_blocks->coupled_blocks() : 0)
                          << " block_fallbacks=" << (local_blocks ? local_blocks->diagonal_fallbacks() : 0)
                          << " preconditioner_reason=" << (preconditioner_reason.empty() ? "none" : "fallback")
                          << " alpha=" << line_alpha
                          << " old_merit=" << line_old_merit
                          << " new_merit=" << line_new_merit
                          << " selected=" << (line_accepted ?
                              (best_choice==0 ? "picard" :
                               (proposals[best_choice]==&newton ? "jfnk" :
                                (proposals[best_choice]==&block_newton ? "jfnk_block" :
                                 (proposals[best_choice]==&v25_newton ? "jfnk_spatial_log_v25" :
 (proposals[best_choice]==&v26_newton ? "jfnk_smooth_v26" :
  (proposals[best_choice]==&v27_newton ? "sst_subspace_trust_v27" :
   (proposals[best_choice]==&pseudo_v29_newton ? "sst_pseudo_transient_v29" :
    (proposals[best_choice]==&pseudo_v30_newton ? "sst_coupled_ptc_v30" :
   (proposals[best_choice]==&v28_newton ? "sst_adaptive_trust_v28" :
    "jfnk_spatial")))))))))
                               : "none")
                          << " accepted=" << line_accepted
                          << " reason=" << newton.reason
                          << " qualified=0\n";
        } else if (controls.sst_anderson_m2_v20 && !inner.converged) {

            // Safeguard: the two-depth extrapolation is only a proposal.
            // Both trial families use the *current* wall omega/G and a fully
            // reassembled SST physical residual. The frozen Krylov residual
            // never controls acceptance.
            const auto raw_k=k,raw_w=omega;
            std::vector<double> xk(n),xw(n),gk(n),gw(n);
            for (std::size_t i=0;i<n;++i) {
                xk[i]=old_k(i);xw[i]=old_w(i);
                gk[i]=raw_k(i);gw[i]=raw_w(i);
            }
            const auto accelerated=anderson_v20.propose(
                xk,xw,gk,gw,controls.k_min,controls.omega_min);
            anderson_v20.record(xk,xw,gk,gw);
            k=old_k;omega=old_w;
            refresh_wall(k,omega);
            line_old_merit=sst_nonlinear_merit(sst_reassembled_returned_state_residual(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,face_values_now(),sources_now()));
            double best_merit=std::numeric_limits<double>::infinity();
            auto best_k=old_k,best_w=old_w;
            double best_alpha=0.,best_theta=0.;
            for (int choice=0;choice<(accelerated.available?2:1);++choice) {
                const auto& tk=choice==0?gk:accelerated.k;
                const auto& tw=choice==0?gw:accelerated.omega;
                for (double alpha:{1.,0.5,0.25,0.125,0.0625,0.03125,0.015625,0.0078125}) {
                    bool admissible=true;
                    for (std::size_t i=0;i<n;++i) {
                        const double ck=old_k(i)+alpha*(tk[i]-old_k(i));
                        const double cw=old_w(i)+alpha*(tw[i]-old_w(i));
                        if (!std::isfinite(ck)||!std::isfinite(cw)||
                            ck<controls.k_min||cw<controls.omega_min) {
                            admissible=false;break;
                        }
                        k(i)=ck;omega(i)=cw;
                    }
                    if (!admissible) continue;
                    refresh_projected(k,omega);
                    const double merit=sst_nonlinear_merit(sst_reassembled_returned_state_residual(
                        mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                        controls,k_bcs,omega_bcs,face_values_now(),sources_now()));
                    if (sst_nonlinear_step_reduces_merit(line_old_merit,merit,alpha)&&
                        merit<best_merit) {
                        best_merit=merit;best_k=k;best_w=omega;
                        best_alpha=alpha;best_theta=choice==0?0.:accelerated.theta;
                    }
                }
            }
            line_accepted=std::isfinite(best_merit);
            if (line_accepted) {
                k=std::move(best_k);omega=std::move(best_w);
                line_alpha=best_alpha;line_new_merit=best_merit;
                consecutive_line_search_stalls=0;
            } else {
                k=old_k;omega=old_w;
                line_alpha=0.;line_new_merit=line_old_merit;
                ++consecutive_line_search_stalls;
            }
            refresh_wall(k,omega);
            if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
                (iter==1||iter%10==0||!line_accepted))
                std::cerr << std::scientific << std::setprecision(10)
                          << "SST_ANDERSON_V20 iteration=" << iter
                          << " theta=" << best_theta
                          << " alpha=" << line_alpha
                          << " old_merit=" << line_old_merit
                          << " new_merit=" << line_new_merit
                          << " accepted=" << line_accepted
                          << " qualified=0\n";
        } else if (controls.sst_nonlinear_line_search && !inner.converged) {
            auto trial_k=k, trial_w=omega;
            k=old_k; omega=old_w;
            refresh_wall(k,omega);
            const auto starting=sst_reassembled_returned_state_residual(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
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
        } else if (!inner.converged && picard_damping_v41 < 1.0) {
            // Historical fixed damping remains unchanged when not using
            // line search, and never certifies an undamped candidate.
            const double alpha=picard_damping_v41;
            for(std::size_t i=0;i<n;++i) {
                k(i)=old_k(i)+alpha*(k(i)-old_k(i));
                omega(i)=old_w(i)+alpha*(omega(i)-old_w(i));
            }
        }
        // Never modify a rejected old state. Any accepted candidate already
        // lies on the owner constraint manifold; refresh both G and omega
        // again to document the exact state used by the residual certificate.
        SSTWallOwnerSyncV21 final_wall_sync;
        if (controls.sst_owner_consistent_projection_v21 && line_accepted)
            final_wall_sync=refresh_projected(k,omega);
        else
            refresh_wall(k,omega);
        if (std::getenv("CFDX_SST_DIAGNOSTICS") &&
            controls.sst_owner_consistent_projection_v21 &&
            (iter==1 || iter%100==0 || iter==max_iterations))
            std::cerr << std::scientific << std::setprecision(12)
                      << "SST_OWNER_SYNC_V21 iteration=" << iter
                      << " owners=" << final_wall_sync.owners
                      << " initial_gap=" << raw_wall_sync.maximum_before
                      << " trial_gap=" << final_wall_sync.maximum_before
                      << " final_gap=" << final_wall_sync.maximum_after
                      << " final_relative=" << final_wall_sync.relative_after
                      << " worst_cell=" << raw_wall_sync.worst_cell
                      << " line_accepted=" << line_accepted
                      << " qualified=0\n";
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
            picard_damping_v41 < 1.0 ||
            controls.sst_nonlinear_line_search ||
            controls.sst_anderson_m2_v20 ||
            controls.sst_jfnk_v22_experimental ||
            !v33_trace_cells.empty() ||
            std::getenv("CFDX_SST_DIAGNOSTICS") != nullptr;
        SSTNonlinearResidualSample returned;
        if (must_audit_returned) {
            returned = sst_reassembled_returned_state_residual(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
            result.k_residual = returned.k_raw_inf;
            result.second_residual = returned.omega_raw_inf;
        }
        // V18: distinguish the exact owner equation from residuals of the
        // unconstrained PDE. The target is refreshed on the SAME trial k/w
        // pair, including corner area weights. Never count this as an
        // OpenFOAM or experimental BFS qualification certificate.
        if (sources_now() && sources_now()->has_owner_omega_constraints() &&
            std::getenv("CFDX_SST_DIAGNOSTICS") &&
            (iter==1 || iter%100==0 || iter==max_iterations)) {
            const auto audit=sst_wall_owner_audit_v18(omega,sources_now());
            std::cerr << std::scientific << std::setprecision(12)
                      << "SST_WALL_OWNER_CONSTRAINT_V18 iteration=" << iter
                      << " owners=" << audit.owners
                      << " max_absolute_mismatch=" << audit.max_absolute_mismatch
                      << " worst_cell=" << audit.worst_cell
                      << " qualified=0\n";
        }
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
        // V33: dense, strictly read-only forensic sampling of selected FV rows.
        // No timestep, source, wall value, Rhie-Chow flux, convergence gate,
        // or accepted state is modified. The sampled state is the RETURNED
        // state, with the same omega/G wall refresh as the residual above.
        if (!v33_trace_cells.empty()) {
            const auto trace_cross=sst_boundary_gradient_dot(
                k,omega,mesh,geometry,k_bcs,omega_bcs,face_values_now(),
                controls.sst_use_wls_boundary_gradients,
                controls.sst_gradient_dimension);
            for (std::size_t trace_i=0; trace_i<v33_trace_cells.size(); ++trace_i) {
                const auto row=v33_trace_cells[trace_i];
                const auto budgets=sst_reassembled_discrete_cell_budgets(
                    mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                    controls,k_bcs,omega_bcs,face_values_now(),row,sources_now());
                const auto [f1,f2]=compute_sst_blending(
                    k(row),omega(row),wall_distance(row),
                    controls.molecular_viscosity,controls.beta_star,
                    trace_cross[row],controls.density,controls.sst_sigma_w2);
                const double wall_g=sources_now()
                    ? sources_now()->for_cell(row) : -1.;
                const auto branches=sst_v33_source_branches(
                    k(row),omega(row),strain_rate(row),f2,wall_g,
                    controls.density,controls.a1,controls.beta_star,
                    controls.sst_production_limiter,controls.k_min,
                    controls.omega_min);
                const auto production=sst_cell_sources(
                    k(row),omega(row),strain_rate(row),f1,f2,
                    trace_cross[row],controls,wall_g);
                const double source_gap=std::abs(
                    production.sk-branches.limited_production);
                const double nut_gap=std::abs(
                    production.nut-branches.turbulent_viscosity);
                // Independent algebraic check. Fail closed if the diagnostic
                // ceases to represent the actual CFDX source equations.
                if (source_gap > 1e-11*std::max(1.,std::abs(production.sk)) ||
                    nut_gap > 1e-11*std::max(1.,std::abs(production.nut)))
                    throw std::runtime_error(
                        "V33 SST trace: source audit differs from production");
                const bool branch_changed=iter>1 && sst_v33_branch_changed(
                    v33_previous_branches[trace_i],branches);
                v33_previous_branches[trace_i]=branches;
                for (const auto& item : {
                    std::pair<const char*,SSTCellDiscreteBudget>{"k",budgets.first},
                    {"omega",budgets.second}}) {
                    const auto& b=item.second;
                    std::cerr << std::scientific << std::setprecision(16)
                        << "SST_ROW_TRACE_V33 iteration=" << iter
                        << " cell=" << row << " field=" << item.first
                        << " k=" << k(row) << " omega=" << omega(row)
                        << " f1=" << f1 << " f2=" << f2
                        << " convection=" << b.convection
                        << " diffusion=" << b.diffusion
                        << " production=" << b.physical_production
                        << " destruction=" << b.physical_destruction
                        << " cross=" << b.physical_cross_diffusion
                        << " full_residual=" << b.full_residual
                        << " local_relative=" << b.relative_local
                        << " closure_error=" << b.closure_error
                        << " source_identity_error=" << b.nonlinear_source_identity_error
                        << " flux_divergence=" << b.continuity_divergence
                        << " owner_constraint=" << b.owner_constraint_row
                        << " omega_viscosity_branch=" << branches.omega_viscosity_branch
                        << " production_limited=" << branches.production_limited
                        << " wall_g_branch=" << branches.wall_production
                        << " k_floor=" << branches.k_on_floor
                        << " omega_floor=" << branches.omega_on_floor
                        << " viscosity_switch_near=" << branches.close_to_viscosity_switch
                        << " production_switch_near=" << branches.close_to_production_switch
                        << " branch_changed=" << branch_changed
                        << " frozen_source_derivative_valid="
                        << branches.frozen_source_derivative_valid
                        << " dRk_dlogk_source="
                        << branches.d_rk_dlog_k_frozen_source
                        << " dRk_dlogomega_source="
                        << branches.d_rk_dlog_omega_frozen_source
                        << " source_units=per_unit_volume"
                        << " cell_volume=" << geometry.cell_volumes[row]
                        << " dRk_dlogk_fv_source="
                        << branches.d_rk_dlog_k_frozen_source*geometry.cell_volumes[row]
                        << " dRk_dlogomega_fv_source="
                        << branches.d_rk_dlog_omega_frozen_source*geometry.cell_volumes[row]
                        << " winner_accepted=" << line_accepted
                        << " stalls=" << consecutive_line_search_stalls
                        << " qualified=0\n";
                }
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
        // V31 forensic audit: the REAL returned-state FV row 341 and its
        // strongest spatial neighbour. Every central difference refreshes
        // omega/G and reassembles the same production operator; no state is
        // committed and no convergence criterion is modified.
        if (std::getenv("CFDX_SST_K341_DERIVATIVES_V31") && n>341 && n<=4096 &&
            must_audit_returned &&
            (iter==1 || iter==10 || iter==20 || iter==25 || iter==27)) {
            const std::size_t probe=341;
            const auto equations=sst_reassembled_equations_v22(
                mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                controls,k_bcs,omega_bcs,face_values_now(),sources_now());
            const auto* offsets=equations.first.matrix.row_offsets_data();
            const auto* cols=equations.first.matrix.columns_data();
            const auto* vals=equations.first.matrix.values_data();
            std::size_t neighbour=probe;
            double strongest=0.;
            for (auto j=offsets[probe];j<offsets[probe+1];++j)
                if (cols[j]!=probe && std::abs(vals[j])>strongest) {
                    neighbour=cols[j];strongest=std::abs(vals[j]);
                }
            if (neighbour==probe)
                throw std::runtime_error("V31 hotspot: no spatial coupling at cell 341");
            std::vector<double> state(2*n);
            for (std::size_t i=0;i<n;++i) {
                state[i]=k(i);state[n+i]=omega(i);
            }
            const SSTV31BudgetEvaluator evaluator=[&](const auto& trial) {
                auto kk=k,ww=omega;
                for (std::size_t i=0;i<n;++i) {
                    kk(i)=trial[i];ww(i)=trial[n+i];
                }
                ScalarBoundaryFaceValues staging_faces;
                SSTWallProductionV17 staging_prod;
                const ScalarBoundaryFaceValues* current_faces=omega_face_values;
                const SSTWallProductionV17* current_prod=nullptr;
                if (wall_refresh) {
                    (*wall_refresh)(kk,ww,staging_faces,staging_prod);
                    staging_prod.validate(n);
                    current_faces=&staging_faces;
                    current_prod=&staging_prod;
                }
                const auto rows=sst_reassembled_discrete_cell_budgets(
                    mesh,geometry,mass_flux,kk,ww,strain_rate,wall_distance,
                    controls,k_bcs,omega_bcs,current_faces,probe,current_prod);
                const auto& b=rows.first;
                SSTV31RowComponents out;
                out.parts={{b.convection,b.diffusion,b.physical_production,
                            b.physical_destruction,b.physical_cross_diffusion,
                            b.full_residual}};
                out.closure_error=b.closure_error+b.nonlinear_source_identity_error;
                out.constrained=b.owner_constraint_row;
                return out;
            };
            const auto baseline=evaluator(state);
            sst_v31_validate_components(baseline);
            std::cerr<<std::scientific<<std::setprecision(11)
                     <<"SST_K341_V31 iteration="<<iter
                     <<" cell="<<probe<<" strongest_neighbour="<<neighbour
                     <<" residual="<<baseline.parts[5]
                     <<" convection="<<baseline.parts[0]
                     <<" diffusion="<<baseline.parts[1]
                     <<" production="<<baseline.parts[2]
                     <<" destruction="<<baseline.parts[3]
                     <<" cross="<<baseline.parts[4]<<" qualified=0\n";
            for (std::size_t variable : {probe,n+probe,neighbour,n+neighbour}) {
                const auto derivative=sst_v31_log_derivative(
                    state,variable,variable<n?controls.k_min:controls.omega_min,
                    1e-4,evaluator);
                std::cerr<<std::scientific<<std::setprecision(11)
                         <<"SST_K341_JACOBIAN_V31 iteration="<<iter
                         <<" row="<<probe<<" variable="
                         <<(variable<n?"k":"omega")
                         <<" variable_cell="<<(variable%n)
                         <<" conv="<<derivative.d_dlog[0]
                         <<" diff="<<derivative.d_dlog[1]
                         <<" prod="<<derivative.d_dlog[2]
                         <<" dest="<<derivative.d_dlog[3]
                         <<" cross="<<derivative.d_dlog[4]
                         <<" dR="<<derivative.d_dlog[5]
                         <<" closure="<<derivative.maximum_relative_closure
                         <<" qualified=0\n";
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
                      << " picard_damping=" << picard_damping_v41
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
        if (two_cycle_v41 && !result.converged) {
            for (std::size_t i=0;i<n;++i) {
                two_cycle_snapshot_k_v41[i]=k(i);
                two_cycle_snapshot_omega_v41[i]=omega(i);
            }
            const auto evidence=two_cycle_v41->observe(
                two_cycle_snapshot_k_v41,two_cycle_snapshot_omega_v41);
            if (evidence.detected && picard_damping_v41>0.5) {
                picard_damping_v41=0.5;
                std::cerr << std::scientific << std::setprecision(9)
                          << "SST_TWO_CYCLE_RESCUE_V41 iteration=" << iter
                          << " cell=" << evidence.cell
                          << " confirmations=" << static_cast<int>(evidence.consecutive)
                          << " one_step_relative=" << evidence.relative_one_step
                          << " two_step_relative=" << evidence.relative_two_step
                          << " next_picard_damping=" << picard_damping_v41
                          << " tolerance_unchanged=" << tolerance
                          << " physical_qualified=0\n";
            }
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
            // V42: resolve the exact failure row automatically, but only
            // on explicit request and after the fail-closed V14 decision.
            // Distinguish max row-scaled (the strict gate) from max cell-local
            // normalized defect (which can be in a different cell).
            if (const char* forensic = std::getenv("CFDX_SST_STALL_FORENSICS_V42")) {
                if (std::string(forensic) == "1") {
                    const auto equations=sst_reassembled_equations_v22(
                        mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                        controls,k_bcs,omega_bcs,face_values_now(),sources_now());
                    cfdx::core::Vector omega_x(n);
                    for (std::size_t i=0; i<n; ++i) omega_x(i)=omega(i);
                    const auto overview=sst_stall_overview_v42(
                        equations.second,omega_x,controls.omega_min);
                    const auto emit=[&](const char* criterion, const SSTStallRowV42& r) {
                        const auto budget=sst_reassembled_discrete_cell_budgets(
                            mesh,geometry,mass_flux,k,omega,strain_rate,wall_distance,
                            controls,k_bcs,omega_bcs,face_values_now(),r.cell,sources_now());
                        const auto& w=budget.second;
                        std::cerr << std::scientific << std::setprecision(17)
                                  << "SST_STALL_ROW_V42 iteration=" << iter
                                  << " criterion=" << criterion
                                  << " cell=" << r.cell
                                  << " rhs=" << r.rhs
                                  << " diag_term=" << r.diagonal_term
                                  << " offdiag_terms=" << r.off_diagonal_terms
                                  << " native_signed=" << r.native_signed
                                  << " compensated_signed=" << r.compensated_signed
                                  << " row_l1=" << r.row_l1
                                  << " local_field_scale=" << r.local_field_scale
                                  << " native_row_scaled=" << r.native_row_scaled
                                  << " compensated_row_scaled=" << r.compensated_row_scaled
                                  << " native_local=" << r.native_local_relative
                                  << " compensated_local=" << r.compensated_local_relative
                                  << " cancellation_ratio=" << r.cancellation_ratio
                                  << " fv_convection=" << w.convection
                                  << " fv_diffusion=" << w.diffusion
                                  << " fv_production=" << w.physical_production
                                  << " fv_destruction=" << w.physical_destruction
                                  << " fv_cross=" << w.physical_cross_diffusion
                                  << " fv_source_identity=" << w.nonlinear_source_identity_error
                                  << " fv_assembled_residual=" << w.full_residual
                                  << " fv_closure_error=" << w.closure_error
                                  << " wall_orthogonal_flux=" << w.wall_orthogonal_flux
                                  << " wall_owner_diffusion=" << w.wall_owner_diffusion_flux
                                  << " wall_molecular_diffusion=" << w.wall_molecular_diffusion_flux
                                  << " wall_dirichlet_faces=" << w.dirichlet_wall_faces
                                  << " owner_constraint=" << (w.owner_constraint_row?1:0)
                                  << " acceptance_unchanged=1 physical_qualified=0\n";
                    };
                    emit("row_scaled_native",overview.native_row_scaled_worst);
                    if (overview.native_local_worst.cell!=overview.native_row_scaled_worst.cell)
                        emit("local_native",overview.native_local_worst);
                    if (overview.compensated_row_scaled_worst.cell!=overview.native_row_scaled_worst.cell &&
                        overview.compensated_row_scaled_worst.cell!=overview.native_local_worst.cell)
                        emit("row_scaled_compensated",overview.compensated_row_scaled_worst);
                }
            }
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
    // Commit the coupled continuation history only after convergence.
    if (controls.sst_coupled_ptc_v30_experimental && result.converged &&
        controls.sst_continuation_state_v30)
        *controls.sst_continuation_state_v30=ptc_v30;
    return result;
}

} // namespace cfdx::physics

