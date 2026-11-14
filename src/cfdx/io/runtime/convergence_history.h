#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::io {

namespace detail {

inline void convergence_json_string(std::ostream& out, const std::string& value)
{
    out << '"';
    for (const char ch : value) {
        switch (ch) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (static_cast<unsigned char>(ch) < 0x20)
                    throw std::invalid_argument("convergence history contains a control character");
                out << ch;
        }
    }
    out << '"';
}

inline void convergence_json_number(std::ostream& out, double value)
{
    if (!std::isfinite(value))
        throw std::invalid_argument("convergence history contains a non-finite value");
    out << value;
}

// The V8 certificate may be disabled or not evaluated after a rejected step.
// JSON null is evidence of an absent number, never a numerical zero.
inline void convergence_json_optional_number(std::ostream& out, double value)
{
    if (std::isfinite(value)) convergence_json_number(out, value);
    else out << "null";
}

inline const char* convergence_algorithm_name(
    cfdx::physics::PressureVelocityAlgorithm algorithm)
{
    using cfdx::physics::PressureVelocityAlgorithm;
    switch (algorithm) {
        case PressureVelocityAlgorithm::SIMPLE: return "simple";
        case PressureVelocityAlgorithm::SIMPLEC: return "simplec";
        case PressureVelocityAlgorithm::PISO: return "piso";
        case PressureVelocityAlgorithm::PIMPLE: return "pimple";
        case PressureVelocityAlgorithm::FRACTIONAL_STEP: return "fractional_step";
        case PressureVelocityAlgorithm::COUPLED: return "coupled";
    }
    return "unknown";
}

inline void convergence_json_linear_plan(
    std::ostream& out, const cfdx::core::LinearSolverPlan& plan, bool resolved)
{
    out << "{\"resolved\": " << (resolved ? "true" : "false");
    if (resolved) {
        out << ", \"krylov\": ";
        convergence_json_string(out, cfdx::core::to_string(plan.krylov));
        out << ", \"preconditioner\": ";
        convergence_json_string(out, cfdx::core::to_string(plan.preconditioner));
        out << ", \"null_space\": ";
        convergence_json_string(out, cfdx::core::to_string(plan.null_space));
        out << ", \"automatic_krylov\": "
            << (plan.automatic_krylov ? "true" : "false");
        out << ", \"automatic_preconditioner\": "
            << (plan.automatic_preconditioner ? "true" : "false");
        out << ", \"reason\": ";
        convergence_json_string(out, plan.reason);
    }
    out << '}';
}

inline void convergence_json_v8(
    std::ostream& out, const cfdx::physics::GlobalNonlinearResidualV8& v)
{
    const bool evaluated = std::isfinite(v.continuity) || v.sst;
    out << "{\"evaluated\": " << (evaluated ? "true" : "false")
        << ", \"algorithm\": ";
    convergence_json_string(out, convergence_algorithm_name(v.algorithm));
    out << ", \"schur\": ";
    if (v.algorithm == cfdx::physics::PressureVelocityAlgorithm::COUPLED)
        convergence_json_string(out, cfdx::physics::to_string(v.schur));
    else out << "null";
    out << ", \"momentum_global_relative\": ";
    convergence_json_optional_number(out, v.momentum.global_relative);
    out << ", \"momentum_max_active_component_relative\": ";
    convergence_json_optional_number(out, v.momentum.max_active_component_relative);
    out << ", \"momentum_component_relative\": [";
    for (std::size_t i=0; i<3; ++i) {
        if (i) out << ',';
        convergence_json_optional_number(out, v.momentum.component_relative[i]);
    }
    out << "], \"continuity\": ";
    convergence_json_optional_number(out, v.continuity);
    out << ", \"k\": "; convergence_json_optional_number(out, v.k);
    out << ", \"omega\": "; convergence_json_optional_number(out, v.omega);
    out << ", \"eddy_viscosity_closure\": ";
    convergence_json_optional_number(out, v.eddy_viscosity_closure);
    out << ", \"merit\": "; convergence_json_optional_number(out, v.merit);
    out << ", \"sst\": " << (v.sst ? "true" : "false")
        << ", \"sst_transport_committed\": "
        << (v.sst_transport_committed ? "true" : "false")
        << ", \"bounded\": " << (v.bounded ? "true" : "false")
        << ", \"all_finite\": " << (v.all_finite ? "true" : "false")
        << ", \"converged\": " << (v.converged ? "true" : "false") << '}';
}

inline void convergence_json_array3(
    std::ostream& out, const std::array<double, 3>& values)
{
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        convergence_json_number(out, values[i]);
    }
    out << ']';
}

inline void convergence_json_size_array3(
    std::ostream& out, const std::array<std::size_t, 3>& values)
{
    out << '[' << values[0] << ',' << values[1] << ',' << values[2] << ']';
}

} // namespace detail

// Persist the authoritative nonlinear convergence history produced by the
// incompressible solver. This is execution evidence, not a validation verdict.
// The artifact is deliberately separate from CFDX-DAT state/checkpoint files.
inline void write_convergence_history(
    const std::string& path,
    const cfdx::physics::IncompressibleSolveResult& result)
{
    using namespace detail;

    if (!std::isfinite(result.reference_momentum_residual))
        throw std::invalid_argument(
            "write_convergence_history: non-finite reference residual");

    const std::string tmp_path = path + ".tmp";
    std::ofstream out(tmp_path, std::ios::trunc);
    if (!out)
        throw std::runtime_error(
            "write_convergence_history: cannot open " + path);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    try {
        out << "{\n";
        out << "  \"format\": \"CFDX-CONVERGENCE\",\n";
        out << "  \"schema_version\": 1,\n";
        out << "  \"converged\": " << (result.converged ? "true" : "false") << ",\n";
        out << "  \"iterations\": " << result.iterations << ",\n";
        out << "  \"reference_momentum_residual\": ";
        convergence_json_number(out, result.reference_momentum_residual);
        out << ",\n";
        out << "  \"convergence_status\": "
            << static_cast<int>(result.convergence_status) << ",\n";
        out << "  \"convergence_reason\": ";
        convergence_json_string(out, result.convergence_reason);
        out << ",\n";
        out << "  \"linear_plan\": {\n";
        out << "    \"schur_resolved\": "
            << (result.coupled_schur_model_resolved ? "true" : "false") << ",\n";
        out << "    \"requested_schur\": ";
        convergence_json_string(out, cfdx::physics::to_string(result.requested_coupled_schur_model));
        out << ",\n    \"resolved_schur\": ";
        if (result.coupled_schur_model_resolved)
            convergence_json_string(out, cfdx::physics::to_string(result.resolved_coupled_schur_model));
        else out << "null";
        out << ",\n    \"pressure\": ";
        convergence_json_linear_plan(out, result.pressure_linear_plan,
                                      result.pressure_linear_plan_resolved);
        out << ",\n    \"coupled\": ";
        convergence_json_linear_plan(out, result.coupled_linear_plan,
                                      result.coupled_linear_plan_resolved);
        out << "\n  },\n";
        out << "  \"history\": [\n";

        for (std::size_t i = 0; i < result.history.size(); ++i) {
            const auto& h = result.history[i];
            out << "    {\n";
            out << "      \"iteration\": " << h.iteration << ",\n";
            out << "      \"global_nonlinear_v8\": ";
            convergence_json_v8(out, h.global_nonlinear_v8);
            out << ",\n";
            out << "      \"high_order_fraction\": "; convergence_json_number(out, h.high_order_fraction); out << ",\n";
            out << "      \"full_order_eligible\": " << (h.full_order_eligible ? "true" : "false") << ",\n";
            out << "      \"momentum_residual\": "; convergence_json_number(out, h.momentum_residual); out << ",\n";
            out << "      \"pressure_residual\": "; convergence_json_number(out, h.pressure_residual); out << ",\n";
            out << "      \"continuity_l1\": "; convergence_json_number(out, h.continuity_l1); out << ",\n";
            out << "      \"continuity_linf\": "; convergence_json_number(out, h.continuity_linf); out << ",\n";
            out << "      \"continuity_normalized\": "; convergence_json_number(out, h.continuity_normalized); out << ",\n";
            out << "      \"momentum_equation_residual\": "; convergence_json_number(out, h.momentum_equation_residual); out << ",\n";
            out << "      \"momentum_rhs_scale\": "; convergence_json_number(out, h.momentum_rhs_scale); out << ",\n";
            out << "      \"momentum_equation_residual_reference\": "; convergence_json_number(out, h.momentum_equation_residual_reference); out << ",\n";
            out << "      \"momentum_equation_residual_rhs_relative\": "; convergence_json_number(out, h.momentum_equation_residual_rhs_relative); out << ",\n";
            out << "      \"momentum_equation_residual_rhs_active_component_max\": "; convergence_json_number(out, h.momentum_equation_residual_rhs_active_component_max); out << ",\n";
            out << "      \"momentum_rhs_scale_components\": "; convergence_json_array3(out, h.momentum_rhs_scale_components); out << ",\n";
            out << "      \"momentum_equation_residual_rhs_relative_components\": "; convergence_json_array3(out, h.momentum_equation_residual_rhs_relative_components); out << ",\n";
            out << "      \"momentum_equation_residual_relative\": "; convergence_json_number(out, h.momentum_equation_residual_relative); out << ",\n";
            out << "      \"velocity_change_inf\": "; convergence_json_number(out, h.velocity_change_inf); out << ",\n";
            out << "      \"pressure_change_inf\": "; convergence_json_number(out, h.pressure_change_inf); out << ",\n";
            out << "      \"effective_alpha_u\": "; convergence_json_number(out, h.effective_alpha_u); out << ",\n";
            out << "      \"effective_alpha_p\": "; convergence_json_number(out, h.effective_alpha_p); out << ",\n";
            out << "      \"nonlinear_convergence_metric\": "; convergence_json_number(out, h.nonlinear_convergence_metric); out << ",\n";
            out << "      \"momentum_linear_iterations\": " << h.momentum_linear_iterations << ",\n";
            out << "      \"pressure_linear_iterations\": " << h.pressure_linear_iterations << ",\n";
            out << "      \"pressure_correctors_used\": " << h.pressure_correctors_used << ",\n";
            out << "      \"linear_tolerance_used\": "; convergence_json_number(out, h.linear_tolerance_used); out << ",\n";
            out << "      \"corrected_flux_continuity_linf\": "; convergence_json_number(out, h.corrected_flux_continuity_linf); out << ",\n";
            out << "      \"reconstructed_velocity_continuity_linf\": "; convergence_json_number(out, h.reconstructed_velocity_continuity_linf); out << ",\n";
            out << "      \"flux_velocity_mismatch_linf\": "; convergence_json_number(out, h.flux_velocity_mismatch_linf); out << ",\n";
            out << "      \"flux_velocity_mismatch_boundary_linf\": "; convergence_json_number(out, h.flux_velocity_mismatch_boundary_linf); out << ",\n";
            out << "      \"flux_velocity_mismatch_normalized\": "; convergence_json_number(out, h.flux_velocity_mismatch_normalized); out << ",\n";
            out << "      \"rhie_chow_identity_linf\": "; convergence_json_number(out, h.rhie_chow_identity_linf); out << ",\n";
            out << "      \"simplec_signed_offdiag_sum_linf\": [" << h.simplec_signed_offdiag_sum_linf[0] << "," << h.simplec_signed_offdiag_sum_linf[1] << "," << h.simplec_signed_offdiag_sum_linf[2] << "],\n";
            out << "      \"simplec_abs_offdiag_sum_linf\": [" << h.simplec_abs_offdiag_sum_linf[0] << "," << h.simplec_abs_offdiag_sum_linf[1] << "," << h.simplec_abs_offdiag_sum_linf[2] << "],\n";
            out << "      \"simplec_denominator_signed_linf\": [" << h.simplec_denominator_signed_linf[0] << "," << h.simplec_denominator_signed_linf[1] << "," << h.simplec_denominator_signed_linf[2] << "],\n";
            out << "      \"simplec_denominator_abs_linf\": [" << h.simplec_denominator_abs_linf[0] << "," << h.simplec_denominator_abs_linf[1] << "," << h.simplec_denominator_abs_linf[2] << "],\n";
            out << "      \"simplec_denominator_gap_linf\": [" << h.simplec_denominator_gap_linf[0] << "," << h.simplec_denominator_gap_linf[1] << "," << h.simplec_denominator_gap_linf[2] << "],\n";
            out << "      \"simplec_denominator_worst_cell\": [" << h.simplec_denominator_worst_cell[0] << "," << h.simplec_denominator_worst_cell[1] << "," << h.simplec_denominator_worst_cell[2] << "],\n";
            out << "      \"momentum_equation_residual_components\": "; convergence_json_array3(out, h.momentum_equation_residual_components); out << ",\n";
            out << "      \"momentum_equation_residual_internal\": "; convergence_json_number(out, h.momentum_equation_residual_internal); out << ",\n";
            out << "      \"momentum_equation_residual_boundary\": "; convergence_json_number(out, h.momentum_equation_residual_boundary); out << ",\n";
            out << "      \"pressure_gradient_linf\": "; convergence_json_number(out, h.pressure_gradient_linf); out << ",\n";
            out << "      \"pressure_gradient_l2\": "; convergence_json_number(out, h.pressure_gradient_l2); out << ",\n";
            out << "      \"momentum_residual_cell\": " << h.momentum_residual_cell << ",\n";
            out << "      \"momentum_residual_no_pressure\": "; convergence_json_number(out, h.momentum_residual_no_pressure); out << ",\n";
            out << "      \"momentum_pressure_contribution\": "; convergence_json_number(out, h.momentum_pressure_contribution); out << ",\n";
            out << "      \"momentum_residual_patch\": "; convergence_json_string(out, h.momentum_residual_patch); out << ",\n";
            out << "      \"mass_boundary_flux\": "; convergence_json_number(out, h.mass_boundary_flux); out << ",\n";
            out << "      \"mass_global_cell_balance\": "; convergence_json_number(out, h.mass_global_cell_balance); out << ",\n";
            out << "      \"mass_local_l1\": "; convergence_json_number(out, h.mass_local_l1); out << ",\n";
            out << "      \"mass_local_linf\": "; convergence_json_number(out, h.mass_local_linf); out << ",\n";
            out << "      \"mass_local_l2\": "; convergence_json_number(out, h.mass_local_l2); out << ",\n";
            out << "      \"mass_normalized_imbalance\": "; convergence_json_number(out, h.mass_normalized_imbalance); out << ",\n";
            out << "      \"mass_worst_cell\": " << h.mass_worst_cell << ",\n";
            out << "      \"mass_nonfinite_faces\": " << h.mass_nonfinite_faces << ",\n";
            out << "      \"boundedness_nonfinite_velocity\": " << h.boundedness_nonfinite_velocity << ",\n";
            out << "      \"momentum_conservation_residual\": "; convergence_json_array3(out, h.momentum_conservation_residual); out << ",\n";
            out << "      \"momentum_conservation_normalized\": "; convergence_json_array3(out, h.momentum_conservation_normalized); out << ",\n";
            out << "      \"momentum_conservation_worst_cell\": "; convergence_json_size_array3(out, h.momentum_conservation_worst_cell); out << "\n";
            out << "    }";
            if (i + 1 != result.history.size()) out << ',';
            out << "\n";
        }

        out << "  ]\n";
        out << "}\n";
    } catch (...) {
        out.close();
        std::remove(tmp_path.c_str());
        throw;
    }

    if (!out) {
        out.close();
        std::remove(tmp_path.c_str());
        throw std::runtime_error(
            "write_convergence_history: write failed for " + path);
    }
    out.close();
    if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
        std::remove(tmp_path.c_str());
        throw std::runtime_error(
            "write_convergence_history: cannot replace " + path);
    }
}

} // namespace cfdx::io
