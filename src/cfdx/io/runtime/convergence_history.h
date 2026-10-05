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
        out << "  \"history\": [\n";

        for (std::size_t i = 0; i < result.history.size(); ++i) {
            const auto& h = result.history[i];
            out << "    {\n";
            out << "      \"iteration\": " << h.iteration << ",\n";
            out << "      \"momentum_residual\": "; convergence_json_number(out, h.momentum_residual); out << ",\n";
            out << "      \"pressure_residual\": "; convergence_json_number(out, h.pressure_residual); out << ",\n";
            out << "      \"continuity_l1\": "; convergence_json_number(out, h.continuity_l1); out << ",\n";
            out << "      \"continuity_linf\": "; convergence_json_number(out, h.continuity_linf); out << ",\n";
            out << "      \"continuity_normalized\": "; convergence_json_number(out, h.continuity_normalized); out << ",\n";
            out << "      \"momentum_equation_residual\": "; convergence_json_number(out, h.momentum_equation_residual); out << ",\n";
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
