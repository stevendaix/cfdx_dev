#include "cfdx/io/runtime/convergence_history.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

int main()
{
    using namespace cfdx::physics;
    using cfdx::io::write_convergence_history;

    IncompressibleSolveResult result;
    result.converged = true;
    result.iterations = 1;
    result.reference_momentum_residual = 2.0;
    result.convergence_reason = "converged \"normally\"";
    IncompressibleIteration h;
    h.iteration = 1;
    h.momentum_residual = 1.0;
    h.pressure_residual = 0.5;
    h.continuity_l1 = 0.1;
    h.continuity_linf = 0.2;
    h.continuity_normalized = 0.3;
    h.momentum_equation_residual = 0.4;
    h.momentum_equation_residual_relative = 0.5;
    h.velocity_change_inf = 0.01;
    h.pressure_change_inf = 0.02;
    h.effective_alpha_u = 0.7;
    h.effective_alpha_p = 0.8;
    h.nonlinear_convergence_metric = 0.5;
    h.momentum_linear_iterations = 11;
    h.pressure_linear_iterations = 12;
    h.pressure_correctors_used = 2;
    h.linear_tolerance_used = 1e-8;
    h.corrected_flux_continuity_linf = 0.01;
    h.reconstructed_velocity_continuity_linf = 0.02;
    h.flux_velocity_mismatch_linf = 0.03;
    h.momentum_equation_residual_components = {0.1, 0.2, 0.3};
    h.momentum_equation_residual_internal = 0.4;
    h.momentum_equation_residual_boundary = 0.5;
    h.pressure_gradient_linf = 0.6;
    h.pressure_gradient_l2 = 0.7;
    h.momentum_residual_cell = 4;
    h.momentum_residual_no_pressure = 0.8;
    h.momentum_pressure_contribution = 0.9;
    h.momentum_residual_patch = "wall";
    h.mass_boundary_flux = 1.0;
    h.mass_global_cell_balance = 1.1;
    h.mass_local_l1 = 1.2;
    h.mass_local_linf = 1.3;
    h.mass_local_l2 = 1.4;
    h.mass_normalized_imbalance = 1.5;
    h.mass_worst_cell = 5;
    h.mass_nonfinite_faces = 0;
    h.boundedness_nonfinite_velocity = 0;
    h.momentum_conservation_residual = {0.2, 0.3, 0.4};
    h.momentum_conservation_normalized = {0.5, 0.6, 0.7};
    h.momentum_conservation_worst_cell = {1, 2, 3};
    result.history.push_back(h);

    const auto path = std::filesystem::temp_directory_path() /
                      "cfdx_convergence_history_test.json";
    std::error_code ec;
    std::filesystem::remove(path, ec);
    write_convergence_history(path.string(), result);

    std::ifstream in(path);
    assert(in.good());
    std::stringstream buffer;
    buffer << in.rdbuf();
    const std::string json = buffer.str();
    assert(json.find("\"format\": \"CFDX-CONVERGENCE\"") != std::string::npos);
    assert(json.find("\"schema_version\": 1") != std::string::npos);
    assert(json.find("\"history\": [") != std::string::npos);
    assert(json.find("\"momentum_residual\": 1") != std::string::npos);
    assert(json.find("\"momentum_conservation_worst_cell\": [1,2,3]") != std::string::npos);
    assert(json.find("\\\"normally\\\"") != std::string::npos);

    IncompressibleSolveResult bad = result;
    bad.history[0].momentum_residual = std::numeric_limits<double>::quiet_NaN();
    bool rejected = false;
    try {
        write_convergence_history((path.string() + ".bad"), bad);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    std::filesystem::remove(path, ec);
    std::filesystem::remove(path.string() + ".bad", ec);
    return 0;
}
