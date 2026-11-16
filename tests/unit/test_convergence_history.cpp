#include "cfdx/io/runtime/convergence_history.h"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include "common/live_assert.hpp"

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
    h.high_order_fraction = 1.0;
    h.full_order_eligible = true;
    h.momentum_equation_residual_rhs_active_component_max = 0.125;
    h.momentum_residual = 1.0;
    h.pressure_residual = 0.5;
    h.continuity_l1 = 0.1;
    h.continuity_linf = 0.2;
    h.continuity_normalized = 0.3;
    h.momentum_equation_residual = 0.4;
    h.momentum_rhs_scale = 0.5;
    h.momentum_equation_residual_reference = 0.25;
    h.momentum_equation_residual_rhs_relative = 0.125;
    h.momentum_rhs_scale_components = {0.5, 0.25, 0.125};
    h.momentum_equation_residual_rhs_relative_components = {0.125, 0.25, 0.5};
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
    h.flux_velocity_mismatch_boundary_linf = 0.3125;
    h.flux_velocity_mismatch_normalized = 0.0625;
    h.rhie_chow_identity_linf = 0.03125;
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
    // Persist the *resolved* Schur/linear plan and independent returned-state
    // certificate, not merely the requested numerical configuration.
    result.requested_coupled_schur_model = CoupledSchurModel::PCD;
    result.resolved_coupled_schur_model = CoupledSchurModel::PCD;
    result.coupled_schur_model_resolved = true;
    result.coupled_linear_plan_resolved = true;
    result.coupled_linear_plan.krylov = cfdx::core::KrylovModel::FGMRES;
    result.coupled_linear_plan.preconditioner = cfdx::core::PreconditionerModel::PCD;
    h.global_nonlinear_v8.algorithm = PressureVelocityAlgorithm::COUPLED;
    h.global_nonlinear_v8.schur = CoupledSchurModel::PCD;
    h.global_nonlinear_v8.momentum.global_relative = 1e-10;
    h.global_nonlinear_v8.momentum.max_active_component_relative = 2e-10;
    h.global_nonlinear_v8.momentum.component_relative = {1e-10, 2e-10, 0};
    h.global_nonlinear_v8.continuity = 3e-10;
    h.global_nonlinear_v8.k = 4e-10;
    h.global_nonlinear_v8.omega = 5e-10;
    h.global_nonlinear_v8.eddy_viscosity_closure = 6e-10;
    h.global_nonlinear_v8.merit = 6e-10;
    h.global_nonlinear_v8.sst = true;
    h.global_nonlinear_v8.sst_transport_committed = true;
    h.global_nonlinear_v8.bounded = true;
    h.global_nonlinear_v8.all_finite = true;
    h.global_nonlinear_v8.converged = true;
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
    assert(json.find("\"requested_schur\": \"pcd\"") != std::string::npos);
    assert(json.find("\"resolved_schur\": \"pcd\"") != std::string::npos);
    assert(json.find("\"preconditioner\": \"pcd\"") != std::string::npos);
    assert(json.find("\"global_nonlinear_v8\": {") != std::string::npos);
    assert(json.find("\"momentum_max_active_component_relative\": 2") != std::string::npos);
    assert(json.find("\"sst_transport_committed\": true") != std::string::npos);
    assert(json.find("\"schema_version\": 1") != std::string::npos);
    assert(json.find("\"history\": [") != std::string::npos);
    assert(json.find("\"high_order_fraction\": 1") != std::string::npos);
    assert(json.find("\"full_order_eligible\": true") != std::string::npos);
    assert(json.find("\"momentum_equation_residual_rhs_active_component_max\": 0.125") != std::string::npos);
    assert(json.find("\"momentum_residual\": 1") != std::string::npos);
    assert(json.find("\"momentum_rhs_scale\": 0.5") != std::string::npos);
    assert(json.find("\"momentum_equation_residual_reference\": 0.25") != std::string::npos);
    assert(json.find("\"momentum_equation_residual_rhs_relative\": 0.125") != std::string::npos);
    assert(json.find("\"momentum_rhs_scale_components\": [0.5,0.25,0.125]") != std::string::npos);
    assert(json.find("\"flux_velocity_mismatch_boundary_linf\": 0.3125") != std::string::npos);
    assert(json.find("\"flux_velocity_mismatch_normalized\": 0.0625") != std::string::npos);
    assert(json.find("\"rhie_chow_identity_linf\": 0.03125") != std::string::npos);
    assert(json.find("\"momentum_equation_residual_rhs_relative_components\": [0.125,0.25,0.5]") != std::string::npos);
    assert(json.find("\"momentum_conservation_worst_cell\": [1,2,3]") != std::string::npos);
    assert(json.find("\\\"normally\\\"") != std::string::npos);

    // If V8 was not evaluated, nulls must stay null rather than becoming
    // a forged 0.0 certificate or invalid NaN JSON.
    IncompressibleSolveResult missing = result;
    missing.history[0].global_nonlinear_v8 = GlobalNonlinearResidualV8{};
    write_convergence_history(path.string() + ".missing", missing);
    std::ifstream missing_stream(path.string() + ".missing");
    std::stringstream missing_content;
    missing_content << missing_stream.rdbuf();
    assert(missing_content.str().find("\"evaluated\": false") != std::string::npos);
    assert(missing_content.str().find("\"k\": null") != std::string::npos);
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
    std::filesystem::remove(path.string() + ".missing", ec);
    return 0;
}
