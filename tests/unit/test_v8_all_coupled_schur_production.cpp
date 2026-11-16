#include "../audit/sst_channel_fixture.h"
#include "common/test_harness.h"
#include <array>
#include <cmath>
#include <iostream>
#include <string>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    const std::array<CoupledSchurModel,6> family{{
        CoupledSchurModel::BlockLocal, CoupledSchurModel::PCD,
        CoupledSchurModel::LSC, CoupledSchurModel::BFBT,
        CoupledSchurModel::SIMPLE, CoupledSchurModel::SIMPLEC}};
    for (const auto model: family) {
        run_case(std::string("v8_production_coupled_sst_")+to_string(model), [model] {
            SSTChannelCase c;
            c.controls.algorithm = PressureVelocityAlgorithm::COUPLED;
            c.controls.coupling.schur_model = model;
            c.controls.require_global_nonlinear_certificate_v8 = true;
            c.runtime.max_transport_iterations = 30;
            c.ubc["wall"].value = {0.0, 0.0, 0.0}; // nonuniform shear challenges all Schur paths
            c.controls.coupling.coupled_linear_tolerance = 1e-12;
            c.controls.coupling.coupled_max_iterations = 200;
            const auto result = solve_steady_incompressible(
                c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls);
            EXPECT_TRUE(result.resolved_coupled_schur_model == model);
            EXPECT_TRUE(result.coupled_linear_plan_resolved);
            std::cout << "selected_plan="
                      << cfdx::core::to_string(result.coupled_linear_plan.krylov)
                      << "/" << cfdx::core::to_string(result.coupled_linear_plan.preconditioner)
                      << " requested_schur=" << to_string(model) << '\n';
            EXPECT_TRUE(result.history.size() == 1);
            const auto& report = result.history.back().global_nonlinear_v8;
            EXPECT_TRUE(report.all_finite);
            EXPECT_TRUE(report.sst);
            EXPECT_TRUE(report.schur == model);
            EXPECT_TRUE(report.continuity >= 0.);
            EXPECT_TRUE(report.k >= 0.);
            EXPECT_TRUE(report.omega >= 0.);
            if (result.converged) EXPECT_TRUE(report.converged);
            std::cout << to_string(model) << " U=" <<
                report.momentum.max_active_component_relative
                      << " p=" << report.continuity
                      << " k=" << report.k
                      << " omega=" << report.omega
                      << " accepted=" << report.converged << '\n';
        });
    }
    // Exercise the Schur+AMG path by selecting its actual linear plan, not
    // merely the Schur label. All six Schur approximations must remain usable.
    for (const auto model : family) {
        run_case(std::string("v8_explicit_block_schur_amg_sst_")+to_string(model), [model] {
            SSTChannelCase c;
            c.controls.algorithm = PressureVelocityAlgorithm::COUPLED;
            c.controls.coupling.schur_model = model;
            c.controls.require_global_nonlinear_certificate_v8 = true;
            c.controls.coupled_linear_solver.krylov = cfdx::core::KrylovModel::FGMRES;
            c.controls.coupled_linear_solver.preconditioner =
                cfdx::core::PreconditionerModel::CoupledBlockSchur;
            c.controls.coupled_linear_solver.gmres_restart = 64;
            c.runtime.max_transport_iterations = 30;
            c.ubc["wall"].value = {0.,0.,0.};
            c.controls.coupling.coupled_max_iterations = 200;
            const auto result=solve_steady_incompressible(
                c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
            EXPECT_TRUE(result.coupled_linear_plan_resolved);
            EXPECT_TRUE(result.coupled_linear_plan.preconditioner ==
                cfdx::core::PreconditionerModel::CoupledBlockSchur);
            EXPECT_TRUE(result.resolved_coupled_schur_model == model);
            EXPECT_TRUE(result.history.size()==1);
            EXPECT_TRUE(result.history.back().global_nonlinear_v8.all_finite);
            if (result.converged) EXPECT_TRUE(result.history.back().global_nonlinear_v8.converged);
            std::cout << "explicit_block_schur=" << to_string(model)
                      << " plan=" << cfdx::core::to_string(result.coupled_linear_plan.preconditioner)
                      << " local=" << result.history.back().global_nonlinear_v8.merit << '\n';
        });
    }
    // A distinct AMG-based block preconditioner, valid for BlockLocal only.
    run_case("v8_explicit_braess_sarazin_block_local_sst", [] {
        SSTChannelCase c;
        c.controls.algorithm = PressureVelocityAlgorithm::COUPLED;
        c.controls.coupling.schur_model = CoupledSchurModel::BlockLocal;
        c.controls.require_global_nonlinear_certificate_v8 = true;
        c.controls.coupled_linear_solver.krylov = cfdx::core::KrylovModel::FGMRES;
        c.controls.coupled_linear_solver.preconditioner =
            cfdx::core::PreconditionerModel::BraessSarazin;
        c.controls.coupled_linear_solver.gmres_restart = 64;
        c.runtime.max_transport_iterations = 30;
        c.ubc["wall"].value={0.,0.,0.};
        c.controls.coupling.coupled_max_iterations=200;
        const auto result=solve_steady_incompressible(
            c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
        EXPECT_TRUE(result.coupled_linear_plan_resolved);
        EXPECT_TRUE(result.coupled_linear_plan.preconditioner ==
            cfdx::core::PreconditionerModel::BraessSarazin);
        EXPECT_TRUE(result.history.size()==1);
        EXPECT_TRUE(result.history.back().global_nonlinear_v8.all_finite);
    });
    // Native MGR and opt-in overlapping Schwarz are separate actual
    // preconditioner families, not aliases for a Schur approximation.
    run_case("v8_explicit_mgr_block_local_sst", [] {
        SSTChannelCase c;
        c.controls.algorithm=PressureVelocityAlgorithm::COUPLED;
        c.controls.coupling.schur_model=CoupledSchurModel::BlockLocal;
        c.controls.require_global_nonlinear_certificate_v8=true;
        c.controls.coupled_linear_solver.krylov=cfdx::core::KrylovModel::FGMRES;
        c.controls.coupled_linear_solver.preconditioner=cfdx::core::PreconditionerModel::MGR;
        c.controls.coupled_linear_solver.gmres_restart=64;
        c.runtime.max_transport_iterations=30;
        c.ubc["wall"].value={0.,0.,0.};
        c.controls.coupling.coupled_max_iterations=200;
        const auto result=solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
        EXPECT_TRUE(result.coupled_linear_plan_resolved);
        EXPECT_TRUE(result.coupled_linear_plan.preconditioner ==
            cfdx::core::PreconditionerModel::MGR);
        EXPECT_TRUE(result.history.size()==1);
        EXPECT_TRUE(result.history.back().global_nonlinear_v8.all_finite);
    });
    run_case("v8_opt_in_schwarz_block_local_sst", [] {
        SSTChannelCase c;
        c.controls.algorithm=PressureVelocityAlgorithm::COUPLED;
        c.controls.coupling.schur_model=CoupledSchurModel::BlockLocal;
        c.controls.require_global_nonlinear_certificate_v8=true;
        c.controls.coupled_linear_solver.krylov=cfdx::core::KrylovModel::FGMRES;
        c.controls.coupled_linear_solver.preconditioner=
            cfdx::core::PreconditionerModel::CoupledBlockSchur;
        c.controls.coupling.coupled_schwarz_experimental=true;
        c.controls.coupled_linear_solver.gmres_restart=64;
        c.runtime.max_transport_iterations=30;
        c.ubc["wall"].value={0.,0.,0.};
        c.controls.coupling.coupled_max_iterations=200;
        const auto result=solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
        EXPECT_TRUE(result.history.size()==1);
        EXPECT_TRUE(result.history.back().global_nonlinear_v8.all_finite);
    });
    return run_all();
}
