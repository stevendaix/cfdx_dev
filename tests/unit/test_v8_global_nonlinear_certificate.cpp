#include "cfdx/physics/global_nonlinear_residual_v8.h"
#include "../audit/sst_channel_fixture.h"
#include "common/test_harness.h"
#include <array>
#include <cmath>
#include <iostream>
#include <string>

using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    const std::array<CoupledSchurModel, 6> all_schur{{
        CoupledSchurModel::BlockLocal, CoupledSchurModel::PCD,
        CoupledSchurModel::LSC, CoupledSchurModel::BFBT,
        CoupledSchurModel::SIMPLE, CoupledSchurModel::SIMPLEC}};
    ConvergenceCriteria criteria;
    criteria.relative_tolerance = 1e-8;
    criteria.continuity_tolerance = 1e-11;
    for (const auto schur : all_schur) {
        const auto name = std::string("v8_true_residual_invariant_to_schur_") + to_string(schur);
        run_case(name, [&,schur] {
            const auto momentum = certify_momentum_residual_components(
                {1e-12, 2e-12, 0.}, {1., 1., 0.}, criteria.relative_tolerance);
            const auto fine = certify_global_nonlinear_residual_v8(
                PressureVelocityAlgorithm::COUPLED, schur, momentum,
                1e-12, criteria, true, 2e-12, 3e-12, 1e-14, true, true);
            EXPECT_TRUE(fine.all_finite);
            EXPECT_TRUE(fine.converged);
            EXPECT_TRUE(fine.schur == schur);
            EXPECT_TRUE(std::abs(fine.merit-3e-12) < 1e-25);
            const auto fail_k = certify_global_nonlinear_residual_v8(
                PressureVelocityAlgorithm::COUPLED, schur, momentum,
                1e-12, criteria, true, 1e-5, 3e-12, 1e-14, true, true);
            EXPECT_TRUE(!fail_k.converged);
            const auto fail_mass = certify_global_nonlinear_residual_v8(
                PressureVelocityAlgorithm::COUPLED, schur, momentum,
                1e-10, criteria, false);
            EXPECT_TRUE(!fail_mass.converged);
            const auto fail_provisional = certify_global_nonlinear_residual_v8(
                PressureVelocityAlgorithm::COUPLED, schur, momentum,
                1e-12, criteria, true, 0., 0., 0., false, true);
            EXPECT_TRUE(!fail_provisional.converged);
            const auto fail_floor = certify_global_nonlinear_residual_v8(
                PressureVelocityAlgorithm::COUPLED, schur, momentum,
                1e-12, criteria, true, 0., 0., 0., true, false);
            EXPECT_TRUE(!fail_floor.converged);
            const auto fail_missing = certify_global_nonlinear_residual_v8(
                PressureVelocityAlgorithm::COUPLED, schur, momentum,
                1e-12, criteria, true);
            EXPECT_TRUE(!fail_missing.all_finite);
        });
    }
    run_case("v8_full_production_sst_operator_on_accepted_state", [] {
        SSTChannelCase c;
        c.controls.require_global_nonlinear_certificate_v8 = true;
        c.runtime.max_transport_iterations = 30;
        const auto result = solve_steady_incompressible(
            c.mesh, c.U, c.p, c.ubc, c.pbc, c.controls);
        EXPECT_TRUE(result.history.size() == 1);
        const auto& certificate = result.history.back().global_nonlinear_v8;
        EXPECT_TRUE(certificate.sst);
        EXPECT_TRUE(certificate.all_finite);
        EXPECT_TRUE(certificate.continuity >= 0.);
        EXPECT_TRUE(certificate.k >= 0.);
        EXPECT_TRUE(certificate.omega >= 0.);
        EXPECT_TRUE(std::isfinite(certificate.eddy_viscosity_closure));
        std::cout << "v8_production_sst U=" <<
            certificate.momentum.max_active_component_relative
                  << " p=" << certificate.continuity
                  << " k=" << certificate.k
                  << " omega=" << certificate.omega
                  << " closure=" << certificate.eddy_viscosity_closure
                  << " converged=" << result.converged << '\n';
        if (result.converged) EXPECT_TRUE(certificate.converged);
    });
    return run_all();
}
