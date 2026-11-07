#include "cfdx/physics/coupled_flux_equivalence_certificate.h"
#include "common/test_harness.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    // Two cells: one internal face of +1, inlet boundary -1, outlet +1.
    // Each cell has exactly zero flux divergence. 8x8 linear system pressure
    // rows explicitly represent the same flux divergence at the chosen state.
    SparseMatrix A(8,8);
    for (std::size_t row=0; row<8; ++row) A.push_back(row,row,1.0);
    A.finalize();
    Vector state(8,0.0), rhs(8,0.0);
    state(6)=0.0; state(7)=0.0;
    const std::vector<std::size_t> owner{0,0,1};
    const std::vector<std::int64_t> neighbour{-1,1,-1};
    const std::vector<double> flux{-1.0,1.0,1.0};
    run_case("coupled_face_flux_certificate_agrees_on_exact_conservative_flux", [&] {
        const auto c=certify_coupled_face_flux_equivalence(
            A,rhs,state,2,owner,neighbour,flux);
        EXPECT_TRUE(c.complete && c.passed);
        EXPECT_TRUE(c.checked_cells == 2);
        EXPECT_NEAR(c.max_identity_error, 0.0, 0.0);
        EXPECT_NEAR(c.max_conservative_flux_divergence, 0.0, 0.0);
    });
    run_case("coupled_face_flux_certificate_catches_wrong_pressure_row", [&] {
        Vector wrong_rhs = rhs; wrong_rhs(6) = 0.4;
        const auto c=certify_coupled_face_flux_equivalence(
            A,wrong_rhs,state,2,owner,neighbour,flux);
        EXPECT_TRUE(c.complete);
        EXPECT_TRUE(!c.passed);
        EXPECT_TRUE(c.max_identity_error > 0.3);
        std::cout << "FACE_FLUX_ORACLE row_vs_flux_error="
                  << c.max_identity_error << '\n';
    });
    run_case("coupled_face_flux_certificate_separates_conservation_gate", [&] {
        const std::vector<double> not_converged{-1.0,1.25,1.0};
        Vector matching_rhs = rhs; matching_rhs(6)=-0.25; matching_rhs(7)=0.25;
        CoupledFluxEquivalenceControls controls;
        controls.continuity_linf_gate=1e-8;
        const auto c=certify_coupled_face_flux_equivalence(
            A,matching_rhs,state,2,owner,neighbour,not_converged,controls);
        EXPECT_TRUE(c.complete);
        EXPECT_TRUE(c.max_identity_tolerance_ratio <= 1.0);
        EXPECT_TRUE(!c.continuity_gate_passed);
        EXPECT_TRUE(!c.passed);
    });
    run_case("coupled_face_flux_relative_only_zero_residual_is_defined", [&] {
        CoupledFluxEquivalenceControls controls;
        controls.absolute_tolerance = 0.0;
        controls.relative_tolerance = 1e-8;
        const auto exact = certify_coupled_face_flux_equivalence(
            A,rhs,state,2,owner,neighbour,flux,controls);
        EXPECT_TRUE(exact.complete && exact.passed);
        EXPECT_NEAR(exact.max_identity_tolerance_ratio, 0.0, 0.0);
        auto wrong = rhs; wrong(6) = 0.1;
        const auto mismatch = certify_coupled_face_flux_equivalence(
            A,wrong,state,2,owner,neighbour,flux,controls);
        EXPECT_TRUE(mismatch.complete && !mismatch.passed);
    });
    run_case("coupled_face_flux_certificate_rejects_invalid_topology", [&] {
        const std::vector<std::int64_t> bad_neighbour{-1,0,-1};
        EXPECT_TRUE(!certify_coupled_face_flux_equivalence(
            A,rhs,state,2,owner,bad_neighbour,flux).passed);
        const std::vector<double> nan_flux{-1.0,
            std::numeric_limits<double>::quiet_NaN(),1.0};
        EXPECT_TRUE(!certify_coupled_face_flux_equivalence(
            A,rhs,state,2,owner,neighbour,nan_flux).passed);
    });
    return run_all();
}
