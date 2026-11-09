#include "cfdx/core/linalg/pressure_conservative_face_operator.h"
#include "common/test_harness.h"
#include <cmath>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("implicit_fv_face_balance_and_dirichlet_anchor", [] {
        ConservativePressureFaceOperator faces(2, {{0, 1, 2.0}, {1, -1, 1.0}});
        SparseMatrix a(2, 2);
        a.push_back(0, 0, 2.0); a.push_back(0, 1, -2.0);
        a.push_back(1, 0, -2.0); a.push_back(1, 1, 3.0); a.finalize();
        Vector high(2, 0.0), low(2, 0.0), rhs(2, 0.0);
        high(0) = 1e8; high(1) = 1e8; low(1) = 2e-9;
        rhs(0) = -4e-9;
        rhs(1) = 1e8 + 4e-9;
        // The boundary term dominates the relative norm here. Verify the
        // face/CSR identity separately with a modest well-scaled field.
        high(0) = 2.0; high(1) = 3.0; low.fill(0.0);
        rhs(0) = -2.0; rhs(1) = 5.0;
        EXPECT_TRUE(faces.active_face_count() == 2);
        EXPECT_TRUE(faces.residual(rhs, high, low).second < 1e-14);
        EXPECT_TRUE(compensated_pressure_residual(a, rhs, high, low).second < 1e-14);
    });
    run_case("face_refinement_rejects_false_csr_certificate", [] {
        ConservativePressureFaceOperator faces(2, {{0, 1, 2.0}, {1, -1, 0.5}});
        SparseMatrix a(2, 2);
        a.push_back(0, 0, 2.0); a.push_back(0, 1, -2.0);
        a.push_back(1, 0, -2.0); a.push_back(1, 1, 2.5);
        a.finalize();
        Vector rhs(2, 0.0), high(2, 0.0), low(2, 0.0);
        rhs(0) = -2.0; rhs(1) = 3.5;
        const auto proof = refine_pressure_on_implicit_faces(
            faces, a, rhs, high, low, 1e-11, 2,
            [](const SparseMatrix&, const Vector&, Vector& delta) {
                delta(0) = 2.0; delta(1) = 3.0; return true;
            });
        EXPECT_TRUE(proof.valid);
        EXPECT_TRUE(proof.implicit_face_converged);
        EXPECT_TRUE(proof.original_csr_converged);
        EXPECT_TRUE(proof.accepted_refinements == 1);
        EXPECT_TRUE(!proof.physical_fv_flux_certified);
    });
    run_case("nonfinite_coefficient_and_wrong_size_fail_closed", [] {
        EXPECT_THROW((ConservativePressureFaceOperator(2, {{0, 1, -2.0}})),
                     std::invalid_argument);
        ConservativePressureFaceOperator faces(2, {{0, 1, 2.0}});
        Vector b(1, 1.0), h(1, 1.0), l(1, 0.0);
        EXPECT_THROW((void)faces.residual(b, h, l), std::invalid_argument);
    });
    run_case("full_face_slots_conservative_bridge", [] {
        // Two pressure-active faces and a Neumann face: native face indices
        // must be preserved, not compacted to only pressure-active faces.
        ConservativePressureFaceOperator faces(2,
            {{0, 1, 2.0}, {1, -1, 1.0}, {0, -1, 0.0}});
        Vector rhs(2, 0.0), h(2, 0.0), l(2, 0.0);
        h(0)=2.0; h(1)=3.0; rhs(0)=-2.0; rhs(1)=5.0;
        // Before the step: divergence [2,-5] = -rhs. Pressure face
        // response then cancels this divergence exactly.
        const std::vector<double> base{2.0,-3.0,0.0};
        const auto audit=faces.audit_implicit_flux_bridge(rhs,h,l,base,1e-11);
        EXPECT_TRUE(audit.valid);
        EXPECT_TRUE(audit.rhs_matches_base_flux);
        EXPECT_TRUE(audit.implicit_linearized_continuity_converged);
        EXPECT_TRUE(audit.implicit_flux_continuity_relative < 1e-14);
        EXPECT_TRUE(!audit.complete_rhie_chow_certified);
        EXPECT_THROW((void)faces.audit_implicit_flux_bridge(
            rhs,h,l,{2.0,-3.0},1e-11),std::invalid_argument);
    });
    run_case("small_true_residual_does_not_hide_wrong_face_rhs", [] {
        ConservativePressureFaceOperator faces(2,{{0,1,2.0},{1,-1,1.0}});
        Vector rhs(2,0.0),h(2,0.0),l(2,0.0);
        h(0)=2.0;h(1)=3.0;rhs(0)=-2.0;rhs(1)=5.0;
        // A_face*x=b is exact, but the flux input is not the one from
        // which b was assembled. This MUST fail the bridge certificate.
        const auto audit=faces.audit_implicit_flux_bridge(
            rhs,h,l,{2.0001,-3.0},1e-11);
        EXPECT_TRUE(audit.implicit_operator_residual_relative < 1e-14);
        EXPECT_TRUE(!audit.rhs_matches_base_flux);
        EXPECT_TRUE(!audit.implicit_linearized_continuity_converged);
        EXPECT_TRUE(audit.implicit_flux_continuity_relative > 1e-6);
        EXPECT_THROW((void)faces.audit_implicit_flux_bridge(
            rhs,h,l,{2.0,std::numeric_limits<double>::quiet_NaN()},1e-11),
            std::domain_error);
    });
    run_case("rounded_csr_is_not_a_face_acceptance_gate", [] {
        ConservativePressureFaceOperator faces(2,{{0,1,2.0},{1,-1,1.0}});
        SparseMatrix a(2,2);
        a.push_back(0,0,2.0000001); a.push_back(0,1,-2.0);
        a.push_back(1,0,-2.0); a.push_back(1,1,3.0); a.finalize();
        Vector rhs(2,0.0),h(2,0.0),l(2,0.0);
        h(0)=2.0;h(1)=3.0;rhs(0)=-2.0;rhs(1)=5.0;
        const auto audit=faces.audit_implicit_flux_bridge(
            rhs,h,l,{2.0,-3.0},1e-11);
        EXPECT_TRUE(audit.implicit_linearized_continuity_converged);
        EXPECT_TRUE(compensated_pressure_residual(a,rhs,h,l).second > 1e-11);
        EXPECT_TRUE(!audit.complete_rhie_chow_certified);
    });
    run_case("full_RC_gap_is_deferred_and_locates_native_face", [] {
        ConservativePressureFaceOperator faces(2,
            {{0,1,2.0},{1,-1,1.0},{0,-1,0.0}});
        Vector rhs(2,0.0),h(2,0.0),l(2,0.0);
        h(0)=2.0;h(1)=3.0;rhs(0)=-2.0;rhs(1)=5.0;
        const std::vector<double> base{2.0,-3.0,0.0};
        // Implicit correction cancels the base divergence. The additional
        // deferred gradient jump on native face 0 causes a physical defect.
        const std::vector<long double> full{0.25L,0.0L,0.0L};
        const auto gap=faces.diagnose_full_flux_gap(rhs,h,l,base,full);
        EXPECT_TRUE(gap.valid);
        EXPECT_TRUE(gap.implicit_continuity_relative<1e-14);
        EXPECT_TRUE(gap.deferred_divergence_relative>1e-2);
        EXPECT_TRUE(gap.full_continuity_relative>1e-2);
        EXPECT_TRUE(gap.worst_deferred_face==0);
        EXPECT_TRUE(gap.worst_divergence_cell==0);
        EXPECT_THROW((void)faces.diagnose_full_flux_gap(
            rhs,h,l,base,{0.0L,0.0L}),std::invalid_argument);
    });
    return run_all();
}
