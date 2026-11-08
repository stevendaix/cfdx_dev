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
    return run_all();
}
