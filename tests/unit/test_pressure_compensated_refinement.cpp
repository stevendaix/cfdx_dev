#include "cfdx/core/linalg/pressure_compensated_refinement.h"
#include "common/test_harness.h"
#include <cmath>
#include <iostream>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("compensated_refinement_restores_unrepresentable_jump", [] {
        SparseMatrix a(2, 2);
        a.push_back(0, 0, 2.0); a.push_back(0, 1, -2.0);
        a.push_back(1, 0, -2.0); a.push_back(1, 1, 2.0);
        a.finalize();
        Vector b(2, 0.0), high(2, 1e8), low(2, 0.0);
        b(0) = -4e-9; b(1) = 4e-9;
        const auto before = compensated_pressure_residual(a, b, high, low);
        EXPECT_TRUE(before.second > 0.99);
        const auto report = refine_pressure_high_low(a, b, high, low, 1e-11, 2,
            [](const SparseMatrix&, const Vector& rhs, Vector& delta) {
                delta(0) = 0.0;
                delta(1) = rhs(1) / 2.0;
                return true;
            });
        EXPECT_TRUE(report.valid);
        EXPECT_TRUE(report.algebraic_converged);
        EXPECT_TRUE(!report.physical_face_flux_certified);
        EXPECT_TRUE(report.passes_accepted == 1);
        EXPECT_TRUE(report.final_relative_residual < 1e-12);
        EXPECT_TRUE(report.final_double_only_residual > 0.99);
        EXPECT_TRUE(std::abs(low(1) - 2e-9) < 1e-24);
    });
    run_case("zero_progress_and_nonfinite_increment_fail_closed", [] {
        SparseMatrix a(2, 2);
        a.push_back(0, 0, 3.0); a.push_back(0, 1, -1.0);
        a.push_back(1, 0, -1.0); a.push_back(1, 1, 3.0);
        a.finalize();
        Vector b(2, 1.0), high(2, 0.0), low(2, 0.0);
        const auto prior = high;
        auto report = refine_pressure_high_low(a, b, high, low, 1e-11, 3,
            [](const SparseMatrix&, const Vector&, Vector& delta) {
                delta(0) = 0.0; delta(1) = 0.0; return true;
            });
        EXPECT_TRUE(!report.algebraic_converged);
        EXPECT_TRUE(report.passes_accepted == 0);
        EXPECT_TRUE(high(0) == prior(0));
        report = refine_pressure_high_low(a, b, high, low, 1e-11, 3,
            [](const SparseMatrix&, const Vector&, Vector& delta) {
                delta(0) = std::numeric_limits<double>::infinity();
                return true;
            });
        EXPECT_TRUE(!report.algebraic_converged);
        EXPECT_TRUE(report.passes_accepted == 0);
        EXPECT_TRUE(!report.physical_face_flux_certified);
    });
    run_case("well_resolved_exact_system_does_not_invoke_callback", [] {
        SparseMatrix a(2, 2);
        a.push_back(0, 0, 3.0); a.push_back(0, 1, -1.0);
        a.push_back(1, 0, -1.0); a.push_back(1, 1, 2.0);
        a.finalize();
        Vector b(2, 0.0), high(2, 1.0), low(2, 0.0);
        b(0) = 2.0; b(1) = 1.0;
        const auto report = refine_pressure_high_low(a, b, high, low, 1e-11, 3,
            [](const SparseMatrix&, const Vector&, Vector&) {
                throw std::runtime_error("must not call"); return false;
            });
        EXPECT_TRUE(report.algebraic_converged);
        EXPECT_TRUE(report.passes_attempted == 0);
    });
    return run_all();
}
