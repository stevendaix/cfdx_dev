#include "cfdx/physics/pcd_fp_equilibration.h"
#include "cfdx/physics/pcd_fp_capture.h"
#include "common/test_harness.h"
#include <cmath>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

int main() {
    run_case("fp_row_equilibration_preserves_original_solution", [] {
        SparseMatrix A(3, 3);
        A.push_back(0,0,1e-9); A.push_back(0,1,-3e-10);
        A.push_back(1,0,-8e3); A.push_back(1,1,2e4); A.push_back(1,2,-2e3);
        A.push_back(2,1,-10); A.push_back(2,2,20);
        A.finalize();
        const auto eq = equilibrate_pcd_fp_rows(A);
        EXPECT_TRUE(eq.smallest_row_max <= 1e-9);
        EXPECT_TRUE(eq.largest_row_max >= 2e4);
        EXPECT_TRUE(eq.row_multipliers.size() == 3);
        Vector x(3, 0.0); x(0) = 1.2; x(1) = -2.3; x(2) = 4.5;
        Vector b(3, 0.0);
        const auto data = A.matvec(x);
        for (std::size_t i=0;i<3;++i) b(i)=data[i];
        const auto scaled_b = pcd_fp_scaled_rhs(b, eq.row_multipliers);
        const auto scaled_ax = eq.matrix.matvec(x);
        for (std::size_t i = 0; i < 3; ++i)
            EXPECT_NEAR(scaled_ax[i], scaled_b(i), 1e-12);
        Vector zero(3, 0.0);
        const auto quality = pcd_fp_true_residual(A, b, zero, x);
        EXPECT_TRUE(quality.relative_rhs < 1e-14);
        const auto bad = pcd_fp_true_residual(A, b, zero, zero);
        EXPECT_TRUE(bad.relative_rhs > 0.999);
    });
    run_case("fp_capture_budget_is_global_and_separates_hard_rhs", [] {
        for (std::size_t n = 1; n <= 3; ++n)
            EXPECT_TRUE(reserve_pcd_fp_capture_slot(false) == n);
        EXPECT_TRUE(reserve_pcd_fp_capture_slot(false) == 0);
        for (std::size_t n = 1; n <= 3; ++n)
            EXPECT_TRUE(reserve_pcd_fp_capture_slot(true) == n);
        EXPECT_TRUE(reserve_pcd_fp_capture_slot(true) == 0);
    });
    run_case("fp_zero_row_fails_closed", [] {
        SparseMatrix A(2,2);
        A.push_back(0,0,1.0); A.finalize();
        bool rejected = false;
        try { (void)equilibrate_pcd_fp_rows(A); }
        catch (const std::invalid_argument&) { rejected = true; }
        EXPECT_TRUE(rejected);
    });
    run_case("fp_amg_action_probe_rejects_catastrophic_amplification", [] {
        struct Amplifier : Preconditioner {
            double factor;
            explicit Amplifier(double f) : factor(f) {}
            bool setup(const SparseMatrix&) override { return true; }
            bool apply(const Vector& rhs, Vector& x) const override {
                if (rhs.size() != x.size()) return false;
                for (std::size_t i = 0; i < rhs.size(); ++i)
                    x(i) = factor * rhs(i);
                return true;
            }
            const char* name() const override { return "AmplifierProbe"; }
        };
        SparseMatrix A(3,3);
        for (std::size_t i=0;i<3;++i) A.push_back(i,i,1.0);
        A.finalize();
        Vector rhs(3,1.0);
        Amplifier good(1.0), bad(1e15);
        EXPECT_TRUE(good.setup(A));
        const auto safe=pcd_fp_action_probe(A,good,rhs);
        EXPECT_TRUE(safe.applied && safe.finite);
        EXPECT_NEAR(safe.amplification,1.0,1e-12);
        const auto unsafe=pcd_fp_action_probe(A,bad,rhs);
        EXPECT_TRUE(unsafe.applied && unsafe.finite);
        EXPECT_TRUE(unsafe.amplification > 1e6);
        Vector zero(3,0.0);
        const auto null_rhs=pcd_fp_action_probe(A,bad,zero);
        EXPECT_TRUE(null_rhs.applied && null_rhs.finite);
        EXPECT_NEAR(null_rhs.amplification,0.0,0.0);
    });
    return run_all();
}
