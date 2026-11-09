#include "cfdx/core/linalg/pressure_roundoff_leakage.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    run_case("near_neumann_row_sum_roundoff_is_not_dirichlet_leakage", [] {
        SparseMatrix a(3, 3);
        a.push_back(0,0,0.3); a.push_back(0,1,-0.1);
        a.push_back(1,0,-0.1); a.push_back(1,1,0.3); a.push_back(1,2,-0.2);
        a.push_back(2,1,-0.2); a.push_back(2,2,0.4);
        a.finalize();
        Vector rhs(3, 1.0), p(3, 0.0);
        p(1) = 1e8;
        auto v = diagnose_pressure_roundoff_leakage(a, rhs, p);
        EXPECT_TRUE(v.valid);
        EXPECT_TRUE(v.rows == 3);
        EXPECT_TRUE(v.pressure_anchor_rows == 2);
        EXPECT_TRUE(v.negligible_roundoff_rows == 1);
        EXPECT_TRUE(v.spurious_leakage_relative_l2 > 1e-11);
        EXPECT_TRUE(v.worst_row == 1);
        EXPECT_TRUE(!v.physical_face_flux_certified);
        p.fill(0.0);
        v = diagnose_pressure_roundoff_leakage(a, rhs, p);
        EXPECT_TRUE(v.valid);
        EXPECT_TRUE(v.spurious_leakage_relative_l2 == 0.0);
    });
    run_case("bad_precision_and_nonfinite_state_fail_closed", [] {
        SparseMatrix a(1, 1);
        a.push_back(0,0,1.0); a.finalize();
        Vector rhs(1,1.0), p(1,1.0);
        EXPECT_TRUE(!diagnose_pressure_roundoff_leakage(a,rhs,p,0.0).valid);
        p(0) = std::numeric_limits<double>::infinity();
        EXPECT_TRUE(!diagnose_pressure_roundoff_leakage(a,rhs,p).valid);
        p(0) = 1.0; rhs(0) = 0.0;
        EXPECT_TRUE(!diagnose_pressure_roundoff_leakage(a,rhs,p).valid);
    });
    return run_all();
}
