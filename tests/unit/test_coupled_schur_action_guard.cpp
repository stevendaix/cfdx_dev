#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <cmath>
#include <iostream>

using namespace cfdx::core;
using namespace cfdx::testing;

int main() {
    // Auu is coupled across cells. Its exact inverse is NOT the local M^-1;
    // choosing C~D*M^-1*G makes the approximate Schur nearly singular,
    // although the physical monolithic A is well-defined. A truncated
    // Schur action then produces a spectacular A*M^-1*r amplification.
    SparseMatrix A(8, 8);
    for (std::size_t i = 0; i < 6; ++i)
        A.push_back(i, i, 2.0);
    A.push_back(0, 1, -1.0);
    A.push_back(1, 0, -1.0);
    for (std::size_t i = 0; i < 2; ++i) {
        A.push_back(i, 6 + i, 1.0);
        A.push_back(6 + i, i, 1.0);
        A.push_back(6 + i, 6 + i, 0.500001);
    }
    A.finalize();
    Vector rhs(8, 0.0);
    rhs(6) = 1.0;
    rhs(7) = -0.3;
    run_case("default_guard_does_not_change_existing_action", [&] {
        CoupledBlockSchurOptions options;
        CoupledBlockSchurAMGPreconditioner pre(2, options);
        EXPECT_TRUE(pre.setup(A));
        Vector z(8, 0.0);
        EXPECT_TRUE(pre.apply(rhs, z));
        EXPECT_TRUE(z.is_valid());
    });
    run_case("catastrophic_image_is_refused_without_fallback", [&] {
        CoupledBlockSchurOptions options;
        options.max_action_amplification = 1e3;
        CoupledBlockSchurAMGPreconditioner pre(2, options);
        EXPECT_TRUE(pre.setup(A));
        Vector z(8, 0.0);
        EXPECT_TRUE(!pre.apply(rhs, z));
        for (std::size_t i = 0; i < z.size(); ++i)
            EXPECT_TRUE(z(i) == 0.0);
    });
    run_case("normal_stable_action_passes_optin_guard", [&] {
        SparseMatrix B(8, 8);
        for (std::size_t i = 0; i < 6; ++i) B.push_back(i, i, 2.0);
        for (std::size_t i = 0; i < 2; ++i) {
            B.push_back(i, 6 + i, 1.0);
            B.push_back(6 + i, i, 1.0);
            B.push_back(6 + i, 6 + i, 2.0);
        }
        B.finalize();
        CoupledBlockSchurOptions options;
        options.max_action_amplification = 1e3;
        CoupledBlockSchurAMGPreconditioner pre(2, options);
        EXPECT_TRUE(pre.setup(B));
        Vector z(8, 0.0);
        EXPECT_TRUE(pre.apply(rhs, z));
        EXPECT_TRUE(z.is_valid());
    });
    run_case("invalid_guard_threshold_rejected_at_setup", [&] {
        CoupledBlockSchurOptions options;
        options.max_action_amplification = -1.0;
        CoupledBlockSchurAMGPreconditioner pre(2, options);
        EXPECT_TRUE(!pre.setup(A));
    });
    std::cout << "Schur catastrophic action guard PASS\n";
    return 0;
}
