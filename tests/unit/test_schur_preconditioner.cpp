#include "cfdx/core/linalg/block_preconditioner.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_block_system() {
    SparseMatrix A(3, 3);
    A.push_back(0,0,4.0); A.push_back(0,1,1.0); A.push_back(0,2,1.0);
    A.push_back(1,0,1.0); A.push_back(1,1,3.0); A.push_back(1,2,2.0);
    A.push_back(2,0,2.0); A.push_back(2,1,1.0); A.push_back(2,2,5.0);
    A.finalize();
    return A;
}

int main() {
    run_case("schur_complement_setup_and_apply", [] {
        const auto A = make_block_system();
        SchurComplementPreconditioner p({0,1}, {2});
        EXPECT_TRUE(p.setup(A));

        Vector r(3);
        r(0) = 9.0;
        r(1) = 10.0;
        r(2) = 19.0;
        Vector z(3);
        EXPECT_TRUE(p.apply(r, z));
        EXPECT_TRUE(z.is_valid());

        // A11^{-1} r1 = (17/11, 21/11).
        EXPECT_NEAR(z(0), 17.0 / 11.0, 1e-12);
        EXPECT_NEAR(z(1), 21.0 / 11.0, 1e-12);

        // S = 5 - [2,1] A11^{-1} [1,2]^T = 46/11.
        // rhs2 = 19 - [2,1] [17/11,21/11]^T = 14.
        EXPECT_NEAR(z(2), 77.0 / 23.0, 1e-12);
        EXPECT_TRUE(std::isfinite(z.norm_inf()));
    });

    run_case("schur_rejects_invalid_partition", [] {
        SparseMatrix A(3,3);
        A.push_back(0,0,1.0); A.push_back(1,1,1.0); A.push_back(2,2,1.0);
        A.finalize();

        SchurComplementPreconditioner overlap({0,1}, {1,2});
        EXPECT_TRUE(!overlap.setup(A));

        SchurComplementPreconditioner incomplete({0}, {2});
        EXPECT_TRUE(!incomplete.setup(A));
    });

    return run_all();
}
