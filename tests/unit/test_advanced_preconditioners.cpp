#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix make_poisson() {
    SparseMatrix A(3,3);
    A.push_back(0,0,2.0); A.push_back(0,1,-1.0);
    A.push_back(1,0,-1.0); A.push_back(1,1,2.0); A.push_back(1,2,-1.0);
    A.push_back(2,1,-1.0); A.push_back(2,2,2.0);
    A.finalize();
    return A;
}

int main() {
    run_case("gauss_seidel_setup_apply", [] {
        auto A = make_poisson();
        GaussSeidelPreconditioner p;
        EXPECT_TRUE(p.setup(A));
        Vector r(3), z(3);
        r(0)=1.0; r(1)=0.0; r(2)=1.0;
        EXPECT_TRUE(p.apply(r,z));
        EXPECT_NEAR(z(0),0.5,1e-14);
        EXPECT_NEAR(z(1),0.25,1e-14);
        EXPECT_NEAR(z(2),0.625,1e-14);
    });

    run_case("ilu0_setup_apply", [] {
        auto A = make_poisson();
        ILU0Preconditioner p;
        EXPECT_TRUE(p.setup(A));
        Vector r(3), z(3);
        r(0)=1.0; r(1)=0.0; r(2)=1.0;
        EXPECT_TRUE(p.apply(r,z));
        EXPECT_TRUE(z.is_valid());
        auto Az = A.matvec(z);
        EXPECT_TRUE(Az.size()==3);
        EXPECT_TRUE(std::abs(Az[0]-r(0)) < 1.0);
    });

    run_case("preconditioner_rejects_missing_diagonal", [] {
        SparseMatrix A(2,2);
        A.push_back(0,1,1.0);
        A.push_back(1,0,1.0);
        A.finalize();
        GaussSeidelPreconditioner gs;
        ILU0Preconditioner ilu;
        EXPECT_TRUE(!gs.setup(A));
        EXPECT_TRUE(!ilu.setup(A));
    });

    return run_all();
}
