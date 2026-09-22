#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "common/test_harness.h"
#include <algorithm>
#include <cmath>
#include <limits>

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

static SparseMatrix make_dense3() {
    SparseMatrix A(3,3);
    A.push_back(0,0,4.0); A.push_back(0,1,1.0); A.push_back(0,2,1.0);
    A.push_back(1,0,2.0); A.push_back(1,1,5.0); A.push_back(1,2,1.0);
    A.push_back(2,0,1.0); A.push_back(2,1,1.0); A.push_back(2,2,3.0);
    A.finalize();
    return A;
}

static double residual_inf(const SparseMatrix& A, const Vector& z, const Vector& r) {
    const auto Az = A.matvec(z);
    double max_err = 0.0;
    for (std::size_t i = 0; i < r.size(); ++i)
        max_err = std::max(max_err, std::abs(Az[i] - r(i)));
    return max_err;
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

    run_case("ilu0_exact_on_zero_fill_complete_pattern", [] {
        auto A = make_dense3();
        ILU0Preconditioner p;
        EXPECT_TRUE(p.setup(A));
        Vector r(3), z(3);
        r(0)=9.0; r(1)=15.0; r(2)=12.0;
        EXPECT_TRUE(p.apply(r,z));
        EXPECT_TRUE(z.is_valid());
        EXPECT_NEAR(z(0),1.0,1e-12);
        EXPECT_NEAR(z(1),2.0,1e-12);
        EXPECT_NEAR(z(2),3.0,1e-12);
        EXPECT_TRUE(residual_inf(A,z,r) < 1e-12);
    });

    run_case("ilu0_sparse_pattern_reduces_residual", [] {
        auto A = make_poisson();
        ILU0Preconditioner p;
        EXPECT_TRUE(p.setup(A));
        Vector r(3), z(3);
        r(0)=1.0; r(1)=0.0; r(2)=1.0;
        EXPECT_TRUE(p.apply(r,z));
        EXPECT_TRUE(z.is_valid());
        EXPECT_TRUE(std::isfinite(residual_inf(A,z,r)));
        EXPECT_TRUE(residual_inf(A,z,r) < r.norm_inf());
    });

    run_case("preconditioner_rejects_missing_or_zero_diagonal", [] {
        SparseMatrix missing(2,2);
        missing.push_back(0,1,1.0); missing.push_back(1,0,1.0); missing.finalize();
        SparseMatrix zero(2,2);
        zero.push_back(0,0,0.0); zero.push_back(1,1,1.0); zero.finalize();
        ILU0Preconditioner p;
        EXPECT_TRUE(!p.setup(missing));
        Vector r(2,1.0), z(2,0.0);
        EXPECT_TRUE(!p.apply(r,z));
        EXPECT_TRUE(!p.setup(zero));
    });

    run_case("ilu0_rejects_nonfinite_input", [] {
        SparseMatrix A = make_poisson();
        A.values_data()[1] = std::numeric_limits<double>::quiet_NaN();
        ILU0Preconditioner p;
        EXPECT_TRUE(!p.setup(A));
    });

    run_case("ilu0_accepts_out_of_order_assembly", [] {
        SparseMatrix A(2,2);
        A.push_back(0,1,1.0);
        A.push_back(0,0,2.0);
        A.push_back(1,1,2.0);
        A.finalize();
        ILU0Preconditioner p;
        EXPECT_TRUE(p.setup(A));
    });

    return run_all();
}
