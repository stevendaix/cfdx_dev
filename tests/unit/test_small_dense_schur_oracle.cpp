#include "cfdx/core/linalg/small_dense_schur_oracle.h"
#include "cfdx/core/linalg/schur_inverse_quality_certificate.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix make(std::size_t n, const std::vector<double>& a) {
    SparseMatrix m(n,n);
    for(std::size_t i=0; i<n; ++i)
        for(std::size_t j=0; j<n; ++j)
            if(a[i*n+j] != 0.0) m.push_back(i,j,a[i*n+j]);
    m.finalize(); return m;
}
}

int main() {
    const auto M = make(2,{4.0,1.0, -0.5,3.0});
    const auto G = make(2,{1.0,0.2,0.3,1.0});
    const auto D = make(2,{1.0,-0.1,0.2,1.0});
    const auto C = make(2,{0.2,0.1,0.0,0.3});
    const BlockOperator blocks(M,G,D,C);
    run_case("small_dense_lu_checks_true_residual_on_nonsymmetric_matrix", [&] {
        SmallDenseCertifiedLU solve;
        EXPECT_TRUE(solve.setup(M));
        Vector rhs(2,0.0),x(2,0.0); rhs(0)=1;rhs(1)=-2;
        EXPECT_TRUE(solve.solve(rhs,x));
        Vector y(2,0.0); M.matvec_into(x,y);
        EXPECT_NEAR(y(0),rhs(0),1e-12);
        EXPECT_NEAR(y(1),rhs(1),1e-12);
    });
    run_case("small_dense_schur_matches_independent_schur_action", [&] {
        SmallDenseExactSchurInverse schur;
        EXPECT_TRUE(schur.setup(blocks));
        Vector rhs(2,0.0), p(2,0.0);
        rhs(0)=1;rhs(1)=-0.6;
        EXPECT_TRUE(schur.apply(rhs,p));
        Vector Sx(2,0.0); schur.schur_matrix().matvec_into(p,Sx);
        EXPECT_NEAR(Sx(0),rhs(0),1e-10);
        EXPECT_NEAR(Sx(1),rhs(1),1e-10);
        const auto probes=schur_inverse_quality_basis(2);
        const auto report=certify_schur_inverse_action(
            blocks, schur,
            [&](const Vector& r,Vector& x) {return schur.solve_velocity(r,x);},probes);
        EXPECT_TRUE(report.complete && report.passed);
        EXPECT_TRUE(report.worst_schur_defect<1e-10);
    });
    run_case("negative_definite_schur_is_not_rejected_as_non_spd", [&] {
        const auto diag=make(2,{2.0,0,0,4.0});
        const auto I=make(2,{1,0,0,1});
        const auto empty=make(2,{0,0,0,0});
        SmallDenseExactSchurInverse schur;
        EXPECT_TRUE(schur.setup(BlockOperator(diag,I,I,empty)));
        Vector b(2,0.0),z(2,0.0);
        b(0)=1;b(1)=-3;
        EXPECT_TRUE(schur.apply(b,z));
        EXPECT_NEAR(z(0),-2.0,1e-12);
        EXPECT_NEAR(z(1),12.0,1e-12);
    });
    run_case("small_dense_schur_rejects_singular_update_without_state_mutation", [&] {
        SmallDenseExactSchurInverse schur;
        EXPECT_TRUE(schur.setup(blocks));
        Vector b(2,1.0),before(2,0.0),after(2,0.0);
        EXPECT_TRUE(schur.apply(b,before));
        const auto empty=make(2,{0,0,0,0});
        const BlockOperator singular(M,empty,empty,empty);
        EXPECT_TRUE(!schur.update_values(singular));
        EXPECT_TRUE(schur.apply(b,after));
        EXPECT_NEAR(before(0),after(0),0.0);
        EXPECT_NEAR(before(1),after(1),0.0);
    });
    run_case("small_dense_reference_rejects_nonfinite_and_oversized", [&] {
        SmallDenseCertifiedLU low_limit(1);
        EXPECT_TRUE(!low_limit.setup(M));
        // The CSR constructor already rejects non-finite coefficients.
        EXPECT_THROW((void)make(2,{1,std::numeric_limits<double>::quiet_NaN(),0,1}), std::invalid_argument);
        SmallDenseExactSchurInverse oracle(1);
        EXPECT_TRUE(!oracle.setup(blocks));
        Vector sentinel(2,7.0);
        EXPECT_TRUE(!oracle.apply(Vector(2,1.0),sentinel));
        EXPECT_NEAR(sentinel(0),7.0,0.0);
    });
    return run_all();
}
