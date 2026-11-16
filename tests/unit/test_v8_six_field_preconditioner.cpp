#include "cfdx/core/linalg/six_field_block_preconditioner_v8.h"
#include "cfdx/core/linalg/globalized_jfnk_v8.h"
#include "common/test_harness.h"
#include <cmath>
#include <string>

using namespace cfdx::core;
using namespace cfdx::testing;

static SparseMatrix sample_coupled_six_field_matrix(std::size_t nc) {
    SparseMatrix A(6*nc,6*nc);
    for(std::size_t i=0;i<6*nc;++i) {
        A.push_back(i,i, i<4*nc ? 2. : (i<5*nc ? 4. : 3.));
        // Real operator may couple velocity/pressure and SST; the block
        // preconditioner is allowed to discard these off-diagonal entries.
        if(i<nc) A.push_back(i,4*nc+i,0.05);
        if(i>=4*nc && i<5*nc) A.push_back(i,i-4*nc,0.01);
    }
    A.finalize();
    return A;
}
int main() {
    run_case("v8_block_adapter_uses_real_preconditioner_setup", [] {
        JacobiPreconditioner up,k,w;
        SixFieldBlockPreconditionerV8 blocks(2,up,k,w);
        const auto A=sample_coupled_six_field_matrix(2);
        EXPECT_TRUE(blocks.setup(A));
        Vector r(12,1.),z(12,0.);
        EXPECT_TRUE(blocks.apply(r,z));
        for(std::size_t i=0;i<8;++i) EXPECT_NEAR(z(i),.5,1e-14);
        for(std::size_t i=8;i<10;++i) EXPECT_NEAR(z(i),.25,1e-14);
        for(std::size_t i=10;i<12;++i) EXPECT_NEAR(z(i),1./3.,1e-14);
    });
    run_case("v8_block_adapter_rejects_wrong_operator_size", [] {
        JacobiPreconditioner up,k,w;
        SixFieldBlockPreconditionerV8 blocks(2,up,k,w);
        EXPECT_TRUE(!blocks.setup(sample_coupled_six_field_matrix(1)));
        Vector r(12,1.),z(12,0.);
        EXPECT_TRUE(!blocks.apply(r,z));
    });
    run_case("v8_jfnk_runs_with_true_six_field_block_adapter", [] {
        JacobiPreconditioner up,k,w;
        SixFieldBlockPreconditionerV8 blocks(1,up,k,w);
        EXPECT_TRUE(blocks.setup(sample_coupled_six_field_matrix(1)));
        Vector state(6,0.);
        state(4)=0.25;state(5)=0.4;
        const auto nonlinear=[](const Vector& x,Vector& f) {
            if(x(4)<=0. || x(5)<=0.) return false;
            f(0)=x(0)-1.+0.15*x(3)+.05*x(4);
            f(1)=x(1)+.5-.1*x(3);
            f(2)=x(2)-.2;
            f(3)=x(3)+.1*x(0);
            f(4)=2.*(x(4)-.8)+.06*x(0);
            f(5)=x(5)-1.2+.08*x(4);
            return true;
        };
        GlobalizedJfnkControlsV8 controls;
        controls.nonlinear_tolerance=1e-9;
        controls.linear_tolerance=1e-5;
        const auto report=globalized_jfnk_v8(state,nonlinear,&blocks,controls);
        Vector r(6,0.);
        EXPECT_TRUE(nonlinear(state,r));
        EXPECT_TRUE(report.converged);
        EXPECT_TRUE(r.norm_inf()<1e-9);
    });
    return run_all();
}
