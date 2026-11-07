#include "cfdx/core/linalg/mgr_preconditioner.h"
#include "common/test_harness.h"

#include <cmath>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix assemble(double fine_diag, double coarse_diag,
                      bool changed_graph=false) {
    SparseMatrix A(8,8);
    for(std::size_t i=0;i<4;++i) {
        A.push_back(i,i,fine_diag);
        A.push_back(i,i+4,-1.0);
        A.push_back(i+4,i,-1.0);
        A.push_back(i+4,i+4,coarse_diag);
    }
    if(changed_graph) A.push_back(0,5,0.125);
    A.finalize();
    return A;
}
}
int main() {
    run_case("failed_MGR_numeric_factor_does_not_mutate_old_factors", [] {
        NativeMGRPreconditioner mgr({0,1,2,3},{4,5,6,7},
            NativeMGRPreconditioner::FineRelaxation::ILU0);
        const auto A=assemble(4.0,3.0);
        EXPECT_TRUE(mgr.setup(A));
        Vector r(8,1.0), previous(8,0.0), current(8,0.0);
        EXPECT_TRUE(mgr.apply(r,previous));
        EXPECT_FALSE(mgr.update_values(assemble(0.0,3.0)));
        EXPECT_TRUE(mgr.numeric_updates()==0);
        EXPECT_TRUE(mgr.apply(r,current));
        for(std::size_t i=0;i<8;++i)
            EXPECT_NEAR(current(i),previous(i),1e-13);
        EXPECT_FALSE(mgr.setup(SparseMatrix{}));
        EXPECT_TRUE(mgr.apply(r,current));
        for(std::size_t i=0;i<8;++i)
            EXPECT_NEAR(current(i),previous(i),1e-13);
    });
    run_case("rejected_MGR_graph_change_preserves_last_solution", [] {
        NativeMGRPreconditioner mgr({0,1,2,3},{4,5,6,7});
        EXPECT_TRUE(mgr.setup(assemble(4.0,3.0)));
        Vector b(8,1.0), prev(8,0.0), after(8,0.0);
        EXPECT_TRUE(mgr.apply(b,prev));
        EXPECT_FALSE(mgr.update_values(assemble(4.0,3.0,true)));
        EXPECT_TRUE(mgr.apply(b,after));
        for(std::size_t i=0;i<8;++i) EXPECT_NEAR(prev(i),after(i),1e-13);
    });
    run_case("accepted_MGR_update_commits_once", [] {
        NativeMGRPreconditioner mgr({0,1,2,3},{4,5,6,7},
            NativeMGRPreconditioner::FineRelaxation::ILU0);
        EXPECT_TRUE(mgr.setup(assemble(4.0,3.0)));
        Vector b(8,1.0), prev(8,0.0), after(8,0.0);
        EXPECT_TRUE(mgr.apply(b,prev));
        EXPECT_TRUE(mgr.update_values(assemble(6.0,2.0)));
        EXPECT_TRUE(mgr.setup_count()==1);
        EXPECT_TRUE(mgr.numeric_updates()==1);
        EXPECT_TRUE(mgr.apply(b,after));
        EXPECT_TRUE(std::abs(prev(0)-after(0))>1e-6);
        EXPECT_TRUE(after.is_valid());
    });
    run_case("MGR_alias_safe_in_place_apply", [] {
        NativeMGRPreconditioner mgr({0,1,2,3},{4,5,6,7},
            NativeMGRPreconditioner::FineRelaxation::ILU0);
        EXPECT_TRUE(mgr.setup(assemble(4.0,3.0)));
        Vector rhs(8,0.0), expected(8,0.0);
        for(std::size_t i=0;i<8;++i) rhs(i)=0.35+0.17*static_cast<double>(i);
        EXPECT_TRUE(mgr.apply(rhs,expected));
        Vector in_place=rhs;
        EXPECT_TRUE(mgr.apply(in_place,in_place));
        for(std::size_t i=0;i<8;++i)
            EXPECT_NEAR(in_place(i),expected(i),1e-12);
    });
    return run_all();
}