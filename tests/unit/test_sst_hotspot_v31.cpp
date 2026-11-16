#include "cfdx/physics/sst_hotspot_derivatives_v31.h"
#include "common/test_harness.h"
#include <cmath>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main() {
    run_case("k_and_omega_cross_derivatives_use_true_row_and_close", []{
        const std::vector<double> state{2.,3.,4.,5.}; // k(0),k(1),omega(0),omega(1)
        const SSTV31BudgetEvaluator sample=[](const auto& x){
            SSTV31RowComponents b;
            // Deliberately nontrivial coupling, including neighbouring fields.
            const double conv=1.5*x[0]-0.8*x[1];
            const double diff=2.*x[0]-0.2*x[1]+0.1*x[3];
            const double prod=-0.4*x[0]*x[2];
            const double dest=0.2*x[0]*x[2]*x[2];
            const double cross=0.03*x[1]*x[3];
            const double all=conv+diff+prod+dest+cross;
            b.parts={{conv,diff,prod,dest,cross,all}};
            return b;
        };
        const auto w=sst_v31_log_derivative(state,2,1e-8,1e-4,sample);
        // dR/d(log(omega(0))) = -0.4*k*omega + 0.4*k*omega^2
        EXPECT_NEAR(w.d_dlog[5],-0.4*2.*4.+0.4*2.*16.,1e-6);
        EXPECT_TRUE(w.maximum_relative_closure<1e-12);
        const auto kn=sst_v31_log_derivative(state,1,1e-8,1e-4,sample);
        EXPECT_NEAR(kn.d_dlog[5],3.*(-.8-.2+.03*5.),1e-6);
        const auto wn=sst_v31_log_derivative(state,3,1e-8,1e-4,sample);
        EXPECT_NEAR(wn.d_dlog[5],5.*(.1+.03*3.),1e-6);
    });
    run_case("hotspot_rejects_false_FV_closure_even_if_row_is_finite", []{
        const SSTV31BudgetEvaluator wrong=[](const auto& x){
            SSTV31RowComponents b;b.parts={{x[0],0.,0.,0.,0.,100.*x[0]}};
            return b;
        };
        EXPECT_THROW(sst_v31_log_derivative({2.},0,1e-8,1e-4,wrong),std::invalid_argument);
    });
    run_case("hotspot_rejects_constraint_rows_and_floor_violations", []{
        const SSTV31BudgetEvaluator constrained=[](const auto& x){
            SSTV31RowComponents b;b.constrained=true;b.parts={{0,0,0,0,0,x[0]}};return b;
        };
        EXPECT_THROW(sst_v31_log_derivative({2.},0,1e-8,1e-4,constrained),std::invalid_argument);
        EXPECT_THROW(sst_v31_log_derivative({1.00001},0,1.,1e-4,constrained),std::invalid_argument);
    });
    return run_all();
}
