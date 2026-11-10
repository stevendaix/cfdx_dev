#include "cfdx/physics/sst_anderson_v20.h"
#include "common/test_harness.h"
#include <limits>
using namespace cfdx::physics;
using namespace cfdx::testing;
int main(){
 run_case("anderson_two_cycle",[]{
   SSTAndersonV20 h;
   EXPECT_TRUE(!h.propose({1.},{3.},{3.},{5.},1.e-8,1.e-8).available);
   h.record({1.},{3.},{3.},{5.});
   auto p=h.propose({3.},{5.},{1.},{3.},1.e-8,1.e-8);
   EXPECT_TRUE(p.available);
   EXPECT_NEAR(p.theta,0.5,1e-12);
   EXPECT_NEAR(p.k[0],2.,1e-12);
   EXPECT_NEAR(p.omega[0],4.,1e-12);
 });
 run_case("anderson_no_fake_improvement",[]{
   SSTAndersonV20 h;
   h.record({1.},{2.},{2.},{3.});
   EXPECT_TRUE(!h.propose({3.},{4.},{4.},{5.},1e-8,1e-8).available);
 });
 run_case("anderson_fail_closed",[]{
   SSTAndersonV20 h;
   EXPECT_THROW(h.record({1.},{2.},{std::numeric_limits<double>::infinity()},{3.}),std::invalid_argument);
   EXPECT_THROW(h.propose({1.,2.},{2.},{3.,4.},{5.,6.},1e-8,1e-8),std::invalid_argument);
 });
 return run_all();
}
