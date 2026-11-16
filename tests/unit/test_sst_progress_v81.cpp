#include "cfdx/physics/sst_progress_v81.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace cfdx::physics;
static void need(bool condition,const char* text) {
    if (!condition) throw std::runtime_error(text);
}
int main() {
    const auto result=sst_progress_distance_v81(1e-2,1e3,1e-3,1e2,1e-8);
    need(!result.all_gates_met,"V81 must not promote unconverged state");
    need(std::abs(result.worst_factor-1e11)<1e-4*1e11,"omega gate ratio");
    const auto yes=sst_progress_distance_v81(1e-9,1e-9,0.,0.,1e-8);
    need(yes.all_gates_met && yes.worst_factor<1.,"successful arithmetic only");
    need(std::abs(sst_viscous_wall_omega_v81(1.56e-5,0.075,
        5.808688154e-6,6.)-3.6987796735795e7)<10.,"wall scaling");
    bool threw=false;
    try { (void)sst_progress_distance_v81(0.,0.,0.,0.,0.); }
    catch(const std::invalid_argument&) {threw=true;}
    need(threw,"zero tolerance must fail");
    threw=false;
    try { (void)sst_progress_distance_v81(
        std::numeric_limits<double>::quiet_NaN(),0.,0.,0.,1e-8); }
    catch(const std::invalid_argument&) {threw=true;}
    need(threw,"nonfinite must fail");
    threw=false;
    try { (void)sst_viscous_wall_omega_v81(1.56e-5,.075,0.,6.); }
    catch(const std::invalid_argument&) {threw=true;}
    need(threw,"zero y must fail");
    std::cout<<"SST_PROGRESS_V81_PASS gate_factor="<<result.worst_factor<<"\n";
}