#include "cfdx/physics/coupled_full_period_v44.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace cfdx::physics;
static void need(bool ok,const char* why) { if (!ok) throw std::runtime_error(why); }
int main() {
    try {
        auto x=make_coupled_sst_phase_candidate_v44(4.0,2000.0);
        need(x.k==4.0 && x.omega==2000.0,"SST phase average changed");
        x=make_coupled_sst_phase_candidate_v44(0.0,1e-10);
        need(x.k==0.0 && x.omega==1e-10,"valid SST low amplitude changed");
        for (int test=0; test<4; ++test) {
            bool rejected=false;
            try {
                (void)make_coupled_sst_phase_candidate_v44(
                    test==0 ? -1.0 : test==1 ? std::numeric_limits<double>::quiet_NaN() : 4.0,
                    test==2 ? 0.0 : test==3 ? std::numeric_limits<double>::infinity() : 2000.0);
            } catch (const std::invalid_argument&) {rejected=true;}
            need(rejected,"invalid phase-mean SST accepted");
        }
        std::cout<<"COUPLED_FULL_PERIOD_V44_PASS positivity=1 nonfinite=1 no_floors=1\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"COUPLED_FULL_PERIOD_V44_FAIL "<<error.what()<<'\n';
        return 1;
    }
}
