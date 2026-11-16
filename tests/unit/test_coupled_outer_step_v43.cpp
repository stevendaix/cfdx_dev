#include "cfdx/physics/coupled_outer_step_v43.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace cfdx::physics;
static void need(bool value,const char* message) {
    if(!value) throw std::runtime_error(message);
}
int main() {
    try {
        need(coupled_outer_beta_v43(1.,1.)==1.,"full candidate changed");
        need(coupled_outer_beta_v43(1.,0.5)==0.5,"retry does not shrink step");
        need(coupled_outer_beta_v43(0.5,0.5)==0.25,"initial damping ignored");
        need(coupled_outer_blend_v43(4.,8.,0.25)==5.,"physical field not interpolated");
        need(coupled_outer_blend_v43(7.,-1.,1.)==-1.,"full candidate not exact");
        bool invalid=false;
        try { (void)coupled_outer_beta_v43(0.0,1.0); }
        catch(const std::invalid_argument&){invalid=true;}
        need(invalid,"zero beta accepted");
        invalid=false;
        try { (void)coupled_outer_blend_v43(std::numeric_limits<double>::infinity(),1.,0.5); }
        catch(const std::invalid_argument&){invalid=true;}
        need(invalid,"nonfinite field accepted");
        invalid=false;
        try { (void)coupled_outer_beta_v43(1.0,1.1); }
        catch(const std::invalid_argument&){invalid=true;}
        need(invalid,"extrapolated retry accepted");
        std::cout<<"COUPLED_OUTER_STEP_V43_PASS full=1 damping=1 finite=1\n";
        return 0;
    }catch(const std::exception& error) {
        std::cerr<<"COUPLED_OUTER_STEP_V43_FAIL "<<error.what()<<'\n';
        return 1;
    }
}
