#include "cfdx/physics/coupled_composite_backtracking_v47.h"
#include "cfdx/physics/coupled_outer_step_v43.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>
using namespace cfdx::physics;
static void need(bool ok) { if (!ok) throw std::runtime_error("V47 composite backtracking contract"); }
template<class Fn> void must_throw(Fn fn) {
    bool rejected=false;
    try { fn(); } catch (const std::invalid_argument&) { rejected=true; }
    need(rejected);
}
int main() {
    validate_coupled_composite_backtracking_v47(10,0.5);
    need(coupled_composite_beta_v47(0,0.5)==1.0);
    need(coupled_composite_beta_v47(1,0.5)==0.5);
    need(coupled_composite_beta_v47(4,0.5)==0.0625);
    // Unlike the existing segregated alpha_u floor, beta can continue below 0.1.
    need(coupled_composite_beta_v47(10,0.5)==1.0/1024.0);
    need(coupled_composite_beta_v47(10,0.5)<0.1);
    need(coupled_outer_blend_v43(4.0,8.0,coupled_composite_beta_v47(2,0.5))==5.0);
    CoupledCompositeSearchV47 search;
    need(search.trial_beta(3,0.5)==0.125);
    search.observe_rejection(1.,1.531782680995);
    search.observe_rejection(0.5,1.207290045632);
    search.observe_rejection(0.25,1.202203470437);
    const double candidate=search.trial_beta(3,0.5);
    need(candidate>0.35 && candidate<0.38);
    need(search.trial_beta(4,0.5)==0.0625);
    CoupledCompositeSearchV47 monotone;
    monotone.observe_rejection(1.0,4.0);
    monotone.observe_rejection(0.5,3.0);
    monotone.observe_rejection(0.25,2.0);
    need(monotone.trial_beta(3,0.5)==0.125);
    must_throw([]{validate_coupled_composite_backtracking_v47(0,0.5);});
    must_throw([]{validate_coupled_composite_backtracking_v47(13,0.5);});
    must_throw([]{validate_coupled_composite_backtracking_v47(10,1.0);});
    must_throw([]{validate_coupled_composite_backtracking_v47(10,0.1);});
    must_throw([]{(void)coupled_composite_beta_v47(13,0.5);});
    must_throw([]{(void)coupled_composite_beta_v47(1,std::numeric_limits<double>::quiet_NaN());});
    std::cout << "COUPLED_COMPOSITE_BACKTRACKING_V47_PASS beta_independent_of_alpha_u_floor=1\n";
}
