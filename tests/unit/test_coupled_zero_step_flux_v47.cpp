#include "cfdx/physics/coupled_zero_step_flux_v47.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <limits>
using namespace cfdx::physics;
static void need(bool v) {if(!v) throw std::runtime_error("V47 zero-step flux test");}
int main() {
    // Two cells, one shared internal face; consistent owner-neighbour signs.
    const std::vector<std::size_t> owner{0,0,1};
    const std::vector<int> neighbour{-1,1,-1};
    const std::vector<double> old{-2,1,1};
    const std::vector<double> remapped{-2,1.25,1};
    auto evaluate=[&](const auto& p){
        return audit_coupled_zero_step_flux_v47(2,3,
            [&](std::size_t f){return owner[f];},
            [&](std::size_t f){return neighbour[f];},
            [&](std::size_t f){return old[f];},
            [&](std::size_t f){return p[f];});
    };
    auto identical=evaluate(old);
    need(identical.divergence_jump_l1==0.0);
    need(identical.relative_face_flux_change==0.0);
    auto change=evaluate(remapped);
    need(std::abs(change.divergence_jump_l1-0.125)<1e-14);
    need(std::abs(change.relative_face_flux_change-0.0625)<1e-14);
    need(change.maximum_face_flux_change==0.25);
    need(change.worst_cell==0);
    bool bad=false;
    try {
      auto broken=remapped;broken[0]=std::numeric_limits<double>::infinity();
      (void)evaluate(broken);
    } catch(const std::invalid_argument&) {bad=true;}
    need(bad);
    std::cout<<"COUPLED_ZERO_STEP_FLUX_V47_PASS owner_neighbour=1 zero_step=1 nonfinite_rejected=1\n";
}
