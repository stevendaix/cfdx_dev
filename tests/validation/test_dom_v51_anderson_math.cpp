// V51 small, sanitizer-friendly numerical kernel test.
// This covers multisecant coefficients, degeneracy and bounds but not the
// full angular transport library (which requires separate instrumentation).
#include "cfdx/physics/dom_holo_global_v51.h"
#include <cmath>
#include <deque>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
namespace {
void require(bool p,const char* message){if(!p)throw std::runtime_error(message);}
}
int main() {
  try {
    using dom_v51_detail::Secant;
    // Exact affine fixed point for gas and reflecting-wall moments.
    std::deque<Secant> history;
    history.push_back({{0.,0.},{1.,2.}});
    history.push_back({{1.,2.},{.2,.4}});
    const auto step=dom_v51_detail::direction(history,1e-12,30.);
    require(step.size()==2,"V51 math: missing joint multisecant direction");
    require(std::abs(step[0]-.25)<1e-10 &&
            std::abs(step[1]-.5)<2e-10,
            "V51 math: incorrect two-variable fixed-point acceleration");
    require(std::abs(1.+step[0]-1.25)<1e-10 &&
            std::abs(2.+step[1]-2.5)<2e-10,
            "V51 math: predicted gas/wall solution mismatch");
    history.push_back(history.back());
    require(dom_v51_detail::direction(history,1e-10,30.).empty() == false,
            "V51 math: regularized rank-deficient least squares failed");
    history.assign(4,Secant{{1.,2.},{.1,.2}});
    require(dom_v51_detail::direction(history,1e-10,30.).empty(),
            "V51 math: zero-rank history not rejected");
    history={{ {0.},{1.} },{{1.},{.1}}};
    require(dom_v51_detail::direction(history,1e-12,1.).empty(),
            "V51 math: oversize multisecant direction not rejected");
    history.push_back({{1.,2.},{.1,.2}});
    require(dom_v51_detail::direction(history,1e-12,30.).empty(),
            "V51 math: inconsistent history lengths not rejected");
    std::cout<<"V51_MATH_ALL_PASS\n";
  } catch(const std::exception& e) {
    std::cerr<<"V51_MATH_FAIL: "<<e.what()<<'\n';return 1;
  }
}
