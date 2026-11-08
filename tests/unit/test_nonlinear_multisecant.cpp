#include "cfdx/physics/nonlinear_multisecant_acceleration.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#define REQUIRE(cond) do { if (!(cond)) throw std::runtime_error("REQUIRE: " #cond); } while(false)
using namespace cfdx::physics;

int main() {
    for (auto method : {NonlinearMultisecantMethod::ANDERSON,
                        NonlinearMultisecantMethod::NGMRES}) {
        NonlinearMultisecantControls c;
        c.method=method; c.mixing=1.0; c.maximum_step_multiplier=20.0;
        NonlinearMultisecantHistory history(c);
        // Scalar fixed-point map G(x)=0.8-0.6x, whose exact fixed point is 0.5.
        // A synthetic conservative two-face flux [q,-q] is included in the
        // same nonlinear state to ensure the combination is a single affine
        // operation on U, p and phi (NOT independent per-face weights).
        auto make_state=[](double u) { return std::vector<double>{u,2*u,1+u,-1-u}; };
        auto g=[](double u) {return 0.8-0.6*u;};
        const std::vector<double> scales{1,2,2,2};
        for (double input : {0.0,0.8,0.32,0.608}) {
            auto s=make_state(input), next=make_state(g(input));
            std::vector<double> physical{g(input)-0.5,2*(g(input)-0.5)};
            history.append({s,next,physical},scales);
        }
        const auto proposal=history.propose();
        REQUIRE(proposal.available);
        REQUIRE(proposal.state.size()==4);
        REQUIRE(std::abs(proposal.state[0]-0.5)<5e-5);
        REQUIRE(std::abs(proposal.state[1]-2*proposal.state[0])<1e-12);
        REQUIRE(std::abs(proposal.state[2]+proposal.state[3])<1e-12);
        double weight_sum=0;
        for (double w : proposal.coefficients) weight_sum+=w;
        REQUIRE(std::abs(weight_sum-1)<1e-12);
    }
    // Near-collinear history: regularization must either stabilize or
    // explicitly refuse the proposal; never return a nonfinite state.
    {
        NonlinearMultisecantControls c;
        c.method=NonlinearMultisecantMethod::ANDERSON;
        NonlinearMultisecantHistory history(c);
        history.append({{1.0},{2.0},{1.0}},{1.0});
        history.append({{1.0},{2.0},{1.0}},{1.0});
        const auto result=history.propose();
        if (result.available)
            for (double v: result.state) REQUIRE(std::isfinite(v));
        history.reset(); REQUIRE(history.size()==0);
    }
    // Central finite-difference Jv must see nonlinear terms and mixed signs;
    // it cannot silently substitute the UPWIND/frozen Jacobian.
    auto f=[](const std::vector<double>& x) {
        return std::vector<double>{x[0]*x[0]+3*x[1],
                                   std::sin(x[1])-x[0]*x[1]};
    };
    const auto jv=nonlinear_directional_jvp(f,{0.7,-0.4},{0.2,-0.3});
    REQUIRE(std::abs(jv[0]-(1.4*0.2-0.9))<1e-8);
    REQUIRE(std::abs(jv[1]-(-(-0.4)*0.2+(std::cos(-0.4)-0.7)*(-0.3)))<1e-8);
    bool rejected=false;
    try { (void)nonlinear_directional_jvp(f,{1.,2.},{0.,0.}); }
    catch (const std::invalid_argument&) {rejected=true;}
    REQUIRE(rejected);
    rejected=false;
    try {NonlinearMultisecantControls bad;bad.depth=11;
         validate_nonlinear_multisecant(bad);}
    catch (const std::invalid_argument&) {rejected=true;}
    REQUIRE(rejected);
    rejected=false;
    try {NonlinearMultisecantControls cfg;
         NonlinearMultisecantHistory hist(cfg);
         hist.append({{1.},{2.},{0.}}, {0.});}
    catch (const std::invalid_argument&) {rejected=true;}
    REQUIRE(rejected);
    std::cout << "PASS Anderson/NGMRES coefficients, conservative affine flux, "
                 "directional nonlinear Jv and fail-closed safeguards\n";
}
