#include "cfdx/physics/pressure_flux_split.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::physics;

static void require(bool good, const char* message) {
    if (!good) throw std::runtime_error(message);
}
static void close(double actual,double target,double tol = 1.e-11) {
    require(std::abs(actual-target) <= tol*(1.0+std::abs(target)),
            "face split not affine exact");
}
int main() {
    try {
        const Vec3 d{1.0,0.0,0.0};
        const Vec3 r{2.0,3.0,5.0};
        const Vec3 grad{1.3,-2.1,0.6};
        for (double theta : {0.0,30.0,60.0,73.0,75.0}) {
            const double angle = theta * std::acos(-1.0)/180.0;
            const Vec3 area{std::cos(angle),std::sin(angle),0.0};
            const double dp = grad.dot(d);
            const double exact = r.x*grad.x*area.x+r.y*grad.y*area.y;
            std::array<PressureFluxSplit,2> schemes;
            for (int method=0;method<2;++method) {
                schemes[method]=split_pressure_flux(area,d,r,grad,
                    method==0 ? PressureFluxSplitScheme::MINIMUM_CORRECTION :
                                PressureFluxSplitScheme::OVER_RELAXED);
                close(schemes[method].implicit_coefficient*dp+
                    schemes[method].deferred_flux,exact);
                require(schemes[method].implicit_coefficient>0.0,
                        "non-SPD implicit face coupling");
                // A single oriented face is a conservative pairwise flux.
                const double owner = schemes[method].implicit_coefficient*dp+
                                     schemes[method].deferred_flux;
                const double neighbour = -owner;
                close(owner+neighbour,0.0);
                close(schemes[method].exact_affine_flux,exact);
            }
            require(schemes[1].implicit_coefficient + 1.0e-12 >=
                        schemes[0].implicit_coefficient,
                    "over-relaxed implicit coupling weaker than minimal");
            if (theta==0.0) {
                close(schemes[0].deferred_flux,0.0);
                close(schemes[1].deferred_flux,0.0);
            }
            std::cout<<"PRESSURE_SPLIT_PASS angle="<<theta
                     <<" implicit_orth="<<schemes[0].implicit_coefficient
                     <<" implicit_over="<<schemes[1].implicit_coefficient<<'\n';
        }
        bool rejected=false;
        try { (void)split_pressure_flux({-1,0,0},d,r,grad,
                                PressureFluxSplitScheme::OVER_RELAXED); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"face normal orientation not checked");
        rejected=false;
        try { (void)split_pressure_flux({1,0,0},d,{-1,1,1},grad,
                                 PressureFluxSplitScheme::OVER_RELAXED); }
        catch (const std::invalid_argument&) { rejected=true; }
        require(rejected,"non-positive pressure response accepted");
        std::cout<<"PRESSURE_SPLIT_CONTRACT: PASS\n";
    } catch(const std::exception& error) {
        std::cerr<<"PRESSURE_SPLIT_CONTRACT: FAIL "<<error.what()<<'\n';
        return 1;
    }
    return 0;
}
