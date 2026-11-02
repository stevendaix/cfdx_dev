#include "cfdx/physics/wall_distance_turbulence.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(what);
}
template<class F> void rejects(F f, const char* what) {
    bool threw=false;
    try { f(); }
    catch(const std::invalid_argument&) { threw=true; }
    catch(const std::domain_error&) { threw=true; }
    catch(const std::overflow_error&) { threw=true; }
    require(threw,what);
}
}

int main() {
    using namespace cfdx::physics;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const auto y = wall_distance_y_plus({0.001}, {0.1},1e-5);
    require(std::abs(y[0]-10) < 1e-12, "nominal y+ changed");
    rejects([&] {(void)wall_distance_y_plus({1e300}, {1e200}, 1.0);},
            "non-representable y+ must fail closed");
    // This intermediate product overflows in double, but the correct final
    // value is finite. The fallback must not reject a representable y+.
    const auto large = wall_distance_y_plus({1e250}, {1e250}, 1e250);
    require(std::isfinite(large[0]) && std::abs(large[0]/1e250-1.0)<1e-12,
            "representable y+ after intermediate overflow");
    const auto tiny = wall_distance_y_plus({1e-200}, {1e-200}, 1e-200);
    require(tiny[0] > 0.0 && std::abs(tiny[0]/1e-200-1.0)<1e-12,
            "representable y+ after intermediate underflow");
    rejects([&] {(void)des_wall_transition_ratio(.1,.2,nan);},
            "DES NaN constant must fail closed");
    rejects([&] {(void)des_wall_transition_ratio(.1,.2,inf);},
            "DES inf constant must fail closed");
    rejects([&] {(void)audit_wall_distance_turbulence_coupling(
        {.01},{.1},1e-5,.2,nan);},
            "audit must reject non-finite DES constant");
    require(std::abs(des_wall_transition_ratio(.13,.1,.65)-2)<1e-12,
            "nominal DES transition changed");

    WallDistanceGrid g;
    g.points.resize(2);
    g.solid.resize(1, 0);
    rejects([&] {require_wall_distance_resolves_surface(g,{1.0,2.0});},
            "malformed solid mask must not be indexed out of bounds");
    g.solid.resize(2,0);
    require_wall_distance_resolves_surface(g,{1.0,2.0});
    std::cout << "WALL_TURBULENCE_CONTRACT_GUARD: PASS\n";
}
