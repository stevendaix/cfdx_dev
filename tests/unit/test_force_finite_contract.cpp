#include "cfdx/physics/forces.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void require(bool ok, const char* case_name) {
    if (!ok) throw std::runtime_error(case_name);
}
template <class F>
void must_reject(F f, const char* case_name) {
    bool rejected = false;
    try { f(); }
    catch (const std::invalid_argument&) { rejected = true; }
    catch (const std::overflow_error&) { rejected = true; }
    catch (const std::domain_error&) { rejected = true; }
    require(rejected, case_name);
}
}

int main() {
    using namespace cfdx::physics::forces;
    SurfaceSample s;
    s.normal = {1.0, 0.0, 0.0};
    s.position = {0.0, 1.0, 0.0};
    s.area = 1.0;
    s.pressure = 2.0;
    Reference ref;
    ref.rho_inf = 1.0;
    ref.u_inf = 1.0;
    ref.area = 1.0;
    ref.length = 1.0;
    const auto nominal = integrate({s}, ref);
    require(std::abs(nominal.cd + 4.0) < 1e-12,
            "nominal aerodynamic drag sign or value changed");

    must_reject([&] {
        auto bad = ref;
        bad.moment_center.x = std::numeric_limits<double>::quiet_NaN();
        (void)integrate({s}, bad);
    }, "non-finite moment origin must be rejected");
    must_reject([&] {
        auto bad = ref;
        bad.rho_inf = 1e-300;
        bad.area = 1e-300;
        (void)integrate({s}, bad);
    }, "force reference q*A underflow must be rejected");
    must_reject([&] {
        auto bad = ref;
        bad.area = 1e308;
        bad.length = 1e308;
        (void)integrate({s}, bad);
    }, "moment reference q*A*L overflow must be rejected");
    must_reject([&] {
        auto bad = s;
        bad.area = 1e308;
        bad.pressure = 1e308;
        (void)integrate({bad}, ref);
    }, "surface pressure force overflow must be rejected");

    AxisymmetricSample axis;
    axis.r = 0.5;
    axis.ds = 0.25;
    axis.nx = -1.0;
    axis.nr = 0.0;
    axis.pressure = 1.0;
    const auto axis_nominal = integrate_axisymmetric({axis},1.0,1.0,1.0);
    require(std::isfinite(axis_nominal.cd_total) && axis_nominal.cd_total > 0.0,
            "nominal axisymmetric force changed");
    must_reject([&] {
        (void)integrate_axisymmetric({axis}, 1e-300, 1.0, 1e-300);
    }, "axisymmetric reference q*A underflow must be rejected");
    must_reject([&] {
        auto bad = axis;
        bad.r = 1e150;
        bad.ds = 1e5;
        bad.pressure = 1e200;
        (void)integrate_axisymmetric({bad},1.0,1.0,1.0);
    }, "axisymmetric pressure force overflow must be rejected");
    std::cout << "FORCE_FINITE_CONTRACT: PASS\n";
}
