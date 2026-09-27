#include "cfdx/physics/forces.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

int main()
{
    using namespace cfdx::physics::forces;

    // Flat wall with outward body normal +x. Pressure acts on the body in -x.
    SurfaceSample s;
    s.position={1.0,0.0,0.0};
    s.normal={1.0,0.0,0.0};
    s.area=2.0;
    s.pressure=3.0;
    s.viscous_traction={0.5,0.0,0.0};

    Reference ref;
    ref.rho_inf=1.0;
    ref.u_inf=2.0;
    ref.area=2.0;
    ref.drag_direction={-1.0,0.0,0.0};
    ref.moment_center={0.0,0.0,0.0};

    const auto r=integrate({s},ref);
    assert(std::abs(r.pressure_force.x+6.0)<1e-12);
    assert(std::abs(r.viscous_force.x-1.0)<1e-12);
    assert(std::abs(r.drag-5.0)<1e-12);
    assert(std::abs(r.cd-5.0/(0.5*1.0*4.0*2.0))<1e-12);

    // Axisymmetric meridional traction: sphere normal points body -> fluid.
    AxisymmetricSample a;
    a.r=0.5; a.ds=0.25; a.nx=-1.0; a.nr=0.0;
    a.pressure=2.0; a.tau_xx=0.4; a.tau_xr=0.0; a.tau_rr=0.0;
    const auto ar=integrate_axisymmetric({a},1.0,1.0,3.14159265358979323846/4.0);
    const double dS=2.0*3.14159265358979323846*0.5*0.25;
    assert(std::abs(ar.pressure_force-2.0*dS)<1e-12);
    assert(std::abs(ar.viscous_force+0.4*dS)<1e-12);
    assert(std::abs(ar.total_force-1.6*dS)<1e-12);
    const auto ar_ref=integrate_axisymmetric({a},1.0,1.0,3.14159265358979323846/4.0,1.0);
    assert(std::abs(ar_ref.pressure_force-1.0*dS)<1e-12);

    // Reject invalid reference data rather than silently producing invalid coefficients.
    bool rejected=false;
    try {
        Reference bad=ref; bad.area=std::numeric_limits<double>::infinity();
        (void)integrate({s},bad);
    } catch (...) { rejected=true; }
    assert(rejected);

    rejected=false;
    try {
        AxisymmetricSample bad=a; bad.ds=std::numeric_limits<double>::quiet_NaN();
        (void)integrate_axisymmetric({bad},1.0,1.0,3.14159265358979323846/4.0);
    } catch (...) { rejected=true; }
    assert(rejected);

    std::cout<<"FORCES_INTEGRATION: PASS\n";
    return 0;
}
