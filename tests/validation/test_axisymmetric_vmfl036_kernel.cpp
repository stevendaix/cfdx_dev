#include "cfdx/physics/axisymmetric_fvm.h"
#include <cassert>
#include <cmath>
#include <iostream>

int main()
{
    using namespace cfdx::physics::axisymmetric;

    const auto re100=make_vmfl036_re100();
    validate_definition(re100);
    assert(std::abs(re100.re()-100.0)<1e-12);
    assert(std::abs(re100.reference_cd()-1.0895)<1e-12);

    const auto fluent=make_vmfl036_fluent_exact();
    validate_definition(fluent);
    assert(std::abs(fluent.re()-50.0)<1e-12);
    assert(std::abs(fluent.reference_cd()-1.0875)<1e-12);

    assert(std::abs(revolution_face_measure(1.0,1.0)-2.0*pi)<1e-12);
    assert(std::abs(revolution_volume(1.0,1.0)-2.0*pi)<1e-12);
    assert(std::abs(cylindrical_divergence(2.0,1.0,3.0,2.0)-4.5)<1e-12);

    const auto visc=cylindrical_momentum_viscous(
        0.01, 1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 2.0, 2.0);
    assert(std::isfinite(visc.axial));
    // d2ur_dx2 + d2ur_dr2 + (1/r)dur_dr - ur/r^2 = 3 + 4 + 2.5 - 0.5 = 9.0\n    // The cylindrical viscous term is therefore 0.01 * 9 = 0.09.\n    assert(std::abs(visc.radial-0.09)<1e-12);

    check_axis_regularity(0.0,0.0,0.0);
    bool rejected=false;
    try { check_axis_regularity(1e-3,0.0,0.0); } catch (...) { rejected=true; }
    assert(rejected);

    PressureFace pf{0,1,2.0,0.5,1.0};
    assert(std::abs(pressure_correction_coefficient(pf,1.0,0.25)-1.0)<1e-12);

    // The generic force evaluator is the sole force-integration implementation.
    using cfdx::physics::forces::AxisymmetricSample;
    std::vector<AxisymmetricSample> samples{
        {0.0,0.25,-1.0,0.0,0.1,1.0,-0.5,0.0,0.0},
        {0.0,0.50,-1.0,0.0,0.1,2.0,-1.0,0.0,0.0}
    };
    const auto drag=cfdx::physics::forces::integrate_axisymmetric(
        samples,1.0,1.0,pi/4.0);
    const double expected_p=2.0*pi*(0.25+1.0);
    const double expected_v=2.0*pi*(0.125+0.5);
    assert(std::abs(drag.pressure_force-expected_p)<1e-12);
    assert(std::abs(drag.viscous_force-expected_v)<1e-12);
    assert(std::abs(drag.total_force-(expected_p+expected_v))<1e-12);

    std::cout << "AXISYMMETRIC_VMFL036_KERNEL: PASS\n";
    std::cout << "VMFL036_RE100 Re=" << re100.re()
              << " Cd_reference=" << re100.reference_cd() << "\n";
    std::cout << "VMFL036_FLUENT_EXACT Re=" << fluent.re()
              << " Cd_reference=" << fluent.reference_cd() << "\n";
    return 0;
}
