// VMFL036 axisymmetric FVM kernel test: verification of the discrete
// operators of the axisymmetric physics against analytic fields.
//
// Fail-closed contract (audit 2026-10-16): every numerical check throws
// std::runtime_error instead of using assert(), because the Release CI
// configuration defines NDEBUG and assert() is compiled out - the previous
// assert-based version was silently inert in Release. A literal-escape
// corruption had additionally swallowed the radial viscous check into a
// comment, so it never executed in any configuration. The stdout markers
// below are preserved byte-for-byte for the evidence pipeline. The test
// remains at kernel level (metrics, discrete operators, closures); the
// SIMPLE loop and the Cd oracles are covered by
// test_vmfl036_axisymmetric.cpp.
#include "cfdx/physics/axisymmetric_fvm.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool ok, const std::string& what)
{
    if(!ok) throw std::runtime_error(what);
}
}

int main()
{
    using namespace cfdx::physics::axisymmetric;
    try {
        const auto re100=make_vmfl036_re100();
        validate_definition(re100);
        require(std::abs(re100.re()-100.0)<1e-12,
                "VMFL036 Re=100 definition reports a different Reynolds number");
        require(std::abs(re100.reference_cd()-1.0895)<1e-12,
                "VMFL036 Re=100 definition reports a different reference Cd (Tabata-Itakura 1.0895)");

        const auto fluent=make_vmfl036_fluent_exact();
        validate_definition(fluent);
        require(std::abs(fluent.re()-50.0)<1e-12,
                "VMFL036 Re=50 definition reports a different Reynolds number");
        require(std::abs(fluent.reference_cd()-1.5785)<1e-12,
                "VMFL036 Re=50 definition reports a different reference Cd (Tabata-Itakura 1.5785)");

        require(std::abs(revolution_face_measure(1.0,1.0)-2.0*pi)<1e-12,
                "VMFL036 revolution_face_measure identity failed");
        require(std::abs(revolution_volume(1.0,1.0)-2.0*pi)<1e-12,
                "VMFL036 revolution_volume identity failed");
        require(std::abs(cylindrical_divergence(2.0,1.0,3.0,2.0)-4.5)<1e-12,
                "VMFL036 cylindrical_divergence identity failed");

        const auto visc=cylindrical_momentum_viscous(
            0.01, 1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 2.0, 2.0);
        require(std::isfinite(visc.axial),
                "VMFL036 viscous axial term is not finite on the analytic field");
        // Restored check: the literal-escape corruption had swallowed this
        // assertion into the preceding comment line, so the radial viscous
        // identity never executed. On the analytic field:
        // d2ur_dx2 + d2ur_dr2 + (1/r)dur_dr - ur/r^2 = 3 + 4 + 2.5 - 0.5 = 9.0
        // and the cylindrical viscous term is 0.01 * 9 = 0.09.
        require(std::abs(visc.radial-0.09)<1e-12,
                "VMFL036 cylindrical viscous radial term != 0.09 on the analytic field");

        check_axis_regularity(0.0,0.0,0.0);
        bool rejected=false;
        try { check_axis_regularity(1e-3,0.0,0.0); } catch (...) { rejected=true; }
        require(rejected,
                "VMFL036 axis regularity did not reject a non-zero radial velocity on the axis");

        PressureFace pf{0,1,2.0,0.5,1.0};
        require(std::abs(pressure_correction_coefficient(pf,1.0,0.25)-1.0)<1e-12,
                "VMFL036 pressure_correction_coefficient identity failed");

        // The generic force evaluator is the sole force-integration implementation.
        using cfdx::physics::forces::AxisymmetricSample;
        std::vector<AxisymmetricSample> samples{
            {0.0,0.25,-1.0,0.0,1.0,1.0,-0.5,0.0,0.0},
            {0.0,0.50,-1.0,0.0,1.0,2.0,-1.0,0.0,0.0}
        };
        const auto drag=cfdx::physics::forces::integrate_axisymmetric(
            samples,1.0,1.0,pi/4.0);
        const double expected_p=2.0*pi*(0.25+1.0);
        const double expected_v=2.0*pi*(0.125+0.5);
        require(std::abs(drag.pressure_force-expected_p)<1e-12,
                "VMFL036 axisymmetric pressure force integration failed");
        require(std::abs(drag.viscous_force-expected_v)<1e-12,
                "VMFL036 axisymmetric viscous force integration failed");
        require(std::abs(drag.total_force-(expected_p+expected_v))<1e-12,
                "VMFL036 axisymmetric total force integration failed");

        std::cout << "AXISYMMETRIC_VMFL036_KERNEL: PASS\n";
        std::cout << "VMFL036_RE100 Re=" << re100.re()
                  << " Cd_reference=" << re100.reference_cd() << "\n";
        std::cout << "VMFL036_FLUENT_EXACT Re=" << fluent.re()
                  << " Cd_reference=" << fluent.reference_cd() << "\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "AXISYMMETRIC_VMFL036_KERNEL: FAIL: " << e.what() << "\n";
        return 1;
    }
}
