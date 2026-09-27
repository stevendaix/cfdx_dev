#include "cfdx/physics/axisymmetric_fvm.h"
#include <cmath>
#include <iostream>

int main()
{
    using namespace cfdx::physics::axisymmetric;
    const auto re100 = make_vmfl036_re100();
    const auto fluent = make_vmfl036_fluent_exact();
    validate_definition(re100);
    validate_definition(fluent);

    if (std::abs(re100.re()-100.0)>1e-12 || std::abs(fluent.re()-50.0)>1e-12)
        return 1;
    const double force = re100.reference_cd()*0.5*re100.rho*re100.velocity*re100.velocity*re100.reference_area();
    if (std::abs(force-0.4278456495)>1e-10)
        return 1;

    std::cout << "VMFL036_REFERENCE: PASS Re100_Re=" << re100.re()
              << " Re100_Cd=1.0895 Re100_F=" << force
              << " FluentExact_Re=" << fluent.re()
              << " FluentExact_Cd=1.0875\n";
    std::cout << "VMFL036_REFERENCE: oracle-only; physical validation is test_vmfl036_axisymmetric\n";
    return 0;
}
