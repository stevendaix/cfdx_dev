#include "cfdx/physics/sst_solver.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

int main() {
    using namespace cfdx::physics;
    using namespace cfdx::core;
    Field<double, Location::CELL> ka(3, "k_old", "m2/s2", 1);
    Field<double, Location::CELL> kb(3, "k_new", "m2/s2", 1);
    Field<double, Location::CELL> wa(3, "omega_old", "1/s", 1);
    Field<double, Location::CELL> wb(3, "omega_new", "1/s", 1);
    for (int i=0;i<3;++i) {
        ka(i) = 1.0e-3;
        kb(i) = 1.001e-3;
        wa(i) = 1.0e5;
        wb(i) = 1.001e5;
    }
    const auto k = sst_field_update_diagnostic(ka,kb,1e-12);
    const auto w = sst_field_update_diagnostic(wa,wb,1e-12);
    if (!(std::abs(k.relative_global_inf-w.relative_global_inf) < 1e-14)) return 1;
    if (!(std::abs(k.absolute_inf-w.absolute_inf) > 1.0)) return 2;
    if (!(k.reference_scale < 0.01 && w.reference_scale > 1e4)) return 3;
    // Rescaling the units of BOTH iterates must not alter the relative change.
    for (int i=0;i<3;++i) { wa(i) *= 60; wb(i) *= 60; }
    if (!(std::abs(sst_field_update_diagnostic(wa,wb,1e-12).relative_global_inf-
                   w.relative_global_inf) < 1e-14)) return 4;
    wb(2) = std::numeric_limits<double>::quiet_NaN();
    bool rejected=false;
    try { (void)sst_field_update_diagnostic(wa,wb,1e-12); }
    catch(const std::invalid_argument&) { rejected=true; }
    if (!rejected) return 5;
    std::cout << "PASS SST V10 dimensional diagnostics (original nonlinear gate retained)\n";
}
