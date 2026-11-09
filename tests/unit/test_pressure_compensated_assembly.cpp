#include "cfdx/physics/pressure_compensated_assembly.h"
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

using cfdx::physics::PressureFaceAssemblyAccumulator;
int main() {
    if (std::numeric_limits<long double>::digits <=
        std::numeric_limits<double>::digits) {
        bool failed = false;
        try { PressureFaceAssemblyAccumulator bad(2); }
        catch (const std::invalid_argument&) { failed = true; }
        assert(failed);
        return 0;
    }
    // Internal owner/neighbour coefficient contributes equally to each
    // diagonal; the off-diagonal contribution is handled by the FV caller.
    PressureFaceAssemblyAccumulator a(2);
    a.add_internal_pair(0, 1, 4.0);
    a.add_pressure_boundary(1, 2.0, 3.0);
    a.add_continuity(0, 1.5);
    a.add_continuity(1, -1.5);
    assert(a.diagonal(0) == 4.0);
    assert(a.diagonal(1) == 6.0);
    assert(a.rhs(0) == -1.5);
    assert(a.rhs(1) == 4.5);

    // Alternating huge/small predictor fluxes are accumulated before one
    // final rounding; the legacy double path loses the middle +1 term.
    PressureFaceAssemblyAccumulator large(2);
    large.add_continuity(0, 1.0e16);
    large.add_continuity(0, 1.0);
    large.add_continuity(0, -1.0e16);
    assert(large.rhs(0) == -1.0);
    double legacy = 0.0;
    legacy += 1.0e16; legacy += 1.0; legacy -= 1.0e16;
    assert(legacy == 0.0);
    // Zero-mean pairs cannot manufacture a global source.
    PressureFaceAssemblyAccumulator pair(2);
    pair.add_continuity(0, 2.25);
    pair.add_continuity(1, -2.25);
    assert(pair.rhs(0) + pair.rhs(1) == 0.0);
    bool bad = false;
    try { pair.add_internal_pair(0, 1, -1.0); }
    catch (const std::invalid_argument&) { bad = true; }
    assert(bad);
    bad = false;
    try { pair.add_pressure_boundary(1, 1.0,
        std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::invalid_argument&) { bad = true; }
    assert(bad);
    return 0;
}
