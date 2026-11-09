#include "cfdx/core/numerics/portable_double_double.h"
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
using cfdx::core::numerics::PortableDoubleDouble;
int main() {
    // 1: Compensate a term that would be lost by binary64 summation.
    PortableDoubleDouble a;
    a.add(1e16); a.add(1); a.add(-1e16);
    assert(a.narrow_to_double() == 1);
    // 2: The very same pattern also works for negative pressure residuals.
    a = {}; a.add(-1e16); a.add(-1); a.add(1e16);
    assert(a.narrow_to_double() == -1);
    // 3: Pair merging preserves increments across cancellation.
    PortableDoubleDouble b;
    b.add(1e16); b.add(1);
    PortableDoubleDouble c;
    c.add(-1e16); c.add(2);
    b.add(c);
    assert(b.narrow_to_double() == 3);
    // 4: An initially zero pair has a real zero representation.
    b = {}; b.add(0.25); b.add(-0.25);
    assert(b.high == 0 && b.low == 0);
    // 5: Invalid fluxes must fail rather than contaminate the accumulator.
    bool caught = false;
    try { b.add(std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::invalid_argument&) { caught=true; }
    assert(caught);
    // 6: An overflow in a diagonal is explicitly rejected.
    b = {}; b.add(std::numeric_limits<double>::max());
    caught = false;
    try { b.add(std::numeric_limits<double>::max()); }
    catch (const std::overflow_error&) { caught=true; }
    assert(caught);
}
