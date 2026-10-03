#include "cfdx/core/numerics/backend_equivalence.h"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

int main()
{
    using cfdx::core::numerics::compare_vectors;
    using cfdx::core::numerics::numerically_equivalent;

    const std::vector<double> reference{1.0, -2.0, 4.0, 8.0};
    const std::vector<double> candidate{1.0 + 1e-14, -2.0, 4.0 - 2e-14, 8.0};
    const auto m = compare_vectors(reference, candidate);

    assert(m.size == reference.size());
    assert(m.l2_abs > 0.0);
    assert(m.linf_abs > 0.0);
    assert(m.l2_relative < 1e-13);
    assert(numerically_equivalent(m, 1e-12, 1e-12));

    const std::vector<double> different{1.0, -2.0, 4.0, 8.1};
    const auto bad = compare_vectors(reference, different);
    assert(!numerically_equivalent(bad, 1e-12, 1e-12));

    bool rejected = false;
    try {
        (void)compare_vectors(reference, {1.0, 2.0});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    rejected = false;
    try {
        (void)numerically_equivalent(m, -1.0, 1e-12);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    return 0;
}
