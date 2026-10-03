#include "cfdx/core/numerics/backend_equivalence.h"
#include "cfdx/core/numerics/temporal.h"
#include "cfdx/core/parallel/restart_mapping.h"

#include <cassert>
#include <cstdio>
#include <vector>

int main()
{
    using namespace cfdx::core;
    using namespace cfdx::core::numerics;
    using namespace cfdx::core::parallel;

    const std::vector<std::uint64_t> ids{10, 20, 30, 40};
    const std::vector<double> direct{1.0, 2.0, 3.0, 4.0};
    require_numerically_equivalent("steady restart", direct, direct, 0.0, 0.0);

    const auto remapped = remap_restart_values(ids, direct, {30, 10, 40, 20});
    require_numerically_equivalent(
        "N-to-M field mapping", {3.0, 1.0, 4.0, 2.0}, remapped, 0.0, 0.0);

    TimeIntegrationContext history(2, 1, "phi");
    Field<double, Location::CELL> initial(2, "phi", "1", 1);
    initial(0) = 10.0;
    initial(1) = 20.0;
    history.initialize(initial);

    Field<double, Location::CELL> current(2, "phi", "1", 1);
    current(0) = 11.0;
    current(1) = 21.0;
    history.shift(current);

    require_numerically_equivalent(
        "transient restart current", {11.0, 21.0},
        {history.phi_curr(0), history.phi_curr(1)}, 0.0, 0.0);
    require_numerically_equivalent(
        "transient restart previous", {10.0, 20.0},
        {history.phi_prev(0), history.phi_prev(1)}, 0.0, 0.0);
    assert(history.has_prev);

    std::printf("N13 restart equivalence PASS\n");
    return 0;
}
