#include "cfdx/core/numerics/backend_equivalence.h"
#include "cfdx/core/parallel/restart_mapping.h"
#include "cfdx/core/numerics/temporal.h"

#include <cassert>
#include <cstdio>
#include <vector>

int main()
{
    using namespace cfdx::core::numerics;
    using namespace cfdx::core::parallel;

    const std::vector<std::uint64_t> ids{10, 20, 30, 40};
    const std::vector<double> direct{1.0, 2.0, 3.0, 4.0};
    const std::vector<double> checkpointed{1.0, 2.0, 3.0, 4.0};
    require_numerically_equivalent(
        "steady restart", direct, checkpointed, 0.0, 0.0);

    const auto remapped = remap_restart_values(ids, direct, {30, 10, 40, 20});
    const std::vector<double> expected{3.0, 1.0, 4.0, 2.0};
    require_numerically_equivalent(
        "N-to-M field mapping", expected, remapped, 0.0, 0.0);

    TemporalHistory history;
    history.push(1.0, {10.0, 20.0});
    history.push(2.0, {11.0, 21.0});
    history.push(3.0, {12.0, 22.0});
    const auto latest = history.current();
    const auto previous = history.previous();
    require_numerically_equivalent(
        "transient restart current", {12.0, 22.0}, latest, 0.0, 0.0);
    require_numerically_equivalent(
        "transient restart previous", {11.0, 21.0}, previous, 0.0, 0.0);

    std::printf("N13 restart equivalence PASS\n");
    return 0;
}
