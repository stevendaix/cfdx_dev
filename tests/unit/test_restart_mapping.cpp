#include "cfdx/core/parallel/restart_mapping.h"

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <vector>

int main() {
    using cfdx::core::parallel::build_restart_permutation;
    using cfdx::core::parallel::remap_restart_values;
    using cfdx::core::parallel::remap_restart_components;

    // Simulate a checkpoint produced with one local ordering and a target
    // partition with a different ordering. Values must follow global IDs,
    // never local cell indices.
    const std::vector<std::uint64_t> source_ids{10, 20, 30, 40};
    const std::vector<std::uint64_t> target_ids{30, 10};
    const auto p = build_restart_permutation(source_ids, target_ids);
    assert((p == std::vector<std::size_t>{2, 0}));

    const std::vector<double> scalar{1.0, 2.0, 3.0, 4.0};
    const auto scalar_target = remap_restart_values(source_ids, scalar, target_ids);
    assert((scalar_target == std::vector<double>{3.0, 1.0}));

    const std::vector<double> vector_values{
        1.0, 10.0, 100.0,
        2.0, 20.0, 200.0,
        3.0, 30.0, 300.0,
        4.0, 40.0, 400.0};
    const auto vector_target = remap_restart_components(
        source_ids, vector_values, 3, target_ids);
    assert((vector_target == std::vector<double>{
        3.0, 30.0, 300.0,
        1.0, 10.0, 100.0}));

    bool rejected = false;
    try {
        (void)build_restart_permutation({10, 10}, {10, 20});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    rejected = false;
    try {
        (void)build_restart_permutation({10, 20}, {10, 30});
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    rejected = false;
    try {
        (void)remap_restart_components(source_ids, scalar, 2, target_ids);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    return 0;
}
