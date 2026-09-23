// M0.11-T05 — Rank-independent restart mapping
#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace cfdx::core::parallel {

// Build the source-position -> target-position permutation from persistent
// global cell identities. IDs are deliberately external to local MPI rank
// ordering, so the mapping remains valid when the same cells are repartitioned
// from N ranks to M ranks.
inline std::vector<std::size_t> build_restart_permutation(
    const std::vector<std::uint64_t>& source_ids,
    const std::vector<std::uint64_t>& target_ids)
{
    if (source_ids.size() != target_ids.size())
        throw std::invalid_argument("restart mapping: source/target cell counts differ");

    std::unordered_map<std::uint64_t, std::size_t> source_index;
    source_index.reserve(source_ids.size());
    for (std::size_t i = 0; i < source_ids.size(); ++i) {
        if (!source_index.emplace(source_ids[i], i).second)
            throw std::invalid_argument("restart mapping: duplicate source cell id");
    }

    std::vector<std::size_t> permutation(target_ids.size());
    for (std::size_t i = 0; i < target_ids.size(); ++i) {
        const auto it = source_index.find(target_ids[i]);
        if (it == source_index.end())
            throw std::invalid_argument("restart mapping: target cell id missing from checkpoint");
        permutation[i] = it->second;
    }
    return permutation;
}

template<class T>
inline std::vector<T> remap_restart_values(
    const std::vector<std::uint64_t>& source_ids,
    const std::vector<T>& source_values,
    const std::vector<std::uint64_t>& target_ids)
{
    if (source_values.size() != source_ids.size())
        throw std::invalid_argument("restart mapping: source value count differs from source ids");
    const auto permutation = build_restart_permutation(source_ids, target_ids);
    std::vector<T> target_values(target_ids.size());
    for (std::size_t i = 0; i < target_values.size(); ++i)
        target_values[i] = source_values[permutation[i]];
    return target_values;
}

template<class T>
inline std::vector<T> remap_restart_components(
    const std::vector<std::uint64_t>& source_ids,
    const std::vector<T>& source_values,
    std::size_t dimension,
    const std::vector<std::uint64_t>& target_ids)
{
    if (dimension == 0)
        throw std::invalid_argument("restart mapping: dimension must be positive");
    if (source_values.size() != source_ids.size() * dimension)
        throw std::invalid_argument("restart mapping: source component count mismatch");

    const auto permutation = build_restart_permutation(source_ids, target_ids);
    std::vector<T> target_values(target_ids.size() * dimension);
    for (std::size_t i = 0; i < target_ids.size(); ++i) {
        for (std::size_t c = 0; c < dimension; ++c)
            target_values[i * dimension + c] =
                source_values[permutation[i] * dimension + c];
    }
    return target_values;
}

} // namespace cfdx::core::parallel
