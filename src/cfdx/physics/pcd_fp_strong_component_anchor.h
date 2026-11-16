#pragma once

// Experimental PCD-only inverse repair for convection-dominated, floating
// pressure-space characteristics. The OUTER coupled FV matrix is never changed.
//
// Identify terminal strongly connected components of the directed graph of
// large Fp row dependencies. If every row of such a component has zero
// normalized row sum, that transport basin has no effective pressure-space
// anchor. An auxiliary zero correction is imposed at one root cell. This
// alters the PCD preconditioner, NOT the physical pressure gauge or the flow.
// The original (unaltered) Fp residual must be reported separately and never
// confused with the auxiliary solve's true residual.

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::physics {

struct PcdFpStrongComponentAnchor {
    cfdx::core::SparseMatrix matrix;
    std::vector<std::size_t> anchored_rows;
    double strong_row_fraction = 1e-2;
    double row_sum_relative_cutoff = 1e-8;

    cfdx::core::Vector rhs(const cfdx::core::Vector& original) const {
        if (original.size() != matrix.n_rows() || !original.is_valid())
            throw std::invalid_argument("PCD Fp anchored RHS dimension/finiteness mismatch");
        auto result = original;
        for (const auto row : anchored_rows) result(row) = 0.0;
        return result;
    }

    void pin_initial(cfdx::core::Vector& x) const {
        if (x.size() != matrix.n_rows() || !x.is_valid())
            throw std::invalid_argument("PCD Fp anchored initial vector invalid");
        for (const auto row : anchored_rows) x(row) = 0.0;
    }
};

inline PcdFpStrongComponentAnchor anchor_pcd_fp_strong_components(
    const cfdx::core::SparseMatrix& A,
    double strong_row_fraction = 1e-2,
    double row_sum_relative_cutoff = 1e-8)
{
    const auto n = A.n_rows();
    if (!n || A.n_cols() != n || !A.is_consistent() ||
        !(strong_row_fraction > 0.0 && strong_row_fraction < 1.0) ||
        !(row_sum_relative_cutoff > 0.0 && row_sum_relative_cutoff < 0.01))
        throw std::invalid_argument("PCD Fp strong-component anchor: invalid matrix/parameters");

    const auto* ro = A.row_offsets_data();
    const auto* co = A.columns_data();
    const auto* va = A.values_data();
    std::vector<std::vector<std::size_t>> graph(n), reverse_graph(n);
    std::vector<double> row_scales(n, 0.0), sums(n, 0.0), diagonals(n, 0.0);
    std::vector<unsigned char> diagonal_found(n, 0);
    for (std::size_t row = 0; row < n; ++row) {
        long double sum = 0.0L;
        for (std::size_t k = ro[row]; k < ro[row + 1]; ++k) {
            const double value = va[k];
            if (!std::isfinite(value))
                throw std::invalid_argument("PCD Fp strong-component anchor: nonfinite matrix");
            row_scales[row] = std::max(row_scales[row], std::abs(value));
            sum += static_cast<long double>(value);
            if (co[k] == row) {
                diagonal_found[row] = 1;
                diagonals[row] = value;
            }
        }
        if (!diagonal_found[row] || !(diagonals[row] > 0.0) ||
            !(row_scales[row] > 0.0) || !std::isfinite(static_cast<double>(sum)))
            throw std::invalid_argument("PCD Fp strong-component anchor: invalid diagonal/row");
        sums[row] = static_cast<double>(sum);
    }
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t k = ro[row]; k < ro[row + 1]; ++k) {
            const auto col = static_cast<std::size_t>(co[k]);
            if (col != row && std::abs(va[k]) >= strong_row_fraction * row_scales[row]) {
                graph[row].push_back(col);
                reverse_graph[col].push_back(row);
            }
        }
    }

    // Iterative Kosaraju: no recursive stack overflow on million-cell meshes.
    std::vector<unsigned char> visited(n, 0);
    std::vector<std::size_t> finishing;
    finishing.reserve(n);
    for (std::size_t start = 0; start < n; ++start) {
        if (visited[start]) continue;
        std::vector<std::pair<std::size_t, std::size_t>> stack{{start, 0}};
        visited[start] = 1;
        while (!stack.empty()) {
            auto& top = stack.back();
            if (top.second < graph[top.first].size()) {
                const auto child = graph[top.first][top.second++];
                if (!visited[child]) {
                    visited[child] = 1;
                    stack.emplace_back(child, 0);
                }
            } else {
                finishing.push_back(top.first);
                stack.pop_back();
            }
        }
    }

    const auto unassigned = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> component(n, unassigned);
    std::size_t component_count = 0;
    for (auto it = finishing.rbegin(); it != finishing.rend(); ++it) {
        if (component[*it] != unassigned) continue;
        std::vector<std::size_t> stack{*it};
        component[*it] = component_count;
        while (!stack.empty()) {
            const auto node = stack.back(); stack.pop_back();
            for (const auto child : reverse_graph[node]) {
                if (component[child] == unassigned) {
                    component[child] = component_count;
                    stack.push_back(child);
                }
            }
        }
        ++component_count;
    }

    std::vector<unsigned char> terminal(component_count, 1);
    std::vector<unsigned char> unanchored(component_count, 1);
    std::vector<std::size_t> first(component_count, unassigned);
    for (std::size_t row = 0; row < n; ++row) {
        const auto cid = component[row];
        first[cid] = std::min(first[cid], row);
        if (std::abs(sums[row]) > row_sum_relative_cutoff * row_scales[row])
            unanchored[cid] = 0;
        for (const auto col : graph[row])
            if (component[col] != cid) terminal[cid] = 0;
    }

    PcdFpStrongComponentAnchor result;
    result.strong_row_fraction = strong_row_fraction;
    result.row_sum_relative_cutoff = row_sum_relative_cutoff;
    for (std::size_t c = 0; c < component_count; ++c)
        if (terminal[c] && unanchored[c])
            result.anchored_rows.push_back(first[c]);
    std::sort(result.anchored_rows.begin(), result.anchored_rows.end());
    if (result.anchored_rows.empty())
        throw std::runtime_error("PCD Fp strong-component anchor requested but no floating terminal basin detected");
    if (result.anchored_rows.size() > std::max<std::size_t>(2, n / 16))
        throw std::runtime_error("PCD Fp strong-component anchor detected too many floating basins");

    // Simultaneous row+column elimination with a diagonal in the native Fp
    // scale. No arbitrary unit pivot or reduction of the physical matrix.
    std::vector<unsigned char> pinned(n, 0);
    for (const auto row : result.anchored_rows) pinned[row] = 1;
    std::vector<cfdx::core::SparseMatrix::Index> row_offsets(ro, ro + n + 1);
    std::vector<cfdx::core::SparseMatrix::Index> cols(co, co + A.nnz());
    std::vector<double> values(va, va + A.nnz());
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t k = ro[row]; k < ro[row + 1]; ++k) {
            if (pinned[row] || pinned[cols[k]])
                values[k] = row == cols[k] ? diagonals[row] : 0.0;
        }
    }
    result.matrix = cfdx::core::SparseMatrix::from_sorted_csr(
        n, n, std::move(row_offsets), std::move(cols), std::move(values));
    return result;
}

}  // namespace cfdx::physics
