#pragma once

#include "cfdx/core/linalg/preconditioner.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace cfdx::core {

// CFDX-native algebraic line-block Jacobi, designed for stretched FV meshes.
// This is an independent implementation of an anisotropy-oriented smoother,
// NOT vendored HYPRE/PETSc code. No geometry or external library is required.
//
// Opposite-sign, nearly symmetric strong CSR edges are assembled into disjoint
// paths (degree <= 2, no cycles). Each tridiagonal principal block is factored
// once, then applied in O(n). Remaining edges are left for the outer Krylov/AMG.
// It is a fixed linear preconditioner while the factorization is unchanged.
class AlgebraicLinePreconditioner final : public Preconditioner {
public:
    explicit AlgebraicLinePreconditioner(double direction_ratio = 32.0,
                                         double strength_fraction = 0.1,
                                         bool multiplicative = false)
        : direction_ratio_(direction_ratio), strength_fraction_(strength_fraction),
          multiplicative_(multiplicative) {}

    bool setup(const SparseMatrix& A) override {
        // Construct and factor the replacement before touching the active
        // preconditioner. A rejected setup must not discard a valid solver.
        if (!valid_matrix(A)) return false;
        if (!(direction_ratio_ > 1.0) || !(strength_fraction_ > 0.0 &&
             strength_fraction_ <= 1.0)) return false;

        const std::size_t n = A.n_rows();
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* value = A.values_data();
        std::vector<double> max_negative(n, 0.0);
        std::vector<double> min_negative(n, std::numeric_limits<double>::infinity());
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] != i && value[k] < 0.0) {
                    const double v = -value[k];
                    max_negative[i] = std::max(max_negative[i], v);
                    min_negative[i] = std::min(min_negative[i], v);
                }
            }
        }

        struct Edge { std::size_t i, j; double magnitude; };
        std::vector<Edge> edges;
        for (std::size_t i = 0; i < n; ++i) {
            if (!(max_negative[i] >= direction_ratio_ * min_negative[i])) continue;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                if (j <= i || value[k] >= 0.0 ||
                    -value[k] < strength_fraction_ * max_negative[i] ||
                    !(max_negative[j] >= direction_ratio_ * min_negative[j]) ||
                    -value[k] < strength_fraction_ * max_negative[j]) continue;
                const double transposed = entry(A, j, i);
                const double scale = std::max(-value[k], std::abs(transposed));
                if (transposed >= 0.0 ||
                    std::abs(transposed - value[k]) > 1e-10 * scale) continue;
                edges.push_back({i, j, -value[k]});
            }
        }
        std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
            if (a.magnitude != b.magnitude) return a.magnitude > b.magnitude;
            if (a.i != b.i) return a.i < b.i;
            return a.j < b.j;
        });

        // Union-find enforces acyclic paths; vertex degrees enforce
        // tridiagonal blocks even on unstructured CSR graphs.
        std::vector<std::size_t> parent(n);
        std::iota(parent.begin(), parent.end(), 0);
        const auto root = [&parent](std::size_t a) {
            while (parent[a] != a) {
                parent[a] = parent[parent[a]];
                a = parent[a];
            }
            return a;
        };
        std::vector<std::vector<std::size_t>> neighbours(n);
        for (const auto& edge : edges) {
            if (neighbours[edge.i].size() == 2 ||
                neighbours[edge.j].size() == 2) continue;
            const std::size_t ri = root(edge.i);
            const std::size_t rj = root(edge.j);
            if (ri == rj) continue;
            parent[rj] = ri;
            neighbours[edge.i].push_back(edge.j);
            neighbours[edge.j].push_back(edge.i);
        }
        for (auto& list : neighbours) std::sort(list.begin(), list.end());

        std::vector<bool> seen(n, false);
        std::vector<Line> candidate;
        candidate.reserve(n);
        std::size_t longest = 0;
        std::size_t covered = 0;
        for (std::size_t start = 0; start < n; ++start) {
            if (seen[start] || neighbours[start].size() > 1) continue;
            Line line;
            std::size_t previous = n;
            std::size_t current = start;
            while (current != n && !seen[current]) {
                seen[current] = true;
                line.indices.push_back(current);
                std::size_t next = n;
                for (const auto adjacent : neighbours[current])
                    if (adjacent != previous && !seen[adjacent]) { next = adjacent; break; }
                previous = current;
                current = next;
            }
            longest = std::max(longest, line.indices.size());
            if (line.indices.size() > 1) covered += line.indices.size();
            candidate.push_back(std::move(line));
        }
        if (candidate.empty()) return false;
        if (!factor_lines(A, candidate)) return false;
        // An algebraic path may capture only one of two equally strong
        // directions (e.g. four strong diagonal neighbours in a rotated
        // diffusion stencil). Line-GS is unsuitable unless the selected
        // tridiagonal blocks contain most of the negative coupling mass.
        long double total_negative = 0.0L, on_line_negative = 0.0L;
        if (multiplicative_) {
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t k = row[i]; k < row[i+1]; ++k)
                    if (col[k] != i && value[k] < 0.0)
                        total_negative += -static_cast<long double>(value[k]);
            for (const auto& line : candidate)
                for (std::size_t k = 1; k < line.indices.size(); ++k) {
                    const auto i = line.indices[k-1], j = line.indices[k];
                    on_line_negative -= static_cast<long double>(entry(A, i, j));
                    on_line_negative -= static_cast<long double>(entry(A, j, i));
                }
        }
        const double coupling_share = total_negative > 0.0L
            ? static_cast<double>(on_line_negative / total_negative) : 0.0;
        // A long line can contain a weaker connection when the strength
        // graph branches. Record each row's *original* best on-line/overall
        // strength ratio to certify later numeric-only matrix updates.
        std::vector<double> dominance(n, 0.0);
        for (const auto& line : candidate) {
            for (std::size_t k = 0; k < line.indices.size(); ++k) {
                const auto i = line.indices[k];
                double best = 0.0;
                if (k) best = std::max(best, -entry(A, i, line.indices[k - 1]));
                if (k + 1 < line.indices.size())
                    best = std::max(best, -entry(A, i, line.indices[k + 1]));
                if (max_negative[i] > 0.0)
                    dominance[i] = best / max_negative[i];
            }
        }
        // Store only the off-block CSR coefficients in compressed form.
        // A multiplicative sweep then visits off-line entries without a
        // full SpMV or per-sweep graph search. Adjacent line unknowns are
        // handled exactly by the factored tridiagonal principal blocks.
        std::vector<std::size_t> off_offsets;
        std::vector<SparseMatrix::Index> off_columns;
        std::vector<std::size_t> off_positions;
        std::vector<double> off_values;
        std::vector<std::size_t> row_line_ids;
        if (multiplicative_)
            build_off_block(A, candidate, off_offsets, off_columns,
                            off_positions, off_values, row_line_ids);
        row_offsets_.assign(row, row + n + 1);
        columns_.assign(col, col + A.nnz());
        off_offsets_ = std::move(off_offsets);
        off_columns_ = std::move(off_columns);
        off_positions_ = std::move(off_positions);
        off_values_ = std::move(off_values);
        row_line_ids_ = std::move(row_line_ids);
        dominant_line_share_ = std::move(dominance);
        coupling_share_ = coupling_share;
        lines_ = std::move(candidate);
        longest_line_ = longest;
        covered_rows_ = covered;
        return true;
    }

    bool update_values(const SparseMatrix& A) override {
        if (lines_.empty() || !valid_matrix(A) ||
            A.n_rows() != row_offsets_.size() - 1 ||
            A.nnz() != columns_.size() ||
            !std::equal(row_offsets_.begin(), row_offsets_.end(), A.row_offsets_data()) ||
            !std::equal(columns_.begin(), columns_.end(), A.columns_data())) return false;
        // Numeric refresh is valid only while the directional strength
        // pattern remains compatible with the existing line graph. Otherwise
        // an explicit full setup() is required; silently reusing stale lines
        // can cause the anisotropic smoother to stall during continuation.
        const std::size_t n = A.n_rows();
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        std::vector<double> strongest(n, 0.0);
        std::vector<double> weakest(n, std::numeric_limits<double>::infinity());
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                if (col[k] != i && val[k] < 0.0) {
                    strongest[i] = std::max(strongest[i], -val[k]);
                    weakest[i] = std::min(weakest[i], -val[k]);
                }
        if (covered_rows_ == 0) {
            for (std::size_t i = 0; i < n; ++i)
                if (strongest[i] >= direction_ratio_ * weakest[i]) return false;
        } else {
            for (const auto& line : lines_) {
                if (line.indices.size() == 1) continue;
                for (std::size_t k = 0; k < line.indices.size(); ++k) {
                    const auto i = line.indices[k];
                    if (!(strongest[i] >= direction_ratio_ * weakest[i]))
                        return false;
                    // min/max anisotropy can remain huge because a tiny
                    // third-direction edge persists even as the dominant
                    // physical coupling rotates. Check the actual old
                    // line's dominance instead of that ratio alone.
                    double best_line = 0.0;
                    if (k) best_line = std::max(best_line,
                        -entry(A, i, line.indices[k - 1]));
                    if (k + 1 < line.indices.size())
                        best_line = std::max(best_line,
                            -entry(A, i, line.indices[k + 1]));
                    constexpr double allowed_dominance_drift = 0.90;
                    if (!(best_line >= allowed_dominance_drift *
                                       dominant_line_share_[i] * strongest[i]))
                        return false;
                }
                for (std::size_t k = 1; k < line.indices.size(); ++k) {
                    const auto i = line.indices[k - 1];
                    const auto j = line.indices[k];
                    if (-entry(A, i, j) < strength_fraction_ * strongest[i] ||
                        -entry(A, j, i) < strength_fraction_ * strongest[j])
                        return false;
                }
            }
        }
        double refreshed_coupling_share = coupling_share_;
        if (multiplicative_) {
            long double total_negative = 0.0L, on_line_negative = 0.0L;
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t p = row[i]; p < row[i+1]; ++p)
                    if (col[p] != i && val[p] < 0.0)
                        total_negative -= static_cast<long double>(val[p]);
            for (const auto& line : lines_)
                for (std::size_t k = 1; k < line.indices.size(); ++k) {
                    const auto i = line.indices[k-1], j = line.indices[k];
                    on_line_negative -= static_cast<long double>(entry(A, i, j));
                    on_line_negative -= static_cast<long double>(entry(A, j, i));
                }
            const double share = total_negative > 0.0L
                ? static_cast<double>(on_line_negative / total_negative) : 0.0;
            if (!(share >= 0.85)) return false;
            refreshed_coupling_share = share;
        }
        auto refreshed = lines_;
        if (!factor_lines(A, refreshed)) return false;
        // Transactionally refresh coefficients of all non-tridiagonal
        // couplings as well as the line factors. No stale off-block edges.
        std::vector<double> off_values(off_positions_.size());
        for (std::size_t k = 0; k < off_positions_.size(); ++k)
            off_values[k] = val[off_positions_[k]];
        lines_ = std::move(refreshed);
        off_values_.swap(off_values);
        coupling_share_ = refreshed_coupling_share;
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (lines_.empty() || r.size() != row_offsets_.size() - 1 || !r.is_valid())
            return false;
        if (z.size() != r.size()) z.resize(r.size());
        // The lines are disjoint. Store the forward Thomas substitution in
        // the output vector itself, then back-substitute in-place. This
        // removes the per-apply heap allocation from every AMG V-cycle,
        // retains O(n) complexity and also supports apply(r, r) aliasing.
        for (const auto& line : lines_) {
            const std::size_t m = line.indices.size();
            z(line.indices[0]) = r(line.indices[0]);
            for (std::size_t k = 1; k < m; ++k)
                z(line.indices[k]) = r(line.indices[k]) -
                    line.lower_factor[k] * z(line.indices[k - 1]);
            z(line.indices[m - 1]) *= line.inverse_pivot[m - 1];
            for (std::size_t k = m - 1; k-- > 0;)
                z(line.indices[k]) =
                    (z(line.indices[k]) - line.upper[k] * z(line.indices[k + 1])) *
                    line.inverse_pivot[k];
        }
        return z.is_valid();
    }

    // A true multiplicative block-line Gauss-Seidel sweep. Unlike block
    // Jacobi, each line consumes updates from previously solved lines.
    // The off-block CSR is cached at setup and refreshed transactionally;
    // no SpMV or heap allocation is required in each sweep if workspace is
    // pre-sized. The operation is generally nonsymmetric: pair with FGMRES.
    bool sweep(const Vector& rhs, Vector& x, bool backward,
               double relaxation, std::vector<double>& workspace) const {
        if (!multiplicative_ || lines_.empty() ||
            off_offsets_.size() != rhs.size() + 1 ||
            rhs.size() != row_line_ids_.size() ||
            x.size() != rhs.size() || &rhs == &x || !rhs.is_valid() ||
            !x.is_valid() || !(relaxation > 0.0 && relaxation <= 1.0) ||
            !std::isfinite(relaxation)) return false;
        if (workspace.size() < longest_line_) workspace.resize(longest_line_);
        for (std::size_t l = 0; l < lines_.size(); ++l) {
            const auto& line = lines_[backward ? lines_.size() - l - 1 : l];
            const auto m = line.indices.size();
            for (std::size_t k = 0; k < m; ++k) {
                const auto i = line.indices[k];
                double value = rhs(i);
                for (std::size_t p = off_offsets_[i]; p < off_offsets_[i+1]; ++p)
                    value -= off_values_[p] * x(off_columns_[p]);
                if (k) value -= line.lower_factor[k] * workspace[k-1];
                if (!std::isfinite(value)) return false;
                workspace[k] = value;
            }
            workspace[m-1] *= line.inverse_pivot[m-1];
            for (std::size_t k = m - 1; k-- > 0;)
                workspace[k] = (workspace[k] - line.upper[k] * workspace[k+1]) *
                               line.inverse_pivot[k];
            for (std::size_t k = 0; k < m; ++k) {
                const auto i = line.indices[k];
                x(i) += relaxation * (workspace[k] - x(i));
                if (!std::isfinite(x(i))) return false;
            }
        }
        return true;
    }

    const std::vector<std::size_t>& row_line_ids() const noexcept {
        return row_line_ids_;
    }
    double coupling_share() const noexcept { return coupling_share_; }
    const char* name() const override { return multiplicative_ ? "CFDXAlgebraicLineGaussSeidel" : "CFDXAlgebraicLineBlockJacobi"; }
    std::size_t line_count() const noexcept { return lines_.size(); }
    std::size_t longest_line() const noexcept { return longest_line_; }
    std::size_t covered_rows() const noexcept { return covered_rows_; }

private:
    struct Line {
        std::vector<std::size_t> indices;
        std::vector<double> lower_factor;
        std::vector<double> upper;
        std::vector<double> inverse_pivot;
    };

    static void build_off_block(
        const SparseMatrix& A, const std::vector<Line>& lines,
        std::vector<std::size_t>& offsets,
        std::vector<SparseMatrix::Index>& columns,
        std::vector<std::size_t>& positions,
        std::vector<double>& values,
        std::vector<std::size_t>& row_line_ids) {
        const auto n = A.n_rows();
        const auto* ro = A.row_offsets_data();
        const auto* ci = A.columns_data();
        const auto* av = A.values_data();
        std::vector<std::size_t> prev(n, n), next(n, n);
        row_line_ids.assign(n, n);
        for (std::size_t l = 0; l < lines.size(); ++l) {
            const auto& nodes = lines[l].indices;
            for (std::size_t k = 0; k < nodes.size(); ++k) {
                const auto i = nodes[k];
                row_line_ids[i] = l;
                if (k) prev[i] = nodes[k - 1];
                if (k + 1 < nodes.size()) next[i] = nodes[k + 1];
            }
        }
        offsets.assign(n + 1, 0);
        columns.clear(); positions.clear(); values.clear();
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t p = ro[i]; p < ro[i+1]; ++p) {
                if (ci[p] == i || ci[p] == prev[i] || ci[p] == next[i])
                    continue;
                columns.push_back(ci[p]);
                positions.push_back(p);
                values.push_back(av[p]);
            }
            offsets[i+1] = columns.size();
        }
    }

    static double entry(const SparseMatrix& A, std::size_t i, std::size_t j) {
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        const auto first = col + row[i];
        const auto last = col + row[i + 1];
        auto it = std::lower_bound(first, last, static_cast<SparseMatrix::Index>(j));
        double sum = 0.0;
        while (it != last && *it == j) {
            sum += val[static_cast<std::size_t>(it - col)];
            ++it;
        }
        return sum;
    }

    static bool valid_matrix(const SparseMatrix& A) {
        if (A.n_rows() == 0 || A.n_rows() != A.n_cols() || !A.is_consistent())
            return false;
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        for (std::size_t i = 0; i < A.n_rows(); ++i)
            for (std::size_t k = row[i]; k < row[i + 1]; ++k)
                if (!std::isfinite(val[k]) || col[k] >= A.n_cols()) return false;
        return true;
    }

    static bool factor_lines(const SparseMatrix& A, std::vector<Line>& lines) {
        constexpr double pivot_guard = 32.0 * std::numeric_limits<double>::epsilon();
        for (auto& line : lines) {
            const std::size_t m = line.indices.size();
            line.upper.assign(m, 0.0);
            line.lower_factor.assign(m, 0.0);
            line.inverse_pivot.assign(m, 0.0);
            for (std::size_t k = 0; k < m; ++k) {
                const std::size_t i = line.indices[k];
                const double diagonal = entry(A, i, i);
                const double lower = k ? entry(A, i, line.indices[k - 1]) : 0.0;
                if (k + 1 < m) line.upper[k] = entry(A, i, line.indices[k + 1]);
                const double scale = std::max({std::abs(diagonal), std::abs(lower),
                                              std::abs(line.upper[k])});
                const double pivot = diagonal -
                    (k ? lower * line.inverse_pivot[k - 1] * line.upper[k - 1] : 0.0);
                if (!(pivot > pivot_guard * scale) || !std::isfinite(pivot))
                    return false;
                line.inverse_pivot[k] = 1.0 / pivot;
                if (!std::isfinite(line.inverse_pivot[k])) return false;
                if (k) line.lower_factor[k] = lower * line.inverse_pivot[k - 1];
            }
        }
        return true;
    }

    double direction_ratio_;
    double strength_fraction_;
    bool multiplicative_ = false;
    double coupling_share_ = 0.0;
    std::vector<Line> lines_;
    std::vector<SparseMatrix::Index> row_offsets_;
    std::vector<SparseMatrix::Index> columns_;
    std::vector<double> dominant_line_share_;
    std::vector<std::size_t> off_offsets_, off_positions_, row_line_ids_;
    std::vector<SparseMatrix::Index> off_columns_;
    std::vector<double> off_values_;
    std::size_t longest_line_ = 0;
    std::size_t covered_rows_ = 0;
};

} // namespace cfdx::core