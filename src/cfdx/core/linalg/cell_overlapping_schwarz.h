#pragma once

// Small-domain, serial reference implementation of cell-aligned overlapping
// additive Schwarz for the CFDX four-field coupled ordering [Ux,Uy,Uz,p].
// Experimental only: no automatic solver dispatch, no MPI communication, no
// coarse level, no claim of scalable or Reynolds-independent convergence.
// Motivated by PETSc PCASM / restricted additive Schwarz (Cai & Sarkis, 1999).

#include "cfdx/core/linalg/preconditioner.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

namespace cfdx::core {

enum class CellSchwarzRestriction {
    Restricted,     // RAS: full residual on overlap; correction only on core.
    WeightedAdditive // Full correction on overlap, normalized by multiplicity.
};

struct CellSchwarzConfig {
    std::size_t cells_per_subdomain{1};
    std::size_t overlap_hops{1};
    std::size_t max_local_unknowns{64};
    CellSchwarzRestriction restriction{CellSchwarzRestriction::Restricted};
    double relative_pivot_threshold{1.0e-12};
};

class CellOverlappingSchwarzPreconditioner final : public Preconditioner {
public:
    CellOverlappingSchwarzPreconditioner(std::size_t n_cells, CellSchwarzConfig config)
        : n_cells_(n_cells), config_(config) {}

    const char* name() const override { return "CellOverlappingSchwarzExperimental"; }
    std::size_t subdomain_count() const noexcept { return domains_.size(); }
    std::size_t largest_subdomain() const noexcept { return largest_; }

    bool setup(const SparseMatrix& A) override {
        // Transactional: a failed setup must not invalidate a working solver.
        if (n_cells_ == 0 || n_cells_ > std::numeric_limits<std::size_t>::max() / 4 ||
            config_.cells_per_subdomain == 0 || config_.max_local_unknowns < 4 ||
            config_.max_local_unknowns % 4 != 0 || config_.max_local_unknowns > 256 ||
            config_.overlap_hops > 8 ||
            (config_.restriction != CellSchwarzRestriction::Restricted &&
             config_.restriction != CellSchwarzRestriction::WeightedAdditive) ||
            !std::isfinite(config_.relative_pivot_threshold) ||
            !(config_.relative_pivot_threshold > 0.0) ||
            !(config_.relative_pivot_threshold < 1.0) ||
            A.n_rows() != 4 * n_cells_ || A.n_cols() != A.n_rows()) return false;
        const auto* offsets = A.row_offsets_data();
        const auto* columns = A.columns_data();
        const auto* values = A.values_data();
        const std::size_t n = A.n_rows();
        if (offsets[0] != 0 || offsets[n] != A.nnz()) return false;

        std::vector<std::vector<std::size_t>> neighbours(n_cells_);
        for (std::size_t i = 0; i < n; ++i) {
            if (offsets[i] > offsets[i + 1] || offsets[i + 1] > A.nnz()) return false;
            const std::size_t cell_i = i % n_cells_;
            for (std::size_t k = offsets[i]; k < offsets[i + 1]; ++k) {
                const auto j = static_cast<std::size_t>(columns[k]);
                if (j >= n || !std::isfinite(values[k])) return false;
                const std::size_t cell_j = j % n_cells_;
                if (cell_i != cell_j) {
                    neighbours[cell_i].push_back(cell_j);
                    neighbours[cell_j].push_back(cell_i);
                }
            }
        }
        for (auto& adj : neighbours) {
            std::sort(adj.begin(), adj.end());
            adj.erase(std::unique(adj.begin(), adj.end()), adj.end());
        }

        std::vector<Domain> prepared;
        std::vector<std::size_t> coverage(n, 0);
        std::size_t new_largest = 0;
        const std::size_t maximum_cells = config_.max_local_unknowns / 4;
        for (std::size_t begin = 0; begin < n_cells_;) {
            const std::size_t end = begin + std::min(config_.cells_per_subdomain, n_cells_ - begin);
            std::vector<std::size_t> cells;
            for (std::size_t c = begin; c < end; ++c) cells.push_back(c);
            if (cells.size() > maximum_cells) return false;
            std::vector<std::size_t> frontier = cells;
            // Expand on the undirected *cell coupling* graph, not CSR index order.
            for (std::size_t hop = 0; hop < config_.overlap_hops; ++hop) {
                std::vector<std::size_t> next = cells;
                for (const auto c : frontier)
                    next.insert(next.end(), neighbours[c].begin(), neighbours[c].end());
                std::sort(next.begin(), next.end());
                next.erase(std::unique(next.begin(), next.end()), next.end());
                if (next.size() > maximum_cells) return false;
                std::vector<std::size_t> added;
                std::set_difference(next.begin(), next.end(), cells.begin(), cells.end(),
                                    std::back_inserter(added));
                cells = std::move(next);
                frontier = std::move(added);
                if (frontier.empty()) break;
            }
            Domain domain;
            domain.global.reserve(4 * cells.size());
            domain.core.reserve(4 * (end - begin));
            for (std::size_t field = 0; field < 4; ++field) {
                for (auto c : cells) domain.global.push_back(field * n_cells_ + c);
                for (std::size_t c = begin; c < end; ++c)
                    domain.core.push_back(field * n_cells_ + c);
            }
            const std::size_t m = domain.global.size();
            new_largest = std::max(new_largest, m);
            domain.lu.assign(m * m, 0.0);
            domain.row_scale.assign(m, 0.0);
            domain.pivots.resize(m);
            for (std::size_t local_i = 0; local_i < m; ++local_i) {
                const std::size_t global_i = domain.global[local_i];
                for (std::size_t k = offsets[global_i]; k < offsets[global_i + 1]; ++k) {
                    const auto it = std::lower_bound(domain.global.begin(), domain.global.end(),
                                                     static_cast<std::size_t>(columns[k]));
                    if (it != domain.global.end() && *it == columns[k])
                        domain.lu[local_i * m + static_cast<std::size_t>(it - domain.global.begin())] += values[k];
                }
            }
            for (std::size_t i = 0; i < m; ++i) {
                double row_norm = 0.0;
                for (std::size_t j = 0; j < m; ++j)
                    row_norm += std::abs(domain.lu[i * m + j]);
                if (!(row_norm > 0.0) || !std::isfinite(row_norm)) return false;
                domain.row_scale[i] = row_norm;
                for (std::size_t j = 0; j < m; ++j)
                    domain.lu[i * m + j] /= row_norm;
            }
            // Row-scaled partial-pivot LU. The exact same row operations are
            // applied to RHS at apply() via the stored pivot indices.
            for (std::size_t k = 0; k < m; ++k) {
                std::size_t pivot = k;
                double best = std::abs(domain.lu[k * m + k]);
                for (std::size_t i = k + 1; i < m; ++i) {
                    const double candidate = std::abs(domain.lu[i * m + k]);
                    if (candidate > best) { pivot = i; best = candidate; }
                }
                if (!std::isfinite(best) || best <= config_.relative_pivot_threshold)
                    return false;
                domain.pivots[k] = pivot;
                if (pivot != k) {
                    for (std::size_t j = 0; j < m; ++j)
                        std::swap(domain.lu[k * m + j], domain.lu[pivot * m + j]);
                }
                for (std::size_t i = k + 1; i < m; ++i) {
                    domain.lu[i * m + k] /= domain.lu[k * m + k];
                    for (std::size_t j = k + 1; j < m; ++j)
                        domain.lu[i * m + j] -= domain.lu[i * m + k] * domain.lu[k * m + j];
                }
            }
            for (const auto global_i : domain.global) ++coverage[global_i];
            prepared.push_back(std::move(domain));
            begin = end;
        }
        if (std::any_of(coverage.begin(), coverage.end(),
                        [](std::size_t count) { return count == 0; })) return false;
        std::vector<SparseMatrix::Index> saved_columns(columns, columns + A.nnz());
        std::vector<SparseMatrix::Index> saved_offsets(offsets, offsets + n + 1);
        domains_.swap(prepared);
        coverage_.swap(coverage);
        csr_columns_.swap(saved_columns);
        csr_offsets_.swap(saved_offsets);
        largest_ = new_largest;
        return true;
    }

    bool update_values(const SparseMatrix& A) override {
        if (domains_.empty() || A.n_rows() != 4 * n_cells_ ||
            A.n_cols() != A.n_rows() || A.nnz() != csr_columns_.size()) return false;
        const auto* rows = A.row_offsets_data();
        const auto* cols = A.columns_data();
        if (!std::equal(csr_offsets_.begin(), csr_offsets_.end(), rows) ||
            !std::equal(csr_columns_.begin(), csr_columns_.end(), cols)) return false;
        return setup(A); // transactional full numeric setup; no stale factors.
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (domains_.empty() || r.size() != 4 * n_cells_ ||
            z.size() != r.size() || !r.is_valid()) return false;
        Vector proposal(r.size(), 0.0);
        for (const auto& d : domains_) {
            const std::size_t m = d.global.size();
            std::vector<double> x(m);
            for (std::size_t i = 0; i < m; ++i)
                x[i] = r(d.global[i]) / d.row_scale[i];
            for (std::size_t k = 0; k < m; ++k)
                std::swap(x[k], x[d.pivots[k]]);
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < i; ++j)
                    x[i] -= d.lu[i * m + j] * x[j];
            for (std::size_t i = m; i-- > 0;) {
                for (std::size_t j = i + 1; j < m; ++j)
                    x[i] -= d.lu[i * m + j] * x[j];
                x[i] /= d.lu[i * m + i];
                if (!std::isfinite(x[i])) return false;
            }
            if (config_.restriction == CellSchwarzRestriction::Restricted) {
                for (const auto gi : d.core) {
                    const auto it = std::lower_bound(d.global.begin(), d.global.end(), gi);
                    proposal(gi) = x[static_cast<std::size_t>(it - d.global.begin())];
                }
            } else {
                for (std::size_t i = 0; i < m; ++i)
                    proposal(d.global[i]) += x[i] / static_cast<double>(coverage_[d.global[i]]);
            }
        }
        if (!proposal.is_valid()) return false;
        z = std::move(proposal);
        return true;
    }

private:
    struct Domain {
        std::vector<std::size_t> global;
        std::vector<std::size_t> core;
        std::vector<double> lu;
        std::vector<double> row_scale;
        std::vector<std::size_t> pivots;
    };
    std::size_t n_cells_;
    CellSchwarzConfig config_;
    std::vector<Domain> domains_;
    std::vector<std::size_t> coverage_;
    std::vector<SparseMatrix::Index> csr_columns_;
    std::vector<SparseMatrix::Index> csr_offsets_;
    std::size_t largest_{0};
};

} // namespace cfdx::core
