#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

// Matrix-free fine-level V-cycle with connectivity-aware pair aggregation.
// The fine operator remains matrix-free during apply(). The SparseMatrix
// passed to setup() is used only for diagonal/graph construction.
// The coarse operator is a diagonal approximation for this foundation stage.
// The mutable workspace makes one instance non-thread-safe across apply calls.
class MatrixFreeVcyclePreconditioner final : public Preconditioner {
public:
    MatrixFreeVcyclePreconditioner(const LinearOperatorBase& op,
                                   double omega = 0.7,
                                   std::size_t pre = 2,
                                   std::size_t post = 2)
        : op_(op), omega_(omega), pre_(pre), post_(post)
    {
        if (op.rows() != op.cols() || !(omega_ > 0.0 && omega_ < 2.0)) {
            throw std::invalid_argument("invalid V-cycle operator");
        }
        setup_diagonal();
    }

    bool setup(const SparseMatrix& A) override
    {
        if (A.n_rows() != op_.rows() || A.n_cols() != op_.cols() ||
            A.n_rows() == 0)
            return false;
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        inv_diag_.assign(A.n_rows(), 0.0);
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            bool found_diagonal = false;
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                if (col[k] >= A.n_cols() || !std::isfinite(val[k]))
                    return false;
                if (col[k] == i) {
                    if (found_diagonal || std::abs(val[k]) <= 1e-30)
                        return false;
                    inv_diag_[i] = 1.0 / val[k];
                    if (!std::isfinite(inv_diag_[i]))
                        return false;
                    found_diagonal = true;
                }
            }
            if (!found_diagonal)
                return false;
        }
        build_connectivity_aware_aggregation(A);
        resize_workspace();
        return true;
    }

    std::size_t coarse_size() const noexcept { return coarse_to_fine_.size(); }

    std::size_t aggregate_of(std::size_t fine_cell) const
    {
        if (fine_cell >= aggregate_of_.size()) {
            throw std::out_of_range("aggregate_of: fine-cell index out of range");
        }
        return aggregate_of_[fine_cell];
    }

    bool apply(const Vector& r, Vector& z) const override
    {
        if (r.size() != op_.rows() || inv_diag_.size() != r.size() ||
            aggregate_of_.size() != r.size() || coarse_to_fine_.empty())
            return false;
        if (!r.is_valid())
            return false;
        if (z.size() != r.size()) z.resize(r.size());
        z.fill(0.0);
        smooth(r, z, pre_);
        if (!z.is_valid())
            return false;

        op_.apply(z, ws_.fine_A);
        for (std::size_t i = 0; i < r.size(); ++i)
            ws_.fine_r(i) = r(i) - ws_.fine_A(i);

        ws_.coarse_r.fill(0.0);
        for (std::size_t i = 0; i < r.size(); ++i)
            ws_.coarse_r(aggregate_of_[i]) += ws_.fine_r(i);

        for (std::size_t c = 0; c < ws_.coarse_x.size(); ++c)
            ws_.coarse_x(c) = ws_.coarse_r(c) * ws_.coarse_inv_diag[c];
            if (!std::isfinite(ws_.coarse_x(c)))
                return false;

        for (std::size_t i = 0; i < r.size(); ++i) {
            z(i) += ws_.coarse_x(aggregate_of_[i]);
            if (!std::isfinite(z(i)))
                return false;
        }

        smooth(r, z, post_);
        return true;
    }

    const char* name() const override { return "matrix-free-agglomerated-vcycle"; }

private:
    struct Workspace {
        Vector fine_r, fine_A, coarse_r, coarse_x;
        std::vector<double> coarse_inv_diag;
    };

    void setup_diagonal()
    {
        Vector d(op_.rows(), 1.0);
        if (op_.has_diagonal()) op_.diagonal(d);
        inv_diag_.resize(d.size());
        for (std::size_t i = 0; i < d.size(); ++i)
            inv_diag_[i] = std::abs(d(i)) > 1e-30 ? 1.0 / d(i) : 1.0;

        aggregate_of_.resize(d.size());
        coarse_to_fine_.resize((d.size() + 1) / 2);
        for (std::size_t i = 0; i < d.size(); ++i) aggregate_of_[i] = i / 2;
        resize_workspace();
    }

    void build_connectivity_aware_aggregation(const SparseMatrix& A)
    {
        const std::size_t n = A.n_rows();
        aggregate_of_.assign(n, 0);
        coarse_to_fine_.clear();
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();
        std::vector<bool> matched(n, false);

        for (std::size_t i = 0; i < n; ++i) {
            if (matched[i]) continue;
            std::size_t best = n;
            double best_strength = -1.0;
            const double di = std::abs(inv_diag_[i]);

            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                if (j == i || j >= n || matched[j]) continue;
                const double dj = std::abs(inv_diag_[j]);
                const double denom = std::sqrt(std::max(
                    std::abs(1.0 / inv_diag_[i]) *
                    std::abs(1.0 / inv_diag_[j]), 1e-60));
                const double strength = std::abs(val[k]) / denom;
                if (strength > best_strength) {
                    best_strength = strength;
                    best = j;
                }
            }

            const std::size_t aggregate = coarse_to_fine_.size();
            coarse_to_fine_.push_back(i);
            aggregate_of_[i] = aggregate;
            matched[i] = true;
            if (best < n) {
                aggregate_of_[best] = aggregate;
                matched[best] = true;
            }
        }
    }

    void resize_workspace()
    {
        const std::size_t nc = coarse_to_fine_.size();
        ws_.fine_r.resize(inv_diag_.size());
        ws_.fine_A.resize(inv_diag_.size());
        ws_.coarse_r.resize(nc);
        ws_.coarse_x.resize(nc);
        ws_.coarse_inv_diag.assign(nc, 0.0);

        for (std::size_t i = 0; i < inv_diag_.size(); ++i) {
            const std::size_t c = aggregate_of_[i];
            const double diagonal =
                1.0 / std::max(std::abs(inv_diag_[i]), 1e-30);
            ws_.coarse_inv_diag[c] += diagonal;
        }
        for (double& diagonal : ws_.coarse_inv_diag)
            diagonal = 1.0 / std::max(diagonal, 1e-30);

        ws_.fine_r.fill(0.0);
        ws_.fine_A.fill(0.0);
        ws_.coarse_r.fill(0.0);
        ws_.coarse_x.fill(0.0);
    }

    void smooth(const Vector& r, Vector& x, std::size_t sweeps) const
    {
        for (std::size_t s = 0; s < sweeps; ++s) {
            op_.apply(x, ws_.fine_A);
            for (std::size_t i = 0; i < r.size(); ++i)
                x(i) += omega_ * inv_diag_[i] * (r(i) - ws_.fine_A(i));
        }
    }

    const LinearOperatorBase& op_;
    double omega_;
    std::size_t pre_, post_;
    std::vector<double> inv_diag_;
    std::vector<std::size_t> aggregate_of_, coarse_to_fine_;
    mutable Workspace ws_;
};

using AgglomeratedAMGPreconditioner = MatrixFreeVcyclePreconditioner;
} // namespace cfdx::core
