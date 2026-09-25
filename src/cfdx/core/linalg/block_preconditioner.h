#pragma once

#include "preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace cfdx::core {

// Dense block-Jacobi preconditioner. Each declared diagonal block is inverted
// during setup and applied independently; off-block couplings are ignored.
class BlockDiagonalPreconditioner final : public Preconditioner {
public:
    explicit BlockDiagonalPreconditioner(std::vector<std::vector<std::size_t>> blocks)
        : blocks_(std::move(blocks)) {}

    bool setup(const SparseMatrix& A) override {
        n_ = 0;
        inverse_blocks_.clear();
        if (A.n_rows() != A.n_cols()) return false;
        n_ = A.n_rows();
        if (!validate()) return false;

        inverse_blocks_.reserve(blocks_.size());
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        for (const auto& block : blocks_) {
            const std::size_t m = block.size();
            std::vector<double> aug(m * 2 * m, 0.0);
            auto at = [m, &aug](std::size_t i, std::size_t j) -> double& {
                return aug[i * (2 * m) + j];
            };

            for (std::size_t i = 0; i < m; ++i) {
                const std::size_t global_row = block[i];
                for (std::size_t k = row[global_row]; k < row[global_row + 1]; ++k) {
                    for (std::size_t j = 0; j < m; ++j) {
                        if (col[k] == block[j]) {
                            at(i, j) += val[k];
                            break;
                        }
                    }
                }
                at(i, m + i) = 1.0;
            }

            for (std::size_t pivot = 0; pivot < m; ++pivot) {
                std::size_t best = pivot;
                double best_abs = std::abs(at(pivot, pivot));
                for (std::size_t i = pivot + 1; i < m; ++i) {
                    const double a = std::abs(at(i, pivot));
                    if (a > best_abs) { best_abs = a; best = i; }
                }
                if (best_abs <= 1e-30) return false;

                if (best != pivot) {
                    for (std::size_t j = 0; j < 2 * m; ++j)
                        std::swap(at(pivot, j), at(best, j));
                }

                const double inv_pivot = 1.0 / at(pivot, pivot);
                for (std::size_t j = 0; j < 2 * m; ++j) at(pivot, j) *= inv_pivot;

                for (std::size_t i = 0; i < m; ++i) {
                    if (i == pivot) continue;
                    const double factor = at(i, pivot);
                    if (factor == 0.0) continue;
                    for (std::size_t j = 0; j < 2 * m; ++j)
                        at(i, j) -= factor * at(pivot, j);
                }
            }

            std::vector<double> inv(m * m);
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t j = 0; j < m; ++j)
                    inv[i * m + j] = at(i, m + j);
            inverse_blocks_.push_back(std::move(inv));
        }
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != n_ || z.size() != n_ || inverse_blocks_.size() != blocks_.size()) return false;
        for (std::size_t b = 0; b < blocks_.size(); ++b) {
            const auto& block = blocks_[b];
            const auto& inv = inverse_blocks_[b];
            for (std::size_t i = 0; i < block.size(); ++i) {
                double value = 0.0;
                for (std::size_t j = 0; j < block.size(); ++j)
                    value += inv[i * block.size() + j] * r(block[j]);
                z(block[i]) = value;
            }
        }
        return true;
    }

    const char* name() const override { return "BlockJacobi"; }

private:
    bool validate() const {
        if (n_ == 0 || blocks_.empty()) return false;
        std::vector<bool> covered(n_, false);
        for (const auto& block : blocks_) {
            if (block.empty()) return false;
            for (std::size_t i : block) {
                if (i >= n_ || covered[i]) return false;
                covered[i] = true;
            }
        }
        return std::all_of(covered.begin(), covered.end(), [](bool v) { return v; });
    }

    std::vector<std::vector<std::size_t>> blocks_;
    std::vector<std::vector<double>> inverse_blocks_;
    std::size_t n_{0};
};


enum class SchurFactorization {
    Diagonal,
    Lower,
    Upper,
    Full
};

// Block Schur preconditioner for a two-field system
// A = [A11 A12; A21 A22]. A11 is inverted exactly as a dense block and
// the Schur complement S = A22 - A21*A11^{-1}*A12 is formed explicitly.
// This is a deterministic baseline for coupled pressure/velocity-like
// systems; it is intentionally limited to two disjoint blocks.
class SchurComplementPreconditioner final : public Preconditioner {
public:
    SchurComplementPreconditioner(std::vector<std::size_t> first,
                                   std::vector<std::size_t> second,
                                   SchurFactorization factorization = SchurFactorization::Lower,
                                   double diagonal_schur_scale = -1.0)
        : first_(std::move(first)), second_(std::move(second)),
          factorization_(factorization),
          diagonal_schur_scale_(diagonal_schur_scale) {}

    bool setup(const SparseMatrix& A) override {
        n_ = 0;
        inv_a11_.clear();
        inv_schur_.clear();
        a12_.clear();
        a21_.clear();
        if (A.n_rows() != A.n_cols()) return false;
        n_ = A.n_rows();
        if (!validate()) return false;

        const std::size_t n1 = first_.size();
        const std::size_t n2 = second_.size();
        const auto* row = A.row_offsets_data();
        const auto* col = A.columns_data();
        const auto* val = A.values_data();

        std::vector<double> a11(n1 * n1, 0.0);
        std::vector<double> a12(n1 * n2, 0.0);
        std::vector<double> a21(n2 * n1, 0.0);
        std::vector<double> a22(n2 * n2, 0.0);
        std::vector<std::size_t> pos1(n_, n1);
        std::vector<std::size_t> pos2(n_, n2);
        for (std::size_t i = 0; i < n1; ++i) pos1[first_[i]] = i;
        for (std::size_t i = 0; i < n2; ++i) pos2[second_[i]] = i;

        for (std::size_t i = 0; i < n_; ++i) {
            const bool r1 = pos1[i] < n1;
            const std::size_t ri = r1 ? pos1[i] : pos2[i];
            for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
                const std::size_t j = col[k];
                if (r1) {
                    if (pos1[j] < n1) a11[ri * n1 + pos1[j]] += val[k];
                    else if (pos2[j] < n2) a12[ri * n2 + pos2[j]] += val[k];
                } else {
                    if (pos1[j] < n1) a21[ri * n1 + pos1[j]] += val[k];
                    else if (pos2[j] < n2) a22[ri * n2 + pos2[j]] += val[k];
                }
            }
        }

        if (!invert_dense(a11, n1, inv_a11_)) return false;
        a12_ = a12;
        a21_ = a21;

        // S = A22 - A21 * inv(A11) * A12.
        std::vector<double> tmp(n2 * n1, 0.0);
        for (std::size_t i = 0; i < n2; ++i)
            for (std::size_t k = 0; k < n1; ++k)
                for (std::size_t j = 0; j < n1; ++j)
                    tmp[i * n1 + k] += a21[i * n1 + j] * inv_a11_[j * n1 + k];

        std::vector<double> schur = a22;
        for (std::size_t i = 0; i < n2; ++i)
            for (std::size_t k = 0; k < n2; ++k)
                for (std::size_t j = 0; j < n1; ++j)
                    schur[i * n2 + k] -= tmp[i * n1 + j] * a12[j * n2 + k];

        if (!invert_dense(schur, n2, inv_schur_)) return false;
        return true;
    }

    bool apply(const Vector& r, Vector& z) const override {
        if (r.size() != n_ || z.size() != n_ ||
            inv_a11_.empty() || inv_schur_.empty()) return false;
        const std::size_t n1 = first_.size();
        const std::size_t n2 = second_.size();

        std::vector<double> y1(n1, 0.0);
        std::vector<double> y2(n2, 0.0);
        const auto solve_first = [&](const std::vector<double>& rhs) {
            std::vector<double> solution(n1, 0.0);
            for (std::size_t i = 0; i < n1; ++i)
                for (std::size_t j = 0; j < n1; ++j)
                    solution[i] += inv_a11_[i * n1 + j] * rhs[j];
            return solution;
        };
        const auto solve_second = [&](const std::vector<double>& rhs) {
            std::vector<double> solution(n2, 0.0);
            for (std::size_t i = 0; i < n2; ++i)
                for (std::size_t j = 0; j < n2; ++j)
                    solution[i] += inv_schur_[i * n2 + j] * rhs[j];
            return solution;
        };

        std::vector<double> rhs1(n1, 0.0);
        std::vector<double> rhs2(n2, 0.0);
        for (std::size_t i = 0; i < n1; ++i) rhs1[i] = r(first_[i]);
        for (std::size_t i = 0; i < n2; ++i) rhs2[i] = r(second_[i]);

        if (factorization_ == SchurFactorization::Upper) {
            y2 = solve_second(rhs2);
            for (std::size_t i = 0; i < n1; ++i)
                for (std::size_t j = 0; j < n2; ++j)
                    rhs1[i] -= a12_[i * n2 + j] * y2[j];
            y1 = solve_first(rhs1);
        } else {
            y1 = solve_first(rhs1);
            if (factorization_ != SchurFactorization::Diagonal)
                rhs2 = second_rhs_(y1, r);
            y2 = solve_second(rhs2);
            if (factorization_ == SchurFactorization::Diagonal)
                for (double& value : y2) value *= diagonal_schur_scale_;
            if (factorization_ == SchurFactorization::Full) {
                std::vector<double> upper_coupling(n1, 0.0);
                for (std::size_t i = 0; i < n1; ++i)
                    for (std::size_t j = 0; j < n2; ++j)
                        upper_coupling[i] += a12_[i * n2 + j] * y2[j];
                const auto correction = solve_first(upper_coupling);
                for (std::size_t i = 0; i < n1; ++i) y1[i] -= correction[i];
            }
        }

        for (std::size_t i = 0; i < n2; ++i) {
            z(second_[i]) = y2[i];
        }
        for (std::size_t i = 0; i < n1; ++i) z(first_[i]) = y1[i];
        return true;
    }

    const char* name() const override { return "SchurComplement"; }

private:
    bool validate() const {
        if (n_ == 0 || first_.empty() || second_.empty()) return false;
        std::vector<bool> covered(n_, false);
        for (const auto i : first_) {
            if (i >= n_ || covered[i]) return false;
            covered[i] = true;
        }
        for (const auto i : second_) {
            if (i >= n_ || covered[i]) return false;
            covered[i] = true;
        }
        return std::all_of(covered.begin(), covered.end(), [](bool v) { return v; });
    }

    static bool invert_dense(const std::vector<double>& input, std::size_t n,
                             std::vector<double>& inverse) {
        if (n == 0 || input.size() != n * n) return false;
        std::vector<double> aug(n * 2 * n, 0.0);
        auto at = [n, &aug](std::size_t i, std::size_t j) -> double& {
            return aug[i * (2 * n) + j];
        };
        double scale = 0.0;
        for (const double v : input) scale = std::max(scale, std::abs(v));
        if (!std::isfinite(scale) || scale == 0.0) return false;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) at(i, j) = input[i * n + j];
            at(i, n + i) = 1.0;
        }
        const double pivot_tol = 1e-14 * scale;
        for (std::size_t p = 0; p < n; ++p) {
            std::size_t best = p;
            for (std::size_t i = p + 1; i < n; ++i)
                if (std::abs(at(i, p)) > std::abs(at(best, p))) best = i;
            if (std::abs(at(best, p)) <= pivot_tol) return false;
            if (best != p)
                for (std::size_t j = 0; j < 2 * n; ++j) std::swap(at(p, j), at(best, j));
            const double pivot = at(p, p);
            for (std::size_t j = 0; j < 2 * n; ++j) at(p, j) /= pivot;
            for (std::size_t i = 0; i < n; ++i) {
                if (i == p) continue;
                const double factor = at(i, p);
                if (factor == 0.0) continue;
                for (std::size_t j = 0; j < 2 * n; ++j) at(i, j) -= factor * at(p, j);
            }
        }
        inverse.assign(n * n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) inverse[i * n + j] = at(i, n + j);
        return true;
    }

    std::vector<double> second_rhs_(const std::vector<double>& y1,
                                    const Vector& r) const {
        // A21 is reconstructed from the compact setup data.
        std::vector<double> rhs(second_.size(), 0.0);
        // The coupling is stored as a dense n2*n1 matrix to keep apply()
        // independent of the original SparseMatrix lifetime.
        for (std::size_t i = 0; i < second_.size(); ++i) {
            rhs[i] = r(second_[i]);
            for (std::size_t j = 0; j < first_.size(); ++j)
                rhs[i] -= a21_[i * first_.size() + j] * y1[j];
        }
        return rhs;
    }

    std::vector<std::size_t> first_;
    std::vector<std::size_t> second_;
    std::vector<double> inv_a11_;
    std::vector<double> inv_schur_;
    std::vector<double> a12_;
    std::vector<double> a21_;
    SchurFactorization factorization_;
    double diagonal_schur_scale_;
    std::size_t n_{0};
};

} // namespace cfdx::core
