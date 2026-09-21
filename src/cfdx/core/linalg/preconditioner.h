#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <unordered_map>

namespace cfdx {
namespace core {

/** Jacobi, Gauss-Seidel and ILU(0) preconditioners. */
class JacobiPreconditioner {
public:
    bool setup(const SparseMatrix& A) {
        if (A.n_rows() != A.n_cols()) return false;
        diag_.assign(A.n_rows(), 0.0);
        for (std::size_t i = 0; i < A.n_rows(); ++i) {
            diag_[i] = A(i, i);
            if (diag_[i] == 0.0 || !std::isfinite(diag_[i])) return false;
        }
        return true;
    }
    bool apply(const Vector& r, Vector& z) const {
        if (r.size() != diag_.size()) return false;
        z.resize(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) z(i) = r(i) / diag_[i];
        return true;
    }
private:
    std::vector<double> diag_;
};

class GaussSeidelPreconditioner {
public:
    bool setup(const SparseMatrix& A) { A_ = &A; return A.n_rows() == A.n_cols(); }
    bool apply(const Vector& r, Vector& z) const {
        if (!A_ || r.size() != A_->n_rows()) return false;
        z.resize(r.size());
        for (std::size_t i = 0; i < A_->n_rows(); ++i) {
            double rhs = r(i);
            double diag = 0.0;
            const auto* v = A_->values_data();
            const auto* c = A_->columns_data();
            const auto* ro = A_->row_offsets_data();
            for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
                if (c[k] == i) diag = v[k];
                else if (c[k] < i) rhs -= v[k] * z(c[k]);
            }
            if (diag == 0.0) return false;
            z(i) = rhs / diag;
        }
        return true;
    }
private:
    const SparseMatrix* A_ = nullptr;
};

/** ILU(0) using the sparsity pattern of A. */
class ILU0Preconditioner {
public:
    bool setup(const SparseMatrix& A) {
        if (A.n_rows() != A.n_cols()) return false;
        n_ = A.n_rows();
        L_.assign(n_, {});
        U_.assign(n_, {});
        const auto* v = A.values_data();
        const auto* c = A.columns_data();
        const auto* ro = A.row_offsets_data();
        for (std::size_t i = 0; i < n_; ++i) {
            for (std::size_t k = ro[i]; k < ro[i+1]; ++k) {
                if (c[k] < i) L_[i][c[k]] = v[k];
                else U_[i][c[k]] = v[k];
            }
        }
        for (std::size_t i = 0; i < n_; ++i) {
            for (auto& [j, lij] : L_[i]) {
                auto it = U_[j].find(j);
                if (it == U_[j].end() || it->second == 0.0) return false;
                lij /= it->second;
                for (const auto& [k, ujk] : U_[j]) {
                    if (k > j) {
                        auto uit = U_[i].find(k);
                        if (uit != U_[i].end()) uit->second -= lij * ujk;
                    }
                }
            }
            auto di = U_[i].find(i);
            if (di == U_[i].end() || std::abs(di->second) <= 1e-30) return false;
        }
        return true;
    }
    bool apply(const Vector& r, Vector& z) const {
        if (r.size() != n_) return false;
        z.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            double s = r(i);
            for (const auto& [j, lij] : L_[i]) s -= lij * z(j);
            z(i) = s;
        }
        for (std::size_t ii = n_; ii-- > 0;) {
            const std::size_t i = ii;
            double s = z(i);
            for (const auto& [j, uij] : U_[i]) if (j > i) s -= uij * z(j);
            z(i) = s / U_[i].at(i);
        }
        return true;
    }
private:
    std::size_t n_ = 0;
    std::vector<std::unordered_map<std::size_t,double>> L_, U_;
};

} // namespace core
} // namespace cfdx
