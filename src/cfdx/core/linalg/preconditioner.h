#pragma once

#include "sparse_matrix.h"
#include "vector.h"
#include "mixed_precision.h"
#include <cstddef>
#include <vector>

namespace cfdx::core {

class Preconditioner {
public:
    virtual ~Preconditioner() = default;
    virtual bool setup(const SparseMatrix& A) = 0;
    virtual bool apply(const Vector& r, Vector& z) const = 0;
    virtual const char* name() const = 0;
};

class JacobiPreconditioner final : public Preconditioner {
public:
    bool setup(const SparseMatrix& A) override;
    bool apply(const Vector& r, Vector& z) const override;
    const char* name() const override { return "Jacobi"; }
private:
    std::vector<double> inv_diag_;
};

inline bool JacobiPreconditioner::setup(const SparseMatrix& A) {
    if (A.n_rows() != A.n_cols()) return false;
    const auto* row = A.row_offsets_data();
    const auto* col = A.columns_data();
    const auto* val = A.values_data();
    inv_diag_.assign(A.n_rows(), 0.0);
    for (std::size_t i = 0; i < A.n_rows(); ++i) {
        bool found = false;
        for (std::size_t k = row[i]; k < row[i + 1]; ++k) {
            if (col[k] == i) {
                if (val[k] == 0.0) { inv_diag_.clear(); return false; }
                inv_diag_[i] = 1.0 / val[k];
                found = true;
                break;
            }
        }
        if (!found) { inv_diag_.clear(); return false; }
    }
    return true;
}

inline bool JacobiPreconditioner::apply(const Vector& r, Vector& z) const {
    if (r.size() != inv_diag_.size() || z.size() != r.size()) return false;
    for (std::size_t i = 0; i < r.size(); ++i) z(i) = inv_diag_[i] * r(i);
    return true;
}

class MixedPrecisionJacobiPreconditioner final : public Preconditioner {
public:
    explicit MixedPrecisionJacobiPreconditioner(SolverPrecision precision=SolverPrecision::FP32)
        : precision_(precision) {}
    bool setup(const SparseMatrix& A) override {
        if (A.n_rows()!=A.n_cols()) return false;
        const auto* row=A.row_offsets_data(); const auto* col=A.columns_data(); const auto* val=A.values_data();
        inv_diag_.assign(A.n_rows(),0.0);
        for(std::size_t i=0;i<A.n_rows();++i){
            bool found=false;
            for(std::size_t k=row[i];k<row[i+1];++k) if(col[k]==i){
                if(!std::isfinite(val[k]) || val[k]==0.0) { inv_diag_.clear(); return false; }
                inv_diag_[i]=1.0/val[k]; found=true; break;
            }
            if(!found) { inv_diag_.clear(); return false; }
        }
        return true;
    }
    bool apply(const Vector& r, Vector& z) const override {
        if(r.size()!=inv_diag_.size() || z.size()!=r.size()) return false;
        for(std::size_t i=0;i<r.size();++i){
            if(precision_==SolverPrecision::FP32)
                z(i)=static_cast<double>(static_cast<float>(r(i))*static_cast<float>(inv_diag_[i]));
            else
                z(i)=static_cast<double>(r(i))*static_cast<double>(inv_diag_[i]);
        }
        return true;
    }
    const char* name() const override { return "MixedPrecisionJacobi"; }
private:
    SolverPrecision precision_;
    std::vector<double> inv_diag_;
};

} // namespace cfdx::core
