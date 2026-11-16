#pragma once

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/physics/sst_block_preconditioner_v23.h"
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Experimental, FROZEN, linear right preconditioner in the V22 LOG variables.
// A scalar FV coefficient A_{ij} contributes
//   d (R_i/scale_i/left_i) / d log(x_j) ~= A_{ij}*x_j/(scale_i*left_i)
// for j != i. The local 2x2 diagonal is supplied by V23 and includes an
// approximation of cross-source derivatives. Thus we retain both the actual
// off-cell convection/diffusion stencil and the local k-omega coupling.
//
// One forward plus one backward BLOCK Gauss-Seidel sweep from ZERO is a
// fixed linear map P^-1, suitable for right-preconditioned GMRES. It does
// NOT replace the fully reassembled nonlinear residual, actual Jv, or
// physical acceptance test. No PETSc/HYPRE dependency is introduced.
class SSTSpatialBlockPreconditionerV24 {
public:
    SSTSpatialBlockPreconditionerV24(
        const cfdx::core::SparseMatrix& k_matrix,
        const cfdx::core::SparseMatrix& omega_matrix,
        const std::vector<double>& frozen_physical_state,
        const std::vector<double>& frozen_residual_scales,
        const std::vector<double>& frozen_left_diagonal,
        const SSTBlockPreconditionerV23& local_blocks)
        : n_(local_blocks.cells()), inverses_(local_blocks.inverse_blocks()) {
        if (!n_ || frozen_physical_state.size()!=2*n_ ||
            frozen_residual_scales.size()!=2*n_ ||
            frozen_left_diagonal.size()!=2*n_)
            throw std::invalid_argument("V24 spatial SST: inconsistent vector dimensions");
        build_off_diagonal(0,k_matrix,frozen_physical_state,
                           frozen_residual_scales,frozen_left_diagonal);
        build_off_diagonal(1,omega_matrix,frozen_physical_state,
                           frozen_residual_scales,frozen_left_diagonal);
    }

    std::size_t cells() const noexcept { return n_; }
    std::size_t off_cell_coefficients() const noexcept {
        return off_[0].values.size()+off_[1].values.size();
    }

    std::vector<double> apply(const std::vector<double>& rhs) const {
        if (rhs.size()!=2*n_)
            throw std::invalid_argument("V24 spatial SST: RHS dimension mismatch");
        for (double r:rhs) if (!std::isfinite(r))
            throw std::invalid_argument("V24 spatial SST: nonfinite RHS");
        std::vector<double> x(2*n_,0.0);
        // Both passes are linear in rhs because the initial state is exactly
        // zero and the sparse operator / local inverses are frozen.
        for (std::size_t i=0;i<n_;++i) update(i,rhs,x);
        for (std::size_t rev=0;rev<n_;++rev) update(n_-1-rev,rhs,x);
        return x;
    }

private:
    struct OffCell {
        std::vector<std::size_t> offsets;
        std::vector<std::size_t> columns;
        std::vector<double> values;
    };
    std::size_t n_;
    std::vector<SSTBlockPreconditionerV23::Block> inverses_;
    std::array<OffCell,2> off_;

    void build_off_diagonal(
        std::size_t part,const cfdx::core::SparseMatrix& matrix,
        const std::vector<double>& state,const std::vector<double>& scales,
        const std::vector<double>& left) {
        if (matrix.n_rows()!=n_ || matrix.n_cols()!=n_ ||
            !matrix.is_consistent())
            throw std::invalid_argument("V24 spatial SST: invalid FV CSR matrix");
        auto& out=off_[part];out.offsets.resize(n_+1,0);
        const auto* offsets=matrix.row_offsets_data();
        const auto* cols=matrix.columns_data();
        const auto* vals=matrix.values_data();
        for (std::size_t i=0;i<n_;++i) {
            const auto row=part*n_+i;
            const double denominator=scales[row]*left[row];
            if (!(denominator>0.) || !std::isfinite(denominator) ||
                !(state[row]>0.) || !std::isfinite(state[row]))
                throw std::invalid_argument("V24 spatial SST: invalid frozen scaling");
            for (std::size_t e=offsets[i];e<offsets[i+1];++e) {
                const std::size_t j=cols[e];
                if (j>=n_ || !std::isfinite(vals[e]))
                    throw std::invalid_argument("V24 spatial SST: invalid off-cell coefficient");
                if (i==j || vals[e]==0.) continue;
                const double coefficient=vals[e]*state[part*n_+j]/denominator;
                if (!std::isfinite(coefficient))
                    throw std::domain_error("V24 spatial SST: transformed coefficient overflow");
                out.columns.push_back(j);
                out.values.push_back(coefficient);
            }
            out.offsets[i+1]=out.values.size();
        }
    }

    void update(std::size_t i,const std::vector<double>& rhs,
                std::vector<double>& x) const {
        double rk=rhs[i],rw=rhs[n_+i];
        for (std::size_t part=0;part<2;++part) {
            const auto& data=off_[part];
            double* r=part==0?&rk:&rw;
            for (std::size_t e=data.offsets[i];e<data.offsets[i+1];++e)
                *r -= data.values[e]*x[part*n_+data.columns[e]];
        }
        const auto& M=inverses_[i];
        const double dk=M[0]*rk+M[1]*rw;
        const double dw=M[2]*rk+M[3]*rw;
        if (!std::isfinite(dk)||!std::isfinite(dw))
            throw std::domain_error("V24 spatial SST: preconditioner overflow");
        x[i]=dk;x[n_+i]=dw;
    }
};


// Numerical audit of the *true* two-by-two diagonal Jacobian block of
// [Rk,Romega] at a chosen cell, including dynamic wall closure, coupled
// gradient/diffusion, and source derivatives. This is a measurement, not a
// replacement for the matrix-free Jacobian or an acceptance certificate.
// The callback must have the SAME physical residual contract as V22 JFNK.
struct SSTLocalJacobianAuditV24 {
    SSTBlockPreconditionerV23::Block approximate{};
    SSTBlockPreconditionerV23::Block actual{};
    double maximum_relative_difference=0.;
    std::size_t residual_evaluations=0;
};

inline SSTLocalJacobianAuditV24 sst_local_jacobian_audit_v24(
    std::size_t cell,
    const std::vector<double>& physical_state,
    const std::vector<double>& frozen_scales,
    const std::vector<double>& frozen_left,
    const SSTBlockPreconditionerV23& approx,
    const std::function<std::vector<double>(const std::vector<double>&)>& full_residual,
    double logarithmic_step=1e-6)
{
    const std::size_t n=approx.cells();
    if (cell>=n || !full_residual || physical_state.size()!=2*n ||
        frozen_scales.size()!=2*n || frozen_left.size()!=2*n ||
        !(logarithmic_step>0.) || !(logarithmic_step<0.01) ||
        !std::isfinite(logarithmic_step))
        throw std::invalid_argument("V24 Jacobian audit: invalid contract");
    SSTLocalJacobianAuditV24 out;
    out.approximate=approx.raw_blocks()[cell];
    for (std::size_t part=0;part<2;++part) {
        const auto column=part*n+cell;
        if (!(physical_state[column]>0.) || !std::isfinite(physical_state[column]))
            throw std::invalid_argument("V24 Jacobian audit: invalid physical state");
        auto plus=physical_state,minus=physical_state;
        plus[column]*=std::exp(logarithmic_step);
        minus[column]*=std::exp(-logarithmic_step);
        const auto fp=full_residual(plus),fm=full_residual(minus);
        out.residual_evaluations+=2;
        if (fp.size()!=2*n || fm.size()!=2*n)
            throw std::invalid_argument("V24 Jacobian audit: incomplete residual");
        for(std::size_t row_part=0;row_part<2;++row_part) {
            const auto row=row_part*n+cell;
            const double divisor=2.*logarithmic_step*frozen_scales[row]*frozen_left[row];
            if (!(divisor>0.) || !std::isfinite(divisor) ||
                !std::isfinite(fp[row]) || !std::isfinite(fm[row]))
                throw std::invalid_argument("V24 Jacobian audit: nonfinite residual or scales");
            out.actual[row_part*2+part]=(fp[row]-fm[row])/divisor;
        }
    }
    double max_actual=0.,max_delta=0.;
    for (std::size_t j=0;j<4;++j) {
        if (!std::isfinite(out.actual[j]) || !std::isfinite(out.approximate[j]))
            throw std::domain_error("V24 Jacobian audit: invalid derivative");
        max_actual=std::max(max_actual,std::abs(out.actual[j]));
        max_delta=std::max(max_delta,std::abs(out.actual[j]-out.approximate[j]));
    }
    out.maximum_relative_difference=max_delta/std::max(max_actual,1e-30);
    return out;
}

} // namespace cfdx::physics
