#pragma once

// Independent high-precision true residual of the EXACT V14 Arnoldi operator
// Fp + W Z^T. The ILU0/AMG anchored matrix is only a preconditioner; its
// residual must never be substituted for the Krylov operator true residual.
// Physical/original Fp and left-compatibility residuals are separately audited.
#include "cfdx/physics/pcd_fp_bisided_basis_v13.h"
#include "cfdx/physics/pcd_fp_equilibration.h"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

inline PcdFpTrueResidual pcd_fp_deflated_true_residual_v14(
    const cfdx::core::SparseMatrix& original,
    const PcdFpBiNullspaceV13& basis,
    const cfdx::core::Vector& compatible_rhs,
    const cfdx::core::Vector& initial,
    const cfdx::core::Vector& solution)
{
    const auto n=original.n_rows();
    if (pcd_fp_v13_hash_matrix(original)!=basis.operator_hash ||
        original.n_cols()!=n || n!=basis.size ||
        compatible_rhs.size()!=n || initial.size()!=n || solution.size()!=n ||
        !compatible_rhs.is_valid() || !initial.is_valid() || !solution.is_valid())
        throw std::invalid_argument("V14 deflated residual: stale basis or invalid operands");
    std::vector<long double> projection_initial, projection_final;
    projection_initial.reserve(basis.right.size());
    projection_final.reserve(basis.right.size());
    for (const auto& z:basis.right) {
        projection_initial.push_back(static_cast<long double>(basis.dot(z,initial)));
        projection_final.push_back(static_cast<long double>(basis.dot(z,solution)));
    }
    PcdFpTrueResidual result;
    const auto* offsets=original.row_offsets_data();
    const auto* columns=original.columns_data();
    const auto* values=original.values_data();
    for (std::size_t row=0;row<n;++row) {
        long double initial_ax=0.0L, final_ax=0.0L;
        for (std::size_t k=offsets[row];k<offsets[row+1];++k) {
            initial_ax += static_cast<long double>(values[k])*initial(columns[k]);
            final_ax += static_cast<long double>(values[k])*solution(columns[k]);
        }
        for (std::size_t m=0;m<basis.left.size();++m) {
            initial_ax += static_cast<long double>(basis.left[m](row))*projection_initial[m];
            final_ax += static_cast<long double>(basis.left[m](row))*projection_final[m];
        }
        const double start=static_cast<double>(static_cast<long double>(compatible_rhs(row))-initial_ax);
        const double finish=static_cast<double>(static_cast<long double>(compatible_rhs(row))-final_ax);
        result.initial_l2=std::hypot(result.initial_l2,start);
        result.final_l2=std::hypot(result.final_l2,finish);
        result.rhs_l2=std::hypot(result.rhs_l2,compatible_rhs(row));
        result.final_linf=std::max(result.final_linf,std::abs(finish));
    }
    result.relative_initial=result.final_l2/(result.initial_l2>0.?result.initial_l2:1.);
    result.relative_rhs=result.final_l2/(result.rhs_l2>0.?result.rhs_l2:1.);
    return result;
}

} // namespace cfdx::physics
