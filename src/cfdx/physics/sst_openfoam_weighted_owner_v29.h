#pragma once
// CFDX V29: independent algebraic reference for the weighted owner-cell
// constraint of OpenFOAM Foundation v11 omegaWallFunction/manipulateMatrix.
// Not a claim of full equivalence of the wall turbulence model. CFDX already
// supplies omega/G at the wall; this addresses only the fvMatrix operation.
#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct SSTV29WallMatrixAudit {
    std::size_t full_rows=0, partial_rows=0, skipped_rows=0;
    std::size_t eliminated_offdiagonals=0;
};
// OpenFOAM Foundation v11 weights w = (coverage - 0.1)/(1-0.1), with
// coverage <= 0.1 treated as zero. Use a separate source-level test for this
// convention; do not confuse raw geometric area with patch coverage ratio.
inline double sst_v29_openfoam_wall_weight(double coverage) {
    if (!std::isfinite(coverage) || coverage<0. || coverage>1.)
        throw std::invalid_argument("V29 wall fraction: coverage outside [0,1]");
    return coverage<=0.1 ? 0. : (coverage-0.1)/0.9;
}

inline SSTV29WallMatrixAudit sst_apply_weighted_owner_matrix_v29(
    ScalarEquation& equation, const std::vector<unsigned char>& mask,
    const std::vector<double>& target,
    const std::vector<double>& raw_coverage,
    bool full_coverage_only=true)
{
    const auto n=equation.rhs.size();
    if (!n || equation.matrix.n_rows()!=n || equation.matrix.n_cols()!=n ||
        equation.diagonal.size()!=n || mask.size()!=n || target.size()!=n ||
        raw_coverage.size()!=n)
        throw std::invalid_argument("V29 wall matrix: invalid shape");
    std::vector<double> w(n,0.);
    std::vector<std::size_t> diag(n);
    SSTV29WallMatrixAudit audit;
    const auto* offsets=equation.matrix.row_offsets_data();
    const auto* columns=equation.matrix.columns_data();
    // Validate *all* specified rows and coefficients before modifying data.
    for (std::size_t row=0;row<n;++row) {
        if (!mask[row]) continue;
        if (!(target[row]>0.) || !std::isfinite(target[row]))
            throw std::invalid_argument("V29 wall matrix: invalid owner target");
        w[row]=sst_v29_openfoam_wall_weight(raw_coverage[row]);
        if (w[row]==0.) {++audit.skipped_rows;continue;}
        if (full_coverage_only && w[row]<1.-1e-12)
            throw std::invalid_argument("V29 wall matrix: partial patch not certified for production");
        bool found=false;
        for (auto j=offsets[row];j<offsets[row+1];++j)
            if (columns[j]==row) {diag[row]=j;found=true;break;}
        if (!found || !(equation.diagonal[row]>0.) ||
            !std::isfinite(equation.diagonal[row]))
            throw std::invalid_argument("V29 wall matrix: invalid diagonal");
        if (w[row]>=1.-1e-12) ++audit.full_rows;
        else ++audit.partial_rows;
    }
    // A single staging copy ensures invalid data never half-fix a matrix.
    auto staged=equation;
    auto* a=staged.matrix.values_data();
    // For every fully constrained column, move its contribution to all
    // UNCONSTRAINED neighboring RHS rows; it is not sufficient to zero only
    // the owner row (V18 legacy behaviour).
    for (std::size_t row=0;row<n;++row) {
        if (mask[row] && w[row]>=1.-1e-12) continue;
        for (auto j=offsets[row];j<offsets[row+1];++j) {
            const auto col=columns[j];
            if (col<n && mask[col] && w[col]>=1.-1e-12 && col!=row) {
                staged.rhs(row)-=a[j]*target[col];
                if (!std::isfinite(staged.rhs(row)))
                    throw std::invalid_argument("V29 wall matrix: RHS overflow");
                a[j]=0.;
                ++audit.eliminated_offdiagonals;
            }
        }
    }
    for (std::size_t row=0;row<n;++row) {
        if (!mask[row] || w[row]==0.) continue;
        if (w[row]>=1.-1e-12) {
            for (auto j=offsets[row];j<offsets[row+1];++j)
                if (columns[j]!=row) a[j]=0.;
            // Preserve the assembled diagonal as OpenFOAM setValue does.
            staged.rhs(row)=staged.diagonal[row]*target[row];
        } else {
            // Partial setValue in a STEADY equation without extra relaxation
            // or a ddt diagonal: coeff=f/(1-f)*diag, rhs+=coeff*target.
            const double delta=w[row]/(1.-w[row])*staged.diagonal[row];
            if (!(delta>0.) || !std::isfinite(delta) ||
                !std::isfinite(staged.rhs(row)+delta*target[row]))
                throw std::invalid_argument("V29 wall matrix: fractional overflow");
            a[diag[row]]+=delta;
            staged.diagonal[row]+=delta;
            staged.rhs(row)+=delta*target[row];
        }
    }
    equation=std::move(staged);
    return audit;
}
} // namespace cfdx::physics
