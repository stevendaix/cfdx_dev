#pragma once
// Opt-in, read-only audit of the exact post-SST reassembled momentum CSR.
// The row is selected with the SAME FP64 summation order as
// scalar_equation_residual_inf(), i.e. the actual nonlinear V7/V8 gate.
// Extended precision below is a second diagnostic, NEVER an acceptance oracle.
#include "cfdx/physics/sst_stall_forensics_v42.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::physics {

struct PostSSTMomentumHotspotV45 {
    std::size_t component = 0;
    std::size_t cell = 0;
    double rhs_scale = 0.0;
    double original_signed = 0.0;
    double original_abs = 0.0;
    double original_rhs_relative = 0.0;
    SSTStallRowV42 independent_row;
};

inline PostSSTMomentumHotspotV45 post_sst_momentum_hotspot_v45(
    const ScalarEquation& equation,
    const cfdx::core::Vector& solution,
    std::size_t component,
    double rhs_scale)
{
    if (solution.empty() || solution.size()!=equation.rhs.size() ||
        equation.matrix.n_rows()!=solution.size() ||
        equation.matrix.n_cols()!=solution.size() ||
        !equation.matrix.is_consistent() ||
        !std::isfinite(rhs_scale) || !(rhs_scale>0.0) || component>=3)
        throw std::invalid_argument("V45 post-SST momentum: invalid FV operator");
    PostSSTMomentumHotspotV45 best;
    best.component=component;
    best.rhs_scale=rhs_scale;
    const auto* offsets=equation.matrix.row_offsets_data();
    const auto* columns=equation.matrix.columns_data();
    const auto* values=equation.matrix.values_data();
    for (std::size_t i=0;i<solution.size();++i) {
        double sum=-equation.rhs(i);
        if (!std::isfinite(sum))
            throw std::invalid_argument("V45 post-SST momentum: nonfinite rhs");
        for (std::size_t j=offsets[i];j<offsets[i+1];++j) {
            if (columns[j]>=solution.size() || !std::isfinite(values[j]) ||
                !std::isfinite(solution(columns[j])))
                throw std::invalid_argument("V45 post-SST momentum: nonfinite CSR entry");
            sum+=values[j]*solution(columns[j]);
            if (!std::isfinite(sum))
                throw std::invalid_argument("V45 post-SST momentum: nonfinite residual");
        }
        const double magnitude=std::abs(sum);
        if (i==0 || magnitude>best.original_abs) {
            best.cell=i;
            best.original_signed=sum;
            best.original_abs=magnitude;
        }
    }
    best.original_rhs_relative=best.original_abs/rhs_scale;
    // One independent row, not O(n*nnz) re-validation. This may differ
    // slightly from the native summation because it uses long double.
    best.independent_row=sst_stall_row_v42(
        equation,solution,best.cell,1e-300,true);
    if (!std::isfinite(best.original_rhs_relative) ||
        best.independent_row.native_signed!=best.original_signed)
        throw std::domain_error("V45 post-SST row did not reproduce true FP64 gate");
    return best;
}

} // namespace cfdx::physics
