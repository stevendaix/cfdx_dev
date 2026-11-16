#pragma once
// V34: experimental branch-aware, physically reassembled diagonal corrections
// to the V23/V24 SST right preconditioner. NEVER changes F, wall equations,
// or the physical convergence test. All vectors are [k cells, omega cells].
#include "cfdx/physics/sst_row_trace_v33.h"
#include "cfdx/physics/sst_spatial_block_preconditioner_v24.h"
#include "cfdx/physics/sst_jfnk_v22.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Stable tie-breaking by cell ID. The two independent residual fields select
// the worst k row and worst omega row, rather than assuming that a particular
// mesh will always call them cells 341 and 329.
inline std::vector<std::size_t> sst_v34_hotspot_cells(
    const std::vector<double>& physical_residual,
    const std::vector<double>& positive_scales,
    std::size_t n, std::size_t per_field=1)
{
    if (!n || physical_residual.size()!=2*n || positive_scales.size()!=2*n ||
        per_field>4)
        throw std::invalid_argument("V34 hotspot selection: invalid dimensions");
    std::vector<std::size_t> chosen;
    if (!per_field) return chosen;
    for (std::size_t field=0;field<2;++field) {
        std::vector<std::size_t> order(n);
        for (std::size_t i=0;i<n;++i) {
            order[i]=i;
            const auto row=field*n+i;
            if (!std::isfinite(physical_residual[row]) ||
                !(positive_scales[row]>0.) || !std::isfinite(positive_scales[row]))
                throw std::invalid_argument("V34 hotspot selection: nonfinite FV row");
        }
        std::stable_sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b) {
            const double aa=std::abs(physical_residual[field*n+a])/positive_scales[field*n+a];
            const double bb=std::abs(physical_residual[field*n+b])/positive_scales[field*n+b];
            return aa>bb;
        });
        for (std::size_t j=0;j<std::min(per_field,n);++j) {
            if (std::find(chosen.begin(),chosen.end(),order[j])==chosen.end())
                chosen.push_back(order[j]);
        }
    }
    return chosen;
}

struct SSTV34BranchCorrection {
    explicit SSTV34BranchCorrection(const SSTBlockPreconditionerV23& base)
        : preconditioner(base) {}
    SSTBlockPreconditionerV23 preconditioner;
    std::vector<std::size_t> sampled_cells;
    std::vector<std::size_t> corrected_cells;
    std::size_t residual_evaluations=0;
    std::size_t skipped_nonsmooth=0;
    std::size_t skipped_near_switch=0;
    std::size_t skipped_singular=0;
    std::size_t skipped_small_change=0;
    double maximum_relative_difference=0.;
    double maximum_one_sided_disagreement=0.;
};

// Replacing only the local 2x2 blocks leaves the spatial off-diagonal FV
// coefficients frozen exactly as in V24. True local derivatives are computed
// in LOG coordinates from four fresh whole-FV residual assemblies per cell.
// A one-sided consistency gate detects hidden min/max switches; the V33 branch
// metadata is an additional safeguard, NOT a substitute for that check.
inline SSTV34BranchCorrection sst_v34_correct_hotspot_blocks(
    const SSTBlockPreconditionerV23& approximate,
    const std::vector<std::size_t>& candidates,
    const std::vector<SSTV33SourceBranches>& frozen_branches,
    const std::vector<double>& physical_state,
    const std::vector<double>& fixed_scales,
    const std::vector<double>& frozen_left,
    const SSTPhysicalResidualV22& reassembled_residual,
    double log_step=1e-6, double maximum_one_sided_gap=0.03,
    double minimum_relative_change=0.01)
{
    const std::size_t n=approximate.cells();
    if (!n || n>4096 || candidates.size()>8 ||
        candidates.size()!=frozen_branches.size() ||
        physical_state.size()!=2*n || fixed_scales.size()!=2*n ||
        frozen_left.size()!=2*n || !reassembled_residual ||
        !(log_step>0.) || !(log_step<1e-3) ||
        !(maximum_one_sided_gap>0.) || maximum_one_sided_gap>=1. ||
        !(minimum_relative_change>=0.) || minimum_relative_change>=1.)
        throw std::invalid_argument("V34 branch preconditioner: invalid contract");
    for (std::size_t i=0;i<2*n;++i) {
        if (!(physical_state[i]>0.) || !std::isfinite(physical_state[i]) ||
            !(fixed_scales[i]>0.) || !std::isfinite(fixed_scales[i]) ||
            !(frozen_left[i]>0.) || !std::isfinite(frozen_left[i]))
            throw std::invalid_argument("V34 branch preconditioner: invalid scale/state");
    }
    SSTV34BranchCorrection report(approximate);
    auto blocks=approximate.raw_blocks();
    std::vector<unsigned char> seen(n,0);
    for (std::size_t slot=0;slot<candidates.size();++slot) {
        const std::size_t cell=candidates[slot];
        if (cell>=n || seen[cell]++)
            throw std::invalid_argument("V34 branch preconditioner: duplicate/invalid cell");
        report.sampled_cells.push_back(cell);
        const auto& branch=frozen_branches[slot];
        if (!branch.frozen_source_derivative_valid ||
            branch.close_to_viscosity_switch ||
            branch.close_to_production_switch) {
            ++report.skipped_near_switch;
            continue;
        }
        SSTBlockPreconditionerV23::Block exact{};
        bool smooth=true;
        const auto center=reassembled_residual(physical_state);
        ++report.residual_evaluations;
        if (center.size()!=2*n)
            throw std::invalid_argument("V34 branch preconditioner: invalid center residual");
        for (std::size_t part=0;part<2;++part) {
            const std::size_t col=part*n+cell;
            if (!(physical_state[col]>0.) || !std::isfinite(physical_state[col]))
                throw std::invalid_argument("V34 branch preconditioner: invalid positive state");
            auto plus=physical_state,minus=physical_state;
            plus[col]*=std::exp(log_step);
            minus[col]*=std::exp(-log_step);
            const auto rp=reassembled_residual(plus);
            const auto rm=reassembled_residual(minus);
            report.residual_evaluations+=2;
            if (rp.size()!=2*n || rm.size()!=2*n)
                throw std::invalid_argument("V34 branch preconditioner: residual dimension changed");
            for (std::size_t row_part=0;row_part<2;++row_part) {
                const auto row=row_part*n+cell;
                const double scale=fixed_scales[row]*frozen_left[row];
                if (!(scale>0.) || !std::isfinite(scale) ||
                    !std::isfinite(rp[row]) || !std::isfinite(rm[row]) ||
                    !std::isfinite(center[row]))
                    throw std::domain_error("V34 branch preconditioner: invalid sampled row");
                const double forward=(rp[row]-center[row])/(log_step*scale);
                const double backward=(center[row]-rm[row])/(log_step*scale);
                const double derivative=0.5*(forward+backward);
                exact[row_part*2+part]=derivative;
                // Compare unilateral derivatives of the SAME frozen rows.
                // The row-level relative gap is evaluated only when the row
                // derivative is resolvable; floating-point noise around a
                // zero derivative must not masquerade as a limiter switch.
                const double scale_direction=std::max({std::abs(forward),
                    std::abs(backward),1e-8});
                const double disagreement=std::abs(forward-backward)/scale_direction;
                if (!std::isfinite(disagreement))
                    throw std::domain_error("V34 branch preconditioner: nonfinite directional difference");
                report.maximum_one_sided_disagreement=std::max(
                    report.maximum_one_sided_disagreement,disagreement);
                if (disagreement>maximum_one_sided_gap) smooth=false;
            }
        }
        if (!smooth) {++report.skipped_nonsmooth; continue;}
        double norm=0.,difference=0.;
        for (std::size_t j=0;j<4;++j) {
            if (!std::isfinite(exact[j]))
                throw std::domain_error("V34 branch preconditioner: nonfinite Jacobian");
            norm=std::max(norm,std::abs(exact[j]));
            difference=std::max(difference,std::abs(exact[j]-blocks[cell][j]));
        }
        const double relative=difference/std::max(norm,1e-30);
        report.maximum_relative_difference=std::max(report.maximum_relative_difference,relative);
        if (!(relative>minimum_relative_change)) {
            ++report.skipped_small_change;
            continue;
        }
        // Only insert locally invertible factors. Never claim that this
        // protects the nonlinear solver; its physical gate remains separate.
        SSTBlockPreconditionerV23 trial({exact});
        if (trial.diagonal_fallbacks()!=0) {
            ++report.skipped_singular;
            continue;
        }
        blocks[cell]=exact;
        report.corrected_cells.push_back(cell);
    }
    report.preconditioner=SSTBlockPreconditionerV23(blocks);
    return report;
}

// Read-only oracle for an OFF-DIAGONAL, cross-cell true FV Jacobian block.
// Unlike the V24 frozen transport CSR, the full Jv captures how viscosity,
// blending, gradients, and refreshed wall values depend on the neighbor.
struct SSTV34CrossCellAudit {
    SSTBlockPreconditionerV23::Block approximate{};
    SSTBlockPreconditionerV23::Block actual{};
    double maximum_relative_difference=0.;
    std::size_t residual_evaluations=0;
};

inline std::array<std::size_t,2> sst_v34_parse_cross_pair(
    const char* specification,std::size_t cell_count)
{
    if (!specification)
        throw std::invalid_argument("V34 cross-cell pair is missing");
    std::string raw(specification);
    if (std::count(raw.begin(),raw.end(),':')!=1)
        throw std::invalid_argument("V34 cross-cell pair requires row:column");
    std::replace(raw.begin(),raw.end(),':',',');
    const auto cells=sst_v33_trace_cells(raw.c_str(),cell_count);
    if (cells.size()!=2)
        throw std::invalid_argument("V34 cross-cell pair requires two distinct cells");
    return {cells[0],cells[1]};
}

inline double sst_v34_frozen_offcell_coefficient(
    const cfdx::core::SparseMatrix& matrix,
    std::size_t row,std::size_t column)
{
    if (!matrix.is_consistent() || matrix.n_rows()!=matrix.n_cols() ||
        row>=matrix.n_rows() || column>=matrix.n_cols() || row==column)
        throw std::invalid_argument("V34 cross-cell: invalid FV matrix/indices");
    double value=0.;
    const auto* offsets=matrix.row_offsets_data();
    const auto* cols=matrix.columns_data();
    const auto* vals=matrix.values_data();
    for (std::size_t e=offsets[row];e<offsets[row+1];++e)
        if (cols[e]==column) value+=vals[e];
    if (!std::isfinite(value))
        throw std::domain_error("V34 cross-cell: nonfinite FV coefficient");
    return value;
}

inline SSTV34CrossCellAudit sst_v34_cross_cell_jacobian_audit(
    std::size_t row_cell,std::size_t column_cell,
    const std::vector<double>& state,
    const std::vector<double>& fixed_scales,
    const std::vector<double>& frozen_left,
    const SSTPhysicalResidualV22& full_residual,
    double frozen_kk,double frozen_ww,
    double log_step=1e-6)
{
    const std::size_t n=state.size()/2;
    if (!n || n>4096 || state.size()!=2*n ||
        fixed_scales.size()!=2*n || frozen_left.size()!=2*n ||
        row_cell>=n || column_cell>=n || row_cell==column_cell ||
        !full_residual || !std::isfinite(frozen_kk) || !std::isfinite(frozen_ww) ||
        !(log_step>0.) || !(log_step<1e-3))
        throw std::invalid_argument("V34 cross-cell: invalid operator contract");
    SSTV34CrossCellAudit audit;
    for (std::size_t field=0;field<2;++field) {
        const std::size_t col=field*n+column_cell;
        const std::size_t row=field*n+row_cell;
        const double scale=fixed_scales[row]*frozen_left[row];
        if (!(state[col]>0.) || !std::isfinite(state[col]) ||
            !(scale>0.) || !std::isfinite(scale))
            throw std::invalid_argument("V34 cross-cell: invalid log variable/scale");
    }
    audit.approximate={
        frozen_kk*state[column_cell]/(fixed_scales[row_cell]*frozen_left[row_cell]),
        0.,0.,
        frozen_ww*state[n+column_cell]/(fixed_scales[n+row_cell]*frozen_left[n+row_cell])};
    for (std::size_t field=0;field<2;++field) {
        const auto col=field*n+column_cell;
        auto plus=state,minus=state;
        plus[col]*=std::exp(log_step);
        minus[col]*=std::exp(-log_step);
        const auto rp=full_residual(plus),rm=full_residual(minus);
        audit.residual_evaluations+=2;
        if (rp.size()!=2*n || rm.size()!=2*n)
            throw std::invalid_argument("V34 cross-cell: residual dimension changed");
        for (std::size_t row_field=0;row_field<2;++row_field) {
            const auto row=row_field*n+row_cell;
            const double divisor=2.*log_step*fixed_scales[row]*frozen_left[row];
            if (!(divisor>0.) || !std::isfinite(divisor) ||
                !std::isfinite(rp[row]) || !std::isfinite(rm[row]))
                throw std::domain_error("V34 cross-cell: nonfinite sampled residual");
            audit.actual[row_field*2+field]=(rp[row]-rm[row])/divisor;
        }
    }
    double maximum=0.,difference=0.;
    for (std::size_t component=0;component<4;++component) {
        maximum=std::max(maximum,std::abs(audit.actual[component]));
        difference=std::max(difference,std::abs(
            audit.actual[component]-audit.approximate[component]));
    }
    audit.maximum_relative_difference=difference/std::max(maximum,1e-30);
    return audit;
}

} // namespace cfdx::physics
