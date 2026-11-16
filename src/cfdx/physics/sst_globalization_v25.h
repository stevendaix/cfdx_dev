#pragma once
#include "cfdx/physics/sst_spatial_block_preconditioner_v24.h"
#include "cfdx/physics/sst_jfnk_v22.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

// The user-facing physical merit is recomputed from actual FV equations.
// A trial callback MAY project omega owner cells and must return the projected
// fields and their true (dynamic-row-normalization) SST defects together.
// This is NOT the fixed-scaled Krylov merit or an OpenFOAM validation.
struct SSTV25PhysicalTrial {
    std::vector<double> fields;
    double merit=std::numeric_limits<double>::infinity();
    double smooth_l2=std::numeric_limits<double>::infinity(); // optional V26; V25 ignores it
    double k_local=std::numeric_limits<double>::infinity();
    double omega_local=std::numeric_limits<double>::infinity();
    std::size_t k_worst=0, omega_worst=0;
};
using SSTV25PhysicalEvaluator =
    std::function<SSTV25PhysicalTrial(const std::vector<double>&)>;

struct SSTV25GlobalizationResult {
    bool accepted=false;
    std::vector<double> fields;
    double alpha=0.;
    double merit=std::numeric_limits<double>::infinity();
    double best_trial_merit=std::numeric_limits<double>::infinity();
    double full_step_merit=std::numeric_limits<double>::infinity();
    double maximum_absolute_log_direction=0.;
    std::size_t attempted=0, valid=0, improving=0;
    std::size_t k_worst=0, omega_worst=0;
    std::string reason;
};

// This search is deliberately *separate* from the frozen merit of V22.
// It explores log-space directions that have passed the TRUE GMRES residual
// check, even when V22's internal fixed-scaled Linf test did not accept.
// Large components are initially limited by a fixed log trust radius;
// reducing alpha does not clip, change the equation, or loosen SST gates.
inline SSTV25GlobalizationResult sst_globalize_log_direction_v25(
    const std::vector<double>& initial,
    const std::vector<double>& log_direction,
    std::size_t n_cells, double k_floor, double omega_floor,
    double baseline_merit, const SSTV25PhysicalEvaluator& evaluate,
    std::size_t maximum_backtracks=12, double log_trust_radius=1.5)
{
    if (!n_cells || initial.size()!=2*n_cells ||
        log_direction.size()!=initial.size() || !evaluate ||
        !(k_floor>0.) || !(omega_floor>0.) ||
        !(baseline_merit>0.) || !std::isfinite(baseline_merit) ||
        !(log_trust_radius>0.) || !std::isfinite(log_trust_radius) ||
        maximum_backtracks>24)
        throw std::invalid_argument("V25 SST log globalizer: invalid contract");
    SSTV25GlobalizationResult result;
    for (std::size_t i=0;i<initial.size();++i) {
        const double floor=i<n_cells?k_floor:omega_floor;
        if (!(initial[i]>floor) || !std::isfinite(initial[i]) ||
            !std::isfinite(log_direction[i]))
            throw std::invalid_argument("V25 SST log globalizer: invalid physical/log state");
        result.maximum_absolute_log_direction=std::max(
            result.maximum_absolute_log_direction,std::abs(log_direction[i]));
    }
    if (!(result.maximum_absolute_log_direction>0.)) {
        result.reason="zero Newton log direction";
        return result;
    }
    const double initial_alpha=std::min(
        1.0,log_trust_radius/result.maximum_absolute_log_direction);
    for (std::size_t step=0;step<=maximum_backtracks;++step) {
        const double alpha=std::ldexp(initial_alpha,-static_cast<int>(step));
        ++result.attempted;
        std::vector<double> candidate(initial.size());
        bool feasible=true;
        for (std::size_t i=0;i<initial.size();++i) {
            const double floor=i<n_cells?k_floor:omega_floor;
            candidate[i]=initial[i]*std::exp(alpha*log_direction[i]);
            if (!std::isfinite(candidate[i]) || !(candidate[i]>=floor)) {
                feasible=false;break;
            }
        }
        if (!feasible) continue;
        try {
            auto measured=evaluate(candidate);
            if (measured.fields.size()!=initial.size() ||
                !std::isfinite(measured.merit) || measured.merit<0. ||
                !std::isfinite(measured.k_local) ||
                !std::isfinite(measured.omega_local) ||
                measured.k_local<0. || measured.omega_local<0. ||
                std::abs(measured.merit-std::max(measured.k_local,measured.omega_local))>
                    1e-12*std::max(1.,measured.merit) ||
                measured.k_worst>=n_cells || measured.omega_worst>=n_cells)
                throw std::domain_error("V25 SST log globalizer: incomplete physical certificate");
            for (std::size_t i=0;i<measured.fields.size();++i) {
                const double floor=i<n_cells?k_floor:omega_floor;
                if (!std::isfinite(measured.fields[i]) || measured.fields[i]<floor)
                    throw std::domain_error("V25 SST log globalizer: invalid projected field");
            }
            ++result.valid;
            if (step==0) result.full_step_merit=measured.merit;
            result.best_trial_merit=std::min(result.best_trial_merit,measured.merit);
            if (measured.merit<baseline_merit*(1.0-1e-9*alpha)) {
                ++result.improving;
                if (!result.accepted || measured.merit<result.merit) {
                    result.accepted=true;
                    result.fields=std::move(measured.fields);
                    result.alpha=alpha;
                    result.merit=measured.merit;
                    result.k_worst=measured.k_worst;
                    result.omega_worst=measured.omega_worst;
                }
            }
        } catch (const std::exception&) {
            // Invalid wall closures or FV trials are not accepted or committed.
        }
    }
    result.reason=result.accepted ? "strict physical merit decreased" :
        "no projected physical step passed strict decrease";
    return result;
}

// Replace a BOUNDED number of wall-owner source-only blocks using the actual
// local 2x2 derivatives of reassembled FV residual. All off-cell couplings,
// wall target derivatives and changes of G/nu_t are included in the measured
// diagonal block. It is *only* a right preconditioner correction; Jv, F and
// boundary conditions are untouched. If measurement fails, retain V23.
struct SSTV25WallBlocks {
    SSTBlockPreconditionerV23 preconditioner;
    std::size_t corrected=0;
    double max_relative_error=0.;
    std::vector<std::size_t> cells;
};

inline SSTV25WallBlocks sst_correct_wall_blocks_v25(
    const SSTBlockPreconditionerV23& old,
    const std::vector<std::size_t>& wall_owner_ids,
    const std::vector<double>& state,
    const std::vector<double>& scales,
    const std::vector<double>& left,
    const SSTPhysicalResidualV22& physical_residual,
    std::size_t maximum_cells=4)
{
    if (!physical_residual || maximum_cells>16)
        throw std::invalid_argument("V25 SST wall blocks: invalid contract");
    auto blocks=old.raw_blocks();
    SSTV25WallBlocks output{old,0,0.,{}};
    if (maximum_cells==0 || wall_owner_ids.empty()) return output;
    std::vector<bool> unique(old.cells(),false);
    for (const auto cell:wall_owner_ids) {
        if (cell>=old.cells() || unique[cell])
            throw std::invalid_argument("V25 SST wall blocks: malformed owner ids");
        unique[cell]=true;
    }
    for (std::size_t i=0;i<std::min(maximum_cells,wall_owner_ids.size());++i) {
        const auto cell=wall_owner_ids[i];
        try {
            const auto audit=sst_local_jacobian_audit_v24(
                cell,state,scales,left,old,physical_residual);
            if (audit.maximum_relative_difference>0.01) {
                blocks[cell]=audit.actual;
                // Do not replace with an indefinite/near-singular local factor.
                SSTBlockPreconditionerV23 checked({blocks[cell]});
                if (checked.diagonal_fallbacks()!=0) {
                    blocks[cell]=audit.approximate;
                    continue;
                }
                ++output.corrected;
                output.max_relative_error=std::max(output.max_relative_error,
                                                   audit.maximum_relative_difference);
                output.cells.push_back(cell);
            }
        } catch (const std::exception&) {
            // The unchanged original local block is the safe fallback.
            blocks[cell]=old.raw_blocks()[cell];
        }
    }
    output.preconditioner=SSTBlockPreconditionerV23(blocks);
    return output;
}

} // namespace cfdx::physics
