#pragma once
// V61: three-scale candidate-owned V60 remainder, RC and SST audit.
// No finite-difference extrapolation is substituted into accepted Fphys:
// the test is strictly diagnostic and ALL three original fail-closed gates
// must pass independently before it can ever support a Newton direction.
#include "cfdx/physics/physical_remainder_audit_v60.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

struct PhysicalRemainderRatesV61 {
    std::array<PhysicalRemainderAdditivityV60,3> scales;
    std::size_t evaluations=0;
    std::array<double,3> maximum_block_relative{};
    std::array<double,3> maximum_physical_component{};
    // True means all three independently reassembled directional tests pass.
    bool all_scales_additive=false;
    // Informational only; not mathematical proofs of non-differentiability.
    bool persistent_defect_suspected=false;
    bool roundoff_suspected=false;
    std::size_t worst_block=0;
    double final_failure_ratio=0.;
};

template<class Evaluate>
PhysicalRemainderRatesV61 audit_physical_remainder_rates_v61(
    const PhysicalCandidateV51& base,const std::vector<double>& v,
    const std::vector<double>& w,double h,Evaluate&& evaluate,
    double unchanged_consistency_limit=1e-6,double k_min=0.,
    double omega_min=0.) {
    if (!(h>0.) || !std::isfinite(h) ||
        !(unchanged_consistency_limit>0.) ||
        !std::isfinite(unchanged_consistency_limit))
        throw std::invalid_argument("V61 invalid multi-rate physical step/guard");
    PhysicalRemainderRatesV61 out;
    out.all_scales_additive=true;
    for(std::size_t scale=0;scale<3;++scale){
        const double step=std::ldexp(h,-static_cast<int>(scale));
        auto& a=out.scales[scale];
        a=audit_physical_remainder_additivity_v60(base,v,w,step,evaluate,
                                                   k_min,omega_min);
        out.evaluations+=a.evaluations;
        out.maximum_block_relative[scale]=*std::max_element(
            a.remainder_blocks.begin(),a.remainder_blocks.end());
        out.maximum_physical_component[scale]=std::max({
            a.global_remainder,out.maximum_block_relative[scale],
            a.convection,a.rhie_chow_flux});
        out.all_scales_additive &=a.audited &&
            out.maximum_physical_component[scale]<=unchanged_consistency_limit;
    }
    const auto& recent=out.scales.back();
    out.worst_block=static_cast<std::size_t>(std::distance(
        recent.remainder_blocks.begin(),std::max_element(
            recent.remainder_blocks.begin(),recent.remainder_blocks.end())));
    const double coarse=out.maximum_physical_component[0];
    const double middle=out.maximum_physical_component[1];
    const double fine=out.maximum_physical_component[2];
    out.final_failure_ratio=fine/std::max(middle,1e-300);
    // Neither criterion can turn a rejected physical direction into PASS.
    // A persistent kink can be distinguished from O(h^2) truncation when
    // this ratio remains near unity, but roundoff must be investigated too.
    out.persistent_defect_suspected=fine>unchanged_consistency_limit &&
        middle>unchanged_consistency_limit && coarse>unchanged_consistency_limit &&
        fine/middle>0.65 && fine/middle<1.5 &&
        middle/coarse>0.65 && middle/coarse<1.5;
    out.roundoff_suspected=fine>unchanged_consistency_limit &&
        fine>2.*middle && fine>2.*coarse;
    return out;
}
} // namespace cfdx::physics