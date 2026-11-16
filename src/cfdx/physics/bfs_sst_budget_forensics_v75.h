#pragma once
// V75: read-only, common-state forensic ledger for k and omega.
// No turbulence source, wall law, residual tolerance, or accepted state is changed.
// The dominant DELTA is descriptive, not proof of physical causality.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

struct BfsSstTermsV75 {
    double convection=0., diffusion=0., production=0.;
    double destruction=0., cross_diffusion=0., full_residual=0.;
};
struct BfsSstCellSnapshotV75 {
    std::size_t cell=0;
    double k=0., omega=0., nut=0., raux=0., phi3=0.;
    std::array<BfsSstTermsV75,2> equation{}; // [0]=k, [1]=omega
    std::array<double,2> raw{};               // independently reassembled Fphys
};
struct BfsSstEquationDeltaV75 {
    double convection=0., diffusion=0., production=0.;
    double destruction=0., cross_diffusion=0.;
    double raw=0., sum=0., closure_error=0.;
    const char* largest_absolute_contribution="none";
};
struct BfsSstCellDeltaV75 {
    std::size_t cell=0;
    std::array<BfsSstEquationDeltaV75,2> equation{};
    double k_change=0., omega_change=0., nut_change=0.;
    double raux_change=0., phi3_change=0.;
    bool independently_closed=false;
};

inline bool bfs_budget_near_v75(double a,double b,
                                double atol,double rtol,double scale=0.) {
    return std::abs(a-b)<=atol+rtol*std::max({std::abs(a),std::abs(b),scale});
}
inline BfsSstCellDeltaV75 bfs_compare_sst_budget_v75(
    const BfsSstCellSnapshotV75& before,const BfsSstCellSnapshotV75& after,
    double atol=1e-12,double rtol=1e-8) {
    if(before.cell!=after.cell || !(atol>0.) || !(rtol>0.) ||
       !std::isfinite(atol) || !std::isfinite(rtol))
        throw std::invalid_argument("V75 invalid SST budget comparison domain");
    const auto validate=[&](const BfsSstCellSnapshotV75& s) {
        const std::array<double,7> state{s.k,s.omega,s.nut,s.raux,s.phi3,
                                         s.raw[0],s.raw[1]};
        for(const double x:state)
            if(!std::isfinite(x))throw std::invalid_argument("V75 nonfinite state");
        if(s.k<0. || !(s.omega>0.) || s.nut<0. || !(s.raux>0.))
            throw std::invalid_argument("V75 invalid physical SST state");
        for(std::size_t q=0;q<2;++q){
            const auto& t=s.equation[q];
            const std::array<double,6> values{t.convection,t.diffusion,t.production,
                         t.destruction,t.cross_diffusion,t.full_residual};
            for(const double x:values)
                if(!std::isfinite(x))
                    throw std::invalid_argument("V75 nonfinite SST budget term");
            // This checks the existing V16 FVM oracle vs the independently
            // re-evaluated production Fphys row, not the operator's own sum.
            if(!bfs_budget_near_v75(t.full_residual,s.raw[q],atol,rtol))
                throw std::runtime_error("V75 V16/Fphys native-row mismatch");
        }
    };
    validate(before);validate(after);
    BfsSstCellDeltaV75 out;
    out.cell=before.cell;
    out.k_change=after.k-before.k;
    out.omega_change=after.omega-before.omega;
    out.nut_change=after.nut-before.nut;
    out.raux_change=after.raux-before.raux;
    out.phi3_change=after.phi3-before.phi3;
    static constexpr std::array<const char*,5> labels{
        "convection","diffusion","production","destruction","cross_diffusion"};
    for(std::size_t q=0;q<2;++q) {
        const auto& a=before.equation[q];
        const auto& b=after.equation[q];
        auto& d=out.equation[q];
        d.convection=b.convection-a.convection;
        d.diffusion=b.diffusion-a.diffusion;
        d.production=b.production-a.production;
        d.destruction=b.destruction-a.destruction;
        d.cross_diffusion=b.cross_diffusion-a.cross_diffusion;
        d.raw=after.raw[q]-before.raw[q];
        const std::array<double,5> pieces{d.convection,d.diffusion,d.production,
                                           d.destruction,d.cross_diffusion};
        d.sum=0.;
        double sum_abs=0.,largest=0.;
        for(std::size_t i=0;i<pieces.size();++i){
            d.sum+=pieces[i];
            sum_abs+=std::abs(pieces[i]);
            if(std::abs(pieces[i])>largest){
                largest=std::abs(pieces[i]);
                d.largest_absolute_contribution=labels[i];
            }
        }
        d.closure_error=d.sum-d.raw;
        // Include cancellation scale, not a hard 1.0 floor that could mask
        // mismatches when physical SST residuals are << 1.
        if(!bfs_budget_near_v75(d.sum,d.raw,atol,rtol,sum_abs))
            throw std::runtime_error("V75 physical SST delta does not close");
    }
    out.independently_closed=true;
    return out;
}
} // namespace cfdx::physics