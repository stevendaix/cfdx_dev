#pragma once
// V60: forensic partition of a candidate-owned six-field physical Jacobian.
// NEVER evaluate a lagged Picard matrix or accepted-state Rhie-Chow flux.
// At every displaced state, recompute the complete FV F and the V57 stage
// tape; subtract only the *same* candidate's genuine primal convection.
// The reported additivity is a diagnostic, NOT a construction of a Newton
// matrix and NOT a replacement for the physical merit or true Krylov gate.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct PhysicalRemainderAdditivityV60 {
    bool audited=false;
    std::size_t evaluations=0;
    std::array<double,6> remainder_blocks{};
    double momentum_remainder=0.,continuity=0.,sst_k=0.,sst_omega=0.;
    double convection=0.,rhie_chow_flux=0.,inverse_diagonal=0.;
    double viscosity=0.,strain_rate=0.,sst_blend=0.;
    double wall_viscosity=0.,wall_production=0.;
    // Scale-free gate for the exact *remainder* (not a physics tolerance).
    double global_remainder=0.;
};

struct PhysicalRemainderDerivativeV60 {
    std::vector<double> remainder,convection,flux,rau,mu,strain,f2,wall_mu,wall_production;
};

inline std::vector<double> physical_remainder_rows_v60(
    const PhysicalCandidateV51& x,const PhysicalResidualV51& r) {
    validate_physical_result_v51(x,r);
    if(!r.stages_v57)
        throw std::invalid_argument("V60 requires a freshly reassembled V57 stage tape");
    auto rows=r.signed_raw;
    for(std::size_t j=0;j<3*x.cells;++j)
        rows[j]-=r.stages_v57->momentum_convection[j];
    return rows;
}

inline double physical_remainder_norm_v60(const std::vector<double>& a,
                                          std::size_t begin=0,
                                          std::size_t end=std::numeric_limits<std::size_t>::max()) {
    long double norm=0.;
    const auto stop=std::min(end,a.size());
    for(std::size_t i=begin;i<stop;++i) {
        if(!std::isfinite(a[i]))
            throw std::invalid_argument("V60 nonfinite residual trace");
        norm+=static_cast<long double>(a[i])*a[i];
    }
    return std::sqrt(static_cast<double>(norm));
}

inline std::vector<double> physical_remainder_difference_v60(
    const std::vector<double>& plus,const std::vector<double>& minus,double h) {
    if(plus.size()!=minus.size() || !(h>0.) || !std::isfinite(h))
        throw std::invalid_argument("V60 unequal trace sizes or invalid difference step");
    std::vector<double> result(plus.size());
    for(std::size_t i=0;i<result.size();++i){
        result[i]=static_cast<double>((static_cast<long double>(plus[i])-
                     static_cast<long double>(minus[i]))/(2.L*h));
        if(!std::isfinite(result[i]))
            throw std::runtime_error("V60 nonfinite reassembled derivative");
    }
    return result;
}

inline PhysicalRemainderDerivativeV60 physical_remainder_probe_v60(
    const PhysicalCandidateV51& x,const std::vector<double>& v,double h,
    const PhysicalResidualV51& plus,const PhysicalResidualV51& minus) {
    if(v.size()!=x.fields.size() || !plus.stages_v57 || !minus.stages_v57)
        throw std::invalid_argument("V60 incomplete physical probe or missing stage");
    const auto& p=*plus.stages_v57;
    const auto& m=*minus.stages_v57;
    PhysicalRemainderDerivativeV60 out;
    out.remainder=physical_remainder_difference_v60(
        physical_remainder_rows_v60(x,plus),physical_remainder_rows_v60(x,minus),h);
    out.convection=physical_remainder_difference_v60(p.momentum_convection,m.momentum_convection,h);
    out.flux=physical_remainder_difference_v60(plus.face_flux,minus.face_flux,h);
    out.rau=physical_remainder_difference_v60(p.inverse_momentum_diagonal,m.inverse_momentum_diagonal,h);
    out.mu=physical_remainder_difference_v60(p.effective_dynamic_viscosity,m.effective_dynamic_viscosity,h);
    out.strain=physical_remainder_difference_v60(p.strain_rate,m.strain_rate,h);
    out.f2=physical_remainder_difference_v60(p.sst_blend_f2,m.sst_blend_f2,h);
    out.wall_mu=physical_remainder_difference_v60(p.wall_momentum_face_mu,m.wall_momentum_face_mu,h);
    out.wall_production=physical_remainder_difference_v60(p.wall_production,m.wall_production,h);
    return out;
}

// On three ACTUAL Krylov directions v, w, v+w, compare complete fresh
// physical states, including newly recomputed RC/rAU/SST. Six residual
// blocks are kept separate to locate the offender; the decomposition uses
// exactly the already-assembled V57 primal convection, not a second oracle.
template<class Evaluate>
PhysicalRemainderAdditivityV60 audit_physical_remainder_additivity_v60(
    const PhysicalCandidateV51& base,const std::vector<double>& v,
    const std::vector<double>& w,double h,Evaluate&& evaluate,
    double k_min=0.,double omega_min=0.) {
    validate_physical_candidate_v51(base,k_min,omega_min);
    if(v.size()!=base.fields.size() || w.size()!=v.size() ||
       !(h>0.) || !std::isfinite(h))
        throw std::invalid_argument("V60 invalid Krylov directions or step");
    PhysicalRemainderAdditivityV60 out;
    std::array<PhysicalRemainderDerivativeV60,3> d;
    for(std::size_t part=0;part<3;++part) {
        std::vector<double> v_part(v.size());
        for(std::size_t i=0;i<v.size();++i) {
            if(!std::isfinite(v[i]) || !std::isfinite(w[i]))
                throw std::invalid_argument("V60 nonfinite Krylov vector");
            v_part[i]=part==0?v[i]:part==1?w[i]:v[i]+w[i];
        }
        auto xp=base,xm=base;
        for(std::size_t i=0;i<v.size();++i) {
            xp.fields[i]+=h*v_part[i];
            xm.fields[i]-=h*v_part[i];
        }
        const auto rp=evaluate_physical_candidate_v51(xp,evaluate,k_min,omega_min);
        const auto rm=evaluate_physical_candidate_v51(xm,evaluate,k_min,omega_min);
        out.evaluations+=2;
        d[part]=physical_remainder_probe_v60(base,v_part,h,rp,rm);
    }
    const auto defect=[&](const std::vector<double>& a,
                          const std::vector<double>& b,
                          const std::vector<double>& sum,
                          std::size_t begin=0,
                          std::size_t end=std::numeric_limits<std::size_t>::max()) {
        if(a.size()!=b.size() || b.size()!=sum.size())
            throw std::runtime_error("V60 changed physical trace topology");
        std::vector<double> mismatch(a.size(),0.);
        for(std::size_t i=begin;i<std::min(end,a.size());++i)
            mismatch[i]=sum[i]-(a[i]+b[i]);
        const double denominator=std::max({physical_remainder_norm_v60(a,begin,end),
            physical_remainder_norm_v60(b,begin,end),
            physical_remainder_norm_v60(sum,begin,end),1e-12});
        return physical_remainder_norm_v60(mismatch,begin,end)/denominator;
    };
    const auto& a=d[0];const auto& b=d[1];const auto& s=d[2];
    for(std::size_t block=0;block<6;++block)
        out.remainder_blocks[block]=defect(a.remainder,b.remainder,s.remainder,
                                          block*base.cells,(block+1)*base.cells);
    out.momentum_remainder=std::max({out.remainder_blocks[0],out.remainder_blocks[1],
                                     out.remainder_blocks[2]});
    out.continuity=out.remainder_blocks[3];
    out.sst_k=out.remainder_blocks[4];out.sst_omega=out.remainder_blocks[5];
    out.global_remainder=defect(a.remainder,b.remainder,s.remainder);
    out.convection=defect(a.convection,b.convection,s.convection);
    out.rhie_chow_flux=defect(a.flux,b.flux,s.flux);
    out.inverse_diagonal=defect(a.rau,b.rau,s.rau);
    out.viscosity=defect(a.mu,b.mu,s.mu);
    out.strain_rate=defect(a.strain,b.strain,s.strain);
    out.sst_blend=defect(a.f2,b.f2,s.f2);
    out.wall_viscosity=defect(a.wall_mu,b.wall_mu,s.wall_mu);
    out.wall_production=defect(a.wall_production,b.wall_production,s.wall_production);
    out.audited=true;
    return out;
}
} // namespace cfdx::physics
