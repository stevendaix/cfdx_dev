#pragma once
// V53 directional audit of the SAME candidate-owned nonlinear FV function.
// Measures right/left/central and halved-step derivatives on signed raw rows;
// does not infer linear convergence from a derivative agreement.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct JVFidelityAuditV53 {
    bool both_sides_valid=false;
    std::size_t evaluations=0;
    double central_halving_relative=std::numeric_limits<double>::infinity();
    double forward_halving_relative=std::numeric_limits<double>::infinity();
    double sided_nonsmooth_relative=std::numeric_limits<double>::infinity();
    double central_norm=0.;
};

template<class Evaluate>
JVFidelityAuditV53 audit_physical_jv_v53(
    const PhysicalCandidateV51& base,const std::vector<double>& direction,
    Evaluate&& callback,double h,double k_min=0.0,double omega_min=0.0) {
    validate_physical_candidate_v51(base,k_min,omega_min);
    if(direction.size()!=base.fields.size() || !std::isfinite(h) ||
       !(h>0.0) || !(h<=1e-2))
        throw std::invalid_argument("V53 invalid direction or finite-difference step");
    double vnorm=0.;
    for(double v:direction) {
        if(!std::isfinite(v))throw std::invalid_argument("V53 nonfinite direction");
        vnorm=std::max(vnorm,std::abs(v));
    }
    if(vnorm==0.)throw std::invalid_argument("V53 zero Jv direction");
    const auto f0=evaluate_physical_candidate_v51(base,callback,k_min,omega_min);
    JVFidelityAuditV53 audit;audit.evaluations=1;
    auto diff=[&](double step) {
        auto x=base;
        for(std::size_t i=0;i<direction.size();++i)
            x.fields[i]+=step*direction[i];
        const auto r=evaluate_physical_candidate_v51(x,callback,k_min,omega_min);
        ++audit.evaluations;
        return r.signed_raw;
    };
    const auto p=diff(h),ph=diff(0.5*h);
    std::vector<double> backward,bh;
    try {backward=diff(-h);bh=diff(-0.5*h);audit.both_sides_valid=true;}
    catch(const std::invalid_argument&) {audit.both_sides_valid=false;}
    double center_half=0.,forward_half=0.,sided=0.,norm=0.;
    for(std::size_t i=0;i<p.size();++i) {
        const double fwd=(p[i]-f0.signed_raw[i])/h;
        const double fwd_half=(ph[i]-f0.signed_raw[i])/(0.5*h);
        norm=std::max(norm,std::abs(fwd_half));
        forward_half=std::max(forward_half,std::abs(fwd_half-fwd));
        if(audit.both_sides_valid) {
            const double backward_d=(f0.signed_raw[i]-backward[i])/h;
            const double central_d=(p[i]-backward[i])/(2*h);
            const double half_d=(ph[i]-bh[i])/h;
            norm=std::max(norm,std::abs(half_d));
            center_half=std::max(center_half,std::abs(half_d-central_d));
            sided=std::max(sided,std::abs(fwd-backward_d));
        }
    }
    const double divisor=std::max(norm,1e-12);
    audit.central_norm=norm;
    audit.forward_halving_relative=forward_half/divisor;
    if(audit.both_sides_valid) {
        audit.central_halving_relative=center_half/divisor;
        audit.sided_nonsmooth_relative=sided/divisor;
    }
    return audit;
}
} // namespace cfdx::physics
