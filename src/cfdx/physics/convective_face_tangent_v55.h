#pragma once
// V55 forensic, opt-in building block.  This is a conservative, frozen-active-
// set tangent for the first-order donor-cell FACE flux.  It deliberately does
// not pretend to differentiate CFDX's full MUSCL/TVD limiter or Rhie-Chow.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct FaceBranchEventV55 {
    std::size_t face=0;
    double base=0., plus=0., minus=0.;
    double absolute_base=0., relative_base=0.;
    bool sign_crossing=false, near_zero=false;
};
struct FaceBranchSummaryV55 {
    std::size_t crossings=0, near_zero_crossings=0;
    double max_absolute_flux=0.;
    std::vector<FaceBranchEventV55> events;
};
inline FaceBranchSummaryV55 audit_face_branches_v55(
        const std::vector<double>& base,const std::vector<double>& plus,
        const std::vector<double>& minus,double relative_zero=1e-6,
        std::size_t max_events=64) {
    if(base.size()!=plus.size() || base.size()!=minus.size() ||
       !std::isfinite(relative_zero) || relative_zero<0. ||
       !max_events)
        throw std::invalid_argument("V55 invalid face branch audit input");
    FaceBranchSummaryV55 out;
    for(std::size_t i=0;i<base.size();++i){
        if(!std::isfinite(base[i]) || !std::isfinite(plus[i]) ||
           !std::isfinite(minus[i]))
            throw std::invalid_argument("V55 nonfinite face flux");
        out.max_absolute_flux=std::max(out.max_absolute_flux,
            std::max({std::abs(base[i]),std::abs(plus[i]),std::abs(minus[i])}));
    }
    const double scale=out.max_absolute_flux;
    for(std::size_t i=0;i<base.size();++i){
        const bool crossing=(plus[i]>0. && minus[i]<0.) ||
                            (plus[i]<0. && minus[i]>0.);
        if(!crossing)continue;
        ++out.crossings;
        const bool near=scale>0. && std::abs(base[i])<=relative_zero*scale;
        if(near)++out.near_zero_crossings;
        if(out.events.size()<max_events)
            out.events.push_back({i,base[i],plus[i],minus[i],
                std::abs(base[i]),scale>0.?std::abs(base[i])/scale:0.,
                true,near});
    }
    return out;
}
// Positive phi denotes owner-to-neighbour flow.  Exact same face flux must
// be added to one cell and subtracted from the other (conservation).
inline double donor_face_flux_v55(double phi,double owner,double neighbour){
    if(!std::isfinite(phi)||!std::isfinite(owner)||!std::isfinite(neighbour))
        throw std::invalid_argument("V55 nonfinite donor flux input");
    return phi*(phi>=0.?owner:neighbour);
}
inline double donor_face_frozen_tangent_v55(double phi,double owner,
        double neighbour,double dphi,double downer,double dneighbour){
    for(double x:{phi,owner,neighbour,dphi,downer,dneighbour})
        if(!std::isfinite(x))throw std::invalid_argument("V55 nonfinite tangent input");
    return dphi*(phi>=0.?owner:neighbour)+
        phi*(phi>=0.?downer:dneighbour);
}
// Strictly local conservative scatter, independent of the physical solver.
inline void scatter_conservative_face_v55(std::vector<double>& residual,
        std::size_t owner,std::size_t neighbour,double face_flux){
    if(owner>=residual.size() || neighbour>=residual.size() ||
       owner==neighbour || !std::isfinite(face_flux))
        throw std::invalid_argument("V55 invalid face scatter");
    residual[owner]+=face_flux;
    residual[neighbour]-=face_flux;
}
} // namespace cfdx::physics
