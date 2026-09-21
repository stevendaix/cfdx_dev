#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
namespace cfdx::physics {
struct AdaptiveCflControls{double target_cfl=1.0,min_dt=1e-12,max_dt=1e12,growth_limit=1.25,shrink_limit=0.5;};
inline double adaptive_time_step(double dt,double measured_cfl,const AdaptiveCflControls&c={}){
    if(!(dt>0.0)||!(measured_cfl>0.0)||!std::isfinite(dt)||!std::isfinite(measured_cfl))throw std::invalid_argument("invalid CFL state");
    double f=std::sqrt(c.target_cfl/measured_cfl);f=std::clamp(f,c.shrink_limit,c.growth_limit);return std::clamp(dt*f,c.min_dt,c.max_dt);
}
inline double pseudo_transient_cfl(std::size_t it,double cfl0=0.5,double growth=1.2,double cfl_max=100.0){if(cfl0<=0.0||growth<=1.0||cfl_max<cfl0)throw std::invalid_argument("invalid pseudo-CFL");return std::min(cfl_max,cfl0*std::pow(growth,static_cast<double>(it)));}
}