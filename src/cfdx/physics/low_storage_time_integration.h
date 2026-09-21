#pragma once
#include "cfdx/core/field/field.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
namespace cfdx::physics {
template<class RHS> inline void ssprk3_step(cfdx::core::Field<double,cfdx::core::Location::CELL>&u,double dt,RHS&&rhs){
    if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument("invalid dt");const std::size_t n=u.size();cfdx::core::Field<double,cfdx::core::Location::CELL> k(n,"rhs"),u0=u;
    rhs(u,k);for(std::size_t i=0;i<n;++i)u(i)+=dt*k(i);
    rhs(u,k);for(std::size_t i=0;i<n;++i)u(i)=0.75*u0(i)+0.25*(u(i)+dt*k(i));
    rhs(u,k);for(std::size_t i=0;i<n;++i)u(i)=(1.0/3.0)*u0(i)+(2.0/3.0)*(u(i)+dt*k(i));
}
template<class RHS> inline void low_storage_rk2_step(cfdx::core::Field<double,cfdx::core::Location::CELL>&u,double dt,RHS&&rhs){
    if(!(dt>0.0)||!std::isfinite(dt))throw std::invalid_argument("invalid dt");const std::size_t n=u.size();cfdx::core::Field<double,cfdx::core::Location::CELL> k(n,"rhs");rhs(u,k);for(std::size_t i=0;i<n;++i)u(i)+=0.5*dt*k(i);rhs(u,k);for(std::size_t i=0;i<n;++i)u(i)+=0.5*dt*k(i);
}
}