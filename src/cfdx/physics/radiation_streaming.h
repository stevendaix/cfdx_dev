#pragma once
#include "cfdx/core/field/field.h"
#include <cstddef>
#include <stdexcept>
#include <vector>
namespace cfdx::physics {
template<class DirectionSolver> inline void solve_radiation_streaming(std::size_t n,const std::vector<double>&weights,DirectionSolver&&solve,cfdx::core::Field<double,cfdx::core::Location::CELL>&G){
    if(weights.empty()||G.size()!=n)throw std::invalid_argument("invalid radiation streaming dimensions");cfdx::core::Field<double,cfdx::core::Location::CELL>I(n,"I");G.fill(0.0);
    for(std::size_t m=0;m<weights.size();++m){solve(m,I);for(std::size_t c=0;c<n;++c)G(c)+=weights[m]*I(c);}
}
}