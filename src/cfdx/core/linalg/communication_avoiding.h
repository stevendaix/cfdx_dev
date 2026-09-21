#pragma once
#include "cfdx/core/linalg/vector.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
namespace cfdx::core {
struct ReductionPacket { double dot=0.0; double norm2=0.0; double max_abs=0.0; };
inline ReductionPacket fused_reduction(const Vector& a,const Vector& b){
    ReductionPacket r; const std::size_t n=a.size(); if(b.size()!=n)return r;
    for(std::size_t i=0;i<n;++i){r.dot+=a(i)*b(i);r.norm2+=a(i)*a(i);r.max_abs=std::max(r.max_abs,std::abs(a(i)));} return r;
}
} // namespace cfdx::core
