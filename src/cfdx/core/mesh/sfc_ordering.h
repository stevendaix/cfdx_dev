#pragma once
#include "cfdx/core/field/field.h"
#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>
namespace cfdx::core {
inline std::uint64_t morton3d(std::uint32_t x,std::uint32_t y,std::uint32_t z){
 auto split=[](std::uint32_t v){std::uint64_t r=0;for(unsigned i=0;i<21;++i)r|=((static_cast<std::uint64_t>(v)>>i)&1ULL)<<(3*i);return r;};
 return split(x)|(split(y)<<1)|(split(z)<<2);
}
inline std::vector<std::uint32_t> morton_order(const std::vector<Vec3>& centres){
 std::vector<std::uint32_t> idx(centres.size());for(std::uint32_t i=0;i<idx.size();++i)idx[i]=i;if(centres.empty())return idx;
 Vec3 lo=centres[0],hi=centres[0];for(const auto&p:centres){lo.x=std::min(lo.x,p.x);lo.y=std::min(lo.y,p.y);lo.z=std::min(lo.z,p.z);hi.x=std::max(hi.x,p.x);hi.y=std::max(hi.y,p.y);hi.z=std::max(hi.z,p.z);}
 auto q=[&](double v,double a,double b){if(b<=a)return 0U;const double t=std::clamp((v-a)/(b-a),0.0,1.0);return static_cast<std::uint32_t>(t*static_cast<double>((1U<<21)-1U));};
 std::sort(idx.begin(),idx.end(),[&](std::uint32_t a,std::uint32_t b){return morton3d(q(centres[a].x,lo.x,hi.x),q(centres[a].y,lo.y,hi.y),q(centres[a].z,lo.z,hi.z))<morton3d(q(centres[b].x,lo.x,hi.x),q(centres[b].y,lo.y,hi.y),q(centres[b].z,lo.z,hi.z));});return idx;
}
} // namespace cfdx::core
