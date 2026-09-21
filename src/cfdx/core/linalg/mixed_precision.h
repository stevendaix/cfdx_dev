#pragma once
#include "cfdx/core/linalg/vector.h"
#include <cmath>
#include <cstddef>
#include <limits>
namespace cfdx::core {
enum class SolverPrecision{FP64,FP32};
struct PrecisionPolicy{SolverPrecision operator_precision=SolverPrecision::FP32;SolverPrecision reduction_precision=SolverPrecision::FP64;double residual_refresh_factor=10.0;std::size_t refresh_interval=20;bool enabled=true;};
inline double mixed_precision_dot(const Vector&a,const Vector&b,SolverPrecision p=SolverPrecision::FP64){if(a.size()!=b.size())return std::numeric_limits<double>::quiet_NaN();if(p==SolverPrecision::FP32){float s=0.0f;for(std::size_t i=0;i<a.size();++i)s+=static_cast<float>(a(i))*static_cast<float>(b(i));return s;}double s=0.0;for(std::size_t i=0;i<a.size();++i)s+=a(i)*b(i);return s;}
inline double mixed_precision_norm2(const Vector&a,SolverPrecision p){return p==SolverPrecision::FP32?std::sqrt(std::max(0.0,mixed_precision_dot(a,a,p))):a.norm2();}
}