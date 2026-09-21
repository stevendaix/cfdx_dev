// Transport Models - Mass Diffusivity
// 
// Architecture réorganisée P2:
//   diffusivity.h - Mass diffusivity models

#pragma once

#include <cmath>

namespace cfdx {
namespace transport {
namespace diffusivity {

inline double constant(double D) {
    return D;
}

inline double schmidt(double mu, double rho, double Sc) {
    return mu / (rho * Sc);
}

}  // namespace diffusivity
}  // namespace transport
}  // namespace cfdx