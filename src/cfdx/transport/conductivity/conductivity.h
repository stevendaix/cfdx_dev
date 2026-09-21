// Transport Models - Thermal Conductivity
// 
// Architecture réorganisée P2:
//   conductivity.h - Thermal conductivity models

#pragma once

#include <cmath>

namespace cfdx {
namespace transport {
namespace conductivity {

inline double constant(double k) {
    return k;
}

inline double prandtl(double mu, double Cp, double Pr) {
    return mu * Cp / Pr;
}

}  // namespace conductivity
}  // namespace transport
}  // namespace cfdx