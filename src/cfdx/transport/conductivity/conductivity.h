#pragma once

#include <cmath>
#include <stdexcept>

namespace cfdx {
namespace transport {
namespace conductivity {

inline double constant(double k)
{
    if (!std::isfinite(k) || k < 0.0)
        throw std::invalid_argument("conductivity::constant: k must be finite and >= 0");
    return k;
}

inline double prandtl(double mu, double Cp, double Pr)
{
    if (!std::isfinite(mu) || mu <= 0.0 ||
        !std::isfinite(Cp) || Cp <= 0.0 ||
        !std::isfinite(Pr) || Pr <= 0.0)
        throw std::invalid_argument(
            "conductivity::prandtl: mu, Cp and Pr must be finite and > 0");
    const double k = mu * Cp / Pr;
    if (!std::isfinite(k) || k < 0.0)
        throw std::runtime_error("conductivity::prandtl: non-finite result");
    return k;
}

}  // namespace conductivity
}  // namespace transport
}  // namespace cfdx
