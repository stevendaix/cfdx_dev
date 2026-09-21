#pragma once

#include <cmath>
#include <stdexcept>

namespace cfdx {
namespace transport {
namespace diffusivity {

inline double constant(double D)
{
    if (!std::isfinite(D) || D < 0.0)
        throw std::invalid_argument("diffusivity::constant: D must be finite and >= 0");
    return D;
}

inline double schmidt(double mu, double rho, double Sc)
{
    if (!std::isfinite(mu) || mu <= 0.0 ||
        !std::isfinite(rho) || rho <= 0.0 ||
        !std::isfinite(Sc) || Sc <= 0.0)
        throw std::invalid_argument(
            "diffusivity::schmidt: mu, rho and Sc must be finite and > 0");
    const double D = mu / (rho * Sc);
    if (!std::isfinite(D) || D < 0.0)
        throw std::runtime_error("diffusivity::schmidt: non-finite result");
    return D;
}

}  // namespace diffusivity
}  // namespace transport
}  // namespace cfdx
