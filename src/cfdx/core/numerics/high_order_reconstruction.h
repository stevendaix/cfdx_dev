// N16.1 — Experimental high-order finite-volume reconstruction.
//
// This header starts with a self-contained WENO-Z kernel operating on five
// consecutive cell averages on a uniform 1-D stencil. It reconstructs the
// value at the right face of the central cell.
//
// This is an N16 verification primitive, not production WENO support on
// arbitrary meshes.

#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace cfdx::core {

struct WENO5ZDiagnostics {
    double beta0{0.0}, beta1{0.0}, beta2{0.0}, tau5{0.0};
    double weight0{0.0}, weight1{0.0}, weight2{0.0};
    bool finite{true};
    bool bounded{false};
};

namespace detail {
inline void weno5z_weights(const std::array<double, 5>& u, double epsilon,
                           double power, std::array<double, 3>& weights,
                           WENO5ZDiagnostics* diagnostics)
{
    const double u0=u[0], u1=u[1], u2=u[2], u3=u[3], u4=u[4];
    const double beta0=(13.0/12.0)*std::pow(u0-2*u1+u2,2)
        +0.25*std::pow(u0-4*u1+3*u2,2);
    const double beta1=(13.0/12.0)*std::pow(u1-2*u2+u3,2)
        +0.25*std::pow(u1-u3,2);
    const double beta2=(13.0/12.0)*std::pow(u2-2*u3+u4,2)
        +0.25*std::pow(3*u2-4*u3+u4,2);
    const double tau5=std::abs(beta0-beta2);
    constexpr std::array<double,3> d{{0.1,0.6,0.3}};
    if (!(epsilon>0.0) || !std::isfinite(epsilon))
        throw std::invalid_argument("weno5z: epsilon must be finite and > 0");
    if (!(power>0.0) || !std::isfinite(power))
        throw std::invalid_argument("weno5z: power must be finite and > 0");
    const std::array<double,3> beta{{beta0,beta1,beta2}};
    double sum=0.0;
    std::array<double,3> alpha{};
    for (std::size_t k=0;k<3;++k) {
        const double ratio=tau5/(beta[k]+epsilon);
        alpha[k]=d[k]*(1.0+std::pow(ratio,power));
        sum+=alpha[k];
    }
    if (!(sum>0.0) || !std::isfinite(sum))
        throw std::runtime_error("weno5z: invalid nonlinear weight normalization");
    for (std::size_t k=0;k<3;++k) weights[k]=alpha[k]/sum;
    if (diagnostics) {
        diagnostics->beta0=beta0; diagnostics->beta1=beta1; diagnostics->beta2=beta2;
        diagnostics->tau5=tau5; diagnostics->weight0=weights[0];
        diagnostics->weight1=weights[1]; diagnostics->weight2=weights[2];
    }
}
} // namespace detail

// WENO5-Z finite-volume reconstruction at x_{i+1/2}. Input values are cell
// averages and the uniform-grid assumption is part of the mathematical contract.
inline double weno5z_reconstruct_right(const std::array<double,5>& u,
                                       double epsilon=1.0e-40,
                                       double power=2.0,
                                       WENO5ZDiagnostics* diagnostics=nullptr)
{
    for (double value:u)
        if (!std::isfinite(value))
            throw std::invalid_argument("weno5z: stencil contains a non-finite value");
    std::array<double,3> w{};
    detail::weno5z_weights(u,epsilon,power,w,diagnostics);
    const double p0=(1.0/3.0)*u[0]-(7.0/6.0)*u[1]+(11.0/6.0)*u[2];
    const double p1=-(1.0/6.0)*u[1]+(5.0/6.0)*u[2]+(1.0/3.0)*u[3];
    const double p2=(1.0/3.0)*u[2]+(5.0/6.0)*u[3]-(1.0/6.0)*u[4];
    const double value=w[0]*p0+w[1]*p1+w[2]*p2;
    if (!std::isfinite(value))
        throw std::runtime_error("weno5z: non-finite reconstructed value");
    if (diagnostics) diagnostics->finite=true;
    return value;
}

// Explicit bounded policy. Bounding is deliberately separate from WENO
// accuracy so the two properties cannot be conflated.
inline double weno5z_reconstruct_right_bounded(
    const std::array<double,5>& u, double epsilon=1.0e-40, double power=2.0,
    WENO5ZDiagnostics* diagnostics=nullptr)
{
    const double raw=weno5z_reconstruct_right(u,epsilon,power,diagnostics);
    const auto [lo,hi]=std::minmax_element(u.begin(),u.end());
    const double bounded=std::clamp(raw,*lo,*hi);
    if (diagnostics) diagnostics->bounded=(bounded==raw);
    return bounded;
}

} // namespace cfdx::core
