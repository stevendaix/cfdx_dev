#pragma once
// V43 nonlinear relaxation of a COUPLED Picard candidate. This policy
// changes neither the certified full-step linear solution nor its tolerance.
// The damped state must have its authoritative Rhie-Chow flux reassembled and
// be assessed against the SAME physical residual gates after SST transport.
#include <cmath>
#include <stdexcept>

namespace cfdx::physics {

inline double coupled_outer_beta_v43(double requested, double retry_scale) {
    if (!std::isfinite(requested) || !(requested>0.0 && requested<=1.0) ||
        !std::isfinite(retry_scale) || !(retry_scale>0.0 && retry_scale<=1.0))
        throw std::invalid_argument("V43 nonlinear COUPLED carrier step out of range");
    const double value=requested*retry_scale;
    if (!(value>0.0 && value<=1.0) || !std::isfinite(value))
        throw std::invalid_argument("V43 nonlinear COUPLED retry step invalid");
    return value;
}

inline double coupled_outer_blend_v43(double previous, double linear_candidate,
                                     double beta) {
    if (!std::isfinite(previous) || !std::isfinite(linear_candidate) ||
        !(beta>0.0 && beta<=1.0) || !std::isfinite(beta))
        throw std::invalid_argument("V43 nonlinear COUPLED trial contains invalid field");
    const double result=previous+beta*(linear_candidate-previous);
    if (!std::isfinite(result))
        throw std::invalid_argument("V43 nonlinear COUPLED field overflow");
    return result;
}

} // namespace cfdx::physics
