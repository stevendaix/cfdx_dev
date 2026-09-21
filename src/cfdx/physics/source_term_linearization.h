// M0.7-T06 — Source Term Linearization (Physics-specific)
//
// Spécification CFDX v0.7 §32 :
//   Terme source linéarisé : S(φ) = Su + Sp * φ
//
// Note: Core source term utilities are in cfdx/core/numerics/source_term.h
// Transport property models are in cfdx/physics/transport_models.h
// Sutherland parameters are also defined there.

#pragma once

#include <vector>
#include <stdexcept>
#include <cmath>

namespace cfdx {
namespace physics {

inline void validate_mass_fractions(const std::vector<double>& Y, double tol = 1e-10) {
    double sum = 0.0;
    for (size_t i = 0; i < Y.size(); ++i) {
        if (Y[i] < -tol) {
            throw std::runtime_error("Mass fraction Y[" + std::to_string(i) + "] must be >= 0");
        }
        sum += Y[i];
    }
    if (std::abs(sum - 1.0) > tol) {
        throw std::runtime_error("Sum of mass fractions must equal 1, got " + std::to_string(sum));
    }
}

}  // namespace physics
}  // namespace cfdx