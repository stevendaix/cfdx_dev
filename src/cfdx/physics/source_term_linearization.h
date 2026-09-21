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
#include <functional>
#include <string>
#include <utility>
#include <stdexcept>
#include <cmath>
#include "cfdx/core/numerics/source_term.h"

namespace cfdx {
namespace physics {

using LinearizedSourceTerm = cfdx::core::SourceTerm;

// Physics-facing adapter for the core Su + Sp*phi linearization API.
inline LinearizedSourceTerm linearize_source_term(
    const cfdx::core::Field<double, cfdx::core::Location::CELL>& phi,
    std::function<double(double)> source,
    std::function<double(double)> derivative,
    const std::string& name = "physics_source")
{
    return cfdx::core::linearize_source(phi, std::move(source), std::move(derivative), name);
}

inline void validate_source_linearization(const LinearizedSourceTerm& source) {
    if (source.Su.size() != source.Sp.size())
        throw std::runtime_error("source linearization: Su/Sp size mismatch");
    for (std::size_t i = 0; i < source.Su.size(); ++i) {
        if (!std::isfinite(source.Su(i)) || !std::isfinite(source.Sp(i)))
            throw std::runtime_error("source linearization: non-finite coefficient");
        if (source.Sp(i) > 0.0)
            throw std::runtime_error("source linearization: Sp must be <= 0");
    }
}

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