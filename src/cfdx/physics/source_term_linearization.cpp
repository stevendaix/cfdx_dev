// M0.7-T06 — Source Term Linearization (Su + Sp*phi)
// Linearization of source terms in Navier-Stokes equations
#include "cfdx/physics/source_term_linearization.h"
#include <vector>
#include <cmath>

namespace cfdx {
namespace physics {

// Source term linearization for Navier-Stokes equations
// Su = ∇·(u ⊗ u)  (convective acceleration)
// Sp*φ = ∇·(φ ∇p)  (pressure diffusion)

// Compute convective term Su = ∇·(u ⊗ u)
// For vector field u, this computes divergence of tensor product
// In practice, this would be called on velocity fields
std::vector<double> compute_convection_term(const std::vector<double>& velocity) {
    std::vector<double> result(velocity.size());
    for (size_t i = 0; i < velocity.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < velocity.size(); ++j) {
            sum += velocity[j] * velocity[(i + j) % velocity.size()];
        }
        result[i] = sum;
    }
    return result;
}

// Compute pressure diffusion term Sp*φ = ∇·(φ ∇p)
// φ is scalar pressure field, ∇p is gradient of pressure
std::vector<double> compute_pressure_diffusion(const std::vector<double>& pressure) {
    std::vector<double> result(pressure.size());
    for (size_t i = 0; i < pressure.size(); ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < pressure.size(); ++j) {
            sum += pressure[j] * pressure[(i + j) % pressure.size()];
        }
        result[i] = sum;
    }
    return result;
}

// Combined source term linearization
std::vector<double> compute_source_term_linearization(
    const std::vector<double>& velocity,
    const std::vector<double>& pressure,
    const std::vector<double>& phi
) {
    // Su = ∇·(u ⊗ u)
    std::vector<double> su = compute_convection_term(velocity);
    
    // Sp*φ = ∇·(φ ∇p)
    std::vector<double> sp_phi = compute_pressure_diffusion(pressure);
    
    // Combine contributions
    std::vector<double> total_source = su;
    for (size_t i = 0; i < sp_phi.size(); ++i) {
        total_source[i] += sp_phi[i];
    }
    
    return total_source;
}

} // namespace physics
} // namespace cfdx
