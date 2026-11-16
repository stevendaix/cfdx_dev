#pragma once
// V49: read-only native Rhie-Chow closure for a supplied physical state.
//
// The *input seed* must be determined solely by the requested U and boundary
// conditions (not by the last accepted face flux). The two callbacks must
// compute their results without mutating the solver's accepted state.
// A(phi) is the momentum operator assembled with current physical properties
// and phi; R(U,p,diag(A(phi))) is the production face-flux reconstruction.
// Solve phi = R(U,p,diag(A(phi))) to obtain a history-independent *candidate*
// flux. Report failure if the fixed-point check does not close; never infer
// physical convergence merely from this subsystem's convergence.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct StateOnlyRhieChowV49 {
    bool converged=false;
    std::size_t iterations=0;
    double relative_linf=std::numeric_limits<double>::infinity();
    std::vector<double> flux;
    std::array<std::vector<double>,3> inverse_diagonal;
};

inline double normalized_flux_fixed_point_defect_v49(
    const std::vector<double>& a, const std::vector<double>& b) {
    if (a.empty() || a.size()!=b.size())
        throw std::invalid_argument("V49 incompatible or empty face flux");
    double delta=0.0, scale=0.0;
    for (std::size_t f=0;f<a.size();++f) {
        if (!std::isfinite(a[f]) || !std::isfinite(b[f]))
            throw std::invalid_argument("V49 nonfinite face flux");
        delta=std::max(delta,std::abs(a[f]-b[f]));
        scale=std::max({scale,std::abs(a[f]),std::abs(b[f])});
    }
    // Relative to the actual face-flux scale, not max(1,|phi|): BFS face
    // mass flows can be far below 1 kg/s. Using 1 here would certify a
    // misleadingly small *absolute* defect as a relative closure.
    return delta/std::max(scale,1e-300);
}

inline void validate_native_rau_v49(
    const std::array<std::vector<double>,3>& rAU,std::size_t n) {
    if (n==0) throw std::invalid_argument("V49 empty momentum state");
    for (const auto& component:rAU) {
        if (component.size()!=n)
            throw std::invalid_argument("V49 momentum diagonal size mismatch");
        for (double value:component)
            if (!(value>0.) || !std::isfinite(value))
                throw std::invalid_argument("V49 invalid native inverse momentum diagonal");
    }
}

template<class AssembleInverseDiagonal,class ReconstructFlux>
StateOnlyRhieChowV49 evaluate_state_only_rhie_chow_v49(
    const std::vector<double>& deterministic_seed,
    std::size_t cells,
    AssembleInverseDiagonal&& assemble,
    ReconstructFlux&& reconstruct,
    std::size_t max_iterations=32,
    double tolerance=1e-10,
    double damping=0.5) {
    if (deterministic_seed.empty() || !cells || max_iterations==0 ||
        max_iterations>4096 || !std::isfinite(tolerance) || tolerance<=0.0 ||
        tolerance>=1.0 || !std::isfinite(damping) || damping<=0.0 || damping>1.0)
        throw std::invalid_argument("V49 invalid fixed-point controls");
    auto phi=deterministic_seed;
    for (double f:phi)
        if (!std::isfinite(f))
            throw std::invalid_argument("V49 nonfinite deterministic seed");
    StateOnlyRhieChowV49 result;
    for (std::size_t it=1;it<=max_iterations;++it) {
        const auto rAU=assemble(phi);
        validate_native_rau_v49(rAU,cells);
        const auto reconstructed=reconstruct(rAU);
        const double defect=normalized_flux_fixed_point_defect_v49(phi,reconstructed);
        result.iterations=it;
        result.relative_linf=defect;
        if (defect<=tolerance) {
            // Certify the operator actually returned to the caller. Never
            // return the stale diagonal from before the final phi change.
            const auto final_rAU=assemble(reconstructed);
            validate_native_rau_v49(final_rAU,cells);
            const auto check_flux=reconstruct(final_rAU);
            const double checked=normalized_flux_fixed_point_defect_v49(
                reconstructed,check_flux);
            result.relative_linf=checked;
            if (checked<=tolerance) {
                result.converged=true;
                result.flux=reconstructed;
                result.inverse_diagonal=final_rAU;
                return result;
            }
        }
        // The blending is only a globalization of the fixed-point iteration;
        // the accepted returned flux, if any, always satisfies the UNBLENDED
        // production reconstruction with its own freshly assembled diagonal.
        for (std::size_t f=0;f<phi.size();++f) {
            phi[f]=(1.0-damping)*phi[f]+damping*reconstructed[f];
            if (!std::isfinite(phi[f]))
                throw std::invalid_argument("V49 nonfinite damped flux");
        }
    }
    // No misleading approximate value for nonconverged native flux. The
    // caller may log the defect but must not use it for physical qualification.
    return result;
}

} // namespace cfdx::physics
