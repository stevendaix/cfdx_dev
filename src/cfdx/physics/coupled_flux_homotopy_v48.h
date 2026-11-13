#pragma once
// Experimental Rhie-Chow coefficient homotopy for a COUPLED trial. Keep the
// previous committed inverse momentum diagonal alongside the previous face
// flux. For beta in [0,1], use rAU(beta)=(1-beta)*rAU_old+beta*rAU_new
// when reconstructing phi from the SAME trial U,p. Never directly interpolate
// fluxes: the production Rhie-Chow reconstruction remains authoritative.
//
// This is a path-dependent globalization device, NOT a state-only F(X) and
// NOT a proof of physical convergence. In particular, a homotopy trial with
// beta<1 must never produce the final global convergence certificate.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

using CoupledInverseMomentumDiagonalV48 =
    std::array<std::vector<double>, 3>;

inline void validate_coupled_flux_homotopy_v48(
    bool enabled, bool composite_backtracking_enabled) {
    if (enabled && !composite_backtracking_enabled)
        throw std::invalid_argument(
            "V48 flux homotopy requires the V47 composite COUPLED search");
}

inline CoupledInverseMomentumDiagonalV48 blend_inverse_momentum_v48(
    const CoupledInverseMomentumDiagonalV48& previous,
    const CoupledInverseMomentumDiagonalV48& current,
    double beta) {
    if (!std::isfinite(beta) || beta<0.0 || beta>1.0)
        throw std::invalid_argument("V48 invalid trial beta");
    const std::size_t cells=current[0].size();
    if (cells==0)
        throw std::invalid_argument("V48 empty momentum diagonal");
    CoupledInverseMomentumDiagonalV48 blended;
    for (std::size_t d=0;d<3;++d) {
        if (current[d].size()!=cells || previous[d].size()!=cells)
            throw std::invalid_argument("V48 momentum diagonal shape mismatch");
        blended[d].resize(cells);
        for (std::size_t c=0;c<cells;++c) {
            const double a=previous[d][c],b=current[d][c];
            if (!(a>0.0) || !(b>0.0) || !std::isfinite(a) || !std::isfinite(b))
                throw std::invalid_argument("V48 invalid inverse momentum diagonal");
            // Respect exact endpoints. An affine interpolation of positive
            // inverse diagonals is positive; avoid overflow when magnitudes
            // are very different by using convex weighting.
            const double value=beta==0.0 ? a : beta==1.0 ? b
                : (1.0-beta)*a+beta*b;
            if (!(value>0.0) || !std::isfinite(value))
                throw std::invalid_argument("V48 nonfinite blended diagonal");
            blended[d][c]=value;
        }
    }
    return blended;
}

struct CoupledFluxIdentityV48 {
    double relative_l1=0.0;
    double absolute_linf=0.0;
    double reference_linf=0.0;
    std::size_t worst_face=0;
};

template<class OldFace,class ReconstructedFace>
CoupledFluxIdentityV48 check_coupled_flux_identity_v48(
    std::size_t faces, OldFace accepted, ReconstructedFace restored) {
    if (faces==0)
        throw std::invalid_argument("V48 empty face flux");
    long double reference=0.0L, difference=0.0L;
    CoupledFluxIdentityV48 audit;
    for (std::size_t f=0;f<faces;++f) {
        const double a=accepted(f),b=restored(f);
        if (!std::isfinite(a) || !std::isfinite(b))
            throw std::invalid_argument("V48 nonfinite face flux");
        const double error=std::abs(a-b);
        reference+=std::abs(static_cast<long double>(a));
        difference+=error;
        audit.reference_linf=std::max(audit.reference_linf,std::abs(a));
        if (error>audit.absolute_linf) {
            audit.absolute_linf=error;
            audit.worst_face=f;
        }
    }
    audit.relative_l1=static_cast<double>(difference/std::max(reference,1e-300L));
    if (!std::isfinite(audit.relative_l1))
        throw std::invalid_argument("V48 nonfinite flux identity");
    return audit;
}

// An experimental path-dependent flux must not meet a final, operator-
// independent convergence gate merely because the continuation reduced its
// continuity residual. Full beta=1 candidates use the ordinary production
// Rhie-Chow coefficient and remain eligible for V8/V13/V14 certification.
inline bool coupled_flux_homotopy_final_gate_v48(bool enabled, double beta) {
    if (!std::isfinite(beta) || beta<0.0 || beta>1.0)
        throw std::invalid_argument("V48 invalid gate beta");
    return !enabled || beta==1.0;
}

} // namespace cfdx::physics
