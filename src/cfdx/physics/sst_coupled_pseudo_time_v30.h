#pragma once
// CFDX V30 experimental coupled SST pseudo-transient continuation.
// This is an algebraic surrogate of a local FV spectral radius, NOT a
// physical-time scheme and NOT a complete SU2/Fluent pseudo-time algorithm.
// The added mass changes only the Newton direction. The final certificate
// MUST evaluate the original, steady, wall-refreshed equations.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct SSTV30ContinuationState {
    double cfl = 0.5;
    bool frozen_by_outer_v31 = false;
    // Disabled unless explicitly enabled by the V32 outer controller. At a
    // positive value this caps the k-omega shifted Jacobian at the same
    // relative stationary diagonal fraction used by momentum.
    double v32_max_diagonal_fraction = std::numeric_limits<double>::infinity();
    std::size_t accepted_inner = 0;
    std::size_t rejected_inner = 0;
    std::size_t accepted_outer = 0;
    // Diagnostic only: accept_inner historically updates CFL after ANY
    // physically accepted SST direction. These counters make attribution
    // explicit without silently changing the V30 continuation policy.
    std::size_t selected_own = 0;
    std::size_t selected_other = 0;
    void record_selected_source(bool own_ptc_direction) noexcept {
        if (own_ptc_direction) ++selected_own;
        else ++selected_other;
    }
    static constexpr double lower_cfl = 1.0/128.0;
    static constexpr double upper_cfl = 16.0;

    void validate() const {
        if (!std::isfinite(cfl) || cfl < lower_cfl || cfl > upper_cfl)
            throw std::invalid_argument("V30 SST continuation: invalid CFL");
    }
    void accept_inner(double before, double after) {
        validate();
        if (!(before > 0.) || !std::isfinite(before) ||
            !(after >= 0.) || !std::isfinite(after) || !(after < before))
            throw std::invalid_argument("V30 SST continuation: no true physical decrease");
        const double fraction = (before-after)/before;
        if (!frozen_by_outer_v31)
            cfl=std::min(upper_cfl,cfl*(fraction > 0.05 ? 1.20 : 1.05));
        ++accepted_inner;
    }
    void reject_inner() {
        validate();
        if (!frozen_by_outer_v31)
            cfl=std::max(lower_cfl,0.5*cfl);
        ++rejected_inner;
    }
    void outer_commit() {
        validate();
        ++accepted_outer;
    }
};

// Both FV equations share one CFL; the local scales are allowed to differ.
// The diagonal is based on the actual assembled matrix, its row coefficient
// norm, and the absolute face mass flux entering a cell. This incorporates
// spatial transport, unlike V29's diagonal-only surrogate. Owner omega rows
// already constrained by the wall function have zero pseudo mass.
inline std::vector<double> sst_v30_coupled_pseudo_mass(
    const std::vector<double>& diagonal,
    const std::vector<double>& row_norm,
    const std::vector<double>& absolute_face_flux,
    const std::vector<unsigned char>& constrained,
    double cfl,
    double maximum_diagonal_fraction=std::numeric_limits<double>::infinity())
{
    if (!(maximum_diagonal_fraction>0.0) ||
        std::isnan(maximum_diagonal_fraction))
        throw std::invalid_argument("V32 SST pseudo mass: invalid relative cap");
    const std::size_t n=diagonal.size();
    if (!n || row_norm.size()!=n || absolute_face_flux.size()!=n ||
        (!constrained.empty() && constrained.size()!=n) ||
        !std::isfinite(cfl) || cfl<SSTV30ContinuationState::lower_cfl ||
        cfl>SSTV30ContinuationState::upper_cfl)
        throw std::invalid_argument("V30 SST pseudo mass: invalid sizes/CFL");
    std::vector<double> masses(n,0.);
    for (std::size_t i=0;i<n;++i) {
        if (!(row_norm[i]>0.) || !std::isfinite(row_norm[i]) ||
            !std::isfinite(diagonal[i]) ||
            !(absolute_face_flux[i]>=0.) || !std::isfinite(absolute_face_flux[i]))
            throw std::invalid_argument("V30 SST pseudo mass: invalid FV coefficients");
        if (!constrained.empty() && constrained[i]) continue;
        const double spectral=std::max({std::abs(diagonal[i]),
                                        0.5*row_norm[i],
                                        0.5*absolute_face_flux[i]});
        masses[i]=spectral/cfl;
        if (std::isfinite(maximum_diagonal_fraction)) {
            const double reference=std::abs(diagonal[i]);
            if (!(reference>0.))
                throw std::invalid_argument("V32 SST pseudo mass: zero steady diagonal");
            masses[i]=std::min(masses[i],maximum_diagonal_fraction*reference);
        }
        if (!(masses[i]>0.) || !std::isfinite(masses[i]))
            throw std::invalid_argument("V30 SST pseudo mass: overflow");
    }
    return masses;
}

// Normalize the augmented, LOG-variable Newton rows with the diagonal of
// J_steady + D_tau. This is LEFT row scaling only, not a substitute for the
// full J*v evaluation. Keep the legacy [0.02,50] numerical bounds unchanged.
inline std::vector<double> sst_v30_augmented_row_scaling(
    const std::vector<double>& steady_diagonal,
    const std::vector<double>& pseudo_mass,
    const std::vector<double>& physical_state,
    const std::vector<double>& residual_scale)
{
    const std::size_t n=steady_diagonal.size();
    if (!n || pseudo_mass.size()!=n || physical_state.size()!=n ||
        residual_scale.size()!=n)
        throw std::invalid_argument("V30 augmented row scaling: dimensions");
    std::vector<double> scaling(n);
    for (std::size_t i=0;i<n;++i) {
        if (!std::isfinite(steady_diagonal[i]) ||
            !(pseudo_mass[i]>=0.) || !std::isfinite(pseudo_mass[i]) ||
            !(physical_state[i]>0.) || !std::isfinite(physical_state[i]) ||
            !(residual_scale[i]>0.) || !std::isfinite(residual_scale[i]))
            throw std::invalid_argument("V30 augmented row scaling: invalid coefficient");
        const double effective=steady_diagonal[i]+pseudo_mass[i];
        const double normalized=std::abs(effective*physical_state[i]/residual_scale[i]);
        if (!std::isfinite(normalized))
            throw std::invalid_argument("V30 augmented row scaling: overflow");
        scaling[i]=std::clamp(normalized,0.02,50.0);
    }
    return scaling;
}

// Return a shifted residual for JFNK: F(x) + D_tau*(x-x_old).
// The reference state has exactly the same residual, regardless of CFL.
// This function is pure: it changes neither the original residual nor state.
inline std::vector<double> sst_v30_shift_residual(
    const std::vector<double>& unshifted,
    const std::vector<double>& candidate,
    const std::vector<double>& reference,
    const std::vector<double>& mass)
{
    const auto n=reference.size();
    if (!n || candidate.size()!=n || unshifted.size()!=n || mass.size()!=n)
        throw std::invalid_argument("V30 SST shifted residual: invalid sizes");
    std::vector<double> shifted(n);
    for (std::size_t i=0;i<n;++i) {
        if (!std::isfinite(unshifted[i]) || !std::isfinite(candidate[i]) ||
            !std::isfinite(reference[i]) || !(mass[i]>=0.) || !std::isfinite(mass[i]))
            throw std::invalid_argument("V30 SST shifted residual: invalid field");
        shifted[i]=unshifted[i]+mass[i]*(candidate[i]-reference[i]);
        if (!std::isfinite(shifted[i]))
            throw std::invalid_argument("V30 SST shifted residual: overflow");
    }
    return shifted;
}
} // namespace cfdx::physics
