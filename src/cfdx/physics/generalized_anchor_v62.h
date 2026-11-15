#pragma once
// V62: one COMMON deterministic candidate-owned base-state perturbation for
// all generalized Jacobian columns. Sampling separate sided branches for each
// column is NOT a Bouligand matrix; all columns here share the SAME anchor.
//
// This is only a bounded CPU forensic selector. The primal F(x), physical
// pressure gauge, accepted candidate and all Newton vetoes are unchanged.
// In particular, the two-sided independent audit AT x still rejects a kink
// that a nearby smooth anchor could hide. No implicit limiter smoothing.
#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

inline double generalized_anchor_coefficient_v62(std::size_t index,
                                                  std::size_t block) {
    // Fixed integer avalanche. In contrast to random, hashed and timer-based
    // seeds, a replay uses exactly the same perturbed *physical state*.
    std::uint64_t z=static_cast<std::uint64_t>(index+1u) ^
        UINT64_C(0x9e3779b97f4a7c15);
    z^=z>>30;
    z*=UINT64_C(0xbf58476d1ce4e5b9);
    z^=z>>27;
    z*=UINT64_C(0x94d049bb133111eb);
    z^=z>>31;
    const double q=static_cast<double>(z & UINT64_C(0xffff))/65535.;
    // k and omega are moved only toward their admissible interiors.
    if(block==4 || block==5) return 0.35+0.65*q;
    return (z & UINT64_C(0x10000)) ? 0.35+0.65*q : -(0.35+0.65*q);
}

inline PhysicalCandidateV51 generalized_anchor_candidate_v62(
    const PhysicalCandidateV51& base,const std::array<double,6>& scales,
    bool pin_pressure,std::size_t pin_cell,double normalized_radius) {
    validate_physical_candidate_v51(base);
    if(!std::isfinite(normalized_radius) || !(normalized_radius>0.) ||
       pin_cell>=base.cells)
        throw std::invalid_argument("V62 invalid anchored radius or pressure gauge");
    for(double s:scales)
        if(!std::isfinite(s) || !(s>0.))
            throw std::invalid_argument("V62 invalid anchored state scale");
    auto center=base;
    // A two-dimensional extruded BFS with Uz==0 must not acquire fictitious
    // transverse velocity in the branch selector. Nonzero-Uz 3D is retained.
    bool extruded_zero_uz=true;
    for(std::size_t c=0;c<base.cells;++c)
        extruded_zero_uz &= (base.fields[2*base.cells+c]==0.);
    for(std::size_t block=0;block<6;++block)
        for(std::size_t c=0;c<base.cells;++c) {
            if((block==3 && pin_pressure && c==pin_cell) ||
               (block==2 && extruded_zero_uz)) continue;
            const std::size_t index=block*base.cells+c;
            center.fields[index]+=normalized_radius*scales[block]*
                generalized_anchor_coefficient_v62(index,block);
        }
    return center;
}
} // namespace cfdx::physics