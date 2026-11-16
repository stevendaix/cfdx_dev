#pragma once
// V50: signed, six-field, dimensionless *diagnostic* residual on a shared
// native Rhie-Chow face flux. This is a storage/measurement contract, not a
// Newton solver or a proof that boundary/closure data are history independent.
//
// Ordering: [R(Ux), R(Uy), R(Uz), div(phi), R(k), R(omega)]. Each block
// contains exactly one entry per cell. Momentum is normalized by the same
// component-wise max RHS magnitude used in the CFDX production certificate;
// pressure is normalized by sum(abs(face mass flux)); SST uses the production
// V22 local matrix-row-norm * local-field-scale normalization.
//
// `closed` must certify the V49 flux/diagonal fixed point *before* calling this
// routine. The inputs must all derive from the SAME face flux and state. It is
// intentionally not exposed as the final physical convergence gate.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct NativeSixFieldResidualV50 {
    std::size_t cells=0;
    std::vector<double> normalized;  // 6 * cells, signed, block-major
    std::array<double,6> linf{{}};
    std::array<bool,3> active_momentum{{}};
    double momentum_active_merit=0.0;
    double momentum_global_relative=0.0;
    std::array<std::size_t,6> worst_cell{{}};
    double continuity_l1=0.0;       // equals the native conservative L1 metric
    double merit=0.0;               // diagnostic only: no viscosity-closure gate
    std::uint64_t fingerprint=0;
};

// Platform-independent byte order for the digest: use the IEEE-754 binary64
// bit pattern, then visit bytes least-significant first. Digest is a
// reproducibility marker, NOT a cryptographic integrity certificate.
inline std::uint64_t fingerprint_native_six_v50(const std::vector<double>& values) {
    static_assert(sizeof(double)==sizeof(std::uint64_t),
                  "V50 expects 64-bit double");
    std::uint64_t hash=UINT64_C(14695981039346656037);
    for (double value:values) {
        if (!std::isfinite(value))
            throw std::invalid_argument("V50 nonfinite native vector entry");
        std::uint64_t bits=0;
        std::memcpy(&bits,&value,sizeof(bits));
        for (unsigned j=0;j<8;++j) {
            hash^=(bits>>(8*j)) & UINT64_C(0xff);
            hash*=UINT64_C(1099511628211);
        }
    }
    return hash;
}

inline NativeSixFieldResidualV50 make_native_six_residual_v50(
    const std::array<std::vector<double>,3>& signed_momentum,
    const std::array<double,3>& momentum_rhs_scale,
    const std::vector<double>& signed_continuity,
    double face_flux_l1,
    const std::vector<double>& signed_sst,
    const std::vector<double>& sst_row_norm,
    const std::vector<double>& sst_field_scale,
    bool closed)
{
    const std::size_t n=signed_continuity.size();
    if (!closed || n==0 || n>std::numeric_limits<std::size_t>::max()/6 ||
        signed_sst.size()!=2*n || sst_row_norm.size()!=2*n ||
        sst_field_scale.size()!=2*n ||
        !std::isfinite(face_flux_l1) || !(face_flux_l1>0.0))
        throw std::invalid_argument("V50 unclosed flux or invalid six-field shape");
    for (std::size_t d=0;d<3;++d)
        if (signed_momentum[d].size()!=n ||
            !std::isfinite(momentum_rhs_scale[d]) ||
            !(momentum_rhs_scale[d]>0.0))
            throw std::invalid_argument("V50 invalid momentum block or RHS scale");

    double global_rhs=0.0;
    for (double component:momentum_rhs_scale)
        global_rhs=std::max(global_rhs,component);
    if (!(global_rhs>0.0))
        throw std::invalid_argument("V50 empty global momentum RHS scale");
    NativeSixFieldResidualV50 out;
    out.cells=n;
    for (std::size_t d=0;d<3;++d)
        out.active_momentum[d]=momentum_rhs_scale[d]>=1e-8*global_rhs;
    out.normalized.resize(6*n);
    long double continuity_sum=0.0L;
    for (std::size_t cell=0;cell<n;++cell) {
        for (std::size_t d=0;d<6;++d) {
            const double signed_raw = d<3 ? signed_momentum[d][cell]
                : d==3 ? signed_continuity[cell]
                : signed_sst[(d-4)*n+cell];
            double divisor = d<3
                ? std::max(momentum_rhs_scale[d],1e-8*global_rhs)
                : d==3 ? face_flux_l1
                : sst_row_norm[(d-4)*n+cell]*sst_field_scale[(d-4)*n+cell];
            if (!std::isfinite(signed_raw) || !std::isfinite(divisor) ||
                !(divisor>0.0))
                throw std::invalid_argument("V50 nonfinite or invalid FV row normalization");
            const double value=signed_raw/divisor;
            if (!std::isfinite(value))
                throw std::domain_error("V50 native residual normalization overflow");
            out.normalized[d*n+cell]=value;
            const double abs_value=std::abs(value);
            if (abs_value>out.linf[d]) {
                out.linf[d]=abs_value;
                out.worst_cell[d]=cell;
            }
            if (d==3) continuity_sum+=abs_value;
        }
    }
    for (std::size_t d=0;d<3;++d)
        if (out.active_momentum[d])
            out.momentum_active_merit=std::max(out.momentum_active_merit,out.linf[d]);
    double max_raw_momentum=0.0;
    for (const auto& block:signed_momentum)
        for (double raw:block)
            max_raw_momentum=std::max(max_raw_momentum,std::abs(raw));
    out.momentum_global_relative=max_raw_momentum/global_rhs;
    if (!std::isfinite(out.momentum_global_relative))
        throw std::domain_error("V50 global momentum normalization overflow");
    out.continuity_l1=static_cast<double>(continuity_sum);
    if (!std::isfinite(out.continuity_l1))
        throw std::domain_error("V50 continuity L1 overflow");
    out.merit=std::max({out.momentum_active_merit,
                        out.continuity_l1,out.linf[4],out.linf[5]});
    out.fingerprint=fingerprint_native_six_v50(out.normalized);
    return out;
}

// Additional OPTIONAL final-stop veto. This can only prevent a premature
// certificate; it must never replace the original production convergence
// test, and never accept an otherwise rejected state. Native flux and the
// carrier's authoritative flux must agree at the SAME tolerance.
inline bool native_six_final_stop_eligible_v50(
    const NativeSixFieldResidualV50& sample,
    double carrier_native_flux_gap_l1,
    double tolerance)
{
    if (!(tolerance>0.0) || !std::isfinite(tolerance))
        throw std::invalid_argument("V50 invalid native stop tolerance");
    if (sample.cells==0 || sample.normalized.size()!=6*sample.cells ||
        !std::isfinite(carrier_native_flux_gap_l1) ||
        carrier_native_flux_gap_l1<0.0 ||
        !std::isfinite(sample.momentum_global_relative) ||
        !std::isfinite(sample.merit) ||
        !std::isfinite(sample.continuity_l1))
        return false;
    return sample.merit<=tolerance &&
           sample.momentum_global_relative<=tolerance &&
           carrier_native_flux_gap_l1<=tolerance;
}

} // namespace cfdx::physics
