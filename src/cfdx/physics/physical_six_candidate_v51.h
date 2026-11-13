#pragma once
// V51 experimental candidate-state nonlinear FUNCTION contract.
// Every evaluation receives ONLY the full six-field candidate (Ux,Uy,Uz,p*,k,w)
// plus immutable mesh, material and boundary data held by the caller.
// The evaluator must rebuild strain, viscosity, wall constraints and production,
// momentum, Rhie-Chow flux and both SST equations. It MUST NOT read accepted
// flux, accepted rAU, SST derived fields, or prior nonlinear trials.
// Residual layout: six cell-major blocks [Ux,Uy,Uz,div(phi),k,w].
// A raw signed residual is used for J*v; a separate scaled diagnostic is
// allowed but its state-dependent scales must never enter the J*v difference.
// This is an OFF-by-default forensic adapter, NOT a Newton solver.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

struct PhysicalCandidateV51 {
    std::size_t cells=0;
    std::vector<double> fields; // [Ux,Uy,Uz,p,k,w], block-major
};

// V57: an optional, candidate-owned staging tape. These are NOT lagged
// coefficients and are never used to alter the primal or its acceptance.
// A trace is absent unless requested explicitly by the diagnostic caller.
struct CandidatePhysicalStagesV57 {
    std::vector<double> inverse_momentum_diagonal; // [Ux,Uy,Uz] 3*N
    std::vector<double> effective_dynamic_viscosity; // N
    std::vector<double> strain_rate; // N
    std::vector<double> sst_blend_f2; // N
    std::vector<double> wall_momentum_face_mu; // Nfaces, -1 for non-wall
    std::vector<double> wall_production; // N, -1 when no production override
    std::vector<double> momentum_convection; // 3*N, exact primal from V56 audit
};

struct PhysicalResidualV51 {
    std::vector<double> signed_raw; // [momentum x3, divergence, SST x2]
    std::vector<double> normalized; // separate, NOT differenced by J*v
    std::vector<double> face_flux;
    std::vector<double> eddy_viscosity;
    double flux_fixed_point_defect=std::numeric_limits<double>::infinity();
    double merit=std::numeric_limits<double>::infinity();
    bool flux_closed=false;
    std::optional<CandidatePhysicalStagesV57> stages_v57;
};

inline void validate_physical_candidate_v51(
    const PhysicalCandidateV51& x, double k_min=0.0, double omega_min=0.0) {
    const std::size_t n=x.cells;
    if (!n || n>std::numeric_limits<std::size_t>::max()/6 ||
        x.fields.size()!=6*n || !std::isfinite(k_min) ||
        !std::isfinite(omega_min) || k_min<0.0 || omega_min<0.0)
        throw std::invalid_argument("V51 six-field candidate shape or bounds invalid");
    for (std::size_t j=0;j<6;++j)
        for (std::size_t i=0;i<n;++i) {
            const double value=x.fields[j*n+i];
            if (!std::isfinite(value) ||
                (j==4 && value<k_min) || (j==5 && !(value>omega_min)))
                throw std::invalid_argument("V51 nonfinite or inadmissible physical candidate");
        }
}

inline void validate_physical_result_v51(const PhysicalCandidateV51& x,
                                        const PhysicalResidualV51& r,
                                        double flux_tolerance=1e-10) {
    const std::size_t n=x.cells;
    if (!r.flux_closed || !std::isfinite(r.flux_fixed_point_defect) ||
        !(r.flux_fixed_point_defect>=0.0) ||
        r.flux_fixed_point_defect>flux_tolerance ||
        r.signed_raw.size()!=6*n || r.normalized.size()!=6*n ||
        r.face_flux.empty() || r.eddy_viscosity.size()!=n ||
        !std::isfinite(r.merit) || r.merit<0.0)
        throw std::runtime_error("V51 incomplete candidate FV operator or unclosed flux");
    for (double v:r.signed_raw)
        if (!std::isfinite(v)) throw std::runtime_error("V51 nonfinite raw residual");
    for (double v:r.normalized)
        if (!std::isfinite(v)) throw std::runtime_error("V51 nonfinite scaled residual");
    for (double v:r.face_flux)
        if (!std::isfinite(v)) throw std::runtime_error("V51 nonfinite face flux");
    for (double v:r.eddy_viscosity)
        if (!(v>=0.0) || !std::isfinite(v))
            throw std::runtime_error("V51 invalid candidate turbulent viscosity");
    if (r.stages_v57) {
        const auto& t=*r.stages_v57;
        if (t.inverse_momentum_diagonal.size()!=3*n ||
            t.effective_dynamic_viscosity.size()!=n ||
            t.strain_rate.size()!=n || t.sst_blend_f2.size()!=n ||
            t.wall_production.size()!=n ||
            t.momentum_convection.size()!=3*n ||
            (!t.wall_momentum_face_mu.empty() &&
             t.wall_momentum_face_mu.size()!=r.face_flux.size()))
            throw std::runtime_error("V57 incomplete candidate-stage trace");
        for (double a:t.inverse_momentum_diagonal)
            if (!(a>0.) || !std::isfinite(a))
                throw std::runtime_error("V57 invalid native rAU trace");
        for (double a:t.effective_dynamic_viscosity)
            if (!(a>0.) || !std::isfinite(a))
                throw std::runtime_error("V57 invalid dynamic-viscosity trace");
        for (double a:t.strain_rate)
            if (!(a>=0.) || !std::isfinite(a))
                throw std::runtime_error("V57 invalid strain-rate trace");
        for (double a:t.sst_blend_f2)
            if (!std::isfinite(a))
                throw std::runtime_error("V57 invalid F2 trace");
        for (double a:t.wall_momentum_face_mu)
            if (!std::isfinite(a) || (a>=0. && !(a>0.)))
                throw std::runtime_error("V57 invalid wall viscosity trace");
        for (double a:t.wall_production)
            if (!std::isfinite(a))
                throw std::runtime_error("V57 nonfinite wall production trace");
        for (double a:t.momentum_convection)
            if (!std::isfinite(a))
                throw std::runtime_error("V57 nonfinite momentum convection trace");
    }
}

template<class Evaluate>
PhysicalResidualV51 evaluate_physical_candidate_v51(
    const PhysicalCandidateV51& x, Evaluate&& callback,
    double k_min=0.0, double omega_min=0.0) {
    validate_physical_candidate_v51(x,k_min,omega_min);
    const auto result=callback(x);
    validate_physical_result_v51(x,result);
    return result;
}

// Guarantees a rejected/perturbed state is not published by this API.
// Directional difference of the SIGNED *raw* physical FV residual, not of
// an independently renormalized residual or a frozen Picard matrix.
template<class Evaluate>
std::vector<double> physical_directional_jv_v51(
    const PhysicalCandidateV51& base, const std::vector<double>& direction,
    Evaluate&& callback, double h, double k_min=0.0,
    double omega_min=0.0) {
    validate_physical_candidate_v51(base,k_min,omega_min);
    if (direction.size()!=base.fields.size() ||
        !(h>0.0) || !std::isfinite(h))
        throw std::invalid_argument("V51 directional Jv shape or step invalid");
    PhysicalCandidateV51 shifted=base;
    bool nonzero=false;
    for (std::size_t i=0;i<direction.size();++i) {
        if (!std::isfinite(direction[i]))
            throw std::invalid_argument("V51 nonfinite Jv direction");
        nonzero|=direction[i]!=0.0;
        shifted.fields[i]+=h*direction[i];
    }
    if (!nonzero)
        throw std::invalid_argument("V51 zero Jv direction");
    validate_physical_candidate_v51(shifted,k_min,omega_min);
    const auto f0=evaluate_physical_candidate_v51(base,callback,k_min,omega_min);
    const auto f1=evaluate_physical_candidate_v51(shifted,callback,k_min,omega_min);
    std::vector<double> jv(f0.signed_raw.size());
    for (std::size_t i=0;i<jv.size();++i) {
        jv[i]=(f1.signed_raw[i]-f0.signed_raw[i])/h;
        if (!std::isfinite(jv[i]))
            throw std::runtime_error("V51 nonfinite directional derivative");
    }
    return jv;
}

inline void verify_candidate_replay_v51(const PhysicalResidualV51& a,
                                        const PhysicalResidualV51& b) {
    if (a.signed_raw!=b.signed_raw || a.normalized!=b.normalized ||
        a.face_flux!=b.face_flux || a.eddy_viscosity!=b.eddy_viscosity ||
        a.flux_fixed_point_defect!=b.flux_fixed_point_defect ||
        a.flux_closed!=b.flux_closed || a.merit!=b.merit ||
        a.stages_v57.has_value()!=b.stages_v57.has_value() ||
        (a.stages_v57 && (a.stages_v57->inverse_momentum_diagonal!=b.stages_v57->inverse_momentum_diagonal ||
        a.stages_v57->effective_dynamic_viscosity!=b.stages_v57->effective_dynamic_viscosity ||
        a.stages_v57->strain_rate!=b.stages_v57->strain_rate ||
        a.stages_v57->sst_blend_f2!=b.stages_v57->sst_blend_f2 ||
        a.stages_v57->wall_momentum_face_mu!=b.stages_v57->wall_momentum_face_mu ||
        a.stages_v57->wall_production!=b.stages_v57->wall_production ||
        a.stages_v57->momentum_convection!=b.stages_v57->momentum_convection)))
        throw std::runtime_error("V51 same-candidate physical operator is history-dependent");
}

} // namespace cfdx::physics
