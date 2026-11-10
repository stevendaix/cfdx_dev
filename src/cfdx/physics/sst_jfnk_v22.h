#pragma once

#include "cfdx/physics/coupled_newton_krylov.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// An optional SST-only Newton candidate. This reuses CFDX's audited GMRES
// kernel; PETSc is NOT installed or emulated. Log unknowns keep every
// differentiated trial strictly positive without clipping the SST equations.
// Caller returns the actual production FV residual [R_k,R_omega] on each
// trial, with current F1/F2, G, nu_t, cross diffusion and owner constraints.
struct SSTJFNKV22Candidate {
    bool available=false;
    std::vector<double> fields; // same [k,omega] ordering as the input
    std::size_t function_evaluations=0;
    double frozen_merit_before=std::numeric_limits<double>::infinity();
    double frozen_merit_after=std::numeric_limits<double>::infinity();
    double linear_true_relative=std::numeric_limits<double>::infinity();
    double step_length=0.0;
    std::string reason;
};

using SSTPhysicalResidualV22 =
    std::function<std::vector<double>(const std::vector<double>&)>;

inline SSTJFNKV22Candidate sst_jfnk_propose_v22(
    const std::vector<double>& physical_state,
    const std::vector<double>& frozen_residual_scales,
    std::size_t n_cells,
    double k_floor,double omega_floor,
    const SSTPhysicalResidualV22& physical_residual,
    const CoupledNewtonKrylovControls& settings = {},
    const std::vector<double>& frozen_left_diagonal = {})
{
    if (!n_cells || physical_state.size()!=2*n_cells ||
        frozen_residual_scales.size()!=2*n_cells || !physical_residual ||
        !(k_floor>0.) || !(omega_floor>0.) ||
        !std::isfinite(k_floor) || !std::isfinite(omega_floor))
        throw std::invalid_argument("V22 SST JFNK: incomplete coupled residual contract");
    if (!frozen_left_diagonal.empty() &&
        frozen_left_diagonal.size()!=2*n_cells)
        throw std::invalid_argument("V22 SST JFNK: incomplete Jacobi preconditioner");
    for (const double d:frozen_left_diagonal)
        if (!(d>0.) || !std::isfinite(d))
            throw std::invalid_argument("V22 SST JFNK: invalid Jacobi preconditioner");
    SSTJFNKV22Candidate out;
    const double k_ref=*std::max_element(physical_state.begin(),
                                        physical_state.begin()+static_cast<std::ptrdiff_t>(n_cells));
    const double w_ref=*std::max_element(physical_state.begin()+static_cast<std::ptrdiff_t>(n_cells),
                                        physical_state.end());
    if (!(k_ref>0.) || !(w_ref>0.) || !std::isfinite(k_ref) || !std::isfinite(w_ref))
        throw std::invalid_argument("V22 SST JFNK: invalid field reference scale");
    std::vector<double> log_state(2*n_cells);
    for (std::size_t i=0;i<2*n_cells;++i) {
        const double floor=i<n_cells?k_floor:omega_floor;
        const double ref=i<n_cells?k_ref:w_ref;
        if (!(physical_state[i]>=floor) || !std::isfinite(physical_state[i]) ||
            !(frozen_residual_scales[i]>0.) || !std::isfinite(frozen_residual_scales[i]))
            throw std::invalid_argument("V22 SST JFNK: nonpositive field or residual scale");
        // Central matrix-free differences must stay inside the admissible
        // domain. An active positivity floor is a variational inequality,
        // not a differentiable interior Newton equation; let Picard handle it.
        if (physical_state[i]<=floor*(1.+1e-4)) {
            out.reason="near active positivity floor; Newton direction skipped";
            return out;
        }
        log_state[i]=std::log(physical_state[i]/ref);
    }
    const auto decode=[&](const std::vector<double>& z) {
        if (z.size()!=2*n_cells)
            throw std::invalid_argument("V22 SST JFNK: state dimension changed");
        std::vector<double> state(z.size());
        for (std::size_t i=0;i<z.size();++i) {
            const double ref=i<n_cells?k_ref:w_ref;
            const double floor=i<n_cells?k_floor:omega_floor;
            if (!std::isfinite(z[i]) || z[i]<-700. || z[i]>700.)
                throw std::domain_error("V22 SST JFNK: exponential overflow");
            const double value=ref*std::exp(z[i]);
            if (!std::isfinite(value) || value<floor)
                throw std::domain_error("V22 SST JFNK: positivity violated");
            state[i]=value;
        }
        return state;
    };
    const CoupledResidual fun=[&](const std::vector<double>& x) {
        auto r=physical_residual(decode(x));
        if (r.size()!=2*n_cells)
            throw std::invalid_argument("V22 SST JFNK: residual is incomplete");
        for (std::size_t i=0;i<r.size();++i)
            r[i]/=frozen_residual_scales[i]*(
                frozen_left_diagonal.empty()?1.:frozen_left_diagonal[i]);
        return r;
    };
    const CoupledStateAdmissible admissible=[&](const auto& z) {
        try { (void)decode(z);return true; }
        catch (const std::exception&) {return false;}
    };
    try {
        // The GMRES true residual gate remains explicit. A small Krylov
        // residual is NEVER used as a nonlinear or BFS convergence certificate.
        const auto r=coupled_newton_krylov_step(fun,log_state,
            std::vector<double>(2*n_cells,1.0),settings,admissible);
        out.function_evaluations=r.residual_evaluations;
        out.frozen_merit_before=r.physical_merit_before;
        out.frozen_merit_after=r.physical_merit_after;
        out.linear_true_relative=r.linear_true_relative_residual;
        out.step_length=r.accepted_step;
        out.reason=r.reason;
        if (r.accepted) {
            out.fields=decode(r.state);
            out.available=true;
        }
    } catch (const std::exception& e) {
        // Opt-in preconditioner candidate fails closed, not an accepted step;
        // the caller decides whether a separately certified Picard candidate
        // is possible. No caller-owned physical field is mutated.
        out.reason=std::string("JFNK rejected: ")+e.what();
    }
    return out;
}

} // namespace cfdx::physics
