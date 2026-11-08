#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::physics {

// Experimental matrix-free Newton-Krylov kernel for the COMPLETE discrete
// residual with state X = [Ux, Uy, Uz, p, face_mass_flux]. The caller MUST
// reassemble momentum, continuity, and the Rhie-Chow flux CONSTRAINT for
// every evaluation, including the perturbed (+/-) states. A scalar TVD-only
// residual does not satisfy this contract. No production solver enables this
// kernel implicitly, and it does not install PETSc.
struct CoupledNewtonKrylovControls {
    std::size_t krylov_dimension = 20;
    double finite_difference_step = 1e-6;
    double linear_relative_tolerance = 2e-2;
    double sufficient_decrease = 1e-4;
    std::size_t maximum_backtracks = 10;
};

struct CoupledNewtonKrylovResult {
    bool accepted = false;
    std::vector<double> state;
    double physical_merit_before = std::numeric_limits<double>::infinity();
    double physical_merit_after = std::numeric_limits<double>::infinity();
    double linear_true_relative_residual = std::numeric_limits<double>::infinity();
    double accepted_step = 0.0;
    std::size_t residual_evaluations = 0;
    std::string reason;
};

inline void validate_coupled_newton_controls(const CoupledNewtonKrylovControls& c) {
    if (c.krylov_dimension == 0 || c.krylov_dimension > 80 ||
        c.maximum_backtracks > 32 ||
        !std::isfinite(c.finite_difference_step) ||
        c.finite_difference_step < 1e-10 || c.finite_difference_step > 1e-2 ||
        !std::isfinite(c.linear_relative_tolerance) ||
        c.linear_relative_tolerance <= 0.0 || c.linear_relative_tolerance >= 1.0 ||
        !std::isfinite(c.sufficient_decrease) ||
        c.sufficient_decrease <= 0.0 || c.sufficient_decrease >= 0.5)
        throw std::invalid_argument("coupled JFNK: invalid numerical controls");
}

using CoupledResidual = std::function<std::vector<double>(const std::vector<double>&)>;
using CoupledStateAdmissible = std::function<bool(const std::vector<double>&)>;

inline double coupled_newton_norm2(const std::vector<double>& v) {
    double sum = 0.0;
    for (const double x : v) {
        if (!std::isfinite(x))
            throw std::runtime_error("coupled JFNK: non-finite physical residual or state");
        sum += x * x;
    }
    if (!std::isfinite(sum))
        throw std::runtime_error("coupled JFNK: non-finite norm");
    return std::sqrt(sum);
}

inline double coupled_newton_fixed_merit(
    const std::vector<double>& v, const std::vector<double>& scales) {
    if (v.empty() || v.size() != scales.size())
        throw std::invalid_argument("coupled JFNK: inconsistent physical scales");
    double merit = 0.0;
    for (std::size_t k=0; k<v.size(); ++k) {
        if (!(scales[k] > 0.0) || !std::isfinite(scales[k]) ||
            !std::isfinite(v[k]))
            throw std::invalid_argument("coupled JFNK: nonfinite residual/scale");
        merit = std::max(merit, std::abs(v[k]) / scales[k]);
    }
    return merit;
}

// Explicitly validate the gauge in the callback. For a pure-Neumann pressure
// problem, one PRESSURE residual row must be p[gauge]-p_reference, replacing
// the redundant continuity row; the FULL continuity must still be audited
// separately before any physical convergence certificate.
inline std::vector<double> coupled_newton_directional_jvp(
    const CoupledResidual& residual, const std::vector<double>& state,
    const std::vector<double>& direction, double step = 1e-6) {
    if (!residual || state.empty() || state.size() != direction.size() ||
        !std::isfinite(step) || step <= 0.0)
        throw std::invalid_argument("coupled JFNK: invalid residual or direction");
    double max_state = 0.0, max_direction = 0.0;
    for (std::size_t i=0; i<state.size(); ++i) {
        if (!std::isfinite(state[i]) || !std::isfinite(direction[i]))
            throw std::invalid_argument("coupled JFNK: non-finite direction");
        max_state = std::max(max_state, std::abs(state[i]));
        max_direction = std::max(max_direction, std::abs(direction[i]));
    }
    if (!(max_direction > 0.0))
        throw std::invalid_argument("coupled JFNK: zero direction");
    const double h = step * (1.0 + max_state) / max_direction;
    auto plus = state, minus = state;
    for (std::size_t i=0; i<state.size(); ++i) {
        plus[i] += h * direction[i];
        minus[i] -= h * direction[i];
    }
    const auto rp = residual(plus);
    const auto rm = residual(minus);
    if (rp.size() != state.size() || rm.size() != state.size())
        throw std::invalid_argument("coupled JFNK: residual dimension changed");
    std::vector<double> jv(state.size());
    for (std::size_t i=0; i<state.size(); ++i)
        jv[i] = (rp[i] - rm[i]) / (2.0 * h);
    (void)coupled_newton_norm2(jv);
    return jv;
}

// A single GMRES (Arnoldi + Givens) linearization at a fixed physical state.
// The independently evaluated true Jv residual, NOT the Hessenberg estimate,
// is returned. Krylov budgets are explicit and failures cannot be promoted.
inline bool coupled_newton_gmres(
    const std::function<std::vector<double>(const std::vector<double>&)>& apply,
    const std::vector<double>& rhs, std::size_t maximum_dimension,
    double relative_tolerance, std::vector<double>& step,
    double& true_relative_residual) {
    const auto n = rhs.size();
    const double beta = coupled_newton_norm2(rhs);
    step.assign(n, 0.0);
    if (!(beta > 0.0)) {true_relative_residual=0.0; return true;}
    const std::size_t m = std::min(maximum_dimension,n);
    std::vector<std::vector<double>> basis;
    basis.reserve(m+1);
    basis.push_back(rhs);
    for (double& v : basis.back()) v /= beta;
    std::vector<std::vector<double>> h(m+1,std::vector<double>(m,0.0));
    std::vector<double> cs(m),sn(m),g(m+1,0.0);
    g[0] = beta;
    std::size_t used=0;
    for (std::size_t j=0;j<m;++j) {
        auto w=apply(basis[j]);
        if (w.size()!=n) throw std::runtime_error("coupled JFNK: malformed operator action");
        // Twice-reorthogonalized MGS for nearly dependent momentum/pressure rows.
        for (int pass=0;pass<2;++pass)
            for (std::size_t i=0;i<=j;++i) {
                double dot=0.0;
                for (std::size_t k=0;k<n;++k) dot+=w[k]*basis[i][k];
                h[i][j]+=dot;
                for (std::size_t k=0;k<n;++k) w[k]-=dot*basis[i][k];
            }
        h[j+1][j]=coupled_newton_norm2(w);
        const bool breakdown=h[j+1][j]<1e-13*std::max(beta,1.0);
        if (!breakdown) {
            for (double& value:w) value/=h[j+1][j];
            basis.push_back(std::move(w));
        }
        for (std::size_t i=0;i<j;++i) {
            const double a=h[i][j], b=h[i+1][j];
            h[i][j]=cs[i]*a+sn[i]*b;
            h[i+1][j]=-sn[i]*a+cs[i]*b;
        }
        const double a=h[j][j],b=h[j+1][j];
        const double r=std::hypot(a,b);
        if (!(r>1e-18)) return false;
        cs[j]=a/r; sn[j]=b/r;
        h[j][j]=r; h[j+1][j]=0.0;
        g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
        used=j+1;
        if (std::abs(g[j+1])<=relative_tolerance*beta || breakdown) break;
    }
    if (!used) return false;
    std::vector<double> y(used,0.0);
    for (std::size_t rev=0;rev<used;++rev) {
        const std::size_t i=used-1-rev;
        double a=g[i];
        for (std::size_t j=i+1;j<used;++j) a-=h[i][j]*y[j];
        if (!(std::abs(h[i][i])>1e-18)) return false;
        y[i]=a/h[i][i];
    }
    for (std::size_t j=0;j<used;++j)
        for (std::size_t k=0;k<n;++k)
            step[k]+=y[j]*basis[j][k];
    const auto computed=apply(step);
    std::vector<double> difference(n);
    for (std::size_t k=0;k<n;++k) difference[k]=rhs[k]-computed[k];
    true_relative_residual=coupled_newton_norm2(difference)/beta;
    return true_relative_residual<=relative_tolerance;
}

inline CoupledNewtonKrylovResult coupled_newton_krylov_step(
    const CoupledResidual& full_residual,
    const std::vector<double>& accepted_state,
    const std::vector<double>& fixed_physical_scales,
    const CoupledNewtonKrylovControls& controls = {},
    const CoupledStateAdmissible& admissible = {}) {
    validate_coupled_newton_controls(controls);
    if (!full_residual || accepted_state.empty() ||
        accepted_state.size()!=fixed_physical_scales.size())
        throw std::invalid_argument("coupled JFNK: missing complete coupled operator");
    CoupledNewtonKrylovResult result;
    result.state=accepted_state; // transactional: never mutate caller's last accepted state
    const auto eval=[&](const std::vector<double>& x) {
        ++result.residual_evaluations;
        auto value=full_residual(x);
        if (value.size()!=accepted_state.size())
            throw std::invalid_argument("coupled JFNK: incomplete coupled residual");
        (void)coupled_newton_norm2(value);
        return value;
    };
    const auto base=eval(accepted_state);
    result.physical_merit_before=coupled_newton_fixed_merit(base,fixed_physical_scales);
    if (result.physical_merit_before==0.0) {
        result.accepted=true;
        result.physical_merit_after=0.0;
        result.linear_true_relative_residual=0.0;
        result.reason="already at the full coupled fixed point";
        return result;
    }
    std::vector<double> rhs=base;
    for (double& value:rhs) value=-value;
    const auto apply=[&](const std::vector<double>& direction) {
        return coupled_newton_directional_jvp(eval,accepted_state,direction,
                                              controls.finite_difference_step);
    };
    std::vector<double> direction;
    if (!coupled_newton_gmres(apply,rhs,controls.krylov_dimension,
                             controls.linear_relative_tolerance,direction,
                             result.linear_true_relative_residual)) {
        result.reason="Krylov true directional residual not satisfied";
        return result;
    }
    for (std::size_t k=0;k<=controls.maximum_backtracks;++k) {
        const double alpha=std::ldexp(1.0,-static_cast<int>(k));
        auto candidate=accepted_state;
        for (std::size_t i=0;i<candidate.size();++i)
            candidate[i]+=alpha*direction[i];
        if (admissible && !admissible(candidate)) continue;
        const auto physical=eval(candidate);
        const double merit=coupled_newton_fixed_merit(physical,fixed_physical_scales);
        if (merit<=result.physical_merit_before*(1.0-
                                                controls.sufficient_decrease*alpha)) {
            result.state=std::move(candidate);
            result.accepted=true;
            result.accepted_step=alpha;
            result.physical_merit_after=merit;
            result.reason="full coupled residual decreased";
            return result;
        }
    }
    result.reason="all coupled Newton candidates rejected; accepted state unchanged";
    return result;
}

} // namespace cfdx::physics
