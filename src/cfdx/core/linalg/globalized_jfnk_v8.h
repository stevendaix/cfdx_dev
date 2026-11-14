#pragma once

#include "cfdx/core/linalg/gmres_solver.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::core {

// Model-independent globalization engine for a *dimensionless*, independently
// evaluated nonlinear operator F(x). A caller owns the physical F, block
// scales, positivity/BC constraints and the Schur/PCD/AMG preconditioner.
// Preconditioner::apply is an APPROXIMATE direction only; only F(x_trial)
// can accept a new state. This class does not silently substitute a Schur
// model, alter FV equations, or claim to assemble a production U-p-k-omega F.
struct GlobalizedJfnkControlsV8 {
    std::size_t max_newton = 12;
    std::size_t max_krylov = 160;
    int restart = 30;
    // This is the inexact-Newton forcing tolerance for a finite-difference
    // directional Jacobian (NOT the physical nonlinear qualification gate).
    // An overly strict 1e-8 request on Jv formed from sqrt(eps) differences
    // can demand accuracy below directional differentiation roundoff.
    double linear_tolerance = 1e-5;
    double nonlinear_tolerance = 1e-8;
    double finite_difference_step = 1.490116119384765625e-8;
    double armijo = 1e-4;
    double minimum_lambda = 1.0 / 1024.0;
};

struct GlobalizedJfnkResultV8 {
    bool converged = false;
    std::size_t newton_steps = 0;
    std::size_t krylov_steps = 0;
    std::size_t rejected_steps = 0;
    double initial_merit = std::numeric_limits<double>::infinity();
    double final_merit = std::numeric_limits<double>::infinity();
    std::string reason;
};

inline GlobalizedJfnkResultV8 globalized_jfnk_v8(
    Vector& state,
    const std::function<bool(const Vector&, Vector&)>& evaluate,
    Preconditioner* user_preconditioner,
    const GlobalizedJfnkControlsV8& controls = {})
{
    if (state.empty() || !state.is_valid() || !evaluate ||
        controls.max_newton == 0 || controls.max_krylov == 0 ||
        controls.restart <= 0 || !(controls.linear_tolerance > 0.) ||
        !(controls.nonlinear_tolerance > 0.) ||
        !(controls.finite_difference_step > 0.) ||
        !(controls.armijo > 0. && controls.armijo < 1.) ||
        !(controls.minimum_lambda > 0. && controls.minimum_lambda <= 1.) ||
        !std::isfinite(controls.linear_tolerance) ||
        !std::isfinite(controls.nonlinear_tolerance) ||
        !std::isfinite(controls.finite_difference_step))
        throw std::invalid_argument("JFNK V8: invalid state or tolerances");
    const std::size_t n = state.size();
    GlobalizedJfnkResultV8 result;
    Vector f(n, 0.);
    auto safe_evaluate = [&](const Vector& x, Vector& r) {
        if (x.size() != n || !x.is_valid()) return false;
        try { return evaluate(x, r) && r.size() == n && r.is_valid(); }
        catch (const std::exception&) { return false; }
    };
    if (!safe_evaluate(state, f)) {
        result.reason = "invalid initial physical residual";
        return result;
    }
    result.initial_merit = f.norm_inf();
    result.final_merit = result.initial_merit;
    for (std::size_t step = 0; step < controls.max_newton; ++step) {
        const double before = f.norm_inf();
        if (before <= controls.nonlinear_tolerance) {
            result.converged = true;
            result.reason = "independent nonlinear residual satisfied";
            return result;
        }
        LinearOperator jacobian;
        jacobian.size = n;
        jacobian.apply = [&](const Vector& v, Vector& jv) {
            const double vnorm = v.norm_inf();
            if (jv.size() != n) jv.resize(n);
            if (vnorm == 0.) { jv.fill(0.); return; }
            const double h = controls.finite_difference_step *
                std::max(1., state.norm_inf()) / vnorm;
            if (!(h > 0.) || !std::isfinite(h))
                throw std::runtime_error("JFNK V8: invalid direction step");
            Vector trial = state;
            Vector sampled(n, 0.);
            for (std::size_t i = 0; i < n; ++i) trial(i) += h * v(i);
            double signed_h = h;
            if (!safe_evaluate(trial, sampled)) {
                trial = state;
                for (std::size_t i = 0; i < n; ++i) trial(i) -= h * v(i);
                signed_h = -h;
                if (!safe_evaluate(trial, sampled))
                    throw std::runtime_error("JFNK V8: direction violates both admissible states");
            }
            for (std::size_t i = 0; i < n; ++i)
                jv(i) = (sampled(i) - f(i)) / signed_h;
            if (!jv.is_valid())
                throw std::runtime_error("JFNK V8: non-finite directional Jacobian");
        };
        Vector rhs(n, 0.), delta(n, 0.);
        for (std::size_t i = 0; i < n; ++i) rhs(i) = -f(i);
        SolverResult linear;
        try {
            linear = solve_fgmres(jacobian, rhs, delta,
                controls.restart, controls.max_krylov,
                controls.linear_tolerance, user_preconditioner);
        } catch (const std::exception& e) {
            result.reason = std::string("JFNK V8: linear operator failed: ") + e.what();
            return result;
        }
        result.krylov_steps += linear.iterations;
        if (linear.status != SolverStatus::CONVERGED || !delta.is_valid()) {
            result.reason = "JFNK V8: true linear residual did not converge";
            return result;
        }
        bool accepted = false;
        for (double lambda = 1.; lambda >= controls.minimum_lambda; lambda *= .5) {
            Vector trial = state;
            for (std::size_t i = 0; i < n; ++i) trial(i) += lambda * delta(i);
            Vector trial_f(n, 0.);
            if (!safe_evaluate(trial, trial_f) ||
                trial_f.norm_inf() > (1. - controls.armijo * lambda) * before) {
                ++result.rejected_steps;
                continue;
            }
            state = std::move(trial);
            f = std::move(trial_f);
            result.final_merit = f.norm_inf();
            ++result.newton_steps;
            accepted = true;
            break;
        }
        if (!accepted) {
            result.reason = "JFNK V8: no globally decreasing admissible step";
            return result;
        }
    }
    result.converged = result.final_merit <= controls.nonlinear_tolerance;
    result.reason = result.converged
        ? "independent nonlinear residual satisfied"
        : "JFNK V8: nonlinear iteration budget exhausted";
    return result;
}

} // namespace cfdx::core
