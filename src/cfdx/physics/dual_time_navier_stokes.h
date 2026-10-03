#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

// Production physical-time transaction around the existing incompressible FVM
// solver. The steady nonlinear loop becomes the pseudo-time iteration by
// adding the BE/BDF2 temporal mass term to the production momentum equations.
// Pressure-velocity coupling, flux reconstruction, conservation diagnostics and
// the selected linear solver remain the production implementations.
struct DualTimeNavierStokesControls {
    double dt_min = 1.0e-6;
    double dt_max = 1.0e-1;
    double dt_growth = 1.5;
    double dt_shrink = 0.5;
    double temporal_absolute_tolerance = 1.0e-6;
    double temporal_relative_tolerance = 1.0e-3;
    int max_retries = 8;
};

struct DualTimeNavierStokesCheckpoint {
    cfdx::core::Field<double, cfdx::core::Location::CELL> U;
    cfdx::core::Field<double, cfdx::core::Location::CELL> p;
    cfdx::core::Field<double, cfdx::core::Location::CELL> U_previous;
    bool history_valid = false;
    double physical_time = 0.0;
    double dt_previous = 0.0;
    std::size_t step = 0;
};

struct DualTimeNavierStokesStepReport {
    IncompressibleSolveResult high_order;
    double physical_time_before = 0.0;
    double physical_time_after = 0.0;
    double dt_requested = 0.0;
    double dt_accepted = 0.0;
    double temporal_error = 0.0;
    double dt_proposed = 0.0;
    int retries = 0;
};

class DualTimeNavierStokesLifecycle {
public:
    DualTimeNavierStokesLifecycle(
        cfdx::core::Field<double, cfdx::core::Location::CELL> U,
        cfdx::core::Field<double, cfdx::core::Location::CELL> p,
        double physical_time = 0.0)
        : U_(std::move(U)), p_(std::move(p)), physical_time_(physical_time)
    {
        validate_state();
        if (!std::isfinite(physical_time_))
            throw std::invalid_argument("dual-time physical time must be finite");
    }

    const auto& U() const noexcept { return U_; }
    const auto& p() const noexcept { return p_; }
    double physical_time() const noexcept { return physical_time_; }
    double dt_previous() const noexcept { return dt_previous_; }
    std::size_t step() const noexcept { return step_; }
    bool history_valid() const noexcept { return history_valid_; }

    DualTimeNavierStokesCheckpoint checkpoint() const {
        return {U_, p_, U_previous_, history_valid_, physical_time_, dt_previous_, step_};
    }

    void restore(const DualTimeNavierStokesCheckpoint& checkpoint) {
        U_ = checkpoint.U;
        p_ = checkpoint.p;
        U_previous_ = checkpoint.U_previous;
        history_valid_ = checkpoint.history_valid;
        physical_time_ = checkpoint.physical_time;
        dt_previous_ = checkpoint.dt_previous;
        step_ = checkpoint.step;
        validate_state();
        if (!std::isfinite(physical_time_))
            throw std::invalid_argument("dual-time checkpoint physical time is invalid");
        if (history_valid_ && (!(dt_previous_ > 0.0) || !std::isfinite(dt_previous_)))
            throw std::invalid_argument("dual-time checkpoint previous dt is invalid");
    }

    DualTimeNavierStokesStepReport advance(
        double dt_requested,
        const IncompressibleSolverControls& base_controls,
        const DualTimeNavierStokesControls& controls)
    {
        validate_controls(controls);
        if (!(dt_requested > 0.0) || !std::isfinite(dt_requested))
            throw std::invalid_argument("dual-time requested dt must be finite and positive");

        const double time_before = physical_time_;
        double dt = std::clamp(dt_requested, controls.dt_min, controls.dt_max);
        int retries = 0;

        for (;;) {
            const auto accepted_U = U_;
            const auto accepted_p = p_;
            const bool accepted_history = history_valid_;
            const double accepted_dt_previous = dt_previous_;

            // First physical step is necessarily BE. Thereafter BDF2 uses the
            // actual previous accepted dt, so variable-step coefficients are
            // applied directly in the production momentum matrix.
            const bool use_bdf2 = history_valid_;
            IncompressibleSolverControls high_controls = base_controls;
            high_controls.transient.enabled = true;
            high_controls.transient.scheme =
                use_bdf2 ? cfdx::core::TimeScheme::BDF2 : cfdx::core::TimeScheme::EULER_IMPLICIT;
            high_controls.transient.dt = dt;
            high_controls.transient.dt_previous = history_valid_ ? dt_previous_ : dt;
            high_controls.transient.previous = &accepted_U;
            high_controls.transient.previous_previous = history_valid_ ? &U_previous_ : &accepted_U;
            high_controls.transient.history_valid = history_valid_;

            IncompressibleSolverControls low_controls = high_controls;
            low_controls.transient.scheme = cfdx::core::TimeScheme::EULER_IMPLICIT;

            auto high_U = accepted_U;
            auto high_p = accepted_p;
            auto low_U = accepted_U;
            auto low_p = accepted_p;

            auto high_result = solve_steady_incompressible(
                mesh_, high_U, high_p, velocity_bcs_, pressure_bcs_,
                high_controls);
            auto low_result = solve_steady_incompressible(
                mesh_, low_U, low_p, velocity_bcs_, pressure_bcs_,
                low_controls);

            const double error = velocity_temporal_error(
                high_U, low_U, controls.temporal_absolute_tolerance,
                controls.temporal_relative_tolerance);
            const bool nonlinear_ok = high_result.converged;
            const bool physical_ok = fields_finite_and_admissible(high_U, high_p);
            const bool temporal_ok = error <= 1.0;

            if (nonlinear_ok && physical_ok && temporal_ok) {
                U_ = std::move(high_U);
                p_ = std::move(high_p);
                U_previous_ = accepted_U;
                history_valid_ = true;
                dt_previous_ = dt;
                physical_time_ = time_before + dt;
                ++step_;
                double factor = controls.dt_growth;
                if (error > 0.0 && std::isfinite(error))
                    factor = std::min(
                        controls.dt_growth,
                        std::max(1.0, 0.9 * std::pow(1.0 / error, 1.0 / 2.0)));
                const double dt_proposed = std::clamp(
                    dt * factor, controls.dt_min, controls.dt_max);
                return {std::move(high_result), time_before, physical_time_,
                        dt_requested, dt, error, dt_proposed, retries};
            }

            U_ = accepted_U;
            p_ = accepted_p;
            U_previous_ = accepted_U;
            history_valid_ = accepted_history;
            dt_previous_ = accepted_dt_previous;

            if (retries >= controls.max_retries)
                throw std::runtime_error(
                    "dual-time Navier-Stokes physical step rejected after bounded retries");

            dt = std::max(controls.dt_min, dt * controls.dt_shrink);
            ++retries;
        }
    }

    void configure_problem(
        const cfdx::core::Mesh& mesh,
        VelocityBoundaryConditions velocity_bcs,
        ScalarBoundaryConditions pressure_bcs)
    {
        mesh_ = mesh;
        velocity_bcs_ = std::move(velocity_bcs);
        pressure_bcs_ = std::move(pressure_bcs);
    }

private:
    cfdx::core::Field<double, cfdx::core::Location::CELL> U_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> p_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> U_previous_;
    bool history_valid_ = false;
    double physical_time_ = 0.0;
    double dt_previous_ = 0.0;
    std::size_t step_ = 0;

    cfdx::core::Mesh mesh_;
    VelocityBoundaryConditions velocity_bcs_;
    ScalarBoundaryConditions pressure_bcs_;

    void validate_state() const {
        if (U_.dimension() != 3 || p_.dimension() != 1 ||
            U_.size() != p_.size() || U_.size() == 0)
            throw std::invalid_argument("dual-time Navier-Stokes state shape is invalid");
    }

    static void validate_controls(const DualTimeNavierStokesControls& c) {
        if (!(c.dt_min > 0.0) || !(c.dt_max >= c.dt_min) ||
            !std::isfinite(c.dt_min) || !std::isfinite(c.dt_max) ||
            !(c.dt_growth > 1.0) || !std::isfinite(c.dt_growth) ||
            !(c.dt_shrink > 0.0 && c.dt_shrink < 1.0) ||
            !std::isfinite(c.dt_shrink) ||
            c.max_retries < 0 ||
            !(c.temporal_absolute_tolerance >= 0.0) ||
            !(c.temporal_relative_tolerance >= 0.0) ||
            (c.temporal_absolute_tolerance == 0.0 &&
             c.temporal_relative_tolerance == 0.0))
            throw std::invalid_argument("invalid dual-time Navier-Stokes controls");
    }

    static double velocity_temporal_error(
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& high,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& low,
        double atol, double rtol)
    {
        double maximum = 0.0;
        for (std::size_t d = 0; d < high.dimension(); ++d) {
            for (std::size_t c = 0; c < high.size(); ++c) {
                const double h = high.component_data(d)[c];
                const double l = low.component_data(d)[c];
                if (!std::isfinite(h) || !std::isfinite(l))
                    return std::numeric_limits<double>::infinity();
                maximum = std::max(
                    maximum,
                    std::abs(h - l) / (atol + rtol * std::max(std::abs(h), std::abs(l))));
            }
        }
        return maximum;
    }

    static bool fields_finite_and_admissible(
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& U,
        const cfdx::core::Field<double, cfdx::core::Location::CELL>& p)
    {
        for (std::size_t d = 0; d < U.dimension(); ++d)
            for (std::size_t c = 0; c < U.size(); ++c)
                if (!std::isfinite(U.component_data(d)[c])) return false;
        for (std::size_t c = 0; c < p.size(); ++c)
            if (!std::isfinite(p(c))) return false;
        return true;
    }
};

} // namespace cfdx::physics
