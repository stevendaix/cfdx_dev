#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/numerics/temporal.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

// Failure class used only for bounded nonlinear continuation. Permanent
// geometry/data/diagnostic errors remain std::runtime_error (or more specific
// standard exceptions) and must never be hidden by the retry controller.
class NonlinearRetryableFailure : public std::runtime_error {
public:
    explicit NonlinearRetryableFailure(const std::string& message)
        : std::runtime_error(message) {}
};

enum class TimeStepChangeReason {
    Initial,
    TargetCfl,
    MinDtBound,
    MaxDtBound,
    NonlinearFailureRollback,
    NoChange
};

inline const char* to_string(TimeStepChangeReason reason)
{
    switch (reason) {
    case TimeStepChangeReason::Initial: return "initial";
    case TimeStepChangeReason::TargetCfl: return "target_cfl";
    case TimeStepChangeReason::MinDtBound: return "min_dt_bound";
    case TimeStepChangeReason::MaxDtBound: return "max_dt_bound";
    case TimeStepChangeReason::NonlinearFailureRollback: return "nonlinear_failure_rollback";
    case TimeStepChangeReason::NoChange: return "no_change";
    }
    return "unknown";
}

struct CflDiagnostics {
    double global_cfl = 0.0;
    double local_cfl_min = 0.0;
    double local_cfl_max = 0.0;
    std::size_t limiting_cell = 0;
};

inline std::vector<double> compute_local_convective_cfl(
    const cfdx::core::Mesh& mesh,
    const cfdx::core::Field<double, cfdx::core::Location::FACE>& volumetric_flux,
    const std::vector<double>& cell_volume,
    double dt)
{
    if (!(dt > 0.0) || !std::isfinite(dt))
        throw std::invalid_argument("CFL dt must be finite and strictly positive");
    if (volumetric_flux.size() != mesh.n_faces() ||
        cell_volume.size() != mesh.n_cells())
        throw std::invalid_argument("CFL inputs do not match mesh");

    const auto& cells = mesh.cells();
    const auto* faces = cells.faces_data();
    const auto* offsets = cells.offsets_data();
    std::vector<double> local(mesh.n_cells(), 0.0);

    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const double volume = cell_volume[cell];
        if (!std::isfinite(volume) || !(volume > 0.0))
            throw std::invalid_argument("CFL requires finite positive cell volumes");

        const auto off = offsets[cell];
        const auto count = offsets[cell + 1] - off;
        double flux_sum = 0.0;
        for (std::size_t k = 0; k < count; ++k) {
            const double flux = volumetric_flux(faces[off + k]);
            if (!std::isfinite(flux))
                throw std::invalid_argument("CFL flux contains a non-finite value");
            flux_sum += std::abs(flux);
        }
        local[cell] = dt * flux_sum / volume;
    }
    return local;
}

inline CflDiagnostics summarize_cfl(const std::vector<double>& local_cfl)
{
    if (local_cfl.empty())
        throw std::invalid_argument("cannot summarize an empty CFL field");

    CflDiagnostics result;
    result.local_cfl_min = std::numeric_limits<double>::infinity();
    result.local_cfl_max = 0.0;
    result.limiting_cell = 0;

    for (std::size_t i = 0; i < local_cfl.size(); ++i) {
        const double value = local_cfl[i];
        if (!std::isfinite(value) || value < 0.0)
            throw std::invalid_argument("CFL values must be finite and non-negative");
        result.local_cfl_min = std::min(result.local_cfl_min, value);
        if (value > result.local_cfl_max) {
            result.local_cfl_max = value;
            result.limiting_cell = i;
        }
    }
    result.global_cfl = result.local_cfl_max;
    return result;
}

inline std::vector<double> compute_cell_characteristic_lengths(
    const cfdx::core::Mesh& mesh)
{
    const std::size_t n_faces = mesh.n_faces();
    std::vector<cfdx::core::Vec3> face_centres(n_faces);
    std::vector<cfdx::core::Vec3> face_Sf(n_faces);

    const auto& points = mesh.points();
    const double* px = points.x_data();
    const double* py = points.y_data();
    const double* pz = points.z_data();
    const auto* face_vertices = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();

    for (std::size_t f = 0; f < n_faces; ++f) {
        const cfdx::core::Offset off = face_offsets[f];
        const cfdx::core::Offset count = face_offsets[f + 1] - off;
        const auto geometry = cfdx::core::compute_face_geometry(
            px, py, pz, face_vertices, off, count);
        face_centres[f] = geometry.centre;
        face_Sf[f] = geometry.Sf;
    }

    const auto& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();
    std::vector<double> lengths(mesh.n_cells(), 0.0);

    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const cfdx::core::Offset off = cell_offsets[cell];
        const cfdx::core::Offset count = cell_offsets[cell + 1] - off;
        const auto geometry = cfdx::core::compute_cell_geometry_oriented(
            face_centres.data(), face_Sf.data(), cell_faces + off, count,
            static_cast<cfdx::core::CellIndex>(cell), mesh.ownership());

        if (!std::isfinite(geometry.volume) || !(geometry.volume > 0.0))
            throw std::invalid_argument("invalid cell volume for CFL characteristic length");

        double surface_area = 0.0;
        for (std::size_t k = 0; k < count; ++k) {
            const std::size_t face = cell_faces[off + k];
            const auto sf = face_Sf[face];
            surface_area += std::sqrt(sf.x * sf.x + sf.y * sf.y + sf.z * sf.z);
        }
        if (!std::isfinite(surface_area) || !(surface_area > 0.0))
            throw std::invalid_argument("invalid cell surface area for CFL characteristic length");

        const double h = 2.0 * geometry.volume / surface_area;
        if (!std::isfinite(h) || !(h > 0.0))
            throw std::invalid_argument("invalid CFL characteristic length");
        lengths[cell] = h;
    }
    return lengths;
}

struct TimeStepDecision {
    double old_dt = 0.0;
    double new_dt = 0.0;
    double measured_cfl = 0.0;
    TimeStepChangeReason reason = TimeStepChangeReason::NoChange;
};

struct TimeStepControllerControls {
    double target_cfl = 1.0;
    double min_dt = 1e-12;
    double max_dt = 1e12;
    double growth_limit = 1.25;
    double shrink_limit = 0.5;
    std::size_t max_retries = 8;
};

inline void validate_time_step_controls(const TimeStepControllerControls& c)
{
    if (!std::isfinite(c.target_cfl) || !(c.target_cfl > 0.0) ||
        !std::isfinite(c.min_dt) || !(c.min_dt > 0.0) ||
        !std::isfinite(c.max_dt) || c.max_dt < c.min_dt ||
        !std::isfinite(c.growth_limit) || !(c.growth_limit >= 1.0) ||
        !std::isfinite(c.shrink_limit) || !(c.shrink_limit > 0.0) ||
        c.shrink_limit > 1.0 || c.max_retries == 0) {
        throw std::invalid_argument("invalid timestep controller controls");
    }
}

class AuditableTimeStepController {
public:
    explicit AuditableTimeStepController(TimeStepControllerControls controls = {})
        : controls_(controls)
    {
        validate_time_step_controls(controls_);
    }

    TimeStepDecision propose(double dt, double measured_cfl)
    {
        if (!std::isfinite(dt) || !(dt > 0.0) ||
            !std::isfinite(measured_cfl) || !(measured_cfl > 0.0))
            throw std::invalid_argument("invalid timestep/CFL state");

        const double factor = std::clamp(
            std::sqrt(controls_.target_cfl / measured_cfl),
            controls_.shrink_limit, controls_.growth_limit);
        const double unconstrained = dt * factor;
        const double candidate = std::clamp(
            unconstrained, controls_.min_dt, controls_.max_dt);

        TimeStepChangeReason reason = TimeStepChangeReason::TargetCfl;
        if (candidate == controls_.min_dt && unconstrained < controls_.min_dt)
            reason = TimeStepChangeReason::MinDtBound;
        else if (candidate == controls_.max_dt && unconstrained > controls_.max_dt)
            reason = TimeStepChangeReason::MaxDtBound;
        else if (candidate == dt)
            reason = TimeStepChangeReason::NoChange;

        const TimeStepDecision decision{dt, candidate, measured_cfl, reason};
        record(decision);
        return decision;
    }

    TimeStepDecision rollback(double dt, std::size_t retry)
    {
        if (!std::isfinite(dt) || !(dt > 0.0))
            throw std::invalid_argument("rollback requires a finite positive dt");
        if (retry == 0 || retry > controls_.max_retries)
            throw std::invalid_argument("invalid deterministic timestep retry index");

        const double candidate = std::max(controls_.min_dt, dt * controls_.shrink_limit);
        const TimeStepDecision decision{
            dt, candidate, std::numeric_limits<double>::quiet_NaN(),
            TimeStepChangeReason::NonlinearFailureRollback};
        record(decision);
        return decision;
    }

    void record(const TimeStepDecision& decision)
    {
        if (!std::isfinite(decision.old_dt) || !(decision.old_dt > 0.0) ||
            !std::isfinite(decision.new_dt) || !(decision.new_dt > 0.0))
            throw std::invalid_argument("cannot record invalid timestep decision");
        history_.push_back(decision);
    }

    const std::vector<TimeStepDecision>& history() const { return history_; }
    void clear_history() { history_.clear(); }
    const TimeStepControllerControls& controls() const { return controls_; }

private:
    TimeStepControllerControls controls_;
    std::vector<TimeStepDecision> history_;
};

class TimeStepRollback {
public:
    explicit TimeStepRollback(cfdx::core::TimeIntegrationContext& context)
        : context_(context) {}

    void begin(const cfdx::core::Field<double, cfdx::core::Location::CELL>& state)
    {
        state_snapshot_ = state;
        prev_snapshot_ = context_.phi_prev;
        curr_snapshot_ = context_.phi_curr;
        has_prev_snapshot_ = context_.has_prev;
        dt_prev_snapshot_ = context_.dt_prev;
        active_ = true;
    }

    void reject(cfdx::core::Field<double, cfdx::core::Location::CELL>& state)
    {
        if (!active_)
            throw std::logic_error("cannot rollback without an active timestep transaction");
        state = state_snapshot_;
        context_.phi_prev = prev_snapshot_;
        context_.phi_curr = curr_snapshot_;
        context_.has_prev = has_prev_snapshot_;
        context_.dt_prev = dt_prev_snapshot_;
        active_ = false;
    }

    void commit() { active_ = false; }
    bool active() const { return active_; }

private:
    cfdx::core::TimeIntegrationContext& context_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> state_snapshot_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> prev_snapshot_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> curr_snapshot_;
    bool has_prev_snapshot_ = false;
    double dt_prev_snapshot_ = 0.0;
    bool active_ = false;
};

class NonlinearStateRollback {
public:
    NonlinearStateRollback(
        cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
        cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure)
        : velocity_(velocity), pressure_(pressure) {}

    void begin()
    {
        velocity_snapshot_ = velocity_;
        pressure_snapshot_ = pressure_;
        active_ = true;
    }

    void reject()
    {
        if (!active_)
            throw std::logic_error("cannot rollback nonlinear state without an active transaction");
        velocity_ = velocity_snapshot_;
        pressure_ = pressure_snapshot_;
        active_ = false;
    }

    void commit() { active_ = false; }
    bool active() const { return active_; }

private:
    cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity_;
    cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> velocity_snapshot_;
    cfdx::core::Field<double, cfdx::core::Location::CELL> pressure_snapshot_;
    bool active_ = false;
};

struct NonlinearRetryControls {
    std::size_t max_retries = 3;
    double relaxation_shrink = 0.5;
    double minimum_alpha_u = 0.1;
    double minimum_alpha_p = 0.05;
};

inline void validate_nonlinear_retry_controls(const NonlinearRetryControls& c)
{
    if (c.max_retries == 0 ||
        !std::isfinite(c.relaxation_shrink) || !(c.relaxation_shrink > 0.0) ||
        c.relaxation_shrink >= 1.0 ||
        !std::isfinite(c.minimum_alpha_u) || !(c.minimum_alpha_u > 0.0) ||
        c.minimum_alpha_u > 1.0 ||
        !std::isfinite(c.minimum_alpha_p) || !(c.minimum_alpha_p > 0.0) ||
        c.minimum_alpha_p > 1.0)
        throw std::invalid_argument("invalid nonlinear retry controls");
}

class NonlinearRetryController {
public:
    explicit NonlinearRetryController(NonlinearRetryControls controls = {})
        : controls_(controls)
    {
        validate_nonlinear_retry_controls(controls_);
    }

    void reset()
    {
        retry_count_ = 0;
        alpha_scale_u_ = 1.0;
        alpha_scale_p_ = 1.0;
    }

    bool can_retry() const { return retry_count_ < controls_.max_retries; }

    void reject(bool shrink_u = true, bool shrink_p = true)
    {
        if (!can_retry())
            throw std::runtime_error("nonlinear retry limit exhausted");
        if (!shrink_u && !shrink_p)
            throw std::invalid_argument("nonlinear retry must shrink at least one relaxation channel");
        ++retry_count_;
        if (shrink_u)
            alpha_scale_u_ *= controls_.relaxation_shrink;
        if (shrink_p)
            alpha_scale_p_ *= controls_.relaxation_shrink;
    }

    double alpha_u(double base) const
    {
        return std::max(controls_.minimum_alpha_u, base * alpha_scale_u_);
    }

    double alpha_p(double base) const
    {
        return std::max(controls_.minimum_alpha_p, base * alpha_scale_p_);
    }

    std::size_t retries() const { return retry_count_; }

private:
    NonlinearRetryControls controls_;
    std::size_t retry_count_ = 0;
    double alpha_scale_u_ = 1.0;
    double alpha_scale_p_ = 1.0;
};

} // namespace cfdx::physics
