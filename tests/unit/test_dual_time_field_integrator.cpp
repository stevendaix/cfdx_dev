#include "cfdx/physics/dual_time_field_integrator.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace cfdx::core;
    using namespace cfdx::physics;

    DualTimeField state(1, "u", "1", 1);
    state(0) = 1.0;

    DualTimeFieldSolveControls controls;
    controls.pseudo_time.physical_scheme = DualTimePhysicalScheme::BDF2;
    controls.pseudo_time.max_pseudo_iterations = 200;
    controls.pseudo_time.pseudo_dt_initial = 0.5;
    controls.pseudo_time.pseudo_dt_min = 1.0e-6;
    controls.pseudo_time.pseudo_dt_max = 2.0;
    controls.pseudo_time.absolute_tolerance = 1.0e-12;
    controls.pseudo_time.relative_tolerance = 1.0e-10;
    controls.physical_dt_min = 0.01;
    controls.physical_dt_max = 0.2;
    controls.temporal_absolute_tolerance = 1.0e-2;
    controls.temporal_relative_tolerance = 1.0e-2;
    controls.physical_admissibility = [](const DualTimeField& field) {
        return std::isfinite(field(0)) && field(0) > 0.0;
    };

    DualTimeFieldIntegrator integrator(state);

    auto rhs = [](const DualTimeField& u, DualTimeField& out) {
        out(0) = -u(0);
    };

    const auto first = integrator.advance(0.05, controls, rhs);
    if (!first.acceptance.accepted())
        throw std::runtime_error("bootstrap dual-time step was not accepted");
    if (!integrator.history.has_previous)
        throw std::runtime_error("accepted physical step did not commit history");
    if (!(first.nonlinear_report.converged &&
          first.nonlinear_report.pseudo_iterations > 0))
        throw std::runtime_error("bootstrap pseudo-time solve did not converge");

    const auto second = integrator.advance(0.05, controls, rhs);
    if (!second.acceptance.accepted())
        throw std::runtime_error("BDF2 dual-time step was not accepted");
    if (!(second.state(0) > 0.0 && second.state(0) < first.state(0)))
        throw std::runtime_error("dual-time decay state is not physically consistent");
    if (!(second.temporal_error.normalized_error >= 0.0))
        throw std::runtime_error("invalid temporal error diagnostic");

    const double exact = std::exp(-0.1);
    if (std::abs(second.state(0) - exact) > 2.0e-2)
        throw std::runtime_error("dual-time BDF2 decay is outside verification bound");

    bool rejected = false;
    controls.physical_admissibility = [](const DualTimeField&) { return false; };
    try {
        (void)integrator.advance(0.05, controls, rhs);
    } catch (const DualTimeConvergenceFailure&) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("inadmissible physical step was not rejected");
    if (std::abs(integrator.history.current(0) - second.state(0)) > 1.0e-14)
        throw std::runtime_error("rejected physical step modified accepted history");

    std::cout << "dual-time field integrator tests: PASS\n";
    return 0;
}
