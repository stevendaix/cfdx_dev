#include "cfdx/physics/dual_time_driver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
}

int main() {
    using namespace cfdx::physics;
    DualTimeStepControls controls;
    controls.physical_scheme = DualTimePhysicalScheme::BDF2;
    controls.max_pseudo_iterations = 200;
    controls.pseudo_dt_initial = 0.25;
    controls.pseudo_dt_min = 1.0e-6;
    controls.pseudo_dt_max = 2.0;
    controls.absolute_tolerance = 1.0e-11;
    controls.relative_tolerance = 1.0e-10;

    const double dt = 0.2;
    const double u_nm1 = 1.0;
    const double u_n = std::exp(-dt);
    const double bdf2_root = (2.0 * u_n - 0.5 * u_nm1) / (1.5 + dt);
    auto rhs = [](const std::vector<double>& u) {
        return std::vector<double>{-u[0]};
    };
    DualTimeUpdate update = [](const std::vector<double>& state,
                               const std::vector<double>& residual,
                               double pseudo_dt) {
        return std::vector<double>{state[0] - pseudo_dt * residual[0]};
    };

    const auto result = run_dual_time_step({u_n}, {u_n}, {u_nm1}, dt,
                                           controls, rhs, update);
    require(result.second.converged, "BDF2 dual-time step did not converge");
    require(result.second.pseudo_iterations > 1, "pseudo-time loop was not exercised");
    require(std::abs(result.first[0] - bdf2_root) < 1.0e-9,
            "BDF2 physical root is incorrect");

    controls.physical_scheme = DualTimePhysicalScheme::BACKWARD_EULER;
    const double be_root = u_n / (1.0 + dt);
    const auto be = run_dual_time_step({u_n}, {u_n}, {}, dt, controls, rhs, update);
    require(be.second.converged, "backward-Euler dual-time step did not converge");
    require(std::abs(be.first[0] - be_root) < 1.0e-9,
            "backward-Euler physical root is incorrect");

    controls.max_pseudo_iterations = 3;
    bool failed = false;
    try {
        (void)run_dual_time_step({u_n}, {u_n}, {}, dt, controls, rhs,
            [](const std::vector<double>& state,
               const std::vector<double>&,
               double) { return state; });
    } catch (const DualTimeConvergenceFailure&) {
        failed = true;
    }
    require(failed, "iteration exhaustion must be reported as failure");

    std::cout << "dual-time stepping tests: PASS\n";
    return 0;
}
