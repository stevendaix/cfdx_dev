#include "cfdx/physics/dual_time_navier_stokes.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    try {
        Field<double, Location::CELL> U(2, "U", "m/s", 3);
        Field<double, Location::CELL> p(2, "p", "Pa", 1);
        U.fill(0.0);
        p.fill(0.0);
        DualTimeNavierStokesLifecycle lifecycle(U, p);
        const auto initial = lifecycle.checkpoint();

        auto rejected = [&lifecycle](const DualTimeNavierStokesCheckpoint& bad) {
            bool threw = false;
            try { lifecycle.restore(bad); }
            catch (const std::invalid_argument&) { threw = true; }
            require(threw, "invalid checkpoint must be rejected");
            const auto actual = lifecycle.checkpoint();
            require(actual.step == 0 && actual.physical_time == 0.0 &&
                    !actual.history_valid && actual.dt_previous == 0.0,
                    "rejected checkpoint must not change lifecycle metadata");
            require(actual.U.component_data(0)[0] == 0.0 &&
                    actual.p(0) == 0.0, "rejected checkpoint must preserve accepted fields");
        };

        auto bad = initial;
        bad.U = Field<double, Location::CELL>(3, "U", "m/s", 3);
        rejected(bad);
        bad = initial;
        bad.p(0) = std::numeric_limits<double>::quiet_NaN();
        rejected(bad);
        bad = initial;
        bad.physical_time = std::numeric_limits<double>::infinity();
        rejected(bad);
        bad = initial;
        bad.history_valid = true;
        bad.step = 1;
        bad.dt_previous = 0.01;
        // Missing previous U is invalid, even when the current fields are valid.
        rejected(bad);
        bad.U_previous = U;
        bad.U_previous.component_data(1)[1] = std::numeric_limits<double>::infinity();
        rejected(bad);

        auto valid = initial;
        valid.U.component_data(0)[0] = 0.5;
        valid.p(1) = 2.0;
        valid.U_previous = U;
        valid.history_valid = true;
        valid.physical_time = 0.1;
        valid.dt_previous = 0.1;
        valid.step = 1;
        lifecycle.restore(valid);
        require(lifecycle.step() == 1 && lifecycle.history_valid() &&
                lifecycle.U().component_data(0)[0] == 0.5 &&
                lifecycle.p()(1) == 2.0, "valid checkpoint must restore correctly");

        bad = valid;
        bad.U_previous.component_data(0)[0] = std::numeric_limits<double>::quiet_NaN();
        bool threw = false;
        try { lifecycle.restore(bad); }
        catch (const std::invalid_argument&) { threw = true; }
        require(threw && lifecycle.step() == 1 &&
                lifecycle.U().component_data(0)[0] == 0.5,
                "bad history must preserve an already accepted state");

        std::cout << "dual-time checkpoint transactional contract PASS\n";
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "dual-time checkpoint transactional contract FAIL: " << exc.what() << '\n';
        return 1;
    }
}
