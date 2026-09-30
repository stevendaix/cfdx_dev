// Issue #461 — N5 measured temporal order matrix.
//
// Canonical ODE du/dt = -u, u(0) = 1, exact u(t) = e^{-t}, integrated to
// t = 1 with a 2:1 sweep (dt = 1/32, 1/64, 1/128, 1/256). For every registered
// time scheme the observed order is measured and asserted:
//   euler_explicit  (~1)   crank_nicolson (~2)   bdf2 (~2)
//   rk2 (~2)               rk3 (~3)

#include "cfdx/core/numerics/temporal.h"
#include "cfdx/physics/low_storage_time_integration.h"
#include "verification_metrics.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::verification;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

Field<double, Location::CELL> ode_field(double u) {
    Field<double, Location::CELL> f(1, "u", "1", 1);
    f(0) = u;
    return f;
}

RhsFunction ode_rhs = [](const Field<double, Location::CELL>& x,
                         Field<double, Location::CELL>& rhs) {
    rhs(0) = -x(0);
};

double integrate(TimeScheme scheme, std::size_t n_steps)
{
    const double dt = 1.0 / static_cast<double>(n_steps);
    Field<double, Location::CELL> u = ode_field(1.0);

    if (scheme == TimeScheme::BDF2) {
        // Two-step method: use the exact first step as history.
        Field<double, Location::CELL> u_prev = ode_field(1.0);
        Field<double, Location::CELL> u_curr = ode_field(std::exp(-dt));
        TimeIntegrationContext ctx;
        ctx.n_cells = 1;
        ctx.field_dim = 1;
        for (std::size_t step = 1; step < n_steps; ++step) {
            ctx.phi_prev = u_prev;
            Field<double, Location::CELL> next =
                advance_time(u_curr, dt, ode_rhs, TimeScheme::BDF2, &ctx);
            u_prev = u_curr;
            u_curr = next;
        }
        return u_curr(0);
    }

    if (scheme == TimeScheme::EULER_EXPLICIT ||
        scheme == TimeScheme::CRANK_NICOLSON) {
        for (std::size_t step = 0; step < n_steps; ++step) {
            u = advance_time(u, dt, ode_rhs, scheme);
        }
        return u(0);
    }
    return u(0);
}

double integrate_rk(std::size_t n_steps, bool rk3)
{
    const double dt = 1.0 / static_cast<double>(n_steps);
    Field<double, Location::CELL> u = ode_field(1.0);
    for (std::size_t i = 0; i < n_steps; ++i) {
        if (rk3) {
            low_storage_rk3_step(u, dt, ode_rhs);
        } else {
            low_storage_rk2_step(u, dt, ode_rhs);
        }
    }
    return u(0);
}

struct SchemeCase {
    const char* name;
    double (*integrate_fn)(std::size_t);
    double expected_order;
    double floor;
};

void check_scheme(const SchemeCase& c)
{
    const std::vector<std::size_t> steps = {32u, 64u, 128u, 256u};
    std::vector<double> errors;
    for (const std::size_t n : steps) {
        errors.push_back(std::abs(c.integrate_fn(n) - std::exp(-1.0)));
    }
    std::cout << "N5_ORDER scheme=" << c.name;
    for (std::size_t i = 0; i < errors.size(); ++i) {
        std::cout << " dt=1/" << steps[i] << " err=" << errors[i];
        if (i) std::cout << " order=" << observed_order(errors[i - 1], errors[i]);
    }
    std::cout << "\n";
    require_order(errors, c.expected_order, c.floor, std::string(c.name) + " temporal order");
}

} // namespace

int main()
{
    try {
        std::cout << std::setprecision(12);

        check_scheme({"euler_explicit",
            [](std::size_t n) { return integrate(TimeScheme::EULER_EXPLICIT, n); },
            1.0, 0.7});
        check_scheme({"crank_nicolson",
            [](std::size_t n) { return integrate(TimeScheme::CRANK_NICOLSON, n); },
            2.0, 1.5});
        check_scheme({"bdf2",
            [](std::size_t n) { return integrate(TimeScheme::BDF2, n); },
            2.0, 1.5});
        check_scheme({"rk2",
            [](std::size_t n) { return integrate_rk(n, false); },
            2.0, 1.5});
        check_scheme({"rk3",
            [](std::size_t n) { return integrate_rk(n, true); },
            3.0, 2.5});

        std::cout << "TEMPORAL_ORDER_MATRIX: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "TEMPORAL_ORDER_MATRIX: FAIL: " << e.what() << "\n";
        return 1;
    }
}