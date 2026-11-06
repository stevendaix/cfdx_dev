#include "cfdx/physics/steady_incompressible_solver.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::physics;
using namespace cfdx::core;

static void expect(bool good, const char* message) {
    if (!good) throw std::runtime_error(message);
}

int main() {
    try {
        SparseMatrix m(2, 2);
        m.push_back(0, 0, 2.0);
        m.push_back(0, 1, -1.0);
        m.push_back(1, 0, -1.0);
        m.push_back(1, 1, 2.0);
        m.finalize();
        Vector rhs(2, 0.0), x(2, 0.0);
        rhs(0) = 1.0;
        x(0) = 2.0/3.0;
        x(1) = 1.0/3.0;
        const std::vector<double> incident{1.0, 1.0};
        const auto exact = diagnose_pressure_linear_cancellation(m, rhs, x, incident, 1.e-12);
        expect(exact.valid, "roundoff-level cancellation must be accepted provisionally");
        expect(exact.max_local_mass_ratio < 1.e-12, "roundoff ratio measured incorrectly");
        x(0) += 1.e-6;
        const auto wrong = diagnose_pressure_linear_cancellation(m, rhs, x, incident, 1.e-12);
        expect(!wrong.valid, "large mass defect accepted");
        x(0) -= 1.e-6;
        const auto no_flux = diagnose_pressure_linear_cancellation(m, rhs, x,
            std::vector<double>{1.e-20, 1.0}, 1.e-12);
        expect(!no_flux.valid, "artificially tiny incident flux accepted");
        x(0) = std::numeric_limits<double>::quiet_NaN();
        expect(!diagnose_pressure_linear_cancellation(m, rhs, x, incident, 1.e-12).valid,
               "nonfinite solution accepted");
        x(0) = 2.0/3.0;
        expect(!diagnose_pressure_linear_cancellation(m, rhs, x, incident, 0.0).valid,
               "zero conservation tolerance accepted");
        std::cout << "PASS pressure roundoff diagnostic / local conservation gates\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
    return 0;
}
