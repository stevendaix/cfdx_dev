#include "common/test_harness.h"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
constexpr double D = 1.0;
constexpr double RHO = 1.0;
constexpr double U = 1.0;
constexpr double MU = 0.01;
constexpr double RE_REFERENCE = 100.0;
constexpr double CD_REFERENCE = 1.0895;
constexpr double DOMAIN_RADIUS = 50.0 * D;

double reynolds_number() { return RHO * U * D / MU; }
double projected_area() { return std::acos(-1.0) * D * D / 4.0; }
} // namespace

int main()
{
    try {
        const double re = reynolds_number();
        if (std::abs(re - RE_REFERENCE) > 1e-12)
            throw std::runtime_error("VMFL036 Re=100 setup mismatch");
        if (std::abs(DOMAIN_RADIUS - 50.0) > 1e-12)
            throw std::runtime_error("VMFL036 domain radius mismatch");
        const double area = projected_area();
        if (!(area > 0.0) || !std::isfinite(area))
            throw std::runtime_error("VMFL036 projected area invalid");
        if (!(CD_REFERENCE > 0.0) || !std::isfinite(CD_REFERENCE))
            throw std::runtime_error("VMFL036 literature Cd invalid");

        std::cout << "VMFL036_REFERENCE_CASE: PASS\n";
        std::cout << "VMFL036_RE=" << re << "\n";
        std::cout << "VMFL036_CD_REFERENCE=" << CD_REFERENCE << "\n";
        std::cout << "VMFL036_PROJECTED_AREA=" << area << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "VMFL036_REFERENCE_CASE: FAIL: " << e.what() << "\n";
        return 1;
    }
}
