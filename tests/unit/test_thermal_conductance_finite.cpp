#include "cfdx/physics/thermal_zones.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using cfdx::physics::thermal_interface_conductance;

int main() {
    try {
        const double normal = thermal_interface_conductance(2.0, 0.75, 0.25,
                                                              2.0, 4.0, 0.1);
        if (std::abs(normal - 3.7209302325581395) > 1e-14)
            throw std::runtime_error("normal asymmetric network changed");
        const auto reject = [](double area, double distance, double k) {
            try {
                (void)thermal_interface_conductance(area, distance, distance,
                                                    k, k, 0.0);
            } catch (const std::invalid_argument&) {
                return;
            }
            throw std::runtime_error("non-representable conductance was silently accepted");
        };
        reject(1.0, 1.0, std::numeric_limits<double>::denorm_min());
        reject(std::numeric_limits<double>::denorm_min(), 1.e100, 1.e-100);
        reject(1.0, 1.0, std::numeric_limits<double>::infinity());
        std::cout << "finite thermal interface conductance PASS\n";
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "finite thermal interface conductance FAIL: " << exc.what() << '\n';
        return 1;
    }
}
