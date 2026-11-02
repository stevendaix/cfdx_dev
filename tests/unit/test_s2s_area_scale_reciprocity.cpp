#include "cfdx/physics/radiation_s2s.h"

#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
void require(bool ok, const char* case_name) {
    if (!ok) throw std::runtime_error(case_name);
}
void must_reject(const std::vector<double>& A, const std::vector<double>& F,
                 const char* case_name) {
    bool rejected = false;
    try {
        cfdx::physics::validate_s2s_view_factors(F, A, 1e-10);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, case_name);
}
}

int main() {
    using cfdx::physics::validate_s2s_view_factors;
    const std::vector<double> bad{0.0, 0.5, 0.5, 0.0};
    // Scaling the enclosure's physical dimensions cannot make a 100%
    // reciprocity defect disappear. The old max(1, A*F) scale did exactly that.
    must_reject({2.0, 4.0}, bad, "macroscopic reciprocity mismatch");
    must_reject({2e-10, 4e-10}, bad, "micro-patch reciprocity mismatch");
    must_reject({2e-200, 4e-200}, bad, "very-small-patch reciprocity mismatch");

    const std::vector<double> correct{0.0, 0.5, 0.25, 0.0};
    for (double factor : {1.0, 1e-10, 1e-200, 1e200}) {
        const std::vector<double> A{2.0 * factor, 4.0 * factor};
        validate_s2s_view_factors(correct, A, 1e-10);
        const std::vector<double> approximately_correct{0.0, 0.5, 0.25 * (1 + 2e-12), 0.0};
        validate_s2s_view_factors(approximately_correct, A, 1e-10);
    }
    must_reject({1e-12, 1e-12}, {0.0, 0.0, 1e-14, 0.0},
                "unilateral small view factor");

    // Physical oracle: a uniform scaling of surface areas must preserve
    // radiosities/flux densities and scale integrated radiative power by the
    // same factor. This also exercises the *production* radiosity path.
    using namespace cfdx::physics;
    const std::vector<double> eps{0.7,0.8}, T{600.0,400.0};
    const auto macroscopic = solve_s2s_radiosity(
        {2.0,4.0},eps,T,correct);
    const auto microscopic = solve_s2s_radiosity(
        {2e-10,4e-10},eps,T,correct);
    require(macroscopic.converged && microscopic.converged,
            "valid scaled enclosures did not converge");
    for (std::size_t i=0; i<2; ++i) {
        require(std::abs(macroscopic.radiosity[i]-microscopic.radiosity[i]) <
                    1e-12*std::max(1.0,std::abs(macroscopic.radiosity[i])),
                "radiosity must be invariant to area scale");
        require(std::abs(macroscopic.net_flux[i]-microscopic.net_flux[i]) <
                    1e-12*std::max(1.0,std::abs(macroscopic.net_flux[i])),
                "flux density must be invariant to area scale");
    }
    require(std::abs(microscopic.total_power/1e-10-macroscopic.total_power) <
                1e-11*std::max(1.0,std::abs(macroscopic.total_power)),
            "integrated radiation power must scale with area");
    std::cout << "S2S_AREA_SCALE_RECIPROCITY: PASS\n";
}
