#include "cfdx/core/numerics/spatial_scheme_contract.h"

#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void expect(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void expect_throw(const std::function<void()>& fn, const std::string& message) {
    bool thrown = false;
    try { fn(); } catch (const std::invalid_argument&) { thrown = true; }
    expect(thrown, message);
}
}

int main() {
    using namespace cfdx::core;

    {
        const auto s = make_spatial_scheme_selection(
            SpatialSchemeRole::Gradient,
            NumericalMethodFamily::Gradient,
            "numerics.gradient.gauss");
        expect(s.method.resolved, "gradient selection must resolve");
        expect(s.method.method_id == "gradient.gauss_cell", "gradient method id mismatch");
    }

    {
        const auto s = make_spatial_scheme_selection(
            SpatialSchemeRole::Convection,
            NumericalMethodFamily::Convection,
            "numerics.convection.upwind");
        expect(s.method.resolved, "convection selection must resolve");
        expect(s.method.method_id == "convection.upwind", "convection method id mismatch");
    }

    expect_throw([&] {
        (void)make_spatial_scheme_selection(
            SpatialSchemeRole::Interpolation,
            NumericalMethodFamily::Convection,
            "numerics.convection.upwind");
    }, "cross-role spatial selection must be rejected");

    {
        SpatialSchemeContract c;
        c.role = SpatialSchemeRole::Convection;
        c.configuration_key = "numerics.convection.upwind";
        c.formulation = "phi_f = phi_upwind";
        c.bounded = true;
        c.monotone = true;
        c.conservation = ConservationContract::LocalFaceConservative;
        validate_spatial_scheme_contract(c);
    }

    expect_throw([&] {
        SpatialSchemeContract c;
        c.role = SpatialSchemeRole::Convection;
        c.configuration_key = "numerics.convection.invalid";
        c.formulation = "invalid";
        c.bounded = true;
        validate_spatial_scheme_contract(c);
    }, "bounded convection without conservation must be rejected");

    std::cout << "spatial scheme contract tests: PASS\\n";
    return 0;
}
