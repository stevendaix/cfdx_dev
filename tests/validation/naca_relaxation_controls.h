#pragma once
// NACA qualification A/B diagnostics. The defaults and physical gates are
// unchanged. Invalid experimental settings fail closed instead of being ignored.
#include "cfdx/physics/steady_incompressible_solver.h"
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>

namespace cfdx::validation {

inline double read_naca_relaxation(const char* name, double fallback) {
    const char* raw = std::getenv(name);
    if (!raw) return fallback;
    errno = 0;
    char* end = nullptr;
    const double value = std::strtod(raw, &end);
    if (raw == end || *end != '\0' || errno == ERANGE ||
        !std::isfinite(value) || !(value > 0.0 && value <= 1.0)) {
        throw std::invalid_argument(std::string("Invalid ") + name +
            ": expected a finite number in (0, 1]");
    }
    return value;
}

inline void apply_naca_relaxation_ab(
    cfdx::physics::IncompressibleSolverControls& controls) {
    controls.coupling.alpha_u = read_naca_relaxation(
        "CFDX_NACA_AB_ALPHA_U", controls.coupling.alpha_u);
    controls.coupling.alpha_p = read_naca_relaxation(
        "CFDX_NACA_AB_ALPHA_P", controls.coupling.alpha_p);
}

} // namespace cfdx::validation
