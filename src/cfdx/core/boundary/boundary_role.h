#pragma once
#include <string>
namespace cfdx::core {
enum class BoundaryRole { NONE, INLET, OUTLET, WALL, SYMMETRY, PERIODIC, INTERFACE };
inline const char* boundary_role_name(BoundaryRole role) noexcept {
    switch (role) {
        case BoundaryRole::INLET: return "inlet";
        case BoundaryRole::OUTLET: return "outlet";
        case BoundaryRole::WALL: return "wall";
        case BoundaryRole::SYMMETRY: return "symmetry";
        case BoundaryRole::PERIODIC: return "periodic";
        case BoundaryRole::INTERFACE: return "interface";
        default: return "none";
    }
}
}