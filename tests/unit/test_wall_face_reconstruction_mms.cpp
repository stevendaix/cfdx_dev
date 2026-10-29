#include "cfdx/physics/wall_face_reconstruction.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
using cfdx::core::Vec3;
using cfdx::physics::reconstruct_wall_face_linear;

int main() {
    const Vec3 xo{0.44, 0.14, 0.5};
    const Vec3 xf{0.65, 0.08, 0.5};
    const Vec3 Sf{0.5, -0.2, 0};
    const Vec3 gp{1,0.3,0};
    const std::array<Vec3,3> G{{{0.1, 1.5, 0}, {0.7, -0.2, 0}, {0, 0, 0}}};
    auto P = [](Vec3 v) { return 0.5 + v.x + 0.3*v.y; };
    auto U = [](Vec3 v) { return Vec3{0.1*v.x+1.5*v.y,0.7*v.x-0.2*v.y,0}; };
    auto r = reconstruct_wall_face_linear(xo,xf,Sf,P(xo),gp,U(xo),U(xf),G,1e-3);
    const auto n = Sf.normalized();
    const Vec3 expected = n*P(xf) - Vec3{1e-3*(0.2*n.x+2.2*n.y),
                                         1e-3*(2.2*n.x-0.4*n.y),0};
    if (std::abs(r.pressure-P(xf)) > 1e-14 ||
        (r.force_per_area-expected).mag() > 1e-14)
        throw std::runtime_error("manufactured affine traction mismatch");
    for (std::size_t i=0; i<3; ++i)
        if ((r.velocity_gradient[i]-G[i]).mag() > 1e-14)
            throw std::runtime_error("manufactured affine gradient mismatch");
    bool rejected = false;
    try {
        (void)reconstruct_wall_face_linear(xo,xo,Sf,P(xo),gp,U(xo),U(xf),G,1e-3);
    } catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) throw std::runtime_error("degenerate cell-face distance accepted");
    std::cout << "WALL_FACE_MMS PASS pressure_error=" << std::abs(r.pressure-P(xf))
              << " traction_error=" << (r.force_per_area-expected).mag() << "\n";
}
