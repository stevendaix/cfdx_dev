#include "cfdx/physics/wall_face_reconstruction.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
using cfdx::core::Vec3;
using cfdx::physics::reconstruct_wall_face_linear;

// Quadratic exact solution: a linear face-pressure reconstruction is expected
// to be O(h^2) in pressure, but the rank-one gradient correction is only O(h)
// for general quadratic velocity. Do not mislabel this first-order traction
// reconstruction as second-order without a further consistent stencil.
static std::array<double, 2> errors(double h) {
    const Vec3 xo{0., 0., 0.}, xf{h, .25*h, 0.}, Sf{1., .25, 0.};
    const Vec3 grad_p{2., -.5, 0.};
    auto p = [](Vec3 x) { return 1. + 2.*x.x - .5*x.y + .8*x.x*x.x + .3*x.y*x.y + .4*x.x*x.y; };
    auto u = [](Vec3 x) { return Vec3{x.x*x.x + 2.*x.x*x.y + x.y*x.y,
                                        .5*x.x*x.x - x.y*x.y, 0.}; };
    const std::array<Vec3, 3> grad_u{{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}};
    const auto out = reconstruct_wall_face_linear(xo, xf, Sf, p(xo), grad_p,
                                                   u(xo), u(xf), grad_u, 1.);
    const auto n = Sf.normalized();
    const double ux = 2.*xf.x + 2.*xf.y, uy = 2.*xf.x + 2.*xf.y;
    const double vx = xf.x, vy = -2.*xf.y;
    const Vec3 traction = n*p(xf) - Vec3{
        2.*ux*n.x + (uy+vx)*n.y,
        (uy+vx)*n.x + 2.*vy*n.y, 0.};
    const double pe = std::abs(out.pressure - p(xf));
    // Remove pressure error to isolate the viscous traction's convergence.
    const auto t_numeric = out.force_per_area - n*out.pressure;
    const auto t_exact = traction - n*p(xf);
    return {pe, (t_numeric - t_exact).mag()};
}
int main() {
    const auto e1 = errors(.5), e2 = errors(.25), e3 = errors(.125);
    for (int j = 0; j < 2; ++j) {
        if (!(std::isfinite(e1[j]) && e1[j] > 1e-10 &&
              e2[j] > 0. && e3[j] > 0.))
            throw std::runtime_error("invalid quadratic MMS errors");
        const double ratio12 = e1[j]/e2[j], ratio23 = e2[j]/e3[j];
        const double expected = j == 0 ? 4.0 : 2.0;
        if (!(std::abs(ratio12 - expected) < 1e-7 &&
              std::abs(ratio23 - expected) < 1e-7))
            throw std::runtime_error("quadratic MMS observed order mismatch");
    }
    std::cout << "WALL_FACE_QUADRATIC_MMS PASS pressure_order=2 traction_order=1"
              << " pressure_errors=" << e1[0] << "," << e2[0] << "," << e3[0]
              << " traction_errors=" << e1[1] << "," << e2[1] << "," << e3[1] << "\n";
}
