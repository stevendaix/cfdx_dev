#include "cfdx/physics/wall_face_traction_reconstruction.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

static void require(bool condition, const char* msg) {
    if (!condition) throw std::runtime_error(msg);
}

int main() {
    using cfdx::core::Vec3;
    using cfdx::physics::reconstruct_wall_face_traction;
    const Vec3 d{0.2,-0.5,0.0}, n{0.0,-1.0,0.0};
    // Exact affine Couette Ux=2y, Uy=0, and pressure p=1+3x.
    // Owner at x=y=0, Dirichlet wall velocity Ux=-1 at Cf.
    const auto exact = reconstruct_wall_face_traction(
        d,n,1.0,{3,0,0},{0,0,0},{-1,0,0},
        {{{0,2,0},{0,0,0},{0,0,0}}},0.5);
    require(std::abs(exact.pressure-1.6) < 1e-12,
            "pressure tangential extrapolation is not affine exact");
    require(std::abs(exact.viscous_on_fluid_normal.x + 1.0) < 1e-12,
            "Couette wall shear is incorrect");
    require(std::abs(exact.viscous_on_fluid_normal.y) < 1e-12,
            "Couette generated spurious normal traction");
    // Independent affine oracle with a NONZERO normal pressure derivative:
    // p(x,y) = 1 + 3*x + 4*y; at d=(0.2,-0.5), p_wall=-0.4.
    // A forced homogeneous-Neumann extrapolation would yield +1.6 instead.
    const auto nonzero_dp_dn = reconstruct_wall_face_traction(
        d,n,1.0,{3,4,0},{0,0,0},{-1,0,0},
        {{{0,2,0},{0,0,0},{0,0,0}}},0.5,-4.0);
    require(std::abs(nonzero_dp_dn.pressure+0.4) < 1e-12,
            "prescribed wall dp/dn is not affine exact");
    require(std::abs(nonzero_dp_dn.viscous_on_fluid_normal.x+1.0) < 1e-12,
            "nonzero wall dp/dn contaminated the shear traction");
    // An inaccurate raw cell gradient must still meet the Dirichlet
    // constraint along d; normal traction is not guaranteed exact without a
    // convergent tangential reconstruction, so test the known exact input.
    bool refused = false;
    try {
        (void)reconstruct_wall_face_traction(
            {0,0,0},n,1.0,{0,0,0},{0,0,0},{0,0,0},
            {{{0,0,0},{0,0,0},{0,0,0}}},0.5);
    } catch (const std::invalid_argument&) { refused = true; }
    require(refused, "accepted degenerate face distance");
    refused = false;
    try {
        (void)reconstruct_wall_face_traction(
            d,n,std::numeric_limits<double>::quiet_NaN(),{0,0,0},
            {0,0,0},{0,0,0},{{{0,0,0},{0,0,0},{0,0,0}}},0.5);
    } catch (const std::invalid_argument&) { refused = true; }
    require(refused, "accepted NaN pressure");
    refused = false;
    try {
        (void)reconstruct_wall_face_traction(
            d,n,1.0,{3,4,0},{0,0,0},{-1,0,0},
            {{{0,2,0},{0,0,0},{0,0,0}}},0.5,
            std::numeric_limits<double>::quiet_NaN());
    } catch (const std::invalid_argument&) { refused = true; }
    require(refused, "accepted nonfinite pressure normal derivative");
    std::cout << "WALL_FACE_TRACTION_MMS: PASS\n";
}
