#include "cfdx/physics/wall_shear_quadratic_fit.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using cfdx::core::Vec3;
using cfdx::physics::WallShearSample;
using cfdx::physics::reconstruct_wall_shear_quadratic;

static void check(bool condition, const char* msg) {
    if (!condition) throw std::runtime_error(msg);
}

static double run_case(double h, double theta, double cubic) {
    const Vec3 t{std::cos(theta), std::sin(theta), 0.0};
    const Vec3 n{std::sin(theta), -std::cos(theta), 0.0};
    const Vec3 inward = n * -1.0;
    const Vec3 face{1.0, -0.5, 0.0};
    constexpr double a = 2.25, b = 0.75, c = -0.9;
    auto sample = [&](double s, double depth) {
        const double scalar = a*depth + b*depth*depth +
            c*s*depth + cubic*depth*depth*depth;
        const Vec3 grad_s = t * (c*depth) + inward *
            (a + 2*b*depth + c*s + 3*cubic*depth*depth);
        return WallShearSample{face + t*s + inward*depth, t*scalar,
            {grad_s*t.x, grad_s*t.y, grad_s*t.z}};
    };
    const auto one = sample(0.35*h, h);
    const auto two = sample(-0.20*h, 2*h);
    const auto result = reconstruct_wall_shear_quadratic(
        face, n, t, {0,0,0}, one, two, 0.4);
    check(std::abs(result.first_distance-h) < 1e-12, "wrong first distance");
    check(std::abs(result.second_distance-2*h) < 1e-12, "wrong second distance");
    return std::abs(result.shear_on_body_tangent - 0.4*a);
}

int main() {
    for (double theta : {0.0, 0.61, -1.05}) {
        const double exact = run_case(0.3, theta, 0.0);
        check(exact < 2e-13, "quadratic/sheared wall must be exact");
        const double e1 = run_case(0.16, theta, 0.8);
        const double e2 = run_case(0.08, theta, 0.8);
        const double e3 = run_case(0.04, theta, 0.8);
        check(e1 > 0 && e2 > 0 && e3 > 0, "degenerate convergence experiment");
        check(e1/e2 > 3.99 && e2/e3 > 3.99,
              "wall shear must converge at second order on cubic MMS");
        std::cout << "WALL_SHEAR_QUADRATIC_MMS theta=" << theta
                  << " h_0.16=" << e1 << " h_0.08=" << e2
                  << " h_0.04=" << e3 << " ratio12=" << e1/e2
                  << " ratio23=" << e2/e3 << '\n';
    }
    const Vec3 n{0,-1,0}, t{1,0,0};
    WallShearSample bad1{{0,0.1,0}, {0.1,0,0}, {{{0,1,0},{0,0,0},{0,0,0}}}};
    WallShearSample bad2 = bad1;
    bool rejected = false;
    try { (void)reconstruct_wall_shear_quadratic({0,0,0},n,t,{0,0,0},bad1,bad2,1); }
    catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "equal-depth degenerate stencil accepted");
    // Almost coincident depths amplify the derivative in floating-point.
    bad2.cell_position.y = 0.12;
    rejected = false;
    try { (void)reconstruct_wall_shear_quadratic({0,0,0},n,t,{0,0,0},bad1,bad2,1); }
    catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "near-coincident wall-normal samples accepted");
    bad2.cell_position.y = 0.3;
    bad2.cell_velocity.x = std::numeric_limits<double>::quiet_NaN();
    rejected = false;
    try { (void)reconstruct_wall_shear_quadratic({0,0,0},n,t,{0,0,0},bad1,bad2,1); }
    catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "NaN cell velocity accepted");
    std::cout << "WALL_SHEAR_QUADRATIC_MMS: PASS\n";
}
