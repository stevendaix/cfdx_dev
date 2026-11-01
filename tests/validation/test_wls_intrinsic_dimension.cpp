// Reproduce 2-D planar WLS contract without weakening the default 3-D contract.
#include "cfdx/core/fvm/least_squares_gradient.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
static bool close(double x, double target, double tolerance = 2e-11) {
    return std::abs(x - target) <= tolerance;
}
int main() {
    try {
        const Vec3 origin{1.0, -2.0, 0.3};
        const Vec3 tangent_a{1.0, 0.0, 1.0};
        const Vec3 tangent_b{0.0, 1.0, 1.0};
        const Vec3 target{2.0, -3.0, -1.0}; // tangent to this tilted plane
        std::vector<Vec3> points;
        std::vector<double> vals;
        for (double a : {-2.0, -0.5, 1.0, 2.0}) {
            for (double b : {-1.0, 0.7, 2.0}) {
                const Vec3 point = origin + tangent_a*a + tangent_b*b;
                points.push_back(point);
                vals.push_back(5.0 + target.dot(point-origin));
            }
        }
        for (auto weighting : {GradientWeighting::UNIFORM,
                               GradientWeighting::INVERSE_DISTANCE,
                               GradientWeighting::INVERSE_DISTANCE_SQUARED}) {
            StencilQuality q;
            const Vec3 g = weighted_least_squares_gradient(
                origin, 5.0, points, vals, weighting, &q, 2);
            require(q.rank == 2 && q.dimension == 2 && q.full_rank,
                    "2-D rank contract incorrect");
            require(q.well_conditioned(1e5), "valid planar stencil rejected");
            require(close(g.x,target.x) && close(g.y,target.y) && close(g.z,target.z),
                    "2-D tilted-plane affine gradient not exact");
            std::cout << "WLS_2D_PASS weighting=" << static_cast<int>(weighting)
                      << " gradient=" << g.x << ',' << g.y << ',' << g.z
                      << " condition=" << q.condition_estimate << '\n';
        }
        // The robust 3-D solver preserves an affine field on a skewed,
        // full-rank stencil without squaring the condition number.
        {
            const Vec3 exact{1.25, -2.5, 0.75};
            std::vector<Vec3> points3 = {
                {1.0, 0.2, 0.1}, {-1.0, -0.3, 0.2}, {0.1, 1.0, 0.2},
                {0.1, -1.0, 0.3}, {-0.2, -0.1, 0.7}, {0.4, 0.2, -1.0}};
            std::vector<double> values3;
            for (const auto& d : points3) {
                values3.push_back(7.0 + exact.dot(d));
            }
            StencilQuality q3;
            const Vec3 computed = weighted_least_squares_gradient_robust(
                {0,0,0}, 7.0, points3, values3,
                GradientWeighting::INVERSE_DISTANCE_SQUARED, 3, &q3);
            require(q3.dimension == 3 && q3.full_rank && q3.rank == 3,
                    "3-D robust rank contract incorrect");
            require(close(computed.x,exact.x) && close(computed.y,exact.y) &&
                    close(computed.z,exact.z), "3-D robust affine gradient inaccurate");
            std::cout << "WLS_3D_SVD_PASS condition="
                      << q3.condition_estimate << '\n';
        }
        bool rejected = false;
        try {
            (void)weighted_least_squares_gradient(
                origin, 5.0, points, vals, GradientWeighting::UNIFORM);
        } catch (const std::runtime_error&) { rejected = true; }
        require(rejected, "3-D default must reject planar stencil");

        std::vector<Vec3> collinear = {
            origin+tangent_a, origin+tangent_a*2, origin-tangent_a};
        std::vector<double> cval = {6.0,7.0,4.0};
        rejected = false;
        try {
            (void)weighted_least_squares_gradient(origin, 5.0, collinear, cval,
                GradientWeighting::INVERSE_DISTANCE_SQUARED, nullptr, 2);
        } catch (const std::runtime_error&) { rejected = true; }
        require(rejected, "intrinsic 2-D must reject collinear stencil");

        vals[0] = std::numeric_limits<double>::quiet_NaN();
        rejected = false;
        try {
            (void)weighted_least_squares_gradient(origin, 5.0, points, vals,
                GradientWeighting::UNIFORM, nullptr, 2);
        } catch (const std::runtime_error&) { rejected = true; }
        require(rejected, "2-D must reject non-finite data");
        std::cout << "WLS_INTRINSIC_CONTRACT: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "WLS_INTRINSIC_CONTRACT: FAIL " << e.what() << '\n';
        return 1;
    }
}
