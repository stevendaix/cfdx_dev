#include "cfdx/core/fvm/least_squares_gradient.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using cfdx::core::StencilQuality;
using cfdx::core::Vec3;

struct StencilCase {
    const char* name;
    std::vector<Vec3> points;
    double tolerance;
};

double max_error(const Vec3& actual, const Vec3& expected) {
    return std::max({std::abs(actual.x - expected.x),
                     std::abs(actual.y - expected.y),
                     std::abs(actual.z - expected.z)});
}

void check_case(const StencilCase& test) {
    const Vec3 centre{0.0, 0.0, 0.0};
    const Vec3 exact_gradient{2.0, -3.0, 0.5};
    std::vector<double> values;
    values.reserve(test.points.size());
    for (const Vec3& point : test.points)
        values.push_back(exact_gradient.x * point.x +
                         exact_gradient.y * point.y +
                         exact_gradient.z * point.z);

    StencilQuality quality;
    const Vec3 recovered = cfdx::core::least_squares_gradient(
        centre, 0.0, test.points, values, &quality);

    if (!quality.finite || quality.rank != 3 || !quality.full_rank)
        throw std::runtime_error(std::string(test.name) +
                                 ": expected finite full-rank 3-D stencil");
    if (!std::isfinite(quality.condition_estimate) ||
        quality.condition_estimate < 1.0)
        throw std::runtime_error(std::string(test.name) +
                                 ": invalid condition estimate");

    const double error = max_error(recovered, exact_gradient);
    if (!(error <= test.tolerance))
        throw std::runtime_error(std::string(test.name) +
                                 ": linear field reconstruction error exceeds tolerance");

    // Re-evaluate after a constant shift. A fixed LSQ reconstruction must
    // annihilate constants even when coordinates are strongly anisotropic.
    std::vector<double> shifted = values;
    for (double& value : shifted) value += 17.25;
    const Vec3 shifted_gradient = cfdx::core::least_squares_gradient(
        centre, 17.25, test.points, shifted, nullptr);
    if (max_error(shifted_gradient, exact_gradient) > test.tolerance)
        throw std::runtime_error(std::string(test.name) +
                                 ": constant-shift invariance failed");

    std::cout << test.name << ": rank=" << quality.rank
              << ", condition=" << quality.condition_estimate
              << ", max_error=" << error << '\n';
}

}  // namespace

int main() {
    // Orthogonal but highly stretched directions. The current inverse-distance
    // weighting normalizes each row direction; this is a targeted regression
    // probe, not proof of robustness for arbitrary mesh topology.
    const StencilCase anisotropic{
        "anisotropic_1e5_aspect",
        {{1000.0, 0.0, 0.0}, {-1000.0, 0.0, 0.0},
         {0.0, 0.01, 0.0}, {0.0, -0.01, 0.0},
         {0.0, 0.0, 10.0}, {0.0, 0.0, -10.0}},
        2.0e-11
    };

    // Three nearly aligned direction families, with a small but resolvable
    // transverse component. This exercises the singular-value ordering and
    // conditioning diagnostics without demanding recovery below rank cutoff.
    const StencilCase skewed{
        "skewed_near_collinear",
        {{1.0, 0.0, 0.0}, {-1.0, 0.0, 0.0},
         {1.0, 1.0e-2, 0.0}, {-1.0, -1.0e-2, 0.0},
         {1.0, 2.0e-2, 1.0e-3}, {-1.0, -2.0e-2, -1.0e-3}},
        2.0e-8
    };

    check_case(anisotropic);
    check_case(skewed);
    return 0;
}
