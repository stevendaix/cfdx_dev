#include "cfdx/physics/wall_distance.h"

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using cfdx::physics::WallDistanceBvh;
using cfdx::physics::WallDistanceVec3;
using cfdx::physics::WallSurface;
using cfdx::physics::WallTriangle;
using cfdx::physics::wd_dot;

namespace {

WallDistanceVec3 cross_product(const WallDistanceVec3& a,
                               const WallDistanceVec3& b) {
    return {a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
}

double brute_ray_distance(const WallSurface& surface,
                          const WallDistanceVec3& origin,
                          const WallDistanceVec3& direction,
                          double max_distance) {
    const double norm = std::sqrt(wd_dot(direction, direction));
    if (!(norm > 1e-30) || !(max_distance > 0.0)) {
        return std::numeric_limits<double>::infinity();
    }
    const WallDistanceVec3 ray = direction * (1.0 / norm);
    const double eps = 1e-12 * std::max(1.0, max_distance);
    double best = std::numeric_limits<double>::infinity();

    for (const auto& triangle : surface.triangles) {
        const auto& v0 = surface.points[triangle.v[0]];
        const auto& v1 = surface.points[triangle.v[1]];
        const auto& v2 = surface.points[triangle.v[2]];
        const auto e1 = v1 - v0;
        const auto e2 = v2 - v0;
        const auto h = cross_product(ray, e2);
        const double a = wd_dot(e1, h);
        if (std::abs(a) < 1e-14) continue;
        const double inverse_a = 1.0 / a;
        const auto offset = origin - v0;
        const double u = inverse_a * wd_dot(offset, h);
        if (u < -1e-12 || u > 1.0 + 1e-12) continue;
        const auto q = cross_product(offset, e1);
        const double v = inverse_a * wd_dot(ray, q);
        if (v < -1e-12 || u + v > 1.0 + 1e-12) continue;
        const double t = inverse_a * wd_dot(e2, q);
        if (t >= -eps && t <= max_distance + eps) {
            best = std::min(best, std::max(0.0, t));
        }
    }
    return best;
}

WallSurface make_planar_surface(std::size_t divisions) {
    WallSurface surface;
    surface.points.reserve((divisions + 1) * (divisions + 1));
    surface.triangles.reserve(2 * divisions * divisions);
    for (std::size_t j = 0; j <= divisions; ++j) {
        for (std::size_t i = 0; i <= divisions; ++i) {
            surface.points.push_back({
                static_cast<double>(i) / static_cast<double>(divisions),
                static_cast<double>(j) / static_cast<double>(divisions),
                0.0
            });
        }
    }
    const auto point_id = [divisions](std::size_t i, std::size_t j) {
        return j * (divisions + 1) + i;
    };
    for (std::size_t j = 0; j < divisions; ++j) {
        for (std::size_t i = 0; i < divisions; ++i) {
            const auto p00 = point_id(i, j);
            const auto p10 = point_id(i + 1, j);
            const auto p01 = point_id(i, j + 1);
            const auto p11 = point_id(i + 1, j + 1);
            surface.triangles.push_back({{p00, p10, p11}});
            surface.triangles.push_back({{p00, p11, p01}});
        }
    }
    return surface;
}

std::size_t parse_positive(const char* value, const char* label,
                           std::size_t maximum) {
    char* end = nullptr;
    const auto parsed = std::strtoull(value, &end, 10);
    if (end == value || *end != '\0' || parsed == 0 || parsed > maximum) {
        throw std::invalid_argument(std::string(label) + " must be in [1, " +
                                    std::to_string(maximum) + "]");
    }
    return static_cast<std::size_t>(parsed);
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const std::size_t divisions =
            argc > 1 ? parse_positive(argv[1], "divisions", 256) : 64;
        const std::size_t ray_count =
            argc > 2 ? parse_positive(argv[2], "ray_count", 1000000) : 4096;
        if (argc > 3) {
            std::cerr << "Usage: benchmark_wall_distance_bvh_rays [divisions=64] [rays=4096]\n";
            return 2;
        }

        const WallSurface surface = make_planar_surface(divisions);
        const WallDistanceBvh bvh(surface);
        std::mt19937_64 generator(0xCFD0B17ULL);
        std::uniform_real_distribution<double> coordinate(
            0.001, 0.999);
        std::vector<WallDistanceVec3> origins;
        origins.reserve(ray_count);
        for (std::size_t i = 0; i < ray_count; ++i) {
            origins.push_back({coordinate(generator), coordinate(generator), 1.0});
        }

        const WallDistanceVec3 direction{0.0, 0.0, -1.0};
        constexpr double max_distance = 2.0;

        // Warm up both paths before timing. Correctness is checked for every
        // ray after timing; no performance result is emitted on mismatch.
        for (std::size_t i = 0; i < std::min<std::size_t>(32, ray_count); ++i) {
            const double expected = brute_ray_distance(
                surface, origins[i], direction, max_distance);
            const double actual = bvh.ray_distance(
                origins[i], direction, max_distance);
            if (!std::isfinite(expected) || !std::isfinite(actual) ||
                std::abs(expected - actual) > 1e-12) {
                throw std::runtime_error("warm-up ray differs from brute-force oracle");
            }
        }

        volatile double checksum = 0.0;
        const auto brute_start = std::chrono::steady_clock::now();
        std::vector<double> reference(ray_count);
        for (std::size_t i = 0; i < ray_count; ++i) {
            reference[i] = brute_ray_distance(
                surface, origins[i], direction, max_distance);
            checksum = checksum + reference[i];
        }
        const auto brute_end = std::chrono::steady_clock::now();

        const auto bvh_start = std::chrono::steady_clock::now();
        std::vector<double> accelerated(ray_count);
        for (std::size_t i = 0; i < ray_count; ++i) {
            accelerated[i] = bvh.ray_distance(
                origins[i], direction, max_distance);
            checksum = checksum + accelerated[i];
        }
        const auto bvh_end = std::chrono::steady_clock::now();

        for (std::size_t i = 0; i < ray_count; ++i) {
            if (!std::isfinite(reference[i]) || !std::isfinite(accelerated[i]) ||
                std::abs(reference[i] - accelerated[i]) > 1e-12) {
                throw std::runtime_error("BVH ray result differs from brute-force oracle");
            }
        }

        const double brute_seconds =
            std::chrono::duration<double>(brute_end - brute_start).count();
        const double bvh_seconds =
            std::chrono::duration<double>(bvh_end - bvh_start).count();
        std::cout << std::setprecision(12)
                  << "triangles,rays,brute_force_seconds,bvh_seconds,speedup,checksum\n"
                  << surface.triangles.size() << ',' << ray_count << ','
                  << brute_seconds << ',' << bvh_seconds << ','
                  << (bvh_seconds > 0.0 ? brute_seconds / bvh_seconds : 0.0) << ','
                  << checksum << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "benchmark_wall_distance_bvh_rays: " << error.what() << '\n';
        return 1;
    }
}
