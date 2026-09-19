// M0.2-T03/T04 — Tests for mesh quality metrics

#include "cfdx/core/geometry/mesh_quality.h"
#include "test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
const double PI = std::acos(-1.0);
}

int main() {
    run_case("orthogonal_face_zero_non_orthogonality", []() {
        Vec3 face_centre{1.0, 0.0, 0.0};
        Vec3 cell_centre{0.0, 0.0, 0.0};
        Vec3 Sf{1.0, 0.0, 0.0};
        auto q = compute_face_quality(face_centre, cell_centre, Sf);
        EXPECT_NEAR(q.non_orthogonality, 0.0, 1e-12);
        EXPECT_NEAR(q.skewness, 0.0, 1e-12);
    });

    run_case("skewed_face_nonzero_skewness", []() {
        Vec3 face_centre{1.0, 0.5, 0.0};
        Vec3 cell_centre{0.0, 0.0, 0.0};
        Vec3 Sf{1.0, 0.0, 0.0};
        auto q = compute_face_quality(face_centre, cell_centre, Sf);
        double expected = 0.5 / std::sqrt(1.25);
        EXPECT_NEAR(q.skewness, expected, 1e-12);
    });

    run_case("45_degree_non_orthogonality", []() {
        Vec3 face_centre{1.0, 0.0, 0.0};
        Vec3 cell_centre{0.0, 0.0, 0.0};
        Vec3 Sf{1.0, 1.0, 0.0};
        auto q = compute_face_quality(face_centre, cell_centre, Sf);
        EXPECT_NEAR(q.non_orthogonality, PI / 4.0, 1e-12);
        EXPECT_NEAR(q.non_orthogonality_deg, 45.0, 1e-12);
    });

    run_case("zero_distance_safe", []() {
        Vec3 face_centre{0.0, 0.0, 0.0};
        Vec3 cell_centre{0.0, 0.0, 0.0};
        Vec3 Sf{1.0, 0.0, 0.0};
        auto q = compute_face_quality(face_centre, cell_centre, Sf);
        EXPECT_NEAR(q.skewness, 0.0, 1e-12);
        EXPECT_NEAR(q.non_orthogonality, 0.0, 1e-12);
    });

    run_case("zero_area_safe", []() {
        Vec3 face_centre{1.0, 0.0, 0.0};
        Vec3 cell_centre{0.0, 0.0, 0.0};
        Vec3 Sf{0.0, 0.0, 0.0};
        auto q = compute_face_quality(face_centre, cell_centre, Sf);
        EXPECT_NEAR(q.skewness, 0.0, 1e-12);
    });

    run_case("skewness_bounded_0_1", []() {
        Vec3 face_centre{1.0, 10.0, 0.0};
        Vec3 cell_centre{0.0, 0.0, 0.0};
        Vec3 Sf{1.0, 0.0, 0.0};
        auto q = compute_face_quality(face_centre, cell_centre, Sf);
        EXPECT_TRUE(q.skewness >= 0.0);
        EXPECT_TRUE(q.skewness <= 1.0);
    });

    return run_all();
}