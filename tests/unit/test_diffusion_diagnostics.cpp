// Issue #461 — N3 diffusion mesh-quality diagnostics.

#include "cfdx/core/numerics/diffusion_diagnostics.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

Mesh make_two_cell_unit_cubes() {
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
        {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, p[i][0], p[i][1], p[i][2]);

    m.faces().push_face({0,6,9,3});
    m.faces().push_face({0,1,7,6});
    m.faces().push_face({3,9,10,4});
    m.faces().push_face({0,3,4,1});
    m.faces().push_face({6,7,10,9});
    m.faces().push_face({7,10,4,1});
    m.faces().push_face({2,5,11,8});
    m.faces().push_face({1,2,8,7});
    m.faces().push_face({4,10,11,5});
    m.faces().push_face({1,4,5,2});
    m.faces().push_face({7,8,11,10});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 5; ++f) {
        m.ownership().set_owner(f, 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(5, 0);
    m.ownership().set_neighbour(5, 1);
    for (std::size_t f = 6; f < 11; ++f) {
        m.ownership().set_owner(f, 1);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({5,6,7,8,9,10});
    return m;
}

} // namespace

int main() {
    run_case("diffusion_quality_orthogonal_baseline", []() {
        const Mesh mesh = make_two_cell_unit_cubes();
        const auto geometry = make_geometry_cache(mesh);
        const auto q = diagnose_diffusion_quality(mesh, geometry);

        EXPECT_TRUE(q.internal_faces == 1);
        EXPECT_TRUE(q.non_orthogonal_faces == 0);
        EXPECT_TRUE(q.faces_above_60_deg == 0);
        EXPECT_TRUE(q.faces_above_75_deg == 0);
        EXPECT_TRUE(q.singular_over_relaxed_faces == 0);
        EXPECT_NEAR(q.max_non_orthogonality_deg, 0.0, 1e-12);
        EXPECT_NEAR(q.max_skewness, 0.0, 1e-12);
    });

    run_case("diffusion_quality_detects_real_mesh_nonorthogonality", []() {
        Mesh mesh = make_two_cell_unit_cubes();

        // Distort the shared face while keeping a conforming topology. The
        // geometry cache recomputes both cell centres and the face normal;
        // this is therefore a real mesh-quality perturbation, not a synthetic
        // change to Sf after geometry construction.
        mesh.points().set(2, 2.0, 0.35, 0.0);
        mesh.points().set(8, 2.0, 0.35, 1.0);

        const auto geometry = make_geometry_cache(mesh);
        const auto q = diagnose_diffusion_quality(mesh, geometry);

        EXPECT_TRUE(q.internal_faces == 1);
        EXPECT_TRUE(q.max_non_orthogonality_deg > 0.0);
        EXPECT_TRUE(std::isfinite(q.mean_non_orthogonality_deg));
        EXPECT_TRUE(std::isfinite(q.max_skewness));
        EXPECT_TRUE(q.singular_over_relaxed_faces == 0);
    });

    run_case("diffusion_quality_counts_over_relaxed_singularity", []() {
        Mesh mesh = make_two_cell_unit_cubes();
        auto geometry = make_geometry_cache(mesh);

        // Deliberately construct the singular algebraic geometry used by the
        // over-relaxed contract: Sf.d = 0. This test is only for diagnostics;
        // compute_laplacian separately rejects the same condition.
        geometry.face_Sf[5] = Vec3{0.0, 1.0, 0.0};
        const auto q = diagnose_diffusion_quality(mesh, geometry);

        EXPECT_TRUE(q.internal_faces == 1);
        EXPECT_TRUE(q.singular_over_relaxed_faces == 1);
    });

    run_case("diffusion_quality_rejects_invalid_threshold", []() {
        const Mesh mesh = make_two_cell_unit_cubes();
        const auto geometry = make_geometry_cache(mesh);
        EXPECT_THROW(
            diagnose_diffusion_quality(mesh, geometry, 0.0),
            std::invalid_argument);
        EXPECT_THROW(
            diagnose_diffusion_quality(mesh, geometry, -1.0),
            std::invalid_argument);
        EXPECT_THROW(
            diagnose_diffusion_quality(mesh, geometry, std::numeric_limits<double>::quiet_NaN()),
            std::invalid_argument);
    });

    return run_all();
}
