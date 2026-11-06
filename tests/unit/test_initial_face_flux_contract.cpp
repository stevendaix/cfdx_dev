#include "cfdx/physics/initial_face_flux_contract.h"
#include "common/test_harness.h"

#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_cycle()
{
    Mesh m;
    m.points().resize(3);
    // Only topology is used by this algebraic flux certificate. Four internal
    // faces form a circulation with zero boundary throughflow.
    for (std::size_t f = 0; f < 4; ++f)
        m.faces().push_face({0, 1, 2});
    m.ownership().resize(4);
    for (std::size_t f = 0; f < 4; ++f) {
        m.ownership().set_owner(f, f);
        m.ownership().set_neighbour(f, static_cast<std::int64_t>((f + 1) % 4));
    }
    for (std::size_t c = 0; c < 4; ++c)
        m.cells().push_cell({c, (c + 3) % 4});
    return m;
}

int main()
{
    run_case("closed_cavity_circulation_has_nonzero_internal_flux_scale", [] {
        const auto mesh = make_cycle();
        const auto report = initial_face_flux_continuity_report(mesh,
            {0.2, 0.2, 0.2, 0.2});
        EXPECT_NEAR(report.face_flux_l1, 0.8, 1e-14);
        EXPECT_NEAR(report.local_l1, 0.0, 0.0);
        EXPECT_NEAR(report.normalized, 0.0, 0.0);
    });
    run_case("closed_cavity_still_rejects_local_imbalance", [] {
        const auto mesh = make_cycle();
        const auto report = initial_face_flux_continuity_report(mesh,
            {0.2, 0.2, 0.3, 0.2});
        EXPECT_NEAR(report.local_l1, 0.2, 1e-14);
        EXPECT_TRUE(report.normalized > 1e-7);
    });
    run_case("motionless_domain_has_exact_zero_continuity", [] {
        const auto mesh = make_cycle();
        const auto report = initial_face_flux_continuity_report(mesh,
            {0.0, 0.0, 0.0, 0.0});
        EXPECT_NEAR(report.face_flux_l1, 0.0, 0.0);
        EXPECT_NEAR(report.normalized, 0.0, 0.0);
    });
    run_case("flow_certificate_is_invariant_to_flux_scale", [] {
        const auto mesh = make_cycle();
        const auto small = initial_face_flux_continuity_report(mesh,
            {1e-40, 1e-40, 1.2e-40, 1e-40});
        const auto large = initial_face_flux_continuity_report(mesh,
            {1e10, 1e10, 1.2e10, 1e10});
        EXPECT_NEAR(small.normalized, large.normalized, 1e-15);
    });
    run_case("open_channel_uses_the_same_local_conservation_certificate", [] {
        Mesh mesh;
        mesh.points().resize(3);
        for (std::size_t f = 0; f < 3; ++f)
            mesh.faces().push_face({0, 1, 2});
        mesh.ownership().resize(3);
        mesh.ownership().set_owner(0, 0);
        mesh.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
        mesh.ownership().set_owner(1, 0);
        mesh.ownership().set_neighbour(1, 1);
        mesh.ownership().set_owner(2, 1);
        mesh.ownership().set_neighbour(2, FaceOwnership::BOUNDARY);
        mesh.cells().push_cell({0, 1});
        mesh.cells().push_cell({1, 2});
        const auto balanced = initial_face_flux_continuity_report(
            mesh, {-2.0, 2.0, 2.0}, 2.0);
        EXPECT_NEAR(balanced.local_l1, 0.0, 0.0);
        EXPECT_NEAR(balanced.normalized, 0.0, 0.0);
        const auto unbalanced = initial_face_flux_continuity_report(
            mesh, {-2.0, 2.0, 1.0}, 2.0);
        EXPECT_NEAR(unbalanced.normalized, 0.5, 1e-15);
        EXPECT_TRUE(unbalanced.normalized > 1e-7);
    });
    run_case("invalid_face_flux_cannot_pass_certificate", [] {
        const auto mesh = make_cycle();
        EXPECT_THROW(initial_face_flux_continuity_report(mesh,{0.0}),
                     std::invalid_argument);
        EXPECT_THROW(initial_face_flux_continuity_report(mesh,
            {0.2, 0.2, std::numeric_limits<double>::quiet_NaN(), 0.2}),
            std::invalid_argument);
        EXPECT_THROW(initial_face_flux_continuity_report(mesh,
            {0.2, 0.2, 0.2, 0.2}, -1.0), std::invalid_argument);
    });
    // run_case() only registers deferred tests: execute all six and propagate
    // their failures to CTest. A plain return 0 was a false-positive pass.
    return run_all();
}
