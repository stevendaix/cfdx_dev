#include "cfdx/core/numerics/conservation.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

static Mesh two_cell_internal_mesh() {
    Mesh m;
    m.points().resize(3);
    m.faces().push_face({0, 1, 2});
    m.ownership().resize(1);
    m.ownership().set_owner(0, 0);
    m.ownership().set_neighbour(0, 1);
    m.cells().push_cell({0});
    m.cells().push_cell({0});
    return m;
}

static Mesh one_cell_boundary_mesh() {
    Mesh m;
    m.points().resize(3);
    m.faces().push_face({0, 1, 2});
    m.ownership().resize(1);
    m.ownership().set_owner(0, 0);
    m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
    m.cells().push_cell({0});
    return m;
}

int main() {
    run_case("independent_owner_neighbour_contributions_cancel", [] {
        const Mesh m = two_cell_internal_mesh();
        const std::vector<double> owner{2.5};
        const std::vector<double> neighbour{-2.5};
        const auto r = audit_face_assembly(m, owner, neighbour);
        EXPECT_TRUE(r.finite());
        EXPECT_TRUE(r.antisymmetric());
        EXPECT_TRUE(r.internal_faces == 1);
        EXPECT_TRUE(r.mismatched_faces == 0);
        EXPECT_NEAR(r.max_abs_mismatch, 0.0, 0.0);
    });

    run_case("assembly_audit_detects_independent_flux_mismatch", [] {
        const Mesh m = two_cell_internal_mesh();
        const std::vector<double> owner{2.5};
        const std::vector<double> neighbour{-2.4};
        const auto r = audit_face_assembly(m, owner, neighbour);
        EXPECT_TRUE(!r.antisymmetric());
        EXPECT_TRUE(r.mismatched_faces == 1);
        EXPECT_TRUE(r.worst_face == 0);
        EXPECT_NEAR(r.max_abs_mismatch, 0.1, 1e-15);
    });

    run_case("assembly_audit_checks_boundary_neighbour_contribution", [] {
        const Mesh m = one_cell_boundary_mesh();
        const std::vector<double> owner{1.0};
        const std::vector<double> neighbour{0.0};
        const auto r = audit_face_assembly(m, owner, neighbour);
        EXPECT_TRUE(r.antisymmetric());
        EXPECT_TRUE(r.boundary_faces == 1);

        const auto bad = audit_face_assembly(m, owner, std::vector<double>{1e-6});
        EXPECT_TRUE(!bad.antisymmetric());
        EXPECT_TRUE(bad.mismatched_faces == 1);
    });

    run_case("assembly_audit_rejects_nonfinite_contributions", [] {
        const Mesh m = two_cell_internal_mesh();
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const auto r = audit_face_assembly(m, {nan}, {-1.0});
        EXPECT_TRUE(!r.finite());
        EXPECT_TRUE(!r.antisymmetric());
        EXPECT_TRUE(r.nonfinite_faces == 1);
    });

    run_case("independent_cell_balance_reconstruction", [] {
        const Mesh m = two_cell_internal_mesh();
        Field<double, Location::FACE> flux(1, "phi", "kg/s", 1);
        flux(0) = 3.0;
        const auto balance = reconstruct_cell_balance(m, flux);
        EXPECT_TRUE(balance.size() == 2);
        EXPECT_NEAR(balance[0], 3.0, 0.0);
        EXPECT_NEAR(balance[1], -3.0, 0.0);
        EXPECT_NEAR(balance[0] + balance[1], 0.0, 0.0);
    });

    return run_all();
}
