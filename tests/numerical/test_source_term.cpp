// M0.7-T06 — Tests for Source Term Linearization

#include "cfdx/core/numerics/source_term.h"
#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/boundary.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::testing;

// Helper: build a simple 1D mesh with 4 cells
Mesh make_test_mesh() {
    Mesh m;

    // Points: 10 points for a thin 3D slab (y=0)
    // Bottom (z=0): 0,1,2,3,4
    // Top (z=0.1): 5,6,7,8,9
    m.points().resize(10);
    m.points().set(0, 0.0, 0.0, 0.0);
    m.points().set(1, 0.25, 0.0, 0.0);
    m.points().set(2, 0.5, 0.0, 0.0);
    m.points().set(3, 0.75, 0.0, 0.0);
    m.points().set(4, 1.0, 0.0, 0.0);
    m.points().set(5, 0.0, 0.0, 0.1);
    m.points().set(6, 0.25, 0.0, 0.1);
    m.points().set(7, 0.5, 0.0, 0.1);
    m.points().set(8, 0.75, 0.0, 0.1);
    m.points().set(9, 1.0, 0.0, 0.1);

    // 6 faces (all quads, CCW when viewed from outside):
    // Internal faces (normal in +x direction)
    m.faces().push_face({0, 5, 6, 1});   // 0: between cell 0-1 at x=0.125
    m.faces().push_face({1, 6, 7, 2});   // 1: between cell 1-2 at x=0.375
    m.faces().push_face({2, 7, 8, 3});   // 2: between cell 2-3 at x=0.625
    m.faces().push_face({3, 8, 9, 4});   // 3: between cell 3-boundary at x=0.875
    // Boundary faces
    m.faces().push_face({0, 1, 6, 5});   // 4: boundary at x=0 (normal -x)
    m.faces().push_face({4, 9, 8, 3});   // 5: boundary at x=1 (normal +x)

    // Owner/neighbour
    m.ownership().resize(6);
    m.ownership().set_owner(0, 0); m.ownership().set_neighbour(0, 1);
    m.ownership().set_owner(1, 1); m.ownership().set_neighbour(1, 2);
    m.ownership().set_owner(2, 2); m.ownership().set_neighbour(2, 3);
    m.ownership().set_owner(3, 3); m.ownership().set_neighbour(3, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(4, 0); m.ownership().set_neighbour(4, FaceOwnership::BOUNDARY);
    m.ownership().set_owner(5, 3); m.ownership().set_neighbour(5, FaceOwnership::BOUNDARY);

    // Cells: 4 hex cells
    // Cell 0: faces 0 (internal right), 4 (boundary left), plus 4 more for top/bottom/front/back
    // But we only have 6 faces total. Let's simplify - the source term only needs
    // cell volumes and centres, which compute_cell_geometry will compute from the faces.
    m.cells().push_cell({0, 4});       // cell 0: face 0 (internal) + face 4 (boundary)
    m.cells().push_cell({0, 1});       // cell 1: faces 0, 1
    m.cells().push_cell({1, 2});       // cell 2: faces 1, 2
    m.cells().push_cell({2, 3, 5});    // cell 3: faces 2, 3 (internal) + face 5 (boundary)

    // Boundary patches
    BoundaryPatches bp;
    bp.add_patch("inlet", PatchType::INLET);
    bp.patch(0).face_ids = {4};
    bp.add_patch("outlet", PatchType::OUTLET);
    bp.patch(1).face_ids = {5};
    m.set_boundary(bp);

    return m;
}

int main() {
    run_case("constant_source", [&]() {
        Mesh mesh = make_test_mesh();
        SourceTerm st = make_constant_source(mesh.n_cells(), 10.0, "const");

        EXPECT_NEAR(st.Su(0), 10.0, 1e-12);
        EXPECT_NEAR(st.Su(1), 10.0, 1e-12);
        EXPECT_NEAR(st.Su(2), 10.0, 1e-12);
        EXPECT_NEAR(st.Su(3), 10.0, 1e-12);
        EXPECT_NEAR(st.Sp(0), 0.0, 1e-12);
        EXPECT_NEAR(st.Sp(1), 0.0, 1e-12);
        EXPECT_NEAR(st.Sp(2), 0.0, 1e-12);
        EXPECT_NEAR(st.Sp(3), 0.0, 1e-12);

        // Apply to system
        std::vector<double> diag(4, 5.0);
        std::vector<double> rhs(4, 0.0);
        st.apply_to_system(diag, rhs);

        EXPECT_NEAR(diag[0], 5.0, 1e-12);  // Sp = 0
        EXPECT_NEAR(rhs[0], 10.0, 1e-12);  // Su = 10
    });

    run_case("linear_source", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(2.0);

        SourceTerm st = make_linear_source(phi, -3.0, "linear");  // S = -3*phi

        EXPECT_NEAR(st.Sp(0), -3.0, 1e-12);
        EXPECT_NEAR(st.Sp(1), -3.0, 1e-12);
        EXPECT_NEAR(st.Su(0), 0.0, 1e-12);

        // Apply to system: aP*phi = ... + Sp*phi + Su
        // New diag = aP - Sp = 5 - (-3) = 8
        // New rhs = Su = 0
        std::vector<double> diag(4, 5.0);
        std::vector<double> rhs(4, 0.0);
        st.apply_to_system(diag, rhs);

        EXPECT_NEAR(diag[0], 8.0, 1e-12);
        EXPECT_NEAR(rhs[0], 0.0, 1e-12);
    });

    run_case("decay_source", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(5.0);

        SourceTerm st = make_decay_source(phi, 2.0, "decay");  // S = -2*phi

        EXPECT_NEAR(st.Sp(0), -2.0, 1e-12);
        EXPECT_NEAR(st.Su(0), 0.0, 1e-12);

        std::vector<double> diag(4, 5.0);
        std::vector<double> rhs(4, 0.0);
        st.apply_to_system(diag, rhs);

        EXPECT_NEAR(diag[0], 7.0, 1e-12);  // 5 - (-2) = 7
        EXPECT_NEAR(rhs[0], 0.0, 1e-12);
    });

    run_case("exponential_source", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi(0) = 1.0; phi(1) = 2.0; phi(2) = 0.5; phi(3) = 0.0;

        SourceTerm st = make_exponential_source(phi, 2.0, 1.0, "exp");  // S = 2*exp(-phi)

        // At phi=1: S = 2*e^-1 ≈ 0.7358, dS/dphi = -2*e^-1 ≈ -0.7358
        // Su = f - Sp*phi = 0.7358 - (-0.7358)*1 = 1.4715
        // Sp = -0.7358
        double f1 = 2.0 * std::exp(-1.0);
        double df1 = -2.0 * std::exp(-1.0);
        double su1 = f1 - df1 * 1.0;
        double sp1 = df1;

        EXPECT_NEAR(st.Su(0), su1, 1e-10);
        EXPECT_NEAR(st.Sp(0), sp1, 1e-10);
        EXPECT_TRUE(st.Sp(0) < 0);  // Sp should be negative
    });

    run_case("power_source", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi(0) = 2.0; phi(1) = 3.0;

        SourceTerm st = make_power_source(phi, 2.0, 2.0, "power");  // S = 2*phi^2

        // At phi=2: S = 8, dS/dphi = 8
        // Su = 8 - 8*2 = -8, Sp = 8 (clamped to 0 since >0)
        EXPECT_NEAR(st.Sp(0), 0.0, 1e-12);  // Clamped to 0
        EXPECT_NEAR(st.Su(0), 8.0, 1e-12);  // f(phi) = 2*4 = 8

        // At phi=3: S = 18, dS/dphi = 12
        // Su = 18 - 12*3 = -18, Sp = 12 (clamped to 0)
        EXPECT_NEAR(st.Sp(1), 0.0, 1e-12);
        EXPECT_NEAR(st.Su(1), 18.0, 1e-12);
    });

    run_case("linearize_source_general", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi(0) = 1.0; phi(1) = 2.0; phi(2) = 3.0; phi(3) = 4.0;

        // f(phi) = phi^3, df/dphi = 3*phi^2
        SourceTerm st = linearize_source(phi,
            [](double p) { return p * p * p; },
            [](double p) { return 3.0 * p * p; },
            "cubic");

        // At phi=1: f=1, df=3 → Sp=3 (clamped to 0), Su = 1 - 0*1 = 1
        EXPECT_NEAR(st.Sp(0), 0.0, 1e-12);
        EXPECT_NEAR(st.Su(0), 1.0, 1e-12);

        // At phi=2: f=8, df=12 → Sp=12 (clamped to 0), Su = 8
        EXPECT_NEAR(st.Sp(1), 0.0, 1e-12);
        EXPECT_NEAR(st.Su(1), 8.0, 1e-12);
    });

    run_case("update_linearization", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(1.0);

        // f(phi) = -phi^2, df/dphi = -2*phi
        SourceTerm st = linearize_source(phi,
            [](double p) { return -p * p; },
            [](double p) { return -2.0 * p; },
            "quad");

        // Initial: phi=1, f=-1, df=-2 → Sp=-2, Su = -1 - (-2)*1 = 1
        EXPECT_NEAR(st.Sp(0), -2.0, 1e-12);
        EXPECT_NEAR(st.Su(0), 1.0, 1e-12);

        // Update to phi=2
        phi.fill(2.0);
        update_linearization(phi, st,
            [](double p) { return -p * p; },
            [](double p) { return -2.0 * p; });

        // At phi=2: f=-4, df=-4 → Sp=-4, Su = -4 - (-4)*2 = 4
        EXPECT_NEAR(st.Sp(0), -4.0, 1e-12);
        EXPECT_NEAR(st.Su(0), 4.0, 1e-12);
    });

    run_case("spatial_source", [&]() {
        Mesh mesh = make_test_mesh();

        // Su = x (varies with position)
        // Sp = -1 (constant negative)
        SourceTerm st = make_spatial_source(mesh,
            [](const Vec3& pos) { return pos.x; },
            [](const Vec3& /*pos*/, double /*phi*/) { return -1.0; },
            "spatial");

        // Just verify the function runs and Sp is correctly clamped
        EXPECT_NEAR(st.Sp(0), -1.0, 1e-12);
        EXPECT_NEAR(st.Sp(1), -1.0, 1e-12);
        EXPECT_NEAR(st.Sp(2), -1.0, 1e-12);
        EXPECT_NEAR(st.Sp(3), -1.0, 1e-12);
        // Su values depend on mesh geometry - just verify they're set
        EXPECT_TRUE(st.Su(0) >= 0.0);
        EXPECT_TRUE(st.Su(1) >= 0.0);
        EXPECT_TRUE(st.Su(2) >= 0.0);
        EXPECT_TRUE(st.Su(3) >= 0.0);
    });

    run_case("apply_to_cell", [&]() {
        Mesh mesh = make_test_mesh();
        Field<double, Location::CELL> phi(4, "phi", "1", 1);
        phi.fill(2.0);

        SourceTerm st = make_decay_source(phi, 1.5, "decay");

        double diag = 10.0;
        double rhs = 0.0;
        st.apply_to_cell(0, diag, rhs);

        EXPECT_NEAR(diag, 11.5, 1e-12);  // 10 - (-1.5) = 11.5
        EXPECT_NEAR(rhs, 0.0, 1e-12);
    });

    run_case("source_term_copy_constructor", [&]() {
        Mesh mesh = make_test_mesh();
        SourceTerm st1 = make_constant_source(mesh.n_cells(), 5.0, "test");

        SourceTerm st2 = st1;  // Copy
        EXPECT_NEAR(st2.Su(0), 5.0, 1e-12);
        EXPECT_NEAR(st2.Sp(0), 0.0, 1e-12);

        // Modify original
        st1.Su.fill(99.0);
        EXPECT_NEAR(st2.Su(0), 5.0, 1e-12);  // Independent copy
    });

    return run_all();
}