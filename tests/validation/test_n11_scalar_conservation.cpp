#include "cfdx/physics/finite_volume_transport.h"
#include "cfdx/physics/conservation_boundedness.h"
#include "common/test_harness.h"
#include <cmath>
#include <limits>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube_with_x_patches()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i = 0; i < 8; ++i)
        m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,3,2,1}); // x=0
    m.faces().push_face({4,5,6,7}); // x=1
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});

    m.ownership().resize(6);
    for (std::size_t f = 0; f < 6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});

    Patch x0;
    x0.name = "x0";
    x0.type = PatchType::WALL;
    x0.face_ids = {0};
    m.boundary().add_patch(x0);

    Patch x1;
    x1.name = "x1";
    x1.type = PatchType::WALL;
    x1.face_ids = {1};
    m.boundary().add_patch(x1);

    Patch walls;
    walls.name = "walls";
    walls.type = PatchType::WALL;
    walls.face_ids = {2,3,4,5};
    m.boundary().add_patch(walls);

    return m;
}

int main()
{
    run_case("n11_scalar_post_solve_balance_and_bounds", [] {
        const Mesh mesh = make_unit_cube_with_x_patches();
        const auto geometry = build_fv_geometry(mesh);

        Field<double,Location::FACE> mass_flux(mesh.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> source(mesh.n_cells(),"source","1/s",1);
        Field<double,Location::CELL> implicit_source(mesh.n_cells(),"Sp","1/s",1);
        Field<double,Location::CELL> scalar(mesh.n_cells(),"scalar","1",1);
        mass_flux.fill(0.0);
        source.fill(0.0);
        implicit_source.fill(0.0);
        scalar.fill(0.0);

        ScalarBoundaryConditions bc;
        bc["x0"] = {ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        bc["x1"] = {ScalarBoundaryType::FIXED_VALUE,1.0,0.0};
        bc["walls"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        const auto equation = assemble_scalar_equation(
            mesh, geometry, mass_flux, 1.0, source, implicit_source, bc,
            true, nullptr, nullptr, nullptr, nullptr,
            ConvectionScheme::UPWIND, &scalar);

        Vector solution(mesh.n_cells(), 0.0);
        const auto linear = solve_scalar_equation(
            equation, solution, ScalarSolveControls{200,1e-12,1.0});
        EXPECT_TRUE(linear.status == SolverStatus::CONVERGED);
        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
            scalar(cell) = solution(cell);
        EXPECT_TRUE(std::isfinite(scalar(0)));
        EXPECT_TRUE(scalar(0) >= 0.0 && scalar(0) <= 1.0);

        // Rebuild the physical scalar flux from the accepted post-solve field.
        // This is deliberately independent of the assembled matrix residual.
        const auto reconstructed = reconstruct_scalar_transport_flux(
            mesh, geometry, mass_flux, scalar, 1.0, bc, true,
            ConvectionScheme::UPWIND);

        const auto balance = audit_transport_balance(
            mesh, reconstructed, source);

        EXPECT_TRUE(balance.finite());
        EXPECT_NEAR(balance.residual, 0.0, 1e-12);
        EXPECT_TRUE(balance.normalized_residual < 1e-12);
        EXPECT_TRUE(balance.max_cell_residual < 1e-12);

        const auto bounds = audit_scalar_bounds(scalar, 0.0, 1.0, 0.0);
        EXPECT_TRUE(bounds.bounded());
        EXPECT_TRUE(bounds.nonfinite == 0);
        EXPECT_TRUE(bounds.below == 0);
        EXPECT_TRUE(bounds.above == 0);

        std::cout << "N11_SCALAR_QUALIFICATION: PASS"
                  << " scalar=" << scalar(0)
                  << " balance_residual=" << balance.residual
                  << " balance_normalized=" << balance.normalized_residual
                  << " min=" << bounds.minimum
                  << " max=" << bounds.maximum
                  << "\n";
    });

    run_case("n11_scalar_nonfinite_and_bound_violation_are_hard_failures", [] {
        Field<double,Location::CELL> scalar(3,"scalar","1",1);
        scalar(0) = -1e-3;
        scalar(1) = 0.5;
        scalar(2) = std::numeric_limits<double>::quiet_NaN();

        const auto bounds = audit_scalar_bounds(scalar, 0.0, 1.0, 0.0);
        EXPECT_TRUE(!bounds.bounded());
        EXPECT_TRUE(bounds.below == 1);
        EXPECT_TRUE(bounds.above == 0);
        EXPECT_TRUE(bounds.nonfinite == 1);
    });
}
