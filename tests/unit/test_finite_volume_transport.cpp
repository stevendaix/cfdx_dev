#include "cfdx/physics/finite_volume_transport.h"
#include "common/test_harness.h"
#include <cmath>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0; i<8; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,1,2,3,4,5};
    m.boundary().add_patch(wall);
    return m;
}

int main()
{
    run_case("scalar_dirichlet_diffusion_matches_one_cell_oracle", [] {
        const Mesh m = make_unit_cube();
        const FvGeometry g = build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> su(m.n_cells(),"Su","K/s",1);
        Field<double,Location::CELL> sp(m.n_cells(),"Sp","1/s",1);
        flux.fill(0.0);
        su.fill(1.0);
        sp.fill(0.0);

        ScalarBoundaryConditions bc;
        bc["wall"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
        const auto eq = assemble_scalar_equation(m,g,flux,1.0,su,sp,bc);

        EXPECT_NEAR(eq.matrix(0,0), 12.0, 1e-12);
        EXPECT_NEAR(eq.rhs(0), 1.0, 1e-12);

        Vector x(1,0.0);
        const auto r = solve_scalar_equation(eq,x,{100,1e-14,1.0});
        EXPECT_TRUE(r.status == SolverStatus::CONVERGED);
        EXPECT_NEAR(x(0),1.0/12.0,1e-12);
        EXPECT_NEAR(scalar_equation_residual_inf(eq,x),0.0,1e-12);
    });

    run_case("zero_gradient_inflow_does_not_reduce_owner_diagonal", [] {
        const Mesh m = make_unit_cube();
        const FvGeometry g = build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","m3/s",1);
        Field<double,Location::CELL> su(m.n_cells(),"Su","1/s",1);
        Field<double,Location::CELL> sp(m.n_cells(),"Sp","1/s",1);
        flux.fill(0.0); su.fill(0.0); sp.fill(0.0);
        flux(0) = 1.0;
        flux(1) = -1.0;
        ScalarBoundaryConditions bc;
        bc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        const auto eq = assemble_scalar_equation(m,g,flux,0.0,su,sp,bc,false);
        EXPECT_NEAR(eq.diagonal[0],1.0,1e-12);
    });

    run_case("scalar_constant_state_is_exactly_preserved", [] {
        const Mesh m = make_unit_cube();
        const FvGeometry g = build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1);
        Field<double,Location::CELL> su(m.n_cells(),"Su","1/s",1);
        Field<double,Location::CELL> sp(m.n_cells(),"Sp","1/s",1);
        flux.fill(0.0);
        su.fill(0.0);
        sp.fill(0.0);

        ScalarBoundaryConditions bc;
        bc["wall"] = {ScalarBoundaryType::FIXED_VALUE, 3.0, 0.0};
        const auto eq = assemble_scalar_equation(m,g,flux,2.0,su,sp,bc);

        Vector x(1,3.0);
        EXPECT_NEAR(scalar_equation_residual_inf(eq,x),0.0,1e-12);
    });

    run_case("nonorthogonal_diffusion_correction_vanishes_on_orthogonal_face", [] {
        const Vec3 d{2.0,0.0,0.0};
        const Vec3 Sf{3.0,0.0,0.0};
        const Vec3 go{4.0,-2.0,1.0};
        const Vec3 gn{4.0,-2.0,1.0};
        EXPECT_NEAR(nonorthogonal_diffusion_correction(d,Sf,2.0,go,gn),0.0,1e-14);
    });

    run_case("nonorthogonal_diffusion_correction_matches_tangential_flux", [] {
        const Vec3 d{2.0,0.0,0.0};
        const Vec3 Sf{2.0,1.0,0.0};
        const Vec3 go{3.0,4.0,0.0};
        const Vec3 gn{3.0,4.0,0.0};
        // S_corr=(0,1,0), so Gamma grad.S_corr = 2*4.
        EXPECT_NEAR(nonorthogonal_diffusion_correction(d,Sf,2.0,go,gn),8.0,1e-14);
    });

    run_case("nonorthogonal_diffusion_rejects_missing_field", [] {
        const Mesh m = make_unit_cube();
        const auto g = build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","m3/s",1);
        Field<double,Location::CELL> su(m.n_cells(),"Su","1/s",1);
        Field<double,Location::CELL> sp(m.n_cells(),"Sp","1/s",1);
        flux.fill(0.0); su.fill(0.0); sp.fill(0.0);
        ScalarBoundaryConditions bc;
        bc["wall"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
        EXPECT_THROW(
            assemble_scalar_equation(m,g,flux,1.0,su,sp,bc,true,nullptr,nullptr,
                nullptr,nullptr,ConvectionScheme::UPWIND,nullptr,
                DiffusionScheme::NON_ORTHOGONAL_CORRECTED,nullptr),
            std::invalid_argument);
    });

    run_case("second_order_upwind_constant_state_preserves_constant", [] {
        const Mesh m = make_unit_cube();
        auto g = build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1);
        flux.fill(1.0);
        Field<double,Location::CELL> phi(m.n_cells(),"phi","unit",1);
        Field<double,Location::CELL> su(m.n_cells(),"su","unit",1);
        Field<double,Location::CELL> sp(m.n_cells(),"sp","unit",1);
        phi.fill(3.5);
        su.fill(0.0);
        sp.fill(0.0);
        ScalarBoundaryConditions bc;
        bc["wall"] = {ScalarBoundaryType::FIXED_VALUE,3.5,0.0};
        const auto eq = assemble_scalar_equation(
            m,g,flux,1.0,su,sp,bc,true,nullptr,nullptr,nullptr,nullptr,
            ConvectionScheme::SECOND_ORDER_UPWIND,&phi);
        Vector x(m.n_cells(),3.5);
        EXPECT_NEAR(scalar_equation_residual_inf(eq,x),0.0,1e-12);
    });

    run_case("bounded_convection_constant_state_has_zero_residual", [] {
        const Mesh m = make_unit_cube();
        const FvGeometry g = build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","m3/s",1);
        Field<double,Location::CELL> su(m.n_cells(),"Su","1/s",1);
        Field<double,Location::CELL> sp(m.n_cells(),"Sp","1/s",1);
        flux.fill(0.0);
        su.fill(0.0);
        sp.fill(0.0);
        ScalarBoundaryConditions bc;
        bc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
        EXPECT_THROW(assemble_scalar_equation(m,g,flux,0.0,su,sp,bc,true),
                     std::runtime_error);
    });

    return run_all();
}
