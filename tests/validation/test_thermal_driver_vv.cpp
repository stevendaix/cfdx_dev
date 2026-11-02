// Thermal-only autonomous driver VV (issue #1150, A2+A3).
//
// Every physical case is verified against an exact discrete oracle on an
// orthogonal mesh (the assembled system is small enough to solve by hand)
// plus the independent energy balance audit. Fail-closed batteries verify
// the explicit diagnostics of the deferred couplings, of the boundary
// mapping/value contract, and of the prescribed-mass-flux face field
// contract. Method sheet: docs/audits/17-thermal-driver.md.
#include "cfdx/physics/thermal_driver.h"
#include "common/test_harness.h"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

Mesh two_cell_1d()
{
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{0,1,0},{0,1,1},{0,0,1},
        {0.5,0,0},{0.5,1,0},{0.5,1,1},{0.5,0,1},
        {1,0,0},{1,1,0},{1,1,1},{1,0,1}
    };
    for (std::size_t i=0; i<12; ++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.faces().push_face({8,9,10,11});
    m.faces().push_face({4,8,9,5});
    m.faces().push_face({7,6,10,11});
    m.faces().push_face({4,7,11,8});
    m.faces().push_face({5,9,10,6});

    m.ownership().resize(11);
    m.ownership().set_owner(0,0); m.ownership().set_neighbour(0,FaceOwnership::BOUNDARY);
    m.ownership().set_owner(1,0); m.ownership().set_neighbour(1,1);
    for (std::size_t f=2; f<=5; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(6,1); m.ownership().set_neighbour(6,FaceOwnership::BOUNDARY);
    for (std::size_t f=7; f<=10; ++f) {
        m.ownership().set_owner(f,1);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }

    m.cells().push_cell({0,1,2,3,4,5});
    m.cells().push_cell({1,6,7,8,9,10});

    Patch left; left.name="left"; left.type=PatchType::WALL; left.face_ids={0};
    Patch right; right.name="right"; right.type=PatchType::WALL; right.face_ids={6};
    Patch walls; walls.name="walls"; walls.type=PatchType::WALL;
    walls.face_ids={2,3,4,5,7,8,9,10};
    m.boundary().add_patch(left);
    m.boundary().add_patch(right);
    m.boundary().add_patch(walls);
    return m;
}

struct CubePatchSpec {
    std::string name;
    std::vector<std::size_t> faces;
    PatchType type = PatchType::WALL;
};

Mesh one_cell_cube(const std::vector<CubePatchSpec>& specs)
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    for (std::size_t i=0;i<8;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3});
    m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for (std::size_t f=0;f<6;++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    for (const auto& spec : specs) {
        Patch patch;
        patch.name = spec.name;
        patch.type = spec.type;
        patch.face_ids.assign(spec.faces.begin(), spec.faces.end());
        m.boundary().add_patch(patch);
    }
    return m;
}

Mesh one_cell_cube()
{
    return one_cell_cube({CubePatchSpec{"walls", {0,1,2,3,4,5}}});
}

// Unit solid material: k=1 W/(m.K), cp=1 J/(kg.K), rho=1 kg/m^3, constant
// laws, admissible range [1,5000] K (ThermalMaterial defaults).
ThermalCase solid_case()
{
    ThermalCase c;
    c.name = "thermal_driver_vv";
    ThermalMaterial material;
    material.name = "unit_solid";
    material.properties.conductivity.reference_value = 1.0;
    material.properties.heat_capacity.reference_value = 1.0;
    material.properties.density.rho0 = 1.0;
    c.materials.push_back(material);
    ThermalRegion region;
    region.name = "solid";
    region.kind = ThermalRegionKind::SOLID;
    region.mesh_id = "mesh_solid";
    region.material = "unit_solid";
    region.initial_temperature = 300.0;
    c.regions.push_back(region);
    return c;
}

ThermalBoundaryCondition fixed_temperature(double value)
{
    ThermalBoundaryCondition bc;
    bc.kind = ThermalBoundaryKind::FIXED_TEMPERATURE;
    bc.value = value;
    return bc;
}

ThermalBoundaryCondition convection(double h, double ambient)
{
    ThermalBoundaryCondition bc;
    bc.kind = ThermalBoundaryKind::CONVECTION;
    bc.heat_transfer_coefficient = h;
    bc.ambient_temperature = ambient;
    return bc;
}

// Unit fluid region driven by a prescribed face mass-flux field (A3):
// same unit material as solid_case(), FLUID + PRESCRIBED_MASS_FLUX, the
// face field named by the caller is supplied at solve time.
ThermalCase fluid_flux_case(const std::string& field_name)
{
    ThermalCase c = solid_case();
    c.regions[0].kind = ThermalRegionKind::FLUID;
    c.regions[0].fluid_model = ThermalFluidModel::PRESCRIBED_MASS_FLUX;
    c.regions[0].mass_flux_field = field_name;
    c.regions[0].name = "water";
    return c;
}

} // namespace

int main()
{
    // Steady 1-D conduction, Dirichlet 400/300 K, adiabatic walls: the
    // discrete solution is the exact linear profile (375, 325) K and the
    // constant-k system needs a single pass.
    run_case("steady_conduction_linear_profile", [] {
        Mesh m = two_cell_1d();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(2,"T","K",1), source(2,"source","W/m3",1);
        T(0)=350.0; T(1)=350.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["left"]   = fixed_temperature(400.0);
        bc["right"]  = fixed_temperature(300.0);
        bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        const auto result =
            solve_thermal_region(c,"solid",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),375.0,1e-10);
        EXPECT_NEAR(T(1),325.0,1e-10);
        EXPECT_TRUE(result.outer_iterations == 1u);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
    });

    // Steady uniform generation, all walls at 310 K: six unit-area faces at
    // d=0.5 give 6*k*(T-310)/0.5 = q'''*V, hence T = 310 + 100/12.
    run_case("steady_uniform_generation", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source(0)=100.0;

        ThermalBoundaryConditions bc;
        bc["walls"] = fixed_temperature(310.0);

        const auto result =
            solve_thermal_region(c,"solid",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),310.0 + 100.0/12.0,1e-10);
        EXPECT_TRUE(result.outer_iterations == 1u);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
        EXPECT_NEAR(result.balance.source_integral,100.0,1e-9);
    });

    // Steady fixed heat flux 150 W/m^2 into one face, five walls at 300 K:
    // 150 = 5*2*(T-300) hence T = 315 K.
    run_case("steady_fixed_heat_flux", [] {
        Mesh m = one_cell_cube({CubePatchSpec{"heated",{0}},
                                CubePatchSpec{"walls",{1,2,3,4,5}}});
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["heated"].kind = ThermalBoundaryKind::FIXED_HEAT_FLUX;
        bc["heated"].value = 150.0;
        bc["walls"] = fixed_temperature(300.0);

        const auto result =
            solve_thermal_region(c,"solid",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),315.0,1e-10);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
    });

    // Steady POWER patch 225 W spread uniformly over a unit-area patch,
    // five walls at 300 K: 225 = 10*(T-300) hence T = 322.5 K.
    run_case("steady_power_patch", [] {
        Mesh m = one_cell_cube({CubePatchSpec{"heated",{0}},
                                CubePatchSpec{"walls",{1,2,3,4,5}}});
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["heated"].kind = ThermalBoundaryKind::POWER;
        bc["heated"].value = 225.0;
        bc["walls"] = fixed_temperature(300.0);

        const auto result =
            solve_thermal_region(c,"solid",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),322.5,1e-10);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
    });

    // Convection equilibrium: h=1, T_amb=350 K, no source, initial state at
    // equilibrium: the film fixed point leaves the region at 350 K.
    run_case("steady_convection_equilibrium", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=350.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["walls"] = convection(1.0,350.0);

        const auto result =
            solve_thermal_region(c,"solid",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),350.0,1e-10);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
    });

    // Robin fixed point with generation: h=1, T_amb=300 K, q'''=100. Each
    // face is a Dirichlet at T_b=(T_P+0.5*300)/1.5, so the discrete fixed
    // point solves 12*T - 2*6*T_b = 100 and gives exactly T = 325 K.
    run_case("steady_convection_robin_fixed_point", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source(0)=100.0;

        ThermalBoundaryConditions bc;
        bc["walls"] = convection(1.0,300.0);

        const auto result =
            solve_thermal_region(c,"solid",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),325.0,1e-9);
        EXPECT_TRUE(result.outer_iterations >= 2u);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
        EXPECT_NEAR(result.balance.source_integral,100.0,1e-9);
    });

    // Transient insulated heating: all faces adiabatic, q'''=100 W/m3,
    // rho*cp=1000, dt=10 s. Each driver call is one backward-Euler step and
    // the independent audit must see the full generation accumulate:
    // 100 W and +1 K per step.
    run_case("transient_insulated_heating", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        c.materials[0].properties.heat_capacity.reference_value = 1000.0;
        c.time_mode = ThermalTimeMode::TRANSIENT;
        c.time_step = 10.0;
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source(0)=100.0;

        ThermalBoundaryConditions bc;
        bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        for (int step=0; step<4; ++step) {
            const auto result =
                solve_thermal_region(c,"solid",m,g,bc,T,source);
            EXPECT_TRUE(result.converged);
            EXPECT_NEAR(result.balance.accumulation,100.0,1e-9);
            EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
        }
        EXPECT_NEAR(T(0),304.0,1e-9);
    });

    // Transient Dirichlet relaxation: six walls at 310 K, rho*cp=24,
    // dt=1 s. The backward-Euler recurrence is
    // T^{n+1} = (24*T^n + 12*310)/36, closed form 310 - 10*(2/3)^n.
    run_case("transient_dirichlet_relaxation", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        c.materials[0].properties.heat_capacity.reference_value = 24.0;
        c.time_mode = ThermalTimeMode::TRANSIENT;
        c.time_step = 1.0;
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["walls"] = fixed_temperature(310.0);

        double expected = 300.0;
        for (int step=0; step<5; ++step) {
            const auto result =
                solve_thermal_region(c,"solid",m,g,bc,T,source);
            EXPECT_TRUE(result.converged);
            expected = (24.0*expected + 12.0*310.0)/36.0;
            EXPECT_NEAR(T(0),expected,1e-10);
            EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
        }
        EXPECT_NEAR(T(0),310.0 - 10.0*std::pow(2.0/3.0,5.0),1e-9);
    });

    // Stagnant fluid diffusion: FLUID/STAGNANT_DIFFUSION region with k=0.6,
    // walls at 310 K, q'''=50 W/m3. No momentum model is initialized: the
    // mass flux is zero by contract. Six faces of conductance k/d = 1.2
    // give T = 310 + 50/7.2.
    run_case("stagnant_fluid_diffusion", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        c.materials[0].properties.conductivity.reference_value = 0.6;
        c.regions[0].kind = ThermalRegionKind::FLUID;
        c.regions[0].fluid_model = ThermalFluidModel::STAGNANT_DIFFUSION;
        c.regions[0].name = "air";
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source(0)=50.0;

        ThermalBoundaryConditions bc;
        bc["walls"] = fixed_temperature(310.0);

        const auto result =
            solve_thermal_region(c,"air",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),310.0 + 50.0/7.2,1e-10);
        EXPECT_TRUE(result.outer_iterations == 1u);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
    });

    // Nonlinear conductivity Picard: k(T) = 1 + 0.001*(T-300), walls at
    // 300 K, q'''=100. The discrete fixed point solves
    // 12*k(T)*(T-300) = 100, i.e. 0.012*u^2 + 12*u - 100 = 0 with u=T-300.
    run_case("conductivity_picard", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        c.materials[0].properties.conductivity.model =
            ScalarPropertyModel::LINEAR;
        c.materials[0].properties.conductivity.reference_value = 1.0;
        c.materials[0].properties.conductivity.reference_temperature = 300.0;
        c.materials[0].properties.conductivity.coefficients = {0.001};
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source(0)=100.0;

        ThermalBoundaryConditions bc;
        bc["walls"] = fixed_temperature(300.0);

        const auto result =
            solve_thermal_region(c,"solid",m,g,bc,T,source);
        EXPECT_TRUE(result.converged);
        EXPECT_TRUE(result.outer_iterations >= 2u);
        const double expected = 300.0 + (std::sqrt(148.8) - 12.0)/0.024;
        EXPECT_NEAR(T(0),expected,1e-9);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
    });

    // Steady advection of a constant state (A3): two cells, phi=-1 on the
    // left boundary (inflow), phi=+1 on the internal face, phi=+1 on the
    // right boundary (outflow), all other faces zero. The left wall is at
    // 320 K, all other walls adiabatic. Advection of a uniform field is
    // exact in a single pass: the inflow advects 320 K and flushes the
    // initial 350 K state, so T0 = T1 = 320 K exactly with no conduction.
    run_case("steady_advection_constant_state", [] {
        Mesh m = two_cell_1d();
        auto g = build_fv_geometry(m);
        ThermalCase c = fluid_flux_case("feed_flux");
        Field<double,Location::CELL> T(2,"T","K",1), source(2,"source","W/m3",1);
        T(0)=350.0; T(1)=350.0; source.fill(0.0);

        ThermalFaceFields fields;
        fields.emplace("feed_flux",
            Field<double,Location::FACE>(m.n_faces(),"feed_flux","kg/s",1));
        auto& phi = fields.at("feed_flux");
        phi.fill(0.0);
        phi(0) = -1.0;
        phi(1) = +1.0;
        phi(6) = +1.0;

        ThermalBoundaryConditions bc;
        bc["left"] = fixed_temperature(320.0);
        bc["right"].kind = ThermalBoundaryKind::ADIABATIC;
        bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        const auto result =
            solve_thermal_region(c,"water",m,g,bc,T,source,{},&fields);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),320.0,1e-10);
        EXPECT_NEAR(T(1),320.0,1e-10);
        EXPECT_TRUE(result.outer_iterations == 1u);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
    });

    // Steady throughflow with generation (A3): one unit cube, inlet at
    // 300 K with phi=-2, outlet phi=+2 (zero-gradient), walls adiabatic,
    // q'''=100 W/m3. With m*cp=2 and outlet conductance D=k*A/d=2, the
    // discrete cell equation is 4*T = 4*300 + 100, hence T = 325 K exactly.
    run_case("steady_throughflow_generation", [] {
        Mesh m = one_cell_cube({CubePatchSpec{"inlet",{0}},
                                CubePatchSpec{"outlet",{1}},
                                CubePatchSpec{"walls",{2,3,4,5}}});
        auto g = build_fv_geometry(m);
        ThermalCase c = fluid_flux_case("feed_flux");
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source(0)=100.0;

        ThermalFaceFields fields;
        fields.emplace("feed_flux",
            Field<double,Location::FACE>(m.n_faces(),"feed_flux","kg/s",1));
        auto& phi = fields.at("feed_flux");
        phi.fill(0.0);
        phi(0) = -2.0;
        phi(1) = +2.0;

        ThermalBoundaryConditions bc;
        bc["inlet"] = fixed_temperature(300.0);
        bc["outlet"].kind = ThermalBoundaryKind::ADIABATIC;
        bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        const auto result =
            solve_thermal_region(c,"water",m,g,bc,T,source,{},&fields);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T(0),325.0,1e-10);
        EXPECT_TRUE(result.outer_iterations == 1u);
        EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
        EXPECT_NEAR(result.balance.source_integral,100.0,1e-9);
    });

    // Transient washout (A3): same cube and flux as the steady throughflow,
    // no source, TRANSIENT dt=1 s, rho=10, cp=1. Backward Euler gives
    // (a+b)*T^{n+1} = a*T^n + b*300 with a=rho*cp*V/dt=10 and b=m*cp+D=4,
    // so T^{n+1}-300 = (5/7)*(T^n-300): closed form 300+50*(5/7)^n from
    // T^0=350, verified over four steps with the audit at every step.
    run_case("transient_washout_exponential", [] {
        Mesh m = one_cell_cube({CubePatchSpec{"inlet",{0}},
                                CubePatchSpec{"outlet",{1}},
                                CubePatchSpec{"walls",{2,3,4,5}}});
        auto g = build_fv_geometry(m);
        ThermalCase c = fluid_flux_case("feed_flux");
        c.materials[0].properties.density.rho0 = 10.0;
        c.time_mode = ThermalTimeMode::TRANSIENT;
        c.time_step = 1.0;
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=350.0; source.fill(0.0);

        ThermalFaceFields fields;
        fields.emplace("feed_flux",
            Field<double,Location::FACE>(m.n_faces(),"feed_flux","kg/s",1));
        auto& phi = fields.at("feed_flux");
        phi.fill(0.0);
        phi(0) = -2.0;
        phi(1) = +2.0;

        ThermalBoundaryConditions bc;
        bc["inlet"] = fixed_temperature(300.0);
        bc["outlet"].kind = ThermalBoundaryKind::ADIABATIC;
        bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        double expected = 350.0;
        for (int step=0; step<4; ++step) {
            const auto result =
                solve_thermal_region(c,"water",m,g,bc,T,source,{},&fields);
            EXPECT_TRUE(result.converged);
            expected = 300.0 + (expected-300.0)*(5.0/7.0);
            EXPECT_NEAR(T(0),expected,1e-10);
            EXPECT_TRUE(result.balance.normalized_residual < 1e-9);
        }
        EXPECT_NEAR(T(0),300.0 + 50.0*std::pow(5.0/7.0,4.0),1e-9);
    });

    // Fail-closed battery: couplings explicitly out of scope.
    run_case("fail_closed_coupled_and_unsupported", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source.fill(0.0);
        ThermalBoundaryConditions bc;
        bc["walls"] = fixed_temperature(310.0);

        EXPECT_THROW_WITH(
            solve_thermal_region(c,"ghost",m,g,bc,T,source),
            std::invalid_argument, "unknown thermal region");

        c.regions[0].material = "ghost_material";
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",m,g,bc,T,source),
            std::invalid_argument, "references unknown material");
        c.regions[0].material = "unit_solid";

        // Conjugate interfaces are A4.
        ThermalCase interface_case = solid_case();
        ThermalRegion air = interface_case.regions[0];
        air.name = "air";
        air.kind = ThermalRegionKind::FLUID;
        air.fluid_model = ThermalFluidModel::STAGNANT_DIFFUSION;
        interface_case.regions.push_back(air);
        ThermalInterface interface;
        interface.name = "wall_air_cht";
        interface.region_a = "solid";
        interface.patch_a = "walls";
        interface.region_b = "air";
        interface.patch_b = "walls";
        interface_case.interfaces.push_back(interface);
        interface_case.validate();
        EXPECT_THROW_WITH(
            solve_thermal_region(interface_case,"solid",m,g,bc,T,source),
            std::invalid_argument, "deferred to A4");

        // Thermal links to lumped nodes are A5.
        ThermalCase link_case = solid_case();
        ThermalNode node;
        node.name = "reservoir";
        node.material = "unit_solid";
        node.volume = 0.1;
        link_case.nodes.push_back(node);
        ThermalLink link;
        link.name = "wall_to_reservoir";
        link.source.kind = ThermalEndpointKind::REGION_PATCH;
        link.source.object = "solid";
        link.source.patch = "walls";
        link.target.kind = ThermalEndpointKind::LUMPED_NODE;
        link.target.object = "reservoir";
        link.conductance = 1.0;
        link_case.links.push_back(link);
        link_case.validate();
        EXPECT_THROW_WITH(
            solve_thermal_region(link_case,"solid",m,g,bc,T,source),
            std::invalid_argument, "deferred to A5");

        // Anisotropic tensors need the tensor diffusion operator.
        ThermalCase aniso_case = solid_case();
        aniso_case.materials[0].anisotropic_conductivity = true;
        aniso_case.materials[0].conductivity_shape =
            {{1.5, 1.0, 0.5, 0.0, 0.0, 0.0}};
        aniso_case.materials[0].validate();
        EXPECT_THROW_WITH(
            solve_thermal_region(aniso_case,"solid",m,g,bc,T,source),
            std::invalid_argument, "anisotropic");
    });

    // Fail-closed battery: boundary mapping must be total and unambiguous.
    run_case("fail_closed_invalid_boundary_mapping", [] {
        Mesh m = one_cell_cube({CubePatchSpec{"heated",{0}},
                                CubePatchSpec{"walls",{1,2,3,4,5}}});
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["walls"] = fixed_temperature(310.0);
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",m,g,bc,T,source),
            std::invalid_argument, "missing thermal boundary condition for patch 'heated'");

        bc["heated"] = fixed_temperature(310.0);
        bc["ghost"] = fixed_temperature(310.0);
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",m,g,bc,T,source),
            std::invalid_argument, "unknown patch");

        Mesh empty_m = one_cell_cube(
            {CubePatchSpec{"openings",{0,1,2,3,4,5},PatchType::EMPTY}});
        auto empty_g = build_fv_geometry(empty_m);
        ThermalBoundaryConditions empty_bc;
        empty_bc["openings"].kind = ThermalBoundaryKind::ADIABATIC;
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",empty_m,empty_g,empty_bc,T,source),
            std::invalid_argument, "EMPTY patch");
    });

    // Fail-closed battery: boundary values must respect the material range
    // and the typed condition contract; the entry state too.
    run_case("fail_closed_invalid_boundary_values", [] {
        Mesh m = one_cell_cube();
        auto g = build_fv_geometry(m);
        ThermalCase c = solid_case();
        c.materials[0].minimum_temperature = 250.0;
        c.materials[0].maximum_temperature = 500.0;
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["walls"] = fixed_temperature(600.0);
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",m,g,bc,T,source),
            std::invalid_argument, "outside the material range");

        bc["walls"] = convection(0.0,300.0);
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",m,g,bc,T,source),
            std::invalid_argument, "positive heat transfer coefficient");

        bc["walls"].kind = ThermalBoundaryKind::FIXED_HEAT_FLUX;
        bc["walls"].value = std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",m,g,bc,T,source),
            std::invalid_argument, "non-finite");

        bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;
        T(0)=600.0;
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"solid",m,g,bc,T,source),
            std::invalid_argument, "initial temperature is outside the material range");
    });

    // Fail-closed battery (A3): the prescribed-mass-flux contract must be
    // supplied, resolved, well-formed and discretely conservative, and an
    // inflow boundary needs a thermal state to advect in.
    run_case("fail_closed_prescribed_mass_flux", [] {
        Mesh m = one_cell_cube({CubePatchSpec{"inlet",{0}},
                                CubePatchSpec{"outlet",{1}},
                                CubePatchSpec{"walls",{2,3,4,5}}});
        auto g = build_fv_geometry(m);
        ThermalCase c = fluid_flux_case("feed_flux");
        Field<double,Location::CELL> T(1,"T","K",1), source(1,"source","W/m3",1);
        T(0)=300.0; source.fill(0.0);

        ThermalBoundaryConditions bc;
        bc["inlet"] = fixed_temperature(300.0);
        bc["outlet"].kind = ThermalBoundaryKind::ADIABATIC;
        bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        // No face field map supplied.
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"water",m,g,bc,T,source),
            std::invalid_argument, "requires a face mass-flux field map");

        // Map supplied but the region's field is not in it.
        ThermalFaceFields wrong_fields;
        wrong_fields.emplace("other_flux",
            Field<double,Location::FACE>(m.n_faces(),"other_flux","kg/s",1));
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"water",m,g,bc,T,source,{},&wrong_fields),
            std::invalid_argument, "unknown face mass-flux field");

        // Field size does not match the mesh.
        ThermalFaceFields short_fields;
        short_fields.emplace("feed_flux",
            Field<double,Location::FACE>(3,"feed_flux","kg/s",1));
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"water",m,g,bc,T,source,{},&short_fields),
            std::invalid_argument, "does not match the mesh");

        // Non-finite flux value.
        ThermalFaceFields nan_fields;
        nan_fields.emplace("feed_flux",
            Field<double,Location::FACE>(m.n_faces(),"feed_flux","kg/s",1));
        nan_fields.at("feed_flux").fill(0.0);
        nan_fields.at("feed_flux")(0) =
            std::numeric_limits<double>::quiet_NaN();
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"water",m,g,bc,T,source,{},&nan_fields),
            std::invalid_argument, "non-finite value");

        // Net imbalance: in -2, out +1.
        ThermalFaceFields bad_fields;
        bad_fields.emplace("feed_flux",
            Field<double,Location::FACE>(m.n_faces(),"feed_flux","kg/s",1));
        auto& bad_phi = bad_fields.at("feed_flux");
        bad_phi.fill(0.0);
        bad_phi(0) = -2.0;
        bad_phi(1) = +1.0;
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"water",m,g,bc,T,source,{},&bad_fields),
            std::invalid_argument, "not discretely conservative");

        // Inflow boundary with a flux-only condition: no thermal state.
        ThermalFaceFields fields;
        fields.emplace("feed_flux",
            Field<double,Location::FACE>(m.n_faces(),"feed_flux","kg/s",1));
        auto& phi = fields.at("feed_flux");
        phi.fill(0.0);
        phi(0) = -2.0;
        phi(1) = +2.0;
        ThermalBoundaryConditions flux_bc;
        flux_bc["inlet"].kind = ThermalBoundaryKind::FIXED_HEAT_FLUX;
        flux_bc["inlet"].value = 0.0;
        flux_bc["outlet"].kind = ThermalBoundaryKind::ADIABATIC;
        flux_bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"water",m,g,flux_bc,T,source,{},&fields),
            std::invalid_argument, "inflow boundary patch");
    });

    return run_all();
}
