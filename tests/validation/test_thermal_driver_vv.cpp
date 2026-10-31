// Thermal-only autonomous driver VV (issue #1150, A2+A3+A4).
//
// Every physical case is verified against an exact discrete oracle on an
// orthogonal mesh (the assembled system is small enough to solve by hand)
// plus the independent energy balance audit. Fail-closed batteries verify
// the explicit diagnostics of the deferred couplings, of the boundary
// mapping/value contract, of the prescribed-mass-flux face field contract,
// and of the conjugate-interface pairing contract. Method sheet:
// docs/audits/17-thermal-driver.md.
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

Mesh two_cell_bar(double x0)
{
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {
        {0,0,0},{0,1,0},{0,1,1},{0,0,1},
        {0.5,0,0},{0.5,1,0},{0.5,1,1},{0.5,0,1},
        {1,0,0},{1,1,0},{1,1,1},{1,0,1}
    };
    for (std::size_t i=0; i<12; ++i) m.points().set(i,p[i][0]+x0,p[i][1],p[i][2]);

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

Mesh two_cell_1d()
{
    return two_cell_bar(0.0);
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

// Two 1-D bars of unit cross-section joined by one conjugate or
// contact-resistance interface (A4): bar a of conductivity ka, bar b of
// conductivity kb, half-cell distance 0.25 m at the interface plane, so two
// bars placed side by side pair exactly (coincident centroids, opposed
// normals, equal areas). A CONJUGATE interface joins a solid and a
// stagnant fluid; a CONTACT_RESISTANCE interface joins two solids.
ThermalCase coupled_case(ThermalInterfaceModel model,
                         double contact_resistance, double ka, double kb)
{
    ThermalCase c;
    c.name = "thermal_driver_vv";
    ThermalMaterial a_material;
    a_material.name = "bar_a_material";
    a_material.properties.conductivity.reference_value = ka;
    a_material.properties.heat_capacity.reference_value = 1.0;
    a_material.properties.density.rho0 = 10.0;
    ThermalMaterial b_material;
    b_material.name = "bar_b_material";
    b_material.properties.conductivity.reference_value = kb;
    b_material.properties.heat_capacity.reference_value = 1.0;
    b_material.properties.density.rho0 = 10.0;
    c.materials.push_back(a_material);
    c.materials.push_back(b_material);
    ThermalRegion a;
    a.name = "bar_a";
    a.mesh_id = "mesh_a";
    a.material = "bar_a_material";
    a.initial_temperature = 350.0;
    ThermalRegion b = a;
    b.name = "bar_b";
    b.mesh_id = "mesh_b";
    b.material = "bar_b_material";
    if (model == ThermalInterfaceModel::CONJUGATE_CONDUCTION) {
        b.kind = ThermalRegionKind::FLUID;
        b.fluid_model = ThermalFluidModel::STAGNANT_DIFFUSION;
    }
    c.regions.push_back(a);
    c.regions.push_back(b);
    ThermalInterface interface;
    interface.name = "junction";
    interface.region_a = "bar_a";
    interface.patch_a = "right";
    interface.region_b = "bar_b";
    interface.patch_b = "left";
    interface.model = model;
    interface.contact_resistance = contact_resistance;
    c.interfaces.push_back(interface);
    c.validate();
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

    // Steady conjugate two-bar conduction (A4): solid bar a (k=2, unit
    // cross-section, 0.5 m cells) held at 400 K, stagnant-fluid bar b
    // (k=0.5) held at 300 K, joined by a CONJUGATE interface. Series
    // conductances per unit area: 8, 4, 1.6, 1, 2; total resistance 2.5,
    // hence q = 40 W and (A1,A2,B1,B2) = (395,385,360,320) K exactly.
    run_case("steady_conjugate_two_bars", [] {
        Mesh mesh_a = two_cell_bar(0.0), mesh_b = two_cell_bar(1.0);
        auto g_a = build_fv_geometry(mesh_a);
        auto g_b = build_fv_geometry(mesh_b);
        ThermalCase c = coupled_case(
            ThermalInterfaceModel::CONJUGATE_CONDUCTION, 0.0, 2.0, 0.5);
        Field<double,Location::CELL>
            T_a(2,"T","K",1), source_a(2,"source","W/m3",1),
            T_b(2,"T","K",1), source_b(2,"source","W/m3",1);
        T_a.fill(350.0); source_a.fill(0.0);
        T_b.fill(350.0); source_b.fill(0.0);

        ThermalBoundaryConditions bc_a, bc_b;
        bc_a["left"] = fixed_temperature(400.0);
        bc_a["walls"].kind = ThermalBoundaryKind::ADIABATIC;
        bc_b["right"] = fixed_temperature(300.0);
        bc_b["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        ThermalRegionInputs in_a, in_b;
        in_a.mesh = &mesh_a; in_a.geometry = &g_a;
        in_a.temperature = &T_a; in_a.source = &source_a;
        in_a.boundary_conditions = bc_a;
        in_b.mesh = &mesh_b; in_b.geometry = &g_b;
        in_b.temperature = &T_b; in_b.source = &source_b;
        in_b.boundary_conditions = bc_b;
        std::map<std::string,ThermalRegionInputs> inputs;
        inputs.emplace("bar_a", in_a);
        inputs.emplace("bar_b", in_b);

        const auto result = solve_thermal_regions_coupled(c, inputs);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T_a(0),395.0,1e-8);
        EXPECT_NEAR(T_a(1),385.0,1e-8);
        EXPECT_NEAR(T_b(0),360.0,1e-8);
        EXPECT_NEAR(T_b(1),320.0,1e-8);
        EXPECT_TRUE(result.regions.at("bar_a").converged);
        EXPECT_TRUE(result.regions.at("bar_b").converged);
        EXPECT_TRUE(result.regions.at("bar_a").balance.normalized_residual < 1e-9);
        EXPECT_TRUE(result.regions.at("bar_b").balance.normalized_residual < 1e-9);
        EXPECT_TRUE(result.interface_flux_imbalance < 1e-9);
    });

    // Steady contact-resistance jump (A4): two solid bars (k=1) joined by
    // a CONTACT_RESISTANCE interface with R_c = 0.5 m^2.K/W between 400 K
    // and 300 K. Series resistance 0.25+0.5+1+0.5+0.25 = 2.5, hence
    // q = 40 W and (A1,A2,B1,B2) = (390,370,330,310) K exactly, with the
    // temperature jump (360 -> 340 K) carried by the contact resistance.
    run_case("steady_contact_resistance_jump", [] {
        Mesh mesh_a = two_cell_bar(0.0), mesh_b = two_cell_bar(1.0);
        auto g_a = build_fv_geometry(mesh_a);
        auto g_b = build_fv_geometry(mesh_b);
        ThermalCase c = coupled_case(
            ThermalInterfaceModel::CONTACT_RESISTANCE, 0.5, 1.0, 1.0);
        Field<double,Location::CELL>
            T_a(2,"T","K",1), source_a(2,"source","W/m3",1),
            T_b(2,"T","K",1), source_b(2,"source","W/m3",1);
        T_a.fill(350.0); source_a.fill(0.0);
        T_b.fill(350.0); source_b.fill(0.0);

        ThermalBoundaryConditions bc_a, bc_b;
        bc_a["left"] = fixed_temperature(400.0);
        bc_a["walls"].kind = ThermalBoundaryKind::ADIABATIC;
        bc_b["right"] = fixed_temperature(300.0);
        bc_b["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        ThermalRegionInputs in_a, in_b;
        in_a.mesh = &mesh_a; in_a.geometry = &g_a;
        in_a.temperature = &T_a; in_a.source = &source_a;
        in_a.boundary_conditions = bc_a;
        in_b.mesh = &mesh_b; in_b.geometry = &g_b;
        in_b.temperature = &T_b; in_b.source = &source_b;
        in_b.boundary_conditions = bc_b;
        std::map<std::string,ThermalRegionInputs> inputs;
        inputs.emplace("bar_a", in_a);
        inputs.emplace("bar_b", in_b);

        const auto result = solve_thermal_regions_coupled(c, inputs);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T_a(0),390.0,1e-8);
        EXPECT_NEAR(T_a(1),370.0,1e-8);
        EXPECT_NEAR(T_b(0),330.0,1e-8);
        EXPECT_NEAR(T_b(1),310.0,1e-8);
        EXPECT_TRUE(result.regions.at("bar_a").converged);
        EXPECT_TRUE(result.regions.at("bar_b").converged);
        EXPECT_TRUE(result.regions.at("bar_a").balance.normalized_residual < 1e-9);
        EXPECT_TRUE(result.regions.at("bar_b").balance.normalized_residual < 1e-9);
        EXPECT_TRUE(result.interface_flux_imbalance < 1e-9);
    });

    // Transient coupled step (A4): the conjugate bars of the first case as
    // a TRANSIENT case, dt = 1 s, rho*cp = 10, V = 0.5 m^3 per cell, one
    // implicit step from T^0 = 350 K. The exact 4x4 backward-Euler system
    // solves to A1 = 508510/1353, A2 = 486830/1353, B1 = 474130/1353 and
    // B2 = 456710/1353: every coupling sweep restarts from T^0, so the
    // coupled solve advances the state by exactly one physical step.
    run_case("transient_coupled_step", [] {
        Mesh mesh_a = two_cell_bar(0.0), mesh_b = two_cell_bar(1.0);
        auto g_a = build_fv_geometry(mesh_a);
        auto g_b = build_fv_geometry(mesh_b);
        ThermalCase c = coupled_case(
            ThermalInterfaceModel::CONJUGATE_CONDUCTION, 0.0, 2.0, 0.5);
        c.time_mode = ThermalTimeMode::TRANSIENT;
        c.time_step = 1.0;
        Field<double,Location::CELL>
            T_a(2,"T","K",1), source_a(2,"source","W/m3",1),
            T_b(2,"T","K",1), source_b(2,"source","W/m3",1);
        T_a.fill(350.0); source_a.fill(0.0);
        T_b.fill(350.0); source_b.fill(0.0);

        ThermalBoundaryConditions bc_a, bc_b;
        bc_a["left"] = fixed_temperature(400.0);
        bc_a["walls"].kind = ThermalBoundaryKind::ADIABATIC;
        bc_b["right"] = fixed_temperature(300.0);
        bc_b["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        ThermalRegionInputs in_a, in_b;
        in_a.mesh = &mesh_a; in_a.geometry = &g_a;
        in_a.temperature = &T_a; in_a.source = &source_a;
        in_a.boundary_conditions = bc_a;
        in_b.mesh = &mesh_b; in_b.geometry = &g_b;
        in_b.temperature = &T_b; in_b.source = &source_b;
        in_b.boundary_conditions = bc_b;
        std::map<std::string,ThermalRegionInputs> inputs;
        inputs.emplace("bar_a", in_a);
        inputs.emplace("bar_b", in_b);

        const auto result = solve_thermal_regions_coupled(c, inputs);
        EXPECT_TRUE(result.converged);
        EXPECT_NEAR(T_a(0),508510.0/1353.0,1e-8);
        EXPECT_NEAR(T_a(1),486830.0/1353.0,1e-8);
        EXPECT_NEAR(T_b(0),474130.0/1353.0,1e-8);
        EXPECT_NEAR(T_b(1),456710.0/1353.0,1e-8);
        EXPECT_TRUE(result.regions.at("bar_a").converged);
        EXPECT_TRUE(result.regions.at("bar_b").converged);
        EXPECT_TRUE(result.regions.at("bar_a").balance.normalized_residual < 1e-9);
        EXPECT_TRUE(result.regions.at("bar_b").balance.normalized_residual < 1e-9);
        EXPECT_TRUE(result.interface_flux_imbalance < 1e-9);
    });

    // Fail-closed battery: unknown lookups and couplings out of scope.
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

    // Fail-closed battery (A4): an interface touching a region must be
    // covered by explicit per-face conditions in the single-region solve,
    // and the coupled solver must find a conformal, exclusively coupled
    // pairing of every interface patch.
    run_case("fail_closed_conjugate_interfaces", [] {
        Mesh mesh_a = two_cell_bar(0.0), mesh_b = two_cell_bar(1.0);
        auto g_a = build_fv_geometry(mesh_a);
        auto g_b = build_fv_geometry(mesh_b);
        ThermalCase c = coupled_case(
            ThermalInterfaceModel::CONJUGATE_CONDUCTION, 0.0, 2.0, 0.5);
        Field<double,Location::CELL>
            T_a(2,"T","K",1), source_a(2,"source","W/m3",1),
            T_b(2,"T","K",1), source_b(2,"source","W/m3",1);
        T_a.fill(350.0); source_a.fill(0.0);
        T_b.fill(350.0); source_b.fill(0.0);

        ThermalBoundaryConditions bc_a, bc_b;
        bc_a["left"] = fixed_temperature(400.0);
        bc_a["walls"].kind = ThermalBoundaryKind::ADIABATIC;
        bc_b["right"] = fixed_temperature(300.0);
        bc_b["walls"].kind = ThermalBoundaryKind::ADIABATIC;

        // Single-region solve: the case interface is not covered.
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"bar_a",mesh_a,g_a,bc_a,T_a,source_a),
            std::invalid_argument, "conjugate interface");

        // Interface condition vector longer than the patch.
        ThermalInterfaceFaceConditions oversized;
        oversized.emplace("right",
            std::vector<ThermalInterfaceFaceCondition>(2));
        EXPECT_THROW_WITH(
            solve_thermal_region(c,"bar_a",mesh_a,g_a,bc_a,T_a,source_a,
                                 {},nullptr,&oversized),
            std::invalid_argument, "does not match the patch");

        // A prescribed mass flux crossing an interface face.
        ThermalCase flux_case = fluid_flux_case("feed_flux");
        ThermalRegion other = flux_case.regions[0];
        other.name = "bar_b";
        other.kind = ThermalRegionKind::SOLID;
        other.fluid_model = ThermalFluidModel::STAGNANT_DIFFUSION;
        other.mass_flux_field.clear();
        flux_case.regions.push_back(other);
        ThermalInterface flux_interface;
        flux_interface.name = "junction";
        flux_interface.region_a = "water";
        flux_interface.patch_a = "right";
        flux_interface.region_b = "bar_b";
        flux_interface.patch_b = "left";
        flux_case.interfaces.push_back(flux_interface);
        flux_case.validate();
        Mesh m = two_cell_1d();
        auto g = build_fv_geometry(m);
        Field<double,Location::CELL> T(2,"T","K",1), source(2,"source","W/m3",1);
        T.fill(350.0); source.fill(0.0);
        ThermalFaceFields fields;
        fields.emplace("feed_flux",
            Field<double,Location::FACE>(m.n_faces(),"feed_flux","kg/s",1));
        auto& phi = fields.at("feed_flux");
        phi.fill(0.0);
        phi(0) = -1.0;
        phi(1) = +1.0;
        phi(6) = +1.0;
        ThermalBoundaryConditions flux_bc;
        flux_bc["left"] = fixed_temperature(300.0);
        flux_bc["walls"].kind = ThermalBoundaryKind::ADIABATIC;
        ThermalInterfaceFaceConditions junction;
        junction.emplace("right",
            std::vector<ThermalInterfaceFaceCondition>(
                1, ThermalInterfaceFaceCondition{1.0, 300.0}));
        EXPECT_THROW_WITH(
            solve_thermal_region(flux_case,"water",m,g,flux_bc,T,source,
                                 {},&fields,&junction),
            std::invalid_argument, "carries a prescribed mass flux");

        // Coupled solver: region inputs must be exact and complete.
        ThermalRegionInputs in_a, in_b;
        in_a.mesh = &mesh_a; in_a.geometry = &g_a;
        in_a.temperature = &T_a; in_a.source = &source_a;
        in_a.boundary_conditions = bc_a;
        in_b.mesh = &mesh_b; in_b.geometry = &g_b;
        in_b.temperature = &T_b; in_b.source = &source_b;
        in_b.boundary_conditions = bc_b;
        std::map<std::string,ThermalRegionInputs> inputs;
        inputs.emplace("bar_a", in_a);
        inputs.emplace("bar_b", in_b);

        std::map<std::string,ThermalRegionInputs> ghost_inputs = inputs;
        ghost_inputs.emplace("ghost", in_a);
        EXPECT_THROW_WITH(
            solve_thermal_regions_coupled(c, ghost_inputs),
            std::invalid_argument, "unknown region input");

        std::map<std::string,ThermalRegionInputs> missing_inputs = inputs;
        missing_inputs.erase("bar_b");
        EXPECT_THROW_WITH(
            solve_thermal_regions_coupled(c, missing_inputs),
            std::invalid_argument, "missing region input");

        std::map<std::string,ThermalRegionInputs> null_inputs = inputs;
        null_inputs.at("bar_a").mesh = nullptr;
        EXPECT_THROW_WITH(
            solve_thermal_regions_coupled(c, null_inputs),
            std::invalid_argument, "requires a mesh, a geometry");

        // Coupled solver: unknown interface patch.
        ThermalCase ghost_patch_case = c;
        ghost_patch_case.interfaces[0].patch_a = "ghost_patch";
        EXPECT_THROW_WITH(
            solve_thermal_regions_coupled(ghost_patch_case, inputs),
            std::invalid_argument, "references unknown patch");

        // Coupled solver: an interface patch must not also carry an
        // external boundary condition.
        ThermalBoundaryConditions bc_a_external = bc_a;
        bc_a_external["right"] = fixed_temperature(300.0);
        std::map<std::string,ThermalRegionInputs> external_inputs = inputs;
        external_inputs.at("bar_a").boundary_conditions = bc_a_external;
        EXPECT_THROW_WITH(
            solve_thermal_regions_coupled(c, external_inputs),
            std::invalid_argument, "also has an external boundary condition");

        // Coupled solver: non-conformal patches (two faces against one).
        Mesh wide_a = one_cell_cube({CubePatchSpec{"iface",{0,1}},
                                     CubePatchSpec{"walls",{2,3,4,5}}});
        auto wide_g = build_fv_geometry(wide_a);
        Mesh narrow_b = one_cell_cube({CubePatchSpec{"iface",{0}},
                                       CubePatchSpec{"walls",{1,2,3,4,5}}});
        auto narrow_g = build_fv_geometry(narrow_b);
        ThermalCase cubes = solid_case();
        cubes.regions[0].name = "cube_a";
        ThermalRegion cube_b_region = cubes.regions[0];
        cube_b_region.name = "cube_b";
        cubes.regions.push_back(cube_b_region);
        ThermalInterface cube_interface;
        cube_interface.name = "cube_junction";
        cube_interface.region_a = "cube_a";
        cube_interface.patch_a = "iface";
        cube_interface.region_b = "cube_b";
        cube_interface.patch_b = "iface";
        cubes.interfaces.push_back(cube_interface);
        Field<double,Location::CELL>
            T_wide(1,"T","K",1), source_wide(1,"source","W/m3",1),
            T_narrow(1,"T","K",1), source_narrow(1,"source","W/m3",1);
        T_wide.fill(300.0); source_wide.fill(0.0);
        T_narrow.fill(300.0); source_narrow.fill(0.0);
        ThermalBoundaryConditions bc_wide, bc_narrow;
        bc_wide["walls"] = fixed_temperature(310.0);
        bc_narrow["walls"] = fixed_temperature(310.0);
        ThermalRegionInputs in_wide, in_narrow;
        in_wide.mesh = &wide_a; in_wide.geometry = &wide_g;
        in_wide.temperature = &T_wide; in_wide.source = &source_wide;
        in_wide.boundary_conditions = bc_wide;
        in_narrow.mesh = &narrow_b; in_narrow.geometry = &narrow_g;
        in_narrow.temperature = &T_narrow; in_narrow.source = &source_narrow;
        in_narrow.boundary_conditions = bc_narrow;
        std::map<std::string,ThermalRegionInputs> cube_inputs;
        cube_inputs.emplace("cube_a", in_wide);
        cube_inputs.emplace("cube_b", in_narrow);
        EXPECT_THROW_WITH(
            solve_thermal_regions_coupled(cubes, cube_inputs),
            std::invalid_argument, "not conformal (different face counts)");

        // Coupled solver: paired faces of mismatched areas (unit cube
        // against the same face squashed to half the cross-section, so the
        // centroids still coincide).
        Mesh squashed_b;
        squashed_b.points().resize(8);
        const double q[8][3] = {{0,0.25,0},{1,0.25,0},{1,0.75,0},{0,0.75,0},
                                {0,0.25,1},{1,0.25,1},{1,0.75,1},{0,0.75,1}};
        for (std::size_t i=0;i<8;++i)
            squashed_b.points().set(i,q[i][0],q[i][1],q[i][2]);
        squashed_b.faces().push_face({0,3,2,1});
        squashed_b.faces().push_face({4,5,6,7});
        squashed_b.faces().push_face({0,1,5,4});
        squashed_b.faces().push_face({3,7,6,2});
        squashed_b.faces().push_face({0,4,7,3});
        squashed_b.faces().push_face({1,2,6,5});
        squashed_b.ownership().resize(6);
        for (std::size_t f=0;f<6;++f) {
            squashed_b.ownership().set_owner(f,0);
            squashed_b.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
        }
        squashed_b.cells().push_cell({0,1,2,3,4,5});
        Patch squash_iface;
        squash_iface.name = "iface";
        squash_iface.type = PatchType::WALL;
        squash_iface.face_ids = {0};
        Patch squash_walls;
        squash_walls.name = "walls";
        squash_walls.type = PatchType::WALL;
        squash_walls.face_ids = {1,2,3,4,5};
        squashed_b.boundary().add_patch(squash_iface);
        squashed_b.boundary().add_patch(squash_walls);
        auto squash_g = build_fv_geometry(squashed_b);
        Field<double,Location::CELL>
            T_sq(1,"T","K",1), source_sq(1,"source","W/m3",1);
        T_sq.fill(300.0); source_sq.fill(0.0);
        ThermalRegionInputs in_sq;
        in_sq.mesh = &squashed_b; in_sq.geometry = &squash_g;
        in_sq.temperature = &T_sq; in_sq.source = &source_sq;
        in_sq.boundary_conditions = bc_narrow;
        std::map<std::string,ThermalRegionInputs> area_inputs;
        area_inputs.emplace("cube_a", in_narrow);
        area_inputs.emplace("cube_b", in_sq);
        EXPECT_THROW_WITH(
            solve_thermal_regions_coupled(cubes, area_inputs),
            std::invalid_argument, "mismatched areas");
    });

    return run_all();
}
