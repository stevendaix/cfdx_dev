#include "cfdx/physics/turbulence_solver.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/radiation_solver.h"
#include "cfdx/physics/cht_solver.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

Mesh cube_with_patches()
{
    Mesh m;
    m.points().resize(8);
    const double p[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
    };
    for (std::size_t i=0; i<8; ++i)
        m.points().set(i,p[i][0],p[i][1],p[i][2]);

    m.faces().push_face({0,3,2,1}); // x-
    m.faces().push_face({4,5,6,7}); // x+
    m.faces().push_face({0,1,5,4}); // y-
    m.faces().push_face({3,7,6,2}); // y+ interface
    m.faces().push_face({0,4,7,3}); // z-
    m.faces().push_face({1,2,6,5}); // z+

    m.ownership().resize(6);
    for (std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});

    const char* names[6] = {"x_min","x_max","y_min","interface","z_min","z_max"};
    for (std::size_t f=0; f<6; ++f) {
        Patch patch;
        patch.name=names[f];
        patch.type=PatchType::WALL;
        patch.face_ids={f};
        m.boundary().add_patch(patch);
    }
    return m;
}

std::vector<DiscreteDirection> isotropic_directions()
{
    return {
        {1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},
        {0,1,0,2*M_PI/3},{0,-1,0,2*M_PI/3},
        {0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}
    };
}

ScalarBoundaryConditions fixed_wall_conditions(
    double value, const std::string& interface_name)
{
    ScalarBoundaryConditions bc;
    for (const char* name : {"x_min","x_max","y_min","z_min","z_max"})
        bc[name] = {ScalarBoundaryType::FIXED_VALUE,value,0.0};
    // The interface value is supplied by the CHT coupling itself.
    (void)interface_name;
    return bc;
}

void require_close(double value, double reference, double tolerance,
                   const char* message)
{
    if (!std::isfinite(value) || std::abs(value-reference)>tolerance)
        throw std::runtime_error(message);
}

} // namespace

int main()
{
    try {
        std::cout<<"LEVEL_C_CASE_01\n";
        // Level C-01: turbulence transport must preserve positivity while
        // receiving a non-zero source on an actual CFDX mesh.
        {
            Mesh m=cube_with_patches();
            auto g=build_fv_geometry(m);
            Field<double,Location::FACE> mass_flux(m.n_faces(),"phi","kg/s",1);
            mass_flux.fill(0.0);
            Field<double,Location::CELL> k(1,"k","m2/s2",1);
            Field<double,Location::CELL> epsilon(1,"epsilon","m2/s3",1);
            Field<double,Location::CELL> production(1,"production","m2/s3",1);
            k(0)=0.09;
            epsilon(0)=0.03;
            production(0)=0.02;

            TurbulenceTransportControls controls;
            controls.model=TurbulenceModel::KEPSILON;
            controls.density=1.0;
            controls.molecular_viscosity=1e-3;
            ScalarBoundaryConditions bc;
            bc["interface"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0.0};
            bc["x_min"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0.0};
            bc["x_max"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0.0};
            bc["y_min"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0.0};
            bc["z_min"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0.0};
            bc["z_max"]={ScalarBoundaryType::FIXED_VALUE,1e-6,0.0};

            const auto r=solve_kepsilon_transport(
                m,g,mass_flux,k,epsilon,production,controls,bc,bc,100,1e-8);
            if (r.iterations==0 || k(0)<controls.k_min ||
                epsilon(0)<controls.epsilon_min)
                throw std::runtime_error("turbulence transport positivity gate failed");
        }

        std::cout<<"LEVEL_C_CASE_02\n";
        // Level C-02: radiation/energy must converge on the same physical
        // source state, and the reported energy-balance residual must close.
        {
            Mesh m=cube_with_patches();
            auto g=build_fv_geometry(m);
            Field<double,Location::FACE> mass_flux(m.n_faces(),"phi","kg/s",1);
            mass_flux.fill(0.0);
            Field<double,Location::CELL> T(1,"T","K",1),source(1,"source","W/m3",1);
            Field<double,Location::CELL> irradiation(1,"G","W/m2",1);
            T(0)=800.0; source(0)=0.0; irradiation(0)=0.0;

            ScalarBoundaryConditions thermal_bc;
            ScalarBoundaryConditions radiation_bc;
            for (const char* name : {"x_min","x_max","y_min","interface","z_min","z_max"}) {
                thermal_bc[name]={ScalarBoundaryType::FIXED_VALUE,800.0,0.0};
                radiation_bc[name]={
                    ScalarBoundaryType::FIXED_VALUE,
                    blackbody_intensity(800.0),0.0};
            }

            RadiationEnergyCouplingControls controls;
            controls.radiation.absorption=0.5;
            controls.radiation.max_iterations=100;
            controls.radiation.tolerance=1e-10;
            controls.energy.conductivity=1.0;
            controls.energy.density=1.0;
            controls.energy.cp=1000.0;
            controls.energy.max_iterations=100;
            controls.energy.tolerance=1e-10;
            controls.max_outer_iterations=50;
            controls.tolerance=1e-8;

            const auto r=solve_radiation_energy_coupled(
                m,g,mass_flux,T,source,irradiation,isotropic_directions(),
                controls,radiation_bc,thermal_bc);
            if (!r.converged || r.energy_balance_residuals.empty())
                throw std::runtime_error("radiation-energy coupling did not converge");
            require_close(r.energy_balance_residuals.back(),0.0,1e-8,
                           "radiation-energy energy-balance gate failed");
        }

        std::cout<<"LEVEL_C_CASE_03\n";
        // Level C-03: two-region CHT uses a single matched interface with
        // independent hot/cold boundaries and verifies flux continuity.
        {
            Mesh m1=cube_with_patches();
            Mesh m2=cube_with_patches();
            auto g1=build_fv_geometry(m1);
            auto g2=build_fv_geometry(m2);

            Field<double,Location::FACE> f1(m1.n_faces(),"f1","kg/s",1);
            Field<double,Location::FACE> f2(m2.n_faces(),"f2","kg/s",1);
            f1.fill(0.0); f2.fill(0.0);

            Field<double,Location::CELL> T1(1,"T1","K",1),T2(1,"T2","K",1);
            Field<double,Location::CELL> s1(1,"s1","W/m3",1),s2(1,"s2","W/m3",1);
            T1(0)=500.0; T2(0)=300.0; s1(0)=0.0; s2(0)=0.0;

            EnergySolverControls ec;
            ec.conductivity=1.0;
            ec.max_iterations=100;
            ec.tolerance=1e-10;

            ChtInterfaceControls cc;
            cc.region1_patch="interface";
            cc.region2_patch="interface";
            cc.conductivity1=1.0;
            cc.conductivity2=1.0;
            cc.tolerance=1e-8;
            cc.max_iterations=100;
            cc.matching_tolerance=1e-12;

            auto bc1=fixed_wall_conditions(500.0,"interface");
            auto bc2=fixed_wall_conditions(300.0,"interface");

            const auto r=solve_two_region_cht(
                m1,g1,m2,g2,f1,f2,T1,T2,s1,s2,ec,ec,cc,bc1,bc2);
            if (!r.converged || r.interface_imbalance>cc.tolerance)
                throw std::runtime_error("CHT interface flux-balance gate failed");
            if (!(T1(0)>T2(0)) || !std::isfinite(T1(0)) ||
                !std::isfinite(T2(0)))
                throw std::runtime_error("CHT temperature-state gate failed");
        }

        std::cout<<"LEVEL_C_COUPLED_VERIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr<<"LEVEL_C_COUPLED_VERIFICATION: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
