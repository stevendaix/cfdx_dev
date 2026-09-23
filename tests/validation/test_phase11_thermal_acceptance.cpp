#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/cht_solver.h"
#include "cfdx/physics/thermal_models.h"
#include "common/test_harness.h"

#include <cmath>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

Mesh make_chain(std::size_t n)
{
    if (n == 0) throw std::invalid_argument("chain requires at least one cell");
    Mesh m;
    const std::size_t planes=n+1;
    m.points().resize(4*planes);
    for (std::size_t i=0; i<planes; ++i) {
        const double x=static_cast<double>(i);
        m.points().set(4*i+0,x,0,0);
        m.points().set(4*i+1,x,1,0);
        m.points().set(4*i+2,x,1,1);
        m.points().set(4*i+3,x,0,1);
    }

    std::vector<std::size_t> left(n),right(n),ym(n),yp(n),zm(n),zp(n);
    for (std::size_t i=0; i<n; ++i) {
        const std::size_t a=4*i, b=4*(i+1);
        if (i==0) {
            m.faces().push_face({a+0,a+3,a+2,a+1});
            left[i]=m.n_faces()-1;
        } else {
            left[i]=right[i-1];
        }

        m.faces().push_face({b+0,b+1,b+2,b+3});
        right[i]=m.n_faces()-1;

        m.faces().push_face({a+0,b+0,b+3,a+3}); ym[i]=m.n_faces()-1;
        m.faces().push_face({a+1,a+2,b+2,b+1}); yp[i]=m.n_faces()-1;
        m.faces().push_face({a+0,a+1,b+1,b+0}); zm[i]=m.n_faces()-1;
        m.faces().push_face({a+3,b+3,b+2,a+2}); zp[i]=m.n_faces()-1;
    }

    m.ownership().resize(m.n_faces());
    for (std::size_t f=0; f<m.n_faces(); ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    for (std::size_t i=0; i<n; ++i) {
        const std::size_t ids[6]={left[i],right[i],ym[i],yp[i],zm[i],zp[i]};
        m.cells().push_cell({ids[0],ids[1],ids[2],ids[3],ids[4],ids[5]});
        m.ownership().set_owner(ym[i],static_cast<CellIndex>(i));
        m.ownership().set_owner(yp[i],static_cast<CellIndex>(i));
        m.ownership().set_owner(zm[i],static_cast<CellIndex>(i));
        m.ownership().set_owner(zp[i],static_cast<CellIndex>(i));
        m.ownership().set_owner(right[i],static_cast<CellIndex>(i));
        if (i+1<n)
            m.ownership().set_neighbour(right[i],static_cast<CellIndex>(i+1));
    }

    Patch pleft; pleft.name="left"; pleft.type=PatchType::WALL; pleft.face_ids={left[0]};
    Patch pright; pright.name="right"; pright.type=PatchType::WALL; pright.face_ids={right[n-1]};
    Patch wall; wall.name="wall"; wall.type=PatchType::WALL;
    for (std::size_t i=0; i<n; ++i) {
        wall.face_ids.push_back(ym[i]);
        wall.face_ids.push_back(yp[i]);
        wall.face_ids.push_back(zm[i]);
        wall.face_ids.push_back(zp[i]);
    }
    m.boundary().add_patch(pleft);
    m.boundary().add_patch(pright);
    m.boundary().add_patch(wall);
    return m;
}

Mesh make_cube_with_interface(const char* interface_name)
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
    Patch interface; interface.name=interface_name; interface.type=PatchType::WALL; interface.face_ids={3};
    Patch wall; wall.name="wall"; wall.type=PatchType::WALL;
    wall.face_ids={0,1,2,4,5};
    m.boundary().add_patch(interface);
    m.boundary().add_patch(wall);
    return m;
}

double l2_error(const Field<double,Location::CELL>& T, const std::vector<double>& ref)
{
    double s=0.0;
    for (std::size_t i=0;i<T.size();++i) {
        const double d=T(i)-ref[i];
        s += d*d;
    }
    return std::sqrt(s/static_cast<double>(T.size()));
}

} // namespace

int main()
{
    run_case("phase11_multicell_1d_conduction_exact_linear_profile", [] {
        Mesh m=make_chain(8);
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
        Field<double,Location::CELL> T(m.n_cells(),"T","K",1), source(m.n_cells(),"Q","W/m3",1);
        T.fill(300.0); source.fill(0.0);

        ScalarBoundaryConditions bc;
        bc["left"]={ScalarBoundaryType::FIXED_VALUE,400.0,0.0};
        bc["right"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        bc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        EnergySolverControls c;
        c.conductivity=2.0; c.max_iterations=100; c.tolerance=1e-11; c.relaxation=1.0;
        auto r=solve_energy(m,g,flux,T,source,c,bc);
        EXPECT_TRUE(r.converged);

        std::vector<double> ref(m.n_cells());
        for (std::size_t i=0;i<m.n_cells();++i)
            ref[i]=400.0-100.0*(static_cast<double>(i)+0.5)/8.0;
        EXPECT_NEAR(l2_error(T,ref),0.0,1e-9);
        EXPECT_NEAR(r.history.back().energy_imbalance,0.0,1e-10);
    });

    run_case("phase11_conduction_mms_mesh_convergence", [] {
        auto solve_mms = [](std::size_t n) {
            Mesh m=make_chain(n);
            auto g=build_fv_geometry(m);
            Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
            Field<double,Location::CELL> T(m.n_cells(),"T","K",1), source(m.n_cells(),"Q","W/m3",1);
            T.fill(0.0);
            const double L=static_cast<double>(n);
            const double pi=std::acos(-1.0);
            for (std::size_t i=0;i<n;++i) {
                const double x=static_cast<double>(i)+0.5;
                source(i)=pi*pi/(L*L)*std::sin(pi*x/L);
            }
            ScalarBoundaryConditions bc;
            bc["left"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
            bc["right"]={ScalarBoundaryType::FIXED_VALUE,0.0,0.0};
            bc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};
            EnergySolverControls c;
            c.conductivity=1.0; c.max_iterations=100; c.tolerance=1e-11; c.relaxation=1.0;
            auto r=solve_energy(m,g,flux,T,source,c,bc);
            EXPECT_TRUE(r.converged);
            double err2=0.0, errinf=0.0;
            for (std::size_t i=0;i<n;++i) {
                const double x=static_cast<double>(i)+0.5;
                const double e=T(i)-std::sin(pi*x/L);
                err2 += e*e;
                errinf=std::max(errinf,std::abs(e));
            }
            return std::pair<double,double>{std::sqrt(err2/static_cast<double>(n)),errinf};
        };
        const auto e8=solve_mms(8);
        const auto e16=solve_mms(16);
        const auto e32=solve_mms(32);
        const auto e64=solve_mms(64);
        const double p1=std::log(e8.first/e16.first)/std::log(2.0);
        const double p2=std::log(e16.first/e32.first)/std::log(2.0);
        const double p3=std::log(e32.first/e64.first)/std::log(2.0);
        EXPECT_TRUE(e64.first < e32.first && e32.first < e16.first && e16.first < e8.first);
        EXPECT_TRUE(e64.second < e32.second && e32.second < e16.second && e16.second < e8.second);
        EXPECT_TRUE(p1 > 1.8 && p2 > 1.8 && p3 > 1.8);
    });

    run_case("phase11_cell_conductivity_series_resistance", [] {
        Mesh m=make_chain(2);
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
        Field<double,Location::CELL> T(m.n_cells(),"T","K",1), source(m.n_cells(),"Q","W/m3",1);
        Field<double,Location::CELL> k(m.n_cells(),"k","W/m/K",1);
        T.fill(300.0); source.fill(0.0); k(0)=1.0; k(1)=10.0;

        ScalarBoundaryConditions bc;
        bc["left"]={ScalarBoundaryType::FIXED_VALUE,400.0,0.0};
        bc["right"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        bc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0.0,0.0};

        EnergySolverControls c;
        c.conductivity=1.0; c.conductivity_field=&k; c.max_iterations=100;
        c.tolerance=1e-11; c.relaxation=1.0;
        auto r=solve_energy(m,g,flux,T,source,c,bc);
        EXPECT_TRUE(r.converged);
        const double q=100.0/(0.5/1.0+0.5/10.0);
        EXPECT_NEAR(T(0),400.0-q*0.5,1e-9);
        EXPECT_NEAR(T(1),300.0+q*0.05,1e-9);
    });

    run_case("phase11_temperature_dependent_conductivity_model", [] {
        ThermalConductivityControls c;
        c.model=ThermalConductivityModel::LINEAR_TEMPERATURE;
        c.k0=10.0; c.reference_temperature=300.0; c.slope=0.1;
        EXPECT_NEAR(evaluate_thermal_conductivity(c,300.0),10.0,1e-12);
        EXPECT_NEAR(evaluate_thermal_conductivity(c,400.0),20.0,1e-12);

        c.model=ThermalConductivityModel::POWER_LAW;
        c.k0=4.0; c.reference_temperature=300.0; c.exponent=1.0;
        EXPECT_NEAR(evaluate_thermal_conductivity(c,600.0),8.0,1e-12);
    });

    run_case("phase11_convective_boundary_and_multistep_transient_state", [] {
        Mesh m=make_cube_with_interface("interface");
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1); flux.fill(0.0);
        Field<double,Location::CELL> T(m.n_cells(),"T","K",1), source(m.n_cells(),"Q","W/m3",1);
        T(0)=300.0; source(0)=1000.0;

        ScalarBoundaryConditions bc;
        for (const char* p : {"interface","wall"})
            bc[p]={ScalarBoundaryType::CONVECTIVE,300.0,10.0};

        EnergySolverControls c;
        c.density=1.0; c.cp=1000.0; c.conductivity=1.0; c.dt=1.0;
        c.max_iterations=50; c.tolerance=1e-12; c.relaxation=1.0;

        auto r1=solve_energy(m,g,flux,T,source,c,bc);
        EXPECT_TRUE(r1.converged);
        const double t1=300.0+1000.0/(1000.0+60.0);
        EXPECT_NEAR(T(0),t1,1e-10);

        auto r2=solve_energy(m,g,flux,T,source,c,bc);
        EXPECT_TRUE(r2.converged);
        const double t2=300.0+(t1-300.0)*1000.0/1060.0+1000.0/1060.0;
        EXPECT_NEAR(T(0),t2,1e-10);
    });

    run_case("phase11_cht_contact_resistance_and_flux_continuity", [] {
        Mesh m1=make_cube_with_interface("interface");
        Mesh m2=make_cube_with_interface("interface");
        auto g1=build_fv_geometry(m1), g2=build_fv_geometry(m2);
        Field<double,Location::FACE> f1(m1.n_faces(),"f1","kg/s",1), f2(m2.n_faces(),"f2","kg/s",1);
        f1.fill(0.0); f2.fill(0.0);
        Field<double,Location::CELL> T1(1,"T1","K",1),T2(1,"T2","K",1);
        Field<double,Location::CELL> s1(1,"s1","W/m3",1),s2(1,"s2","W/m3",1);
        T1(0)=400.0; T2(0)=300.0; s1(0)=0.0; s2(0)=0.0;

        ScalarBoundaryConditions bc1,bc2;
        bc1["wall"]={ScalarBoundaryType::FIXED_VALUE,400.0,0.0};
        bc2["wall"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};

        EnergySolverControls e1,e2;
        e1.conductivity=1.0; e2.conductivity=1.0; e1.max_iterations=100; e2.max_iterations=100;
        e1.tolerance=1e-11; e2.tolerance=1e-11; e1.relaxation=1.0; e2.relaxation=1.0;

        ChtInterfaceControls c;
        c.region1_patch="interface"; c.region2_patch="interface";
        c.conductivity1=1.0; c.conductivity2=1.0;
        c.contact_resistance=0.1; c.max_iterations=100; c.tolerance=1e-10;
        c.matching_tolerance=1e-12;

        auto r=solve_two_region_cht(m1,g1,m2,g2,f1,f2,T1,T2,s1,s2,e1,e2,c,bc1,bc2);
        EXPECT_TRUE(r.converged);
        const double q=100.0/(0.5+0.1+0.5);
        EXPECT_NEAR(T1(0),400.0,1e-9);
        EXPECT_NEAR(T2(0),300.0,1e-9);
        const double tint1=400.0-q*0.5;
        const double tint2=300.0+q*0.5;
        EXPECT_NEAR(tint1-tint2,q*0.1,1e-9);
        EXPECT_NEAR(r.interface_imbalance,0.0,1e-10);
    });

    return run_all();
}
