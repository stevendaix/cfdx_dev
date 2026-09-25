#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/cht_solver.h"
#include "cfdx/physics/radiation_solver.h"
#include "cfdx/physics/radiation.h"
#include "common/test_harness.h"

#include <array>
#include <cmath>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh one_d_mesh(std::size_t n, double x0 = 0.0, double x1 = 1.0,
                       const std::string& left_name = "left",
                       const std::string& right_name = "right")
{
    Mesh m;
    const std::size_t np = 4*(n+1);
    m.points().resize(np);
    for (std::size_t i=0;i<=n;++i) {
        const double x=x0+(x1-x0)*static_cast<double>(i)/static_cast<double>(n);
        const std::size_t p=4*i;
        m.points().set(p,x,0,0); m.points().set(p+1,x,1,0);
        m.points().set(p+2,x,1,1); m.points().set(p+3,x,0,1);
    }

    const std::size_t nf=5*n+2;
    // FaceConnectivity grows through push_face(); no resize API is exposed.
    m.ownership().resize(nf);
    std::vector<std::size_t> left_faces, right_faces, walls_faces;
    std::size_t f=0;

    auto add_face = [&](std::initializer_list<std::size_t> ids,
                        std::size_t owner, long long neighbour,
                        std::vector<std::size_t>* patch) {
        m.faces().push_face(ids);
        m.ownership().set_owner(f,owner);
        m.ownership().set_neighbour(f,neighbour);
        if (patch) patch->push_back(f);
        ++f;
    };

    add_face({0,3,2,1},0,FaceOwnership::BOUNDARY,&left_faces);
    for (std::size_t i=0;i<n;++i) {
        const std::size_t p=4*i;
        m.cells().push_cell({f, f+1, f+2, f+3, f+4, f+5});
        if (i>0) {
            // Replace the cell's left-face entry with the already-created
            // internal face. The first cell is created below by construction.
        }
        (void)p;
    }

    // Rebuild the topology deterministically: each cell has left/right and
    // four transverse faces; transverse faces are boundary faces.
    m = Mesh();
    m.points().resize(np);
    for (std::size_t i=0;i<=n;++i) {
        const double x=x0+(x1-x0)*static_cast<double>(i)/static_cast<double>(n);
        const std::size_t p=4*i;
        m.points().set(p,x,0,0); m.points().set(p+1,x,1,0);
        m.points().set(p+2,x,1,1); m.points().set(p+3,x,0,1);
    }
    // FaceConnectivity grows through push_face(); no resize API is exposed.
    m.ownership().resize(nf);
    f=0;
    left_faces.clear(); right_faces.clear(); walls_faces.clear();

    // Left boundary.
    m.faces().push_face({0,3,2,1});
    m.ownership().set_owner(f,0); m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    left_faces.push_back(f++);

    // Internal x-faces.
    std::vector<std::size_t> xf(n+1);
    xf[0]=0;
    for (std::size_t i=1;i<n;++i) {
        const std::size_t p=4*i;
        m.faces().push_face({p,p+1,p+2,p+3});
        m.ownership().set_owner(f,i-1); m.ownership().set_neighbour(f,i);
        xf[i]=f++;
    }
    m.faces().push_face({4*n,4*n+1,4*n+2,4*n+3});
    m.ownership().set_owner(f,n-1); m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    right_faces.push_back(f++); xf[n]=f-1;

    std::vector<std::array<std::size_t,4>> side(n);
    for (std::size_t i=0;i<n;++i) {
        const std::size_t p=4*i, q=4*(i+1);
        m.faces().push_face({p,p+1,q+1,q}); m.ownership().set_owner(f,i); m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY); side[i][0]=f++;
        m.faces().push_face({p+3,q+3,q+2,p+2}); m.ownership().set_owner(f,i); m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY); side[i][1]=f++;
        m.faces().push_face({p,q,q+3,p+3}); m.ownership().set_owner(f,i); m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY); side[i][2]=f++;
        m.faces().push_face({p+1,p+2,q+2,q+1}); m.ownership().set_owner(f,i); m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY); side[i][3]=f++;
        for (auto sf:side[i]) walls_faces.push_back(sf);
    }
    for (std::size_t i=0;i<n;++i)
        m.cells().push_cell({xf[i],xf[i+1],side[i][0],side[i][1],side[i][2],side[i][3]});

    Patch left; left.name=left_name; left.type=PatchType::WALL; left.face_ids=left_faces;
    Patch right; right.name=right_name; right.type=PatchType::WALL; right.face_ids=right_faces;
    Patch walls; walls.name="walls"; walls.type=PatchType::WALL; walls.face_ids=walls_faces;
    m.boundary().add_patch(left); m.boundary().add_patch(right); m.boundary().add_patch(walls);
    return m;
}

static Mesh cube_mesh(double x0, double x1,
                      const std::string& left, const std::string& right)
{
    return one_d_mesh(1,x0,x1,left,right);
}

static void check_history(const EnergySolveResult& r, const char* name)
{
    EXPECT_TRUE(r.converged);
    EXPECT_TRUE(!r.history.empty());
    std::cout << "THERMAL_CHALLENGE: case=" << name
              << " iterations=" << r.iterations
              << " residual=" << r.history.back().residual
              << " energy_imbalance=" << r.history.back().energy_imbalance << '\n';
}

static std::vector<DiscreteDirection> sn_quadrature(int order)
{
    static const double mu2[] = {-0.5773502691896257,0.5773502691896257};
    static const double w2[] = {1.0,1.0};
    static const double mu4[] = {-0.8611363115940526,-0.3399810435848563,
                                  0.3399810435848563,0.8611363115940526};
    static const double w4[] = {0.3478548451374538,0.6521451548625461,
                                0.6521451548625461,0.3478548451374538};
    static const double mu6[] = {-0.9324695142031521,-0.6612093864662645,-0.2386191860831969,
                                  0.2386191860831969,0.6612093864662645,0.9324695142031521};
    static const double w6[] = {0.1713244923791704,0.3607615730481386,0.4679139345726910,
                                0.4679139345726910,0.3607615730481386,0.1713244923791704};
    static const double mu8[] = {-0.9602898564975363,-0.7966664774136267,-0.5255324099163290,-0.1834346424956498,
                                  0.1834346424956498,0.5255324099163290,0.7966664774136267,0.9602898564975363};
    static const double w8[] = {0.1012285362903763,0.2223810344533745,0.3137066458778873,0.3626837833783620,
                                0.3626837833783620,0.3137066458778873,0.2223810344533745,0.1012285362903763};
    const double* mu=nullptr; const double* w=nullptr;
    switch(order) {
        case 2: mu=mu2; w=w2; break;
        case 4: mu=mu4; w=w4; break;
        case 6: mu=mu6; w=w6; break;
        case 8: mu=mu8; w=w8; break;
        default: throw std::invalid_argument("unsupported S_N order");
    }
    const double dphi=2.0*M_PI/static_cast<double>(2*order);
    std::vector<DiscreteDirection> q;
    for(int j=0;j<order;++j)
        for(int k=0;k<2*order;++k) {
            const double phi=(static_cast<double>(k)+0.5)*dphi;
            const double s=std::sqrt(std::max(0.0,1.0-mu[j]*mu[j]));
            q.push_back({s*std::cos(phi),s*std::sin(phi),mu[j],w[j]*dphi});
        }
    return q;
}

int main()
{
    run_case("thermal_multicell_analytical_conduction", [] {
        const std::size_t n=32;
        Mesh m=one_d_mesh(n);
        auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
        Field<double,Location::CELL> T(n,"T","K",1), source(n,"source","W/m3",1); T.fill(350.0); source.fill(0.0);
        ScalarBoundaryConditions bc{{"left",{ScalarBoundaryType::FIXED_VALUE,400,0}},
                                    {"right",{ScalarBoundaryType::FIXED_VALUE,300,0}},
                                    {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
        EnergySolverControls c; c.conductivity=1; c.relaxation=1; c.max_iterations=100; c.tolerance=1e-11;
        const auto r=solve_energy(m,g,phi,T,source,c,bc);
        check_history(r,"multicell_conduction");
        EXPECT_TRUE(r.history.back().residual <= c.tolerance);
        EXPECT_TRUE(std::isfinite(r.history.back().residual));
        EXPECT_TRUE(r.history.back().energy_imbalance <= c.tolerance);
        for(std::size_t i=0;i<n;++i) {
            const double x=(static_cast<double>(i)+0.5)/static_cast<double>(n);
            EXPECT_NEAR(T(i),400.0-100.0*x,2e-9);
        }
    });

    run_case("thermal_mms_quadratic_refinement", [] {
        double previous_error=0.0;
        for(std::size_t n: {4u,8u,16u,32u}) {
            Mesh m=one_d_mesh(n); auto g=build_fv_geometry(m);
            Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
            Field<double,Location::CELL> T(n,"T","K",1), source(n,"source","W/m3",1);
            T.fill(0.25); source.fill(-2.0);
            ScalarBoundaryConditions bc{{"left",{ScalarBoundaryType::FIXED_VALUE,0,0}},
                                        {"right",{ScalarBoundaryType::FIXED_VALUE,1,0}},
                                        {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
            EnergySolverControls c; c.conductivity=1; c.relaxation=1; c.max_iterations=100; c.tolerance=1e-11;
            const auto r=solve_energy(m,g,phi,T,source,c,bc);
            check_history(r,"mms_quadratic");
            double err=0.0;
            for(std::size_t i=0;i<n;++i) {
                const double x=(static_cast<double>(i)+0.5)/static_cast<double>(n);
                err=std::max(err,std::abs(T(i)-x*x));
            }
            if(previous_error>0.0) {
                const double order=std::log(previous_error/err)/std::log(2.0);
                std::cout << "THERMAL_CHALLENGE: case=mms_quadratic n=" << n
                          << " L_inf=" << err << " observed_order=" << order << '\n';
                EXPECT_TRUE(order>1.8);
            } else {
                std::cout << "THERMAL_CHALLENGE: case=mms_quadratic n=" << n
                          << " L_inf=" << err << '\n';
            }
            previous_error=err;
        }
    });

    run_case("thermal_volumetric_generation_refinement", [] {
        double previous_error=0.0;
        for(std::size_t n: {4u,8u,16u,32u}) {
            Mesh m=one_d_mesh(n); auto g=build_fv_geometry(m);
            Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
            Field<double,Location::CELL> T(n,"T","K",1), source(n,"source","W/m3",1);
            T.fill(300.0); source.fill(100.0);
            ScalarBoundaryConditions bc{{"left",{ScalarBoundaryType::FIXED_VALUE,300,0}},
                                        {"right",{ScalarBoundaryType::FIXED_VALUE,300,0}},
                                        {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
            EnergySolverControls c; c.conductivity=2; c.relaxation=1; c.max_iterations=100; c.tolerance=1e-12;
            const auto r=solve_energy(m,g,phi,T,source,c,bc);
            check_history(r,"volumetric_generation");
            double err=0.0, balance= r.history.back().energy_imbalance;
            for(std::size_t i=0;i<n;++i) {
                const double x=(static_cast<double>(i)+0.5)/static_cast<double>(n);
                const double exact=300.0 + 100.0/(4.0)*(x-x*x);
                err=std::max(err,std::abs(T(i)-exact));
            }
            EXPECT_TRUE(balance<1e-12);
            if(previous_error>0.0) {
                const double order=std::log(previous_error/err)/std::log(2.0);
                std::cout << "THERMAL_CHALLENGE: case=volumetric_generation n=" << n
                          << " L_inf=" << err << " observed_order=" << order
                          << " energy_imbalance=" << balance << '\n';
                EXPECT_TRUE(order>1.8);
            } else {
                std::cout << "THERMAL_CHALLENGE: case=volumetric_generation n=" << n
                          << " L_inf=" << err << " energy_imbalance=" << balance << '\n';
            }
            previous_error=err;
        }
    });

    run_case("thermal_volumetric_power_sweep", [] {
        // Uniform volumetric heating in a unit slab with T(0)=T(L)=T0:
        // T(x)=T0 + q''' x(L-x)/(2k), so DeltaT_max=q''' L^2/(8k).
        // This validates both the W/m3 source convention and its scaling.
        const std::vector<double> powers={1.0e2,1.0e3,1.0e4,1.0e5};
        const double k=2.0;
        const double T0=300.0;
        const double L=1.0;
        const std::size_t n=32;
        double previous_dT = 0.0;
        for(const double qv:powers) {
            Mesh m=one_d_mesh(n); auto g=build_fv_geometry(m);
            Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0.0);
            Field<double,Location::CELL> T(n,"T","K",1), source(n,"source","W/m3",1);
            T.fill(T0); source.fill(qv);
            ScalarBoundaryConditions bc{{"left",{ScalarBoundaryType::FIXED_VALUE,T0,0}},
                                        {"right",{ScalarBoundaryType::FIXED_VALUE,T0,0}},
                                        {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
            EnergySolverControls c; c.conductivity=k; c.relaxation=1;
            c.max_iterations=100; c.tolerance=1e-12;
            const auto r=solve_energy(m,g,phi,T,source,c,bc);
            check_history(r,"volumetric_power_sweep");
            double Tmax=T(0), Tmin=T(0);
            for(std::size_t i=0;i<n;++i) {
                Tmax=std::max(Tmax,T(i));
                Tmin=std::min(Tmin,T(i));
            }
            double exact_dT_cell=0.0;
            for(std::size_t i=0;i<n;++i) {
                const double x=(static_cast<double>(i)+0.5)/static_cast<double>(n);
                exact_dT_cell=std::max(
                    exact_dT_cell,qv*x*(L-x)/(2.0*k));
            }
            const double exact_dT_continuous=qv*L*L/(8.0*k);
            const double expected_power=qv*L; // unit cross-sectional area
            const double generated_power=qv*std::accumulate(
                g.cell_volumes.begin(),g.cell_volumes.end(),0.0);
            const double center_x=0.5;
            const double center_exact=T0+qv*center_x*(L-center_x)/(2.0*k);
            std::cout << "THERMAL_POWER_STUDY: qvol=" << qv
                      << " W/m3 Tmax=" << Tmax
                      << " dTmax=" << (Tmax-T0)
                      << " exact_dTmax_cell=" << exact_dT_cell
                      << " exact_dTmax_continuous=" << exact_dT_continuous
                      << " generated_power=" << generated_power
                      << " expected_power=" << expected_power
                      << " center_exact=" << center_exact
                      << " Tmin=" << Tmin << '\n';
            EXPECT_TRUE(Tmax>T0);
            EXPECT_TRUE(Tmin>=T0);
            // Tmax is a cell-centre value, so the exact discrete oracle must
            // use the cell centre nearest x=L/2 rather than the continuous
            // maximum at x=L/2. The latter differs by O(h^2) on an even mesh.
            EXPECT_NEAR(Tmax-T0, exact_dT_cell,
                         1e-12*std::max(1.0,exact_dT_cell));
            EXPECT_NEAR(generated_power,expected_power,1e-12);
            if(previous_dT > 0.0)
                EXPECT_NEAR((Tmax-T0)/previous_dT, 10.0, 1e-12);
            previous_dT=Tmax-T0;
        }
    });

    run_case("thermal_cht_interface_conservation", [] {
        Mesh m1=cube_mesh(0,1,"outer_hot","interface");
        Mesh m2=cube_mesh(1,2,"interface","outer_cold");
        auto g1=build_fv_geometry(m1), g2=build_fv_geometry(m2);
        Field<double,Location::FACE> f1(m1.n_faces(),"f1","kg/s",1), f2(m2.n_faces(),"f2","kg/s",1);
        f1.fill(0); f2.fill(0);
        Field<double,Location::CELL> T1(1,"T1","K",1), T2(1,"T2","K",1);
        T1(0)=400; T2(0)=300;
        Field<double,Location::CELL> s1(1,"s1","W/m3",1), s2(1,"s2","W/m3",1); s1.fill(0); s2.fill(0);
        EnergySolverControls e1; e1.conductivity=2; e1.relaxation=1; e1.max_iterations=20; e1.tolerance=1e-10;
        EnergySolverControls e2=e1;
        ChtInterfaceControls c; c.region1_patch="interface"; c.region2_patch="interface";
        c.conductivity1=2; c.conductivity2=1; c.tolerance=1e-10; c.temperature_tolerance=1e-10;
        ScalarBoundaryConditions b1{{"outer_hot",{ScalarBoundaryType::FIXED_VALUE,400,0}},
                                    {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
        ScalarBoundaryConditions b2{{"outer_cold",{ScalarBoundaryType::FIXED_VALUE,300,0}},
                                    {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
        const auto r=solve_two_region_cht(m1,g1,m2,g2,f1,f2,T1,T2,s1,s2,e1,e2,c,b1,b2);
        std::cout << "THERMAL_CHALLENGE: case=CHT iterations=" << r.iterations
                  << " interface_imbalance=" << r.interface_imbalance
                  << " interface_dT=" << r.interface_temperature_change << '\n';
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(r.interface_imbalance<1e-10);
        EXPECT_TRUE(r.interface_temperature_change<1e-10);
        // T1/T2 are cell-centre temperatures, not interface temperatures.
        // For k1=2, k2=1 and unit-length layers between 400 K and 300 K,
        // q = 100/(1/2 + 1) = 66.666... W/m2.
        EXPECT_NEAR(T1(0),383.3333333333333,1e-10);
        EXPECT_NEAR(T2(0),333.3333333333333,1e-10);
    });

    run_case("radiation_uniform_thermal_equilibrium", [] {
        // Exact LTE equilibrium for a closed gray participating medium:
        // all walls and the cell are at the same temperature, so every
        // ordinate has I=Ib, G=4*pi*Ib=4*sigma*T^4 and q_rad=0.
        Mesh m=one_d_mesh(1); auto g=build_fv_geometry(m);
        Field<double,Location::CELL> T(1,"T","K",1), G(1,"G","W/m2",1), qrad(1,"qrad","W/m3",1);
        T(0)=300.0; G.fill(0.0); qrad.fill(0.0);
        auto dirs=sn_quadrature(4);
        double weight_sum=0.0;
        for(const auto& d:dirs) weight_sum+=d.weight;
        EXPECT_NEAR(weight_sum,4.0*M_PI,1e-12);

        RadiationTransportControls rc;
        rc.absorption=0.1; rc.scattering=0.0;
        rc.max_iterations=100; rc.tolerance=1e-12;
        ScalarBoundaryConditions rbcs{
            {"left",{ScalarBoundaryType::FIXED_VALUE,blackbody_emissive_power(300.0)/M_PI,0}},
            {"right",{ScalarBoundaryType::FIXED_VALUE,blackbody_emissive_power(300.0)/M_PI,0}},
            {"walls",{ScalarBoundaryType::FIXED_VALUE,blackbody_emissive_power(300.0)/M_PI,0}}
        };
        const auto r=solve_participating_radiation(
            m,g,T,G,qrad,dirs,rc,rbcs);
        EXPECT_TRUE(r.converged);
        const double expected_G=4.0*blackbody_emissive_power(300.0);
        std::cout << "RADIATION_EQUILIBRIUM: iterations=" << r.iterations
                  << " G=" << G(0)
                  << " expected_G=" << expected_G
                  << " qrad=" << qrad(0)
                  << " weight_sum=" << weight_sum << '\n';
        EXPECT_NEAR(G(0),expected_G,1e-10*std::max(1.0,expected_G));
        EXPECT_NEAR(qrad(0),0.0,1e-10);
    });

    run_case("thermal_radiation_coupled_conservation", [] {
        Mesh m=one_d_mesh(1); auto g=build_fv_geometry(m);
        Field<double,Location::FACE> phi(m.n_faces(),"phi","kg/s",1); phi.fill(0);
        Field<double,Location::CELL> T(1,"T","K",1), nonrad(1,"nonrad","W/m3",1);
        // Start from the radiation-wall equilibrium temperature so the
        // coupled nonlinear solve tests physical evolution rather than an
        // artificial large first Newton/fixed-point jump.
        T(0)=300; nonrad.fill(100);
        Field<double,Location::CELL> G(1,"G","W/m2",1); G.fill(0);
        auto dirs=sn_quadrature(4);
        RadiationEnergyCouplingControls c;
        c.radiation.absorption=0.1; c.radiation.scattering=0.0;
        c.radiation.max_iterations=100; c.radiation.tolerance=1e-12;
        c.energy.conductivity=1; c.energy.relaxation=0.5; c.energy.max_iterations=200; c.energy.tolerance=1e-12;
        c.max_outer_iterations=100; c.outer_relaxation=0.5; c.tolerance=1e-8;
        ScalarBoundaryConditions rbcs{{"left",{ScalarBoundaryType::FIXED_VALUE,blackbody_emissive_power(300)/M_PI,0}},
                                      {"right",{ScalarBoundaryType::FIXED_VALUE,blackbody_emissive_power(300)/M_PI,0}},
                                      {"walls",{ScalarBoundaryType::FIXED_VALUE,blackbody_emissive_power(300)/M_PI,0}}};
        ScalarBoundaryConditions tbcs{{"left",{ScalarBoundaryType::FIXED_VALUE,300,0}},
                                      {"right",{ScalarBoundaryType::FIXED_VALUE,300,0}},
                                      {"walls",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
        const auto r=solve_radiation_energy_coupled(m,g,phi,T,nonrad,G,dirs,c,rbcs,tbcs);
        std::cout << "THERMAL_CHALLENGE: case=radiation_energy_coupled iterations=" << r.iterations
                  << " source_residual=" << r.source_residuals.back()
                  << " energy_balance=" << r.energy_balance_residuals.back() << '\n';
        EXPECT_TRUE(r.converged);
        EXPECT_TRUE(r.energy_balance_residuals.back()<1e-8);
        EXPECT_TRUE(T(0)>300.0);
        EXPECT_TRUE(T(0)<500.0);
    });

    run_case("radiation_SN_angular_moment_validation", [] {
        for(int order: {2,4,6,8}) {
            auto q=sn_quadrature(order);
            validate_discrete_directions(q,1e-9);
            std::cout << "RADIATION_RESIDUAL: angular=S" << order
                      << " directions=" << q.size() << " quadrature=PASS\n";
        }
    });

    return run_all();
}
