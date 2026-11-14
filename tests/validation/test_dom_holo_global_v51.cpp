// V51: coupled gas/gray wall multisecant and HO-consistent LO predictor.
// Fail closed: every candidate is screened against independently reassembled
// angle-by-angle FV current, Schur fixed point and global energy.
#include "cfdx/physics/dom_holo_global_v51.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void check(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
Mesh grid(std::size_t nx,std::size_t ny,std::size_t nz,double depth) {
    auto m=cfdx::testing::thermal_v36::make_mesh(nx,ny,nz,.018,0.,depth);
    for(std::size_t i=0;i<m.points().size();++i) {
        const double x=m.points().x_data()[i],y=m.points().y_data()[i],
                     z=m.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        m.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,yr,
                         -std::sin(.33)*xr+std::cos(.33)*z);
    }
    return m;
}
DomWallBoundaryConditions boundary(const Mesh& mesh,bool nonuniform) {
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& name=mesh.boundary().patch(p).name;
        walls[name]={name=="bottom"?0.:name=="west"?.45:.85,
            420.+(nonuniform?(name=="east"?20.:name=="west"?-12.:0.):0.)};
    }
    return walls;
}
void exercise(bool high_contrast,bool small,
              std::size_t nx=3,std::size_t ny=3,std::size_t nz=3,double depth=.03) {
    auto mesh=grid(nx,ny,nz,depth);
    auto geometry=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    using Field=DomMomentSchurSystemV44::CellField;
    Field T(mesh.n_cells(),"T","K",1);
    const auto nc=mesh.n_cells();
    std::vector<double> ka(nc),ks(nc);
    for(std::size_t c=0;c<nc;++c) {
        const double x=geometry.cell_centres[c].x;
        T(c)=420.+5.*x;
        if(high_contrast) {
            // Deliberate non-smooth material interface, not just an affine
            // variation in optical properties on a three-dimensional grid.
            const bool thick=(c%3)==0;
            ka[c]=thick?.12:1.75;
            ks[c]=thick?15.:.4;
        } else {
            ka[c]=1.2+.2*std::abs(x);
            ks[c]=.35+.1*std::abs(geometry.cell_centres[c].z);
        }
    }
    RadiationTransportControls settings;
    settings.absorption=1.2;settings.scattering=.35;
    settings.tolerance=1e-9;settings.linear_tolerance=1e-14;
    settings.linear_max_iterations=2000;
    auto wall_bcs=boundary(mesh,true);
    DomMomentSchurSystemV44 ho(mesh,geometry,T,dirs,wall_bcs,settings,ka,ks);
    std::vector<double> x(ho.dimension(),blackbody_intensity(420.));
    for(const auto& wall:ho.wall_faces())
        if(wall.reflect_index<ho.n_reflecting_faces())
            x[nc+wall.reflect_index]=wall.ib;
    DomV51Options options;
    options.base.max_outer=small?5:80;
    options.base.physical_tolerance=1e-9;
    options.base.linear_tolerance=1e-12;
    options.base.max_backtracks=9;
    options.anderson_depth=6;
    options.anderson_ridge=1e-10;
    options.base.preconditioner=DomSmmV49Preconditioner::JacobiAuxiliary;
    // On very large grids, isolate the full HO multisecant accelerator from
    // the currently unqualified LO AMG solve; this is explicitly NOT DSA.
    if(nc>200)options.use_low_order=false;
    auto result=solve_dom_holo_global_v51(ho,x,options);
    std::cout<<"V51_HOLO high_contrast="<<high_contrast
             <<" converged="<<result.converged
             <<" initial="<<result.initial_true<<" final="<<result.final_true
             <<" steps="<<result.accepted_steps
             <<" lo="<<result.lo_steps<<" anderson="<<result.anderson_steps
             <<" picard="<<result.picard_steps
             <<" wall_response_maps="<<result.wall_response_maps
             <<" rejected="<<result.rejected_trials
             <<" maps="<<result.transport_maps
             <<" seconds="<<result.seconds
             <<" diagnostic="<<result.diagnostic<<'\n';
    if(!small) {
        check(result.converged,"V51: no physical convergence");
        check(result.final_true<=options.base.physical_tolerance,
              "V51: solver published physically failing state");
        check(result.anderson_steps>0,
              "V51: multisecant never accepted a coupled gas/wall step");
        if(nc<=200)check(result.lo_solves>0,"V51: low-order correction not exercised");
        check(result.accepted_defect_history.size()==result.accepted_steps+1,
              "V51: accepted history incomplete");
        for(std::size_t k=1;k<result.accepted_defect_history.size();++k)
            check(result.accepted_defect_history[k]<result.accepted_defect_history[k-1],
                  "V51: noncontracting accepted HO state");
        auto final=evaluate_dom_variable_v50(ho,result.final_moments);
        check(final.true_relative<=options.base.physical_tolerance,
              "V51: independent final HO audit rejected");
        double reflector_current=0.;
        for(const auto& w:ho.wall_faces())
            if(w.emissivity==0.)
                reflector_current=std::max(reflector_current,
                    std::abs(final.face_current[w.face]));
        check(reflector_current<1.e-6,
              "V51: pure reflecting gray face has net radiative leakage");
        std::cout<<"V51_REFLECTOR_MAX_WATTS "<<reflector_current<<'\n';
    }
    if(!high_contrast && !small && nc==27) {
        DomV50Options baseline=options.base;
        const auto historical=solve_dom_holo_v50(ho,x,baseline);
        check(historical.converged,"V51 comparator: historical V50 failed");
        double physical_spread=0.;
        for(std::size_t k=0;k<result.final_moments.size();++k)
            physical_spread=std::max(physical_spread,
                std::abs(result.final_moments[k]-historical.final_moments[k])/
                std::max({1.,std::abs(result.final_moments[k]),
                              std::abs(historical.final_moments[k])}));
        check(physical_spread<1.e-7,"V51 comparator: physical states inconsistent");
        check(result.accepted_steps<historical.accepted_steps,
              "V51 comparator: no outer iteration improvement");
        // Different algorithm on the same frozen transport matrices:
        // solve the full exact (J,R) Schur complement by matrix-free GMRES.
        // This checks the accelerated state against a second solver, not
        // against the low-order discretization used to construct directions.
        const auto n=ho.dimension();
        const auto rhs=ho.constant_rhs();
        cfdx::core::LinearOperator op;
        op.size=n;
        op.apply=[&ho,n](const cfdx::core::Vector& v,cfdx::core::Vector& y) {
            std::vector<double> input(n);
            for(std::size_t i=0;i<n;++i)input[i]=v(i);
            auto response=ho.apply_exact_schur(input);
            if(y.size()!=n)y.resize(n);
            for(std::size_t i=0;i<n;++i)y(i)=response[i];
        };
        cfdx::core::Vector b(n,0.),z(n,0.);
        for(std::size_t i=0;i<n;++i){b(i)=rhs[i];z(i)=x[i];}
        const auto exact=cfdx::core::solve_gmres(op,b,z,25,180,1.e-12);
        check(exact.status==cfdx::core::SolverStatus::CONVERGED,
              "V51: independent full Schur GMRES failed");
        double schur_spread=0.;
        for(std::size_t i=0;i<n;++i)
            schur_spread=std::max(schur_spread,
                std::abs(z(i)-result.final_moments[i])/
                std::max({1.,std::abs(z(i)),std::abs(result.final_moments[i])}));
        check(schur_spread<1.e-7,
              "V51: solution differs from independent exact HO Schur");
        std::cout<<"V51_SCHUR_COMPARISON iterations="<<exact.iterations
                 <<" max_relative_spread="<<schur_spread<<'\n';
        std::cout<<"V51_V50_COMPARISON old_steps="<<historical.accepted_steps
                 <<" new_steps="<<result.accepted_steps
                 <<" old_lo="<<historical.lo_accepted_steps
                 <<" old_picard="<<historical.picard_fallbacks
                 <<" old_seconds="<<historical.seconds
                 <<" new_seconds="<<result.seconds
                 <<" solution_spread="<<physical_spread<<'\n';
    }
    // Impossible tolerance must reject all candidates even if accelerated.
    auto impossible=options;
    impossible.base.physical_tolerance=1e-30;
    impossible.base.max_outer=2;
    auto rejected=solve_dom_holo_global_v51(ho,x,impossible);
    check(!rejected.converged && rejected.final_moments.empty(),
          "V51: impossible physical gate published state");
    // Verify that full wall-response derivative is included in HO fixed point.
    std::vector<double> perturb(x.size());
    for(std::size_t k=0;k<perturb.size();++k)
        perturb[k]=0.3+0.05*std::sin(double(k));
    const auto exact=ho.apply_exact_schur(perturb);
    const auto Tx=ho.project(x,true);
    auto xp=x;
    const double h=1e-4;
    for(std::size_t k=0;k<xp.size();++k)xp[k]+=h*perturb[k];
    const auto Tplus=ho.project(xp,true);
    double worst=0.;
    for(std::size_t k=0;k<x.size();++k) {
        const double approximate=perturb[k]-(Tplus[k]-Tx[k])/h;
        worst=std::max(worst,std::abs(approximate-exact[k])/
            std::max({1.,std::abs(exact[k]),std::abs(approximate)}));
    }
    check(worst<2e-7,"V51: exact HO joint cell/gray-wall derivative mismatch");
    std::cout<<"V51_JOINT_JACOBIAN relative="<<worst<<'\n';
}
void degenerate_anderson() {
    using dom_v51_detail::Secant;
    std::deque<Secant> same(4);
    for(auto& h:same) {h.x={1.,2.,3.};h.f={.1,.2,.3};}
    check(dom_v51_detail::direction(same,1e-10,30.).empty(),
          "V51: zero-rank multisecant history unexpectedly accepted");
    bool invalid=false;
    try {DomV51Options o;o.anderson_depth=13;
         auto mesh=grid(3,3,3,.03);const auto geo=build_fv_geometry(mesh);
         DomMomentSchurSystemV44::CellField temp(mesh.n_cells(),"T","K",1);temp.fill(420.);
         RadiationTransportControls controls;controls.absorption=1.;controls.scattering=.1;
         const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
         auto wall_bcs=boundary(mesh,false);
         DomMomentSchurSystemV44 system(mesh,geo,temp,dirs,wall_bcs,controls);
         std::vector<double> x(system.dimension(),blackbody_intensity(420.));
         (void)solve_dom_holo_global_v51(system,x,o);
    }catch(const std::invalid_argument&){invalid=true;}
    check(invalid,"V51: invalid multisecant controls accepted");
}
}
int main(int argc,char** argv) {
    try {
        const std::string mode=argc==1?"--quick":argc==2?argv[1]:"invalid";
        if(mode=="--quick") {degenerate_anderson();exercise(false,false);exercise(true,false);}
        else if(mode=="--sanitizer-small") {degenerate_anderson();exercise(false,true);}
        else if(mode=="--stress-96") exercise(false,false,4,4,6,.012);
        else if(mode=="--stress-200") exercise(false,false,5,5,8,.006);
        else if(mode=="--stress-640") exercise(false,false,8,8,10,.003);
        else throw std::invalid_argument("V51: --quick/--sanitizer-small/--stress-96/--stress-200/--stress-640");
        std::cout<<"V51_ALL_PASS\n";
    } catch(const std::exception& error) {
        std::cerr<<"V51_FAIL: "<<error.what()<<'\n';return 1;
    }
}
