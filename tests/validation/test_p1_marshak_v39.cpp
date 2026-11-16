// Verify the *production* finite-volume P1 operator with its new Marshak
// boundary closure against analytical gray 1D slab and cell-wise FV balance.
#include "cfdx/physics/p1_marshak_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <numeric>
#include <iomanip>
using namespace cfdx::physics;
using namespace cfdx::core;
namespace {
void check(bool yes,const char* why){if(!yes)throw std::runtime_error(why);}
void close(double a,double b,double atol,double rtol=0){
    if(!(std::abs(a-b)<=atol+rtol*std::max(std::abs(a),std::abs(b))))
        throw std::runtime_error("P1 Marshak oracle mismatch");
}
template<class F> void rejects(F f){
    bool bad=false;try{f();}catch(const std::exception&){bad=true;}
    check(bad,"invalid P1 Marshak case accepted");
}
P1MarshakWalls make_walls(const Mesh& mesh,double west_T,double east_T,double west_eps,double east_eps){
    P1MarshakWalls walls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto& name=mesh.boundary().patch(p).name;
        walls.emplace(name,P1MarshakWall{0.,600.});
    }
    walls.at("west")={west_eps,west_T};
    walls.at("east")={east_eps,east_T};
    return walls;
}
// Independent continuous 1-D slab oracle, D G''=a(G-Gb), x in [0,1].
// q_out(0)=D*G'(0), q_out(1)=-D*G'(1), both = c(Gwallface-Gw).
struct SlabOracle {
    double A,B,base,lam,D,c0,c1,G0,G1;
    double value(double x)const{return base + A*std::cosh(lam*x)+B*std::sinh(lam*x);}
    double left_flux()const{return c0*(value(0)-G0);}
    double right_flux()const{return c1*(value(1)-G1);}
};
SlabOracle oracle(double abs,double scat,double T,double west_T,double east_T,double ep0,double ep1){
    const double D=1./(3.*(abs+scat));
    const double lam=std::sqrt(abs/D);
    const double c0=ep0/(2.*(2.-ep0)),c1=ep1/(2.*(2.-ep1));
    const double Gbase=4*blackbody_emissive_power(T);
    const double Gw0=4*blackbody_emissive_power(west_T);
    const double Gw1=4*blackbody_emissive_power(east_T);
    const double m00=-c0, m01=D*lam;
    const double m10=-D*lam*std::sinh(lam)-c1*std::cosh(lam);
    const double m11=-D*lam*std::cosh(lam)-c1*std::sinh(lam);
    const double rhs0=c0*(Gbase-Gw0),rhs1=c1*(Gbase-Gw1);
    const double det=m00*m11-m01*m10;
    check(std::abs(det)>1e-14,"analytic slab singular");
    return {(rhs0*m11-rhs1*m01)/det,(m00*rhs1-m10*rhs0)/det,
            Gbase,lam,D,c0,c1,Gw0,Gw1};
}
struct Measured{double err=0, fluxerr=0, balance=0;};
Measured slab(std::size_t nx,double warp=0.){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(nx,3,3,warp,0,0.15);
    auto g=build_fv_geometry(mesh);
    Field<double,Location::CELL> T(mesh.n_cells(),"T","K",1);
    Field<double,Location::CELL> G(mesh.n_cells(),"G","W/m2",1);
    Field<double,Location::CELL> Q(mesh.n_cells(),"qr","W/m3",1);
    T.fill(600.);G.fill(4*blackbody_emissive_power(600.));Q.fill(-999.);
    const double abs=0.9,scat=2.5;
    std::vector<double> a(mesh.n_cells(),abs),s(mesh.n_cells(),scat);
    const auto walls=make_walls(mesh,900.,300.,.8,.65);
    const auto r=solve_p1_radiation_marshak(mesh,g,T,G,Q,a,s,walls);
    if(!r.converged){
        std::cerr << "MARSH_FAIL outer="<<r.outer_iterations
            <<" phys="<<r.maximum_normalized_cell_residual
            <<" balance="<<r.relative_global_balance<<"\n";
    }
    check(r.converged,"3D P1 Marshak slab did not converge");
    check(r.maximum_normalized_cell_residual<=1.e-10 &&
          r.relative_global_balance<=1.e-10,"physical P1 gates not closed");
    const auto ref=oracle(abs,scat,600,900,300,.8,.65);
    double error2=0,vol=0;
    for(std::size_t i=0;i<mesh.n_cells();++i){
        const double e=G(i)-ref.value(g.cell_centres[i].x);
        error2+=e*e*g.cell_volumes[i];vol+=g.cell_volumes[i];
        close(Q(i),abs*(4*blackbody_emissive_power(600.)-G(i)),1e-12);
    }
    double error=std::sqrt(error2/vol);
    const double qwest=r.patch_outward_W.at("west")/.15;
    const double qeast=r.patch_outward_W.at("east")/.15;
    const double fluxerror=std::max(std::abs(qwest-ref.left_flux()),
                                     std::abs(qeast-ref.right_flux()));
    close(r.total_boundary_outward_W,r.volume_emission_minus_absorption_W,1e-8,1e-10);
    std::cout<<std::scientific<<std::setprecision(11)
             <<"P1_MARSHAK_SLAB nx="<<nx<<" warp="<<warp<<" L2="<<error
             <<" q_err="<<fluxerror<<" physical="<<r.maximum_normalized_cell_residual
             <<" balance="<<r.relative_global_balance<<" outer="<<r.outer_iterations<<'\n';
    return {error,fluxerror,r.relative_global_balance};
}
void equilibrium_and_rejection(){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(4,3,3,0.014,0,0.06);
    auto g=build_fv_geometry(mesh);
    Field<double,Location::CELL> T(mesh.n_cells(),"T","K",1);
    Field<double,Location::CELL> G(mesh.n_cells(),"G","W/m2",1);
    Field<double,Location::CELL> Q(mesh.n_cells(),"qr","W/m3",1);
    T.fill(700.);G.fill(4*blackbody_emissive_power(700.));Q.fill(-15.);
    std::vector<double> a(mesh.n_cells()),s(mesh.n_cells());
    for(std::size_t i=0;i<mesh.n_cells();++i){a[i]=0.2+0.7*(double(i+1)/double(mesh.n_cells()));s[i]=1.+double(i)*.01;}
    auto walls=make_walls(mesh,700,700,.4,.8);
    const auto result=solve_p1_radiation_marshak(mesh,g,T,G,Q,a,s,walls);
    check(result.converged,"variable opacity isothermal P1 not converged");
    const double expected=4*blackbody_emissive_power(700.);
    for(std::size_t i=0;i<G.size();++i){close(G(i),expected,2e-8);close(Q(i),0,2e-8);}
    check(std::abs(result.total_boundary_outward_W)<1e-7 &&
          result.relative_global_balance<1e-8,"variable opacity P1 isothermal balance");
    auto previously_accepted_G=G, previously_accepted_Q=Q;
    rejects([&]{auto bad=walls;bad.at("west").emissivity=-.1;
        solve_p1_radiation_marshak(mesh,g,T,G,Q,a,s,bad);});
    rejects([&]{auto bad=walls;bad.erase("west");
        solve_p1_radiation_marshak(mesh,g,T,G,Q,a,s,bad);});
    rejects([&]{auto bad=walls;for(auto& [name,bc]:bad)bc.emissivity=0;
        auto zeros=a;std::fill(zeros.begin(),zeros.end(),0.);
        solve_p1_radiation_marshak(mesh,g,T,G,Q,zeros,s,bad);});
    rejects([&]{auto bad=a;bad[0]=-1.;
        solve_p1_radiation_marshak(mesh,g,T,G,Q,bad,s,walls);});
    rejects([&]{auto bad=g;
        const auto f=mesh.boundary().patch(0).face_ids.front();
        bad.face_area_vectors[f]=bad.face_area_vectors[f]*(-1.0);
        solve_p1_radiation_marshak(mesh,bad,T,G,Q,a,s,walls);});
    for(std::size_t i=0;i<G.size();++i){
        check(G(i)==previously_accepted_G(i) && Q(i)==previously_accepted_Q(i),
              "P1 Marshak failed input mutated committed fields");
    }
    std::cout<<"P1_MARSHAK_EQUILIBRIUM variable opacity 3D skew passed\n";
}
void rollback_after_nonconvergence(){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(3,3,3,.022,0,.10);
    auto g=build_fv_geometry(mesh);
    Field<double,Location::CELL> T(mesh.n_cells(),"T","K",1);
    Field<double,Location::CELL> G(mesh.n_cells(),"G","W/m2",1);
    Field<double,Location::CELL> Q(mesh.n_cells(),"q","W/m3",1);
    T.fill(500.);G.fill(4*blackbody_emissive_power(500.));Q.fill(123.);
    auto walls=make_walls(mesh,1000.,300.,.8,.9);
    std::vector<double>a(mesh.n_cells(),.8),s(mesh.n_cells(),2.);
    P1MarshakControls ctrl;ctrl.max_outer_iterations=1;
    ctrl.physical_tolerance=1.e-17;ctrl.balance_tolerance=1.e-17;
    const auto r=solve_p1_radiation_marshak(mesh,g,T,G,Q,a,s,walls,ctrl);
    check(!r.converged,"P1 Marshak overly strict fixed-point should reject");
    for(std::size_t i=0;i<G.size();++i){
        check(G(i)==4*blackbody_emissive_power(500.) && Q(i)==123.,
              "P1 Marshak rejection mutated committed state");
    }
    std::cout<<"P1_MARSHAK_ROLLBACK passes\n";
}
}
int main(){
    const auto x=slab(8),y=slab(16),z=slab(32);
    const double p0=std::log(x.err/y.err)/std::log(2.);
    const double p1=std::log(y.err/z.err)/std::log(2.);
    std::cout<<"P1_MARSHAK_ORDER "<<p0<<" "<<p1<<'\n';
    check(p0>1.75 && p1>1.75,"P1 Marshak did not achieve second order on 3D grid");
    equilibrium_and_rejection();
    rollback_after_nonconvergence();
    return 0;
}
