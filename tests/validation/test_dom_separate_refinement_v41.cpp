// CFDX V41: separate spatial/angle DOM convergence on identical 1D transport
// physics. The exact exponential line-of-sight oracle is evaluated by an
// independent angular integral, NOT by the CFDX finite-volume discretization.
#include "cfdx/physics/dom_quadrature_v39.h"
#include "cfdx/physics/radiation_solver.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
constexpr double PI=3.1415926535897932384626433832795;
constexpr double absorption=0.45, gas=500., hot=850., cold=300.;
void require(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
double oracle_G(double x) {
    constexpr int count=32768;
    const double gas_B=blackbody_intensity(gas);
    long double sum=0.0L;
    for(int m=0;m<count;++m) {
        const double mu=(m+0.5)/count;
        sum+=(blackbody_intensity(hot)-gas_B)*std::exp(-absorption*x/mu)
            +(blackbody_intensity(cold)-gas_B)*std::exp(-absorption*(1.0-x)/mu)
            +2.0*gas_B;
    }
    return static_cast<double>(2.0L*PI*sum/count);
}
struct CaseResult {double total_L2=0, spatial_L2=0, angular_L2=0, wall_pi_error=0;std::size_t iterations=0;};
CaseResult run(std::size_t nx,std::size_t n_polar) {
    auto m=cfdx::testing::thermal_v36::make_mesh(nx,3,3,0.,0.,.2);
    auto g=build_fv_geometry(m);
    Field<double,Location::CELL> T(m.n_cells(),"T","K",1);
    Field<double,Location::CELL> G(m.n_cells(),"G","W/m2",1);
    Field<double,Location::CELL> Q(m.n_cells(),"Q","W/m3",1);
    T.fill(gas);G.fill(4.*blackbody_emissive_power(gas));Q.fill(0.);
    DomWallBoundaryConditions walls;
    for(std::size_t p=0;p<m.boundary().n_patches();++p) {
        auto& patch=m.boundary().patch(p);
        if(patch.name=="west") walls[patch.name]={1.,hot};
        else if(patch.name=="east")walls[patch.name]={1.,cold};
        else {walls[patch.name]={1.,gas};patch.type=PatchType::EMPTY;}
    }
    auto dirs=make_dom_tensor_gauss_legendre_v39(n_polar,8);
    for(auto& d:dirs)std::swap(d.dx,d.dz);
    const auto hemisphere=dom_hemisphere_moments(dirs,{1.,0.,0.});
    RadiationTransportControls c;
    c.absorption=absorption;c.scattering=0;
    c.max_iterations=100;c.tolerance=1.e-9;c.linear_tolerance=1.e-12;
    const auto r=solve_participating_radiation_diffuse_gray_walls(
        m,g,T,G,Q,dirs,c,walls,true,true);
    require(r.converged,"V41 DOM failed to converge");
    require(r.wall_outward_W.count("west")==1 && r.wall_outward_W.count("east")==1,
            "V41 DOM missing independent wall power");
    std::map<double,double> ref_by_x;
    long double error=0.,spatial_error=0.,angular_error=0.,vol=0.;
    for(std::size_t i=0;i<m.n_cells();++i) {
        const double x=g.cell_centres[i].x;
        const auto pos=ref_by_x.emplace(x,0.0);
        if(pos.second)pos.first->second=oracle_G(x);
        const double ref=pos.first->second;
        long double angular_exact=0.0L;
        for(const auto& d:dirs) {
            const double mu=d.dx;
            if(mu>0)angular_exact+=d.weight*(blackbody_intensity(gas) +
                (blackbody_intensity(hot)-blackbody_intensity(gas))*
                std::exp(-absorption*x/mu));
            else if(mu<0)angular_exact+=d.weight*(blackbody_intensity(gas) +
                (blackbody_intensity(cold)-blackbody_intensity(gas))*
                std::exp(-absorption*(1.-x)/(-mu)));
        }
        const double angular_ref=static_cast<double>(angular_exact);
        error+=static_cast<long double>(g.cell_volumes[i])*std::pow((G(i)-ref)/ref,2);
        spatial_error+=static_cast<long double>(g.cell_volumes[i])*
            std::pow((G(i)-angular_ref)/ref,2);
        angular_error+=static_cast<long double>(g.cell_volumes[i])*
            std::pow((angular_ref-ref)/ref,2);
        vol+=g.cell_volumes[i];
        require(std::abs(Q(i)-absorption*(4.*blackbody_emissive_power(gas)-G(i)))<1.e-9,
            "V41 DOM independent radiative source mismatch");
    }
    const double e=std::sqrt(static_cast<double>(error/vol));
    const double es=std::sqrt(static_cast<double>(spatial_error/vol));
    const double ea=std::sqrt(static_cast<double>(angular_error/vol));
    std::cout<<std::scientific<<std::setprecision(10)
        <<"DOM_V41 nx="<<nx<<" polar="<<n_polar
        <<" directions="<<dirs.size()<<" total_L2="<<e
        <<" spatial_L2="<<es<<" angular_L2="<<ea
        <<" hemi_pi_error="<<hemisphere.relative_pi_error
        <<" iterations="<<r.iterations
        <<" west_W="<<r.wall_outward_W.at("west")
        <<" east_W="<<r.wall_outward_W.at("east")<<'\n';
    return {e,es,ea,hemisphere.relative_pi_error,r.iterations};
}
}
int main() {
    try {
        // Spatial campaign: identical quadrature and angular error at each h.
        const auto x12=run(12,16),x24=run(24,16),x48=run(48,16);
        require(x24.spatial_L2<x12.spatial_L2 && x48.spatial_L2<x24.spatial_L2,
                "DOM V41 spatial convergence fails versus discrete-angular exact oracle");
        // Angular campaign: identical FV system at every angular level.
        const auto a8=run(48,8),a12=run(48,12),a24=run(48,24);
        require(a12.angular_L2<a8.angular_L2 && x48.angular_L2<a12.angular_L2 &&
                a24.angular_L2<x48.angular_L2,
                "DOM V41 angular quadrature fails versus continuous exact oracle");
        require(a24.total_L2<0.01,"DOM V41 finest combined error insufficient");
        std::cout<<"DOM_V41_SEPARATE_CONVERGENCE PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<"DOM_V41_SEPARATE_CONVERGENCE FAIL: "<<e.what()<<'\n';
        return 1;
    }
}
