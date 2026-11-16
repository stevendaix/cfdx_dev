// V40: independent affine radiation manufactured solution on a skew 3D grid.
// Optically scattering, nonabsorbing medium; every boundary uses a face-wise
// diffuse-gray Marshak emissivity and prescribed temperature derived from the
// CONTINUOUS affine G and normal flux. No solution values are injected into the
// numerical operator except through physically specified wall temperatures.
#include "cfdx/physics/p1_marshak_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void ensure(bool good,const char* text){if(!good)throw std::runtime_error(text);}
double scalar(Vec3 g,Vec3 x){return g.dot(x);}
void affine(std::size_t n,double warp,bool skew){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(n,3,3,warp,0.,.25);
    auto geom=build_fv_geometry(mesh);
    const Vec3 gradient{120.,-75.,55.};
    const double Gbase=3000.,D=1./(3.*1.7),eps=.7;
    const double c=eps/(2.*(2.-eps));
    Field<double,Location::CELL> T(mesh.n_cells(),"T","K",1);
    Field<double,Location::CELL> G(mesh.n_cells(),"G","W/m2",1);
    Field<double,Location::CELL> Q(mesh.n_cells(),"Q","W/m3",1);
    T.fill(620.);G.fill(Gbase);Q.fill(-123.);
    P1MarshakWalls walls;
    P1MarshakFaceWalls facewalls;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p) {
        const auto& patch=mesh.boundary().patch(p);
        walls[patch.name]={eps,610.};
        for(const auto f:patch.face_ids) {
            const double area=geom.face_area_vectors[f].mag();
            const double Gface=Gbase+scalar(gradient,geom.face_centres[f]);
            const double desired_q=-D*scalar(gradient,geom.face_area_vectors[f]);
            const double Gwall=Gface-desired_q/(area*c);
            ensure(Gwall>0,"manufactured wall must be radiative-positive");
            const double Tw=std::pow(Gwall/(4.*5.670374419e-8),.25);
            facewalls[f]={eps,Tw};
        }
    }
    std::vector<double> abs(mesh.n_cells(),0.),scat(mesh.n_cells(),1.7);
    P1MarshakControls ctl;ctl.skew_marshak_robin_v40=skew;
    ctl.max_outer_iterations=90;ctl.robin_gradient_sweeps=16;
    ctl.physical_tolerance=1.e-10;ctl.balance_tolerance=1.e-10;
    auto result=solve_p1_radiation_marshak(mesh,geom,T,G,Q,abs,scat,walls,ctl,&facewalls);
    ensure(result.converged,"affine P1 nonconverged");
    double max_error=0.,max_face_error=0.;
    for(std::size_t i=0;i<mesh.n_cells();++i){
        max_error=std::max(max_error,std::abs(G(i)-Gbase-scalar(gradient,geom.cell_centres[i])));
        ensure(Q(i)==0.0,"nonabsorbing P1 must supply zero volumetric source");
    }
    for(const auto& [f,watts]:result.face_outward_W){
        const double expected=-D*scalar(gradient,geom.face_area_vectors[f]);
        max_face_error=std::max(max_face_error,std::abs(watts-expected));
    }
    std::cout<<std::scientific<<std::setprecision(11)
      <<"P1_AFFINE_V40 n="<<n<<" warp="<<warp<<" skew="<<skew
      <<" max_temperature_like_G_error="<<max_error
      <<" max_wall_watts_error="<<max_face_error
      <<" balance="<<result.relative_global_balance
      <<" robin_sweeps="<<result.maximum_robin_gradient_sweeps
      <<" outer="<<result.outer_iterations<<"\n";
    if(skew){
        ensure(max_error<1.e-6,"affine P1 skew error too large");
        ensure(max_face_error<2.e-7,"affine P1 skew face-flux error too large");
    }
}
}
int main(){
    affine(8,0.0,true);
    affine(8,0.032,true);
    affine(16,0.03,true);
    affine(8,0.03,false); // A/B diagnostic, not an acceptance oracle
}
