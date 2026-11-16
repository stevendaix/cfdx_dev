#include "cfdx/physics/dom_quadrature_v39.h"
#include "thermal_3d_mesh_v36.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <array>
using namespace cfdx::physics;
namespace {
void check(bool ok,const char* why){if(!ok)throw std::runtime_error(why);}
template<class F> void reject(F fn){bool yes=false;try{fn();}catch(const std::exception&){yes=true;}check(yes,"bad angular quadrature accepted");}
}
int main(){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(4,5,3,.025,0,.10);
    const auto geom=build_fv_geometry(mesh);
    const auto axes=std::vector<DiscreteDirection>{{1,0,0,2*M_PI/3},{-1,0,0,2*M_PI/3},
        {0,1,0,2*M_PI/3},{0,-1,0,2*M_PI/3},{0,0,1,2*M_PI/3},{0,0,-1,2*M_PI/3}};
    const auto base=audit_dom_wall_quadrature(mesh,geom,axes);
    check(base.maximum_relative_pi_error>.25,"low angular resolution physical error unreported");
    double prior=1.;
    for(auto [nmu,nphi]:std::vector<std::pair<int,int>>{{4,16},{8,32},{12,48},{20,80}}){
        const auto dirs=make_dom_tensor_gauss_legendre_v39(static_cast<std::size_t>(nmu),static_cast<std::size_t>(nphi));
        const auto c=audit_dom_wall_quadrature(mesh,geom,dirs);
        const auto axis=dom_hemisphere_moments(dirs,{0,0,1});
        check(c.faces_checked>0,"no real wall geometry sampled");
        check(std::abs(axis.outgoing-axis.incoming)<1e-10,"antipodal quadrature unbalanced");
        check(c.maximum_relative_pi_error<prior,"refined angular quadrature did not improve hemisphere moment");
        prior=c.maximum_relative_pi_error;
        std::cout<<"DOM_V39 nPolar="<<nmu<<" nAzimuth="<<nphi
                 <<" angles="<<dirs.size()<<" wall_max_rel_pi_error="<<c.maximum_relative_pi_error
                 <<" axis_rel_pi_error="<<axis.relative_pi_error<<'\n';
        if(nmu==20){
            check(c.maximum_relative_pi_error<.01,"refined 20x80 failed 1% hemisphere gate");
            require_dom_wall_quadrature_accuracy(mesh,geom,dirs,.01);
        }
    }
    reject([&]{make_dom_tensor_gauss_legendre_v39(3,16);});
    reject([&]{make_dom_tensor_gauss_legendre_v39(4,6);});
    reject([&]{require_dom_wall_quadrature_accuracy(mesh,geom,axes,.01);});
    return 0;
}
