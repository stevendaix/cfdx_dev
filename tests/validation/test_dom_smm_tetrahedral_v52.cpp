#include "cfdx/physics/dom_smm_second_moment_v48.h"
#include "radiation_tetrahedral_mesh_v52.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
namespace {
void check(bool good,const char* message){if(!good)throw std::runtime_error(message);}
void test(std::size_t n){
    auto mesh=cfdx::testing::radiation_v52::tetrahedral_mesh(n);
    const auto geo=build_fv_geometry(mesh);
    const auto nc=mesh.n_cells(),nf=mesh.n_faces();
    DomSmmSecondMomentV48 smm(mesh,geo);
    const double pi=std::acos(-1.);
    const Vec3 grad{.08,-.05,.04};
    const double kt=2.;
    std::vector<double> J(nc),kappa(nc,kt),wall(nf,0.);
    std::vector<DomSmmSecondMomentV48::Tensor> E(nc,
            {1./3.,1./3.,1./3.,0.,0.,0.});
    for(std::size_t c=0;c<nc;++c)J[c]=10.+grad.dot(geo.cell_centres[c]);
    for(std::size_t f=0;f<nf;++f)if(mesh.ownership().neighbour(f)<0)
        wall[f]=-(4.*pi/(3.*kt))*geo.face_area_vectors[f].dot(grad);
    const auto result=smm.evaluate(J,E,kappa,wall);
    double worst=0.;
    for(std::size_t f=0;f<nf;++f){
        const double exact=-(4.*pi/(3.*kt))*geo.face_area_vectors[f].dot(grad);
        worst=std::max(worst,std::abs(result.face_watts[f]-exact)/
                        std::max(1.,std::abs(exact)));
    }
    check(worst<1.e-11,"V52 tetra: affine second-moment face flux not exact");
    check(std::abs(result.internal_cancellation_watts)<1.e-12,
          "V52 tetra: internal face non-conservation");
    // Variable tensor with J=constant: the J*grad(E) term MUST survive.
    const double aa=.012,bb=.008;
    std::fill(J.begin(),J.end(),3.);
    for(std::size_t c=0;c<nc;++c){
        const double x=geo.cell_centres[c].x;
        E[c]={1./3.+aa*x,1./3.-aa*x,1./3.,bb*x,0.,0.};
    }
    const Vec3 divE{aa,bb,0.};
    for(std::size_t f=0;f<nf;++f)if(mesh.ownership().neighbour(f)<0)
        wall[f]=-(4.*pi*3./kt)*geo.face_area_vectors[f].dot(divE);
    const auto variable=smm.evaluate(J,E,kappa,wall);
    double tensor_worst=0.;
    for(std::size_t f=0;f<nf;++f){
        const double exact=-(4.*pi*3./kt)*geo.face_area_vectors[f].dot(divE);
        tensor_worst=std::max(tensor_worst,
            std::abs(variable.face_watts[f]-exact)/std::max(1.,std::abs(exact)));
    }
    check(tensor_worst<1.e-11,"V52 tetra: variable tensor current wrong");
    check(std::abs(variable.internal_cancellation_watts)<1.e-12,
          "V52 tetra: internal variable tensor current not conservative");
    std::cout<<"V52_TET n="<<n<<" cells="<<nc<<" faces="<<nf
             <<" affine_rel="<<worst<<" variable_E_rel="<<tensor_worst
             <<" cancel="<<variable.internal_cancellation_watts<<'\n';
}
}
int main(){
    try{test(3);test(4);std::cout<<"V52_TET_PASS\n";}
    catch(const std::exception& e){std::cerr<<"V52_TET_FAIL: "<<e.what()<<'\n';return 1;}
}
