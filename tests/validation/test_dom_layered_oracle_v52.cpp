// Independent analytic ray oracle for a two-layer 1D gray absorbing slab.
// No CFDX solver is called: this certifies the reference and the first-order
// upwind spatial trend, NOT the full heterogeneous 3D DOM implementation.
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
namespace {
constexpr double pi=3.14159265358979323846;
struct Slab {double k1,k2,B1,B2,left,right;};
void require(bool valid,const char* message){if(!valid)throw std::runtime_error(message);}
// Exact characteristic integration across a layer of constant absorption.
double exact_layer(double incoming,double B,double k,double length,double mu){
    return B+(incoming-B)*std::exp(-k*length/mu);
}
double discrete_layer(double incoming,double B,double k,double length,double mu,
                      unsigned n){
    const double optical=k*length/n;
    for(unsigned i=0;i<n;++i)
        incoming=(mu*incoming+optical*B)/(mu+optical);
    return incoming;
}
double flux(const Slab& s,unsigned angles,unsigned ncell=0){
    require(angles>=4,"V52 oracle: invalid angular resolution");
    constexpr double l1=.45,l2=.55;
    const unsigned n1=static_cast<unsigned>(ncell*l1+.5);
    const unsigned n2=ncell-n1;
    if(ncell!=0)require(n1>0 && n2>0,"V52 oracle: invalid layer counts");
    long double total=0.;
    for(unsigned a=0;a<angles;++a){
        const double mu=(a+.5)/double(angles);
        const auto plus=ncell==0 ? exact_layer(s.left,s.B1,s.k1,l1,mu):
            discrete_layer(s.left,s.B1,s.k1,l1,mu,n1);
        const auto minus=ncell==0 ? exact_layer(s.right,s.B2,s.k2,l2,mu):
            discrete_layer(s.right,s.B2,s.k2,l2,mu,n2);
        total+=mu*(plus-minus);
    }
    return 2.*pi*static_cast<double>(total)/angles;
}
void run(const Slab& s,const char* label){
    const auto oracle=flux(s,32768);
    require(std::isfinite(oracle) && std::abs(oracle)>1.,
            "V52 oracle: degenerate flux");
    // Angular and spatial errors are monitored SEPARATELY; coarse angular
    // quadrature may cancel spatial error and must never certify it.
    const auto angle_reference=flux(s,512);
    const auto quadrature_error=std::abs(angle_reference-oracle)/std::abs(oracle);
    require(quadrature_error<5.e-5,"V52 oracle: reference quadrature inaccurate");
    double last=1.e9;
    for(const unsigned cells:{20u,40u,80u,160u}){
        const auto fv=flux(s,512,cells);
        const auto spatial_error=std::abs(fv-angle_reference)/
                                  std::max(1.,std::abs(angle_reference));
        const auto total_error=std::abs(fv-oracle)/
                                std::max(1.,std::abs(oracle));
        std::cout<<"V52_LAYER name="<<label<<" cells="<<cells
                 <<" spatial_rel="<<spatial_error
                 <<" total_rel="<<total_error<<'\n';
        require(spatial_error<last*.86,
                "V52 oracle: first-order upwind layer flux did not contract");
        last=spatial_error;
    }
    std::cout<<"V52_LAYER_REF name="<<label<<" q_ref="<<oracle
             <<" angular_rel="<<quadrature_error<<'\n';
}
}
int main(){
    try{
        run({.25,3.0,240.,420.,1100.,95.},"moderate_jump");
        run({.03,14.,250.,540.,1200.,85.},"strong_jump");
        std::cout<<"V52_LAYER_ORACLE_PASS\n";
    }catch(const std::exception& e){
        std::cerr<<"V52_LAYER_ORACLE_FAIL: "<<e.what()<<'\n';return 1;
    }
}
