#include "chain_fixture.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
using namespace cfdx::core;
using namespace cfdx::physics;
int main(){
 auto m=make_1d_chain(5);
 for(std::size_t i=0;i<m.n_points();++i){double x=m.points().x(i);m.points().set(i,x*x,m.points().y(i),m.points().z(i));}
 auto g=build_fv_geometry(m);
 Field<double,Location::CELL> p(5,"p","",1),su(5,"Su","",1),sp(5,"Sp","",1);
 for(std::size_t c=0;c<5;++c){p(c)=std::sin(11*g.cell_centres[c].x)+.1*c;su(c)=.2*(c+1);sp(c)=-.1;}
 Field<double,Location::FACE> f(m.n_faces(),"F","",1);
 for(std::size_t i=0;i<f.size();++i)f(i)=g.face_area_vectors[i].x;
 std::vector<double> k{1,4,.2,2,3};
 ScalarBoundaryConditions bc;bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,1.3,0};bc["outlet"]={ScalarBoundaryType::FIXED_VALUE,-.4,0};bc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
 int total=0,failed=0;
 for(bool neumann:{false,true})for(double direction:{-1.0,1.0}){
 if(neumann){bc["inlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,.7};bc["outlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,-.4};}
 else{bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,1.3,0};bc["outlet"]={ScalarBoundaryType::FIXED_VALUE,-.4,0};}
 for(std::size_t i=0;i<f.size();++i)f(i)=direction*g.face_area_vectors[i].x;
 for(bool bounded:{false,true})for(auto scheme:{ConvectionScheme::UPWIND,ConvectionScheme::CENTRAL,ConvectionScheme::BLENDED,ConvectionScheme::QUICK,ConvectionScheme::QUICK_BOUNDED,ConvectionScheme::SECOND_ORDER_UPWIND,ConvectionScheme::TVD}){
  auto eq=assemble_scalar_equation(m,g,f,.7,su,sp,bc,bounded,nullptr,nullptr,nullptr,&k,scheme,&p);
  auto flux=reconstruct_scalar_transport_flux(m,g,f,p,.7,bc,bounded,scheme,nullptr,&k);
  Vector x(5,0);for(size_t c=0;c<5;++c)x(c)=p(c);auto ax=eq.matrix.matvec(x);
  std::vector<double> balance(5,0),divf(5,0);
  for(size_t i=0;i<m.n_faces();++i){auto o=m.ownership().owner(i);auto n=m.ownership().neighbour(i);balance[o]+=flux(i);divf[o]+=f(i);if(n>=0){balance[n]-=flux(i);divf[n]-=f(i);}}
  double e=0;for(size_t c=0;c<5;++c){double r=balance[c]-g.cell_volumes[c]*(su(c)+sp(c)*p(c))-(bounded?p(c)*divf[c]:0.0);double error=std::abs(r-(ax[c]-eq.rhs(c)));e=std::isfinite(error)?std::max(e,error):std::numeric_limits<double>::infinity();}
  bool ok=std::isfinite(e)&&e<1e-11;++total;failed+=!ok;
  std::cout<<"flux-contract\t"<<int(scheme)<<"\tneumann="<<neumann<<"\tdirection="<<direction<<"\tbounded="<<bounded<<"\t"<<(ok?"PASS":"FAIL")<<"\t"<<std::setprecision(17)<<e<<'\n';
 }
 }
 std::cout<<"SUMMARY\tchecks="<<total<<"\tfailures="<<failed<<'\n';return failed?1:0;
}
