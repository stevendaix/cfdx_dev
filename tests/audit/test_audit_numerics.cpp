// Independent regression cases for the 2026-10-24 audit.
// Checks remain active with NDEBUG. Reference constants are analytic, not
// extracted from the assembled matrix or reconstructed by the tested code.
#include "affine_mesh_fixture.h"
#include "chain_fixture.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/linalg/jacobi_solver.h"
#include "cfdx/core/linalg/gauss_seidel_solver.h"
#include "cfdx/core/linalg/matrix_free_fv_operator.h"
#include "cfdx/core/numerics/newton.h"
#include "cfdx/core/numerics/interpolation.h"
#include "cfdx/physics/sst_solver.h"
#include "cfdx/physics/thermal_zones.h"
#include "cfdx/physics/thermophysical_models.h"
#include "cfdx/physics/transport_models.h"
#include "cfdx/physics/m1_m4_models.h"
#include "cfdx/physics/low_storage_time_integration.h"
#include "cfdx/physics/dual_time_field_integrator.h"
#include "cfdx/thermodynamics/equation_of_state.h"
#include <limits>
using namespace cfdx::core;
using namespace cfdx::physics;
using Cell=Field<double,Location::CELL>;
using Face=Field<double,Location::FACE>;
static int checks=0, failures=0;
static double max_error(double a,double b){
    // Fail closed: std::max alone can hide a NaN comparison.
    if(!std::isfinite(a)||!std::isfinite(b))return std::numeric_limits<double>::infinity();
    return std::max(a,b);
}
static void check(const std::string& name,bool ok,double value=0){++checks;failures+=!ok;std::cout<<name<<'\t'<<(ok?"PASS":"FAIL")<<'\t'<<std::setprecision(17)<<value<<'\n';}
static void near(const std::string& name,double a,double b,double tol=1e-11){check(name,std::isfinite(a)&&std::abs(a-b)<=tol*std::max(1.0,std::abs(b)),a);}
template<class F>void rejects(const std::string& name,F f){bool yes=false;try{f();}catch(const std::exception&){yes=true;}check(name,yes);}
static SparseMatrix scalar(double v){SparseMatrix a(1,1);a.push_back(0,0,v);a.finalize();return a;}
int main(){
 try {
 for(double s:{1.0,1e200,1e-200}){Vector v(1,s);near("norm/vector",v.norm2()/s,1);near("norm/krylov",krylov_norm2(v)/s,1);}
 {Vector v(1,std::numeric_limits<double>::quiet_NaN());check("norm/NaN",std::isnan(v.norm_inf()));}
 for(int j:{0,1}){Vector b(1,1e200),x(1,0);auto r=j?solve_gauss_seidel(scalar(1),b,x,10,1e-10):solve_jacobi(scalar(1),b,x,10,1e-10);check("stationary/huge/status",r.status==SolverStatus::CONVERGED);near("stationary/huge/solution",x(0)/b(0),1);}
 for(double s:{1.,1e-200}){Cell u(1,"u","",1);u(0)=0;NewtonFieldControls ctl;ctl.max_iterations=1;auto r=solve_newton_field(u,[s](const Cell& q,Cell& res,SparseMatrix& J){res(0)=s*(q(0)-1);J=scalar(s);},ctl);check("newton/last-iteration/status",r.status==SolverStatus::CONVERGED);near("newton/solution",u(0),1);}
 for(double s:{1.,1e12,1e24,1e-100,1e-200}){Vector b(1,s),x(1,0);auto r=solve_cg(scalar(s),b,x,20,1e-10);check("CG/scale/status",r.status==SolverStatus::CONVERGED);near("CG/scale/solution",x(0),1);}
 for(int which:{0,1,2,3}){PrecisionPolicy p;p.enabled=true;Vector b(1,1),x(1,0);auto a=scalar(1.00000001);SolverResult r;
  if(which==0)r=solve_cg(a,b,x,30,1e-12,p);
  if(which==1)r=solve_gmres(a,b,x,1,30,1e-12,nullptr,p);
  if(which==2)r=solve_fgmres(a,b,x,1,30,1e-12,nullptr,p);
  if(which==3)r=solve_bicgstab(a,b,x,30,1e-12,nullptr,p);
  check("mixed/status",r.status==SolverStatus::CONVERGED);near("mixed/FP64-residual",1-a.matvec(x)[0],0,1e-12);
 }
 {Vector b(1,0),x(1,1);auto r=solve_gmres(scalar(1),b,x,1,5,1e-12);check("GMRES/zero-rhs",r.status==SolverStatus::CONVERGED);near("GMRES/zero-rhs/x",x(0),0);}
 for(double s:{1.,1e-16}){Vector b(1,s),x(1,0);auto r=solve_bicgstab(scalar(s),b,x,20,1e-12);check("BiCG/scale",r.status==SolverStatus::CONVERGED);near("BiCG/scale/x",x(0),1);}
 {auto a=scalar(1),g=scalar(1),d=scalar(-1),c=scalar(0);BlockOperator op(a,g,d,c);ExactSchurApproximation sh([](const Vector& b,Vector& x){x=b;return true;});check("Schur/setup",sh.setup(op));Vector b(1,std::numeric_limits<double>::quiet_NaN()),x(1,0);check("Schur/reject-NaN",!sh.apply(b,x));}
 {auto grid=make_affine_cube(3,0,1);auto& m=grid.mesh;Cell p(m.n_cells(),"p","",1);for(size_t c=0;c<m.n_cells();++c)p(c)=grid.geometry.cell_centres[c].x;
  auto g=compute_gradient_weighted_least_squares(p,m,GradientWeighting::INVERSE_DISTANCE_SQUARED,std::numeric_limits<double>::infinity(),BoundaryGradientPolicy::EXCLUDE_BOUNDARY);double e=0;for(size_t c=0;c<m.n_cells();++c)e=max_error(e,std::abs(g(c,0)-1));near("WLS/exclude/all-cells",e,0);
  auto ghost=compute_gradient_least_squares(p,m,BoundaryGradientPolicy::ZERO_GRADIENT_GHOST);near("LS/ghost/corner",ghost(0,0),0.5);
  FvDiffusionOperator zero(m,grid.geometry,0);Vector diag;zero.diagonal(diag);near("matrixfree/zero-diagonal",diag(0),0);
  rejects("matrixfree/reject-NaN",[&]{FvDiffusionOperator bad(m,grid.geometry,std::numeric_limits<double>::quiet_NaN());});
 }
 for(size_t n:{2u,3u}){auto grid=make_affine_cube(n,0,1);Cell p(grid.mesh.n_cells(),"p","",1);for(size_t c=0;c<p.size();++c)p(c)=grid.geometry.cell_centres[c].x;auto q=compute_quadratic_fit_least_squares(p,grid.mesh);double e=0;for(size_t c=0;c<p.size();++c)e=max_error(e,std::abs(q(c,0)-1));near("QUICK/fallback-affine",e,0);}
 for(double l:{1.,1e-6,1e-8}){auto grid=make_affine_cube(4,0,1);auto& m=grid.mesh;for(size_t i=0;i<m.n_points();++i)m.points().set(i,l*m.points().x(i),l*m.points().y(i),l*m.points().z(i));auto geo=make_geometry_cache(m);Cell p(m.n_cells(),"p","",1);for(size_t c=0;c<p.size();++c)p(c)=std::pow(geo.cell_centres[c].x/l,2);auto q=compute_quadratic_fit_least_squares(p,m);double e=0;for(size_t c=0;c<p.size();++c)e=max_error(e,std::abs(l*l*q(c,3)-1));near("QUICK/scale/curvature",e,0);}
 {TurbulenceTransportControls c;auto s=sst_cell_sources(1,1,100,1,1,0,c);near("SST/high-strain",s.sw,161.29032258064516);auto lo=sst_cell_sources(1,1,0.01,1,1,0,c);near("SST/low-strain",lo.sw,(5.0/9.0)*1e-4);
  auto w=komega2006_cell_sources(1,1,100,0,c);near("Wilcox/limited-viscosity",w.nut,0.0034285714285714286);near("Wilcox/unlimited-diffusivity",w.gamma_k,c.density*(c.molecular_viscosity+c.komega_sigma_k));
 }
 {ThermophysicalProperties p;auto& c=p.heat_capacity;c.model=ScalarPropertyModel::TABLE;c.table.temperature={300,400};c.table.value={1000,1200};
  near("enthalpy/clamp/below",p.enthalpy(200,100),100000);near("enthalpy/clamp/above",p.enthalpy(600,500),120000);near("enthalpy/zero",p.enthalpy(200,200),0);
  c.table.extrapolation=ExtrapolationPolicy::LINEAR;near("enthalpy/linear/below",p.enthalpy(200,100),70000);near("enthalpy/linear/above",p.enthalpy(600,500),150000);near("enthalpy/antisymmetry",p.enthalpy(100,600)+p.enthalpy(600,100),0,1e-8);
  near("enthalpy/additivity",p.enthalpy(600,100)-p.enthalpy(350,100)-p.enthalpy(600,350),0,1e-8);
  c.enforce_bounds=true;c.minimum=700;c.maximum=1400;c.model=ScalarPropertyModel::LINEAR;c.reference_temperature=300;c.reference_value=1000;c.coefficients={5};near("enthalpy/clipped-linear",p.enthalpy(500,200),343000);
  c.model=ScalarPropertyModel::POLYNOMIAL;c.coefficients={1000};c.maximum=800;near("enthalpy/clipped-polynomial",p.enthalpy(400,300),80000);
  c.model=ScalarPropertyModel::CONSTANT;c.reference_value=1000;c.maximum=800;near("enthalpy/clipped-constant",p.enthalpy(400,300),80000);
  rejects("enthalpy/NaN-interval-count",[&]{p.enthalpy(400,300,std::numeric_limits<double>::quiet_NaN());});
 }
 {auto t=compute_transport(300,1.2);near("transport/air-Cp",t.Cp,1.4*AIR_R/0.4);cfdx::thermodynamics::IdealGasEOS eos;near("EOS/dp-dT",eos.dp_dT_rho(101325,300),337.75);near("P1/full-solid-angle",m1m4::p1_source(1,0,1000),226814.97676,1e-10);}
 for(int method:{0,1,2}){
  auto step=[method](Cell& u,double dt,auto&& rhs){if(method==0)low_storage_rk2_step(u,dt,rhs);if(method==1)low_storage_rk3_step(u,dt,rhs);if(method==2)ssprk3_step(u,dt,rhs);};
  Cell u(1,"u","",3);for(size_t d=0;d<3;++d)u(0,d)=d+1;
  step(u,.1,[](const Cell& x,Cell& r){r=Cell(x.size(),"r","",x.dimension());r.fill(1);});
  double e=0;for(size_t d=0;d<3;++d)e=max_error(e,std::abs(u(0,d)-(d+1.1)));near("RK/all-components",e,0);
  double old=0;for(int n:{10,20,40}){Cell q(1,"q","",3);q.fill(1);for(int i=0;i<n;++i)step(q,1.0/n,[](const Cell& x,Cell& r){r=x;});double err=std::abs(q(0,2)-std::exp(1.0));if(old>0){double order=std::log(old/err)/std::log(2.0);check("RK/observed-order",order>(method==0?1.9:2.9),order);}old=err;}
 }
 {Cell u(1,"u","",1);u(0)=1;rejects("implicit/reject-nonfinite",[&]{advance_time(u,1,[](const Cell& q,Cell& r){r(0)=-2*std::sqrt(q(0));},TimeScheme::EULER_IMPLICIT);});}
 {Cell u(1,"u","",1);u.fill(0);DualTimeFieldIntegrator it(u);DualTimeFieldSolveControls c;c.use_embedded_be_estimator=false;c.pseudo_time.pseudo_dt_initial=.5;c.pseudo_time.pseudo_dt_max=.5;c.pseudo_time.absolute_tolerance=1e-13;c.pseudo_time.relative_tolerance=1e-13;
  double t=0;for(double dt:{.1,.2,.05,.08}){auto r=it.advance(dt,c,[](const Cell&,Cell& rhs){rhs(0)=1;});t+=dt;near("dual-time/variable-step",r.state(0),t,1e-11);}
 }
 {auto m=make_1d_chain(2);for(size_t i=0;i<m.n_points();++i){double x=m.points().x(i);x=x<=.5?.5*x:.25+1.5*(x-.5);m.points().set(i,x,m.points().y(i),m.points().z(i));}auto g=build_fv_geometry(m);Cell p(2,"p","",1),z(2,"z","",1);z.fill(0);for(size_t c=0;c<2;++c)p(c)=g.cell_centres[c].x;Face f(m.n_faces(),"f","",1);f.fill(0);auto q=interpolate_cell_to_face(p,m,InterpScheme::LINEAR);size_t inner=0;for(size_t i=0;i<m.n_faces();++i)if(m.ownership().neighbour(i)>=0)inner=i;near("linear/nonuniform",q(inner),.25);
  std::vector<double> k{1,4};auto a=assemble_scalar_equation(m,g,f,1,z,z,{},false,nullptr,nullptr,nullptr,&k);near("diffusion/material-series",-a.matrix(0,1),32.0/7.0);
  auto model=build_thermal_zones(m,{0,1},{{"a",0,1,0,.3,true},{"b",1,4,0,.3,true}});ScalarBoundaryConditions bc;bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,300,0};auto at=assemble_thermal_conduction(m,g,model,bc);near("thermal/contact-series",-at.matrix(0,1),1/(.125+.375/4+.3));
  f(inner)=1;auto b=assemble_scalar_equation(m,g,f,1,z,z,{},false);near("transport/local-imbalance",b.max_imbalance,1);
  k[0]=std::numeric_limits<double>::quiet_NaN();rejects("diffusion/reject-NaN-material",[&]{assemble_scalar_equation(m,g,f,1,z,z,{},false,nullptr,nullptr,nullptr,&k);});
 }
 } catch(const std::exception& e) { ++failures;std::cerr<<"UNEXPECTED_EXCEPTION "<<e.what()<<'\n'; }
 std::cout<<"SUMMARY\tchecks="<<checks<<"\tfailures="<<failures<<'\n';return failures?1:0;
}
