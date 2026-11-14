#include "cfdx/physics/global_nonlinear_contraction_v14.h"
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
static void require(bool yes,const char* why){if(!yes) throw std::runtime_error(why);}
static GlobalNonlinearResidualV8 value(double u, double k, double cont) {
 GlobalNonlinearResidualV8 v;
 v.algorithm=PressureVelocityAlgorithm::COUPLED;
 v.schur=CoupledSchurModel::PCD;
 v.momentum.max_active_component_relative=u;
 v.momentum.global_relative=u;
 v.k=k;v.omega=0.0001;v.continuity=cont;
 v.eddy_viscosity_closure=1e-8;
 v.merit=std::max(u,std::max(k,cont));
 v.sst=true;v.sst_transport_committed=true;v.bounded=true;v.all_finite=true;
 return v;
}
int main(){
 try{
  GlobalNonlinearContractionLedgerV14 ledger;
  auto first=value(0.8,0.1,0.01);
  require(ledger.evaluate(first).baseline,"no baseline");
  ledger.commit(first);
  auto second=value(0.6,0.2,0.006);
  auto check=ledger.evaluate(second);
  require(check.allowed&&!check.baseline,"first contraction failed");
  ledger.commit(second);
  auto plateau=value(0.6,0.1,0.003);
  require(!ledger.evaluate(plateau).allowed,"plateau incorrectly accepted");
  require(ledger.committed_steps()==2,"rejected trial modified ledger");
  auto third=value(0.4,0.06,0.003);
  ledger.commit(third);
  require(ledger.committed_steps()==3&&ledger.contraction_steps()==2,
     "multi-step contraction certificate incorrect");
  auto switched=value(0.5,0.03,0.002);
  require(ledger.evaluate(switched,false).baseline,"operator switch not reset");
  ledger.commit(switched,false);
  require(ledger.contraction_steps()==0,"mixed operator convergence mislabeled");
  auto bad=value(0.3,0.02,0.001);bad.sst_transport_committed=false;
  bool rejected=false;try{ledger.commit(bad);}catch(const std::invalid_argument&){rejected=true;}
  require(rejected,"uncommitted state accepted");
  std::cout << "GLOBAL_CONTRACTION_V14_PASS committed=4 contractions=0"<<"\n";
  return 0;
 }catch(const std::exception& e){std::cerr<<"GLOBAL_CONTRACTION_V14_FAIL: "<<e.what()<<"\n";return 1;}
}
