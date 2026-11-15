#include "cfdx/physics/global_nonlinear_contraction_v14.h"
#include "cfdx/physics/global_nonlinear_acceptance_v13.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::physics;
static void require(bool yes, const char* text) {
    if (!yes) throw std::runtime_error(text);
}
template<class F> static void must_throw(F fn,const char* message) {
    bool threw=false;
    try { fn(); } catch (const std::invalid_argument&) { threw=true; }
    require(threw,message);
}
static GlobalNonlinearResidualV8 state(double merit,
                                       bool converged=false,
                                       bool sst_completed=false) {
    GlobalNonlinearResidualV8 v;
    v.sst=true;
    v.sst_transport_committed=sst_completed;
    v.bounded=true;
    v.all_finite=true;
    v.converged=converged;
    v.continuity=1e-13;
    v.k=0.4*merit;
    v.omega=0.7*merit;
    v.eddy_viscosity_closure=0.;
    v.momentum.global_relative=0.85*merit;
    v.momentum.max_active_component_relative=merit;
    v.merit=merit;
    return v;
}
int main() {
 try {
  auto a=state(1.0), b=state(0.95), worse=state(1.01);
  // A provisional state is not a convergence certificate, even in opt-in mode.
  require(!a.converged && !a.sst_transport_committed,"fake qualification");
  GlobalNonlinearContractionLedgerV14 ledger;
  must_throw([&]{ (void)ledger.evaluate(a); },"V14 default incorrectly permits provisional");
  must_throw([&]{ (void)ledger.commit(a); },"V14 default incorrectly commits provisional");
  require(ledger.evaluate(a,true,true).baseline,"V14 provisional baseline rejected");
  ledger.commit(a,true,true);
  must_throw([&]{ (void)ledger.evaluate(b,true,false); },"V14 fallback silently qualified prior provisional");
  const auto decline=ledger.evaluate(b,true,true);
  require(decline.allowed && !decline.baseline && std::abs(decline.ratio-0.95)<1e-12,
      "V14 rejected a genuinely descending physical residual");
  require(!ledger.evaluate(worse,true,true).allowed,"V14 admitted an increasing residual");
  ledger.commit(b,true,true);
  require(ledger.committed_steps()==2 && ledger.contraction_steps()==1,
      "V14 ledger lost provenance");
  auto fake=state(0.3,true,false);
  must_throw([&]{ (void)ledger.evaluate(fake,true,true); },"V14 accepted fictitious final certificate");
  a.all_finite=false;
  must_throw([&]{ (void)ledger.evaluate(a,true,true); },"V14 admitted nonfinite state");
  a=state(1.0); a.bounded=false;
  must_throw([&]{ (void)ledger.evaluate(a,true,true); },"V14 admitted unclipped state");
  a=state(1.0);
  ConvergenceCriteria physical;
  must_throw([&]{ (void)evaluate_global_nonlinear_trial_v13(a,b,physical); },
      "V13 default admitted provisional");
  auto good=evaluate_global_nonlinear_trial_v13(a,b,physical,1.20,10.0,true);
  require(good.allowed,"V13 physical provisional contraction rejected");
  worse=state(2.0);
  auto bad=evaluate_global_nonlinear_trial_v13(b,worse,physical,1.20,10.0,true);
  require(!bad.allowed,"V13 admitted a large physical residual increase");
  auto certified=state(1e-10,true,true);
  must_throw([&]{ (void)ledger.evaluate(certified,true,false); },
      "V14 incorrectly promoted a provisional baseline");
  require(ledger.evaluate(certified,true,true).allowed,
      "V14 did not accept converged endpoint while tracking provisional progress");
  std::cout<<"PASS V83 provisional guard and FINAL convergence separation\n";
  return 0;
 } catch (const std::exception& e) {
  std::cerr<<"FAIL V83 "<<e.what()<<"\n";return 1;
 }
}
