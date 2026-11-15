#include "cfdx/physics/global_outer_progress_v84.h"
#include "cfdx/physics/global_nonlinear_contraction_v14.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::physics;
static void check(bool ok, const char* why) {
    if (!ok) throw std::runtime_error(why);
}
template<class F> static void throws(F f, const char* why) {
    bool yes=false;
    try { f(); } catch (const std::invalid_argument&) { yes=true; }
    check(yes,why);
}
static GlobalNonlinearResidualV8 state(double worst, double global,
                                       double k, double omega,
                                       double continuity=2e-16) {
    GlobalNonlinearResidualV8 r;
    r.algorithm=PressureVelocityAlgorithm::SIMPLE;
    r.sst=true;
    r.sst_transport_committed=false;
    r.bounded=true;
    r.all_finite=true;
    r.converged=false;
    r.momentum.max_active_component_relative=worst;
    r.momentum.global_relative=global;
    r.k=k;r.omega=omega;r.continuity=continuity;
    r.eddy_viscosity_closure=0.;
    r.merit=std::max({worst,global,k,omega,continuity});
    return r;
}
int main() {
  try {
    ConvergenceCriteria gates;
    gates.relative_tolerance=1e-8;
    gates.continuity_tolerance=1e-11;
    const auto first=state(.9951482239106,.58021183,.281450807985,.0037696064);
    const auto second=state(.9950336766582,.44639442,.86215,.0083,1.4e-16);
    GlobalOuterProgressLedgerV84 ledger;
    GlobalNonlinearContractionLedgerV14 original;
    original.commit(first,true,true);
    check(!original.evaluate(second,true,true).allowed,
          "V83 rejects ratio above 0.995, motivating opt-in V84");
    check(ledger.evaluate(first,gates).allowed,"baseline");
    ledger.commit(first,gates);
    const auto g=ledger.evaluate(second,gates);
    check(g.allowed && !g.baseline,"must permit measured step-2 global carrier progress");
    check(g.k_debt==1,"must account for temporary k explosion");
    check(ledger.committed_steps()==1,"trial must not mutate history");
    ledger.commit(second,gates);
    check(ledger.committed_steps()==2 && ledger.k_debt()==1,"commit only accepted");
    const auto third=state(.9949,.35,.89,.0080,1.2e-16);
    const auto h=ledger.evaluate(third,gates);
    check(h.allowed && h.k_debt==2,"second bounded k debt");
    ledger.commit(third,gates);
    const auto fourth=state(.9947,.30,.88,.0079,1.2e-16);
    const auto debt=ledger.evaluate(fourth,gates);
    check(!debt.allowed && debt.reason=="unrecovered_sst_growth_debt",
          "cannot indefinitely hide worsening SST");
    const auto kblow=state(.9947,.30,8.,.008,1.2e-16);
    check(!ledger.evaluate(kblow,gates).allowed,"k unbounded growth blocked");
    const auto notImproving=state(.9951481,.5800,.862,.008,2e-16);
    GlobalOuterProgressLedgerV84 newLedger;
    newLedger.commit(first,gates);
    check(!newLedger.evaluate(notImproving,gates).allowed,
          "SST x3 without carrier progress blocked");
    auto bad=second;bad.all_finite=false;
    throws([&]{ (void)ledger.evaluate(bad,gates); },"nonfinite accepted");
    bad=second;bad.bounded=false;
    throws([&]{ (void)ledger.evaluate(bad,gates); },"nonbounded accepted");
    bad=second;bad.k=std::numeric_limits<double>::quiet_NaN();
    throws([&]{ (void)ledger.evaluate(bad,gates); },"nan SST accepted");
    bad=second;bad.continuity=1e-5;
    check(!ledger.evaluate(bad,gates).allowed,"pressure physical conservation lost");
    bad=second;bad.momentum.max_active_component_relative=2.;
    check(!ledger.evaluate(bad,gates).allowed,"momentum exploded");
    bad=second;bad.converged=true;
    throws([&]{ (void)ledger.evaluate(bad,gates); },"fake final certificate accepted");
    bad=second;bad.converged=true;bad.sst_transport_committed=true;
    check(!ledger.evaluate(bad,gates).allowed,
          "asserted final status cannot bypass physical convergence criteria");
    bad=second;bad.algorithm=PressureVelocityAlgorithm::SIMPLEC;
    throws([&]{ (void)ledger.evaluate(bad,gates); },"SIMPLEC silently enabled");
    throws([&]{ (void)ledger.evaluate(second,gates,false); },"different physical operator");
    bad=second;bad.momentum.global_relative=0.;bad.momentum.max_active_component_relative=0.;
    bad.continuity=0.;bad.k=0.;bad.omega=0.;bad.merit=0.;
    bad.converged=true;bad.sst_transport_committed=true;
    const auto terminal=ledger.evaluate(bad,gates);
    check(terminal.allowed && terminal.reason=="independent_final_certificate",
          "real final certificate vetoed");
    check(!first.converged && !second.converged,"provisional claim changed");
    std::cout<<"PASS V84 nonmonotone provisional debt and non-certifying guards\n";
    return 0;
  } catch(const std::exception& e) {
    std::cerr<<"FAIL V84 "<<e.what()<<"\n"; return 1;
  }
}
