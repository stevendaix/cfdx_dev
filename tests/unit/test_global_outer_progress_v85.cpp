#include "cfdx/physics/global_outer_progress_v85.h"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using namespace cfdx::physics;
static void require(bool cond,const std::string& msg){if(!cond)throw std::runtime_error(msg);}
template<class F>static void must_throw(F fn,const char* msg){bool ok=false;try{fn();}catch(const std::invalid_argument&){ok=true;}require(ok,msg);}
static GlobalNonlinearResidualV8 state(double worst,double glob,double k,double omega,
                                      double cont=1e-16,
                                      PressureVelocityAlgorithm algorithm=PressureVelocityAlgorithm::SIMPLE,
                                      CoupledSchurModel schur=CoupledSchurModel::BlockLocal,
                                      bool sst_committed=false){
    GlobalNonlinearResidualV8 s;
    s.algorithm=algorithm;s.schur=schur;
    s.sst=true;s.sst_transport_committed=sst_committed;s.bounded=true;
    s.all_finite=true;s.converged=false;
    s.momentum.global_relative=glob;
    s.momentum.max_active_component_relative=worst;
    s.k=k;s.omega=omega;s.continuity=cont;
    s.eddy_viscosity_closure=0.;
    s.merit=std::max({worst,glob,k,omega,cont});
    return s;
}
int main(){try{
    ConvergenceCriteria gates;
    gates.relative_tolerance=1e-8;
    gates.continuity_tolerance=1e-11;
    auto one=state(.9951482239106,.580463,.281451,.003770);
    auto two=state(.9950336766582,.446394,.862147,.010256);
    auto three=state(.982421,.347084,.813471,.023871);
    GlobalOuterProgressLedgerV85 ledger;
    require(ledger.evaluate(one,gates).allowed,"initial conservative provisional baseline");
    ledger.commit(one,gates);
    require(ledger.evaluate(two,gates).allowed,"measured V84 second state must remain permitted");
    ledger.commit(two,gates);
    const auto t3=ledger.evaluate(three,gates);
    require(t3.allowed && t3.omega_debt==2 && t3.reason=="bounded_provisional_sst_excursion",
            "measured V84 third candidate must be diagnosed and optionally accepted");
    require(ledger.committed_steps()==2,"evaluation must be pure");
    ledger.commit(three,gates);
    require(ledger.committed_steps()==3,"commit only accepted state");
    // Recovery is required: persistent omega excursions cannot continue forever.
    auto four=state(.96,.28,.70,.020);
    auto five=state(.94,.24,.65,.019);
    auto six=state(.92,.20,.61,.018);
    require(ledger.evaluate(four,gates).allowed,"bounded fourth state");ledger.commit(four,gates);
    require(ledger.evaluate(five,gates).allowed,"bounded fifth state");ledger.commit(five,gates);
    require(ledger.evaluate(six,gates).allowed,"bounded sixth state");ledger.commit(six,gates);
    auto seven=state(.90,.17,.57,.017);
    const auto exhausted=ledger.evaluate(seven,gates);
    require(!exhausted.allowed && exhausted.reason=="sst_exhausted_recovery_window",
            "no unlimited accepted omega debt");
    auto bad=three;bad.omega=std::numeric_limits<double>::quiet_NaN();
    must_throw([&]{(void)ledger.evaluate(bad,gates);},"NaN rejected");
    bad=three;bad.bounded=false;
    must_throw([&]{(void)ledger.evaluate(bad,gates);},"unbounded rejected");
    bad=three;bad.continuity=1e-3;
    require(!ledger.evaluate(bad,gates).allowed,"continuity guard retained");
    bad=three;bad.momentum.global_relative=2.;
    require(!ledger.evaluate(bad,gates).allowed,"momentum blowup rejected");
    bad=three;bad.converged=true;
    must_throw([&]{(void)ledger.evaluate(bad,gates);},"false converged flag rejected");
    // V85 is available to every pressure-velocity algorithm when SST is
    // committed. Only the separate V83 provisional-transport path is SIMPLE-only.
    require(global_outer_progress_v85_algorithm_compatible(
                PressureVelocityAlgorithm::SIMPLE,true),
            "SIMPLE may use the V83 provisional SST transaction");
    for (const auto algorithm : {PressureVelocityAlgorithm::SIMPLEC,
                                 PressureVelocityAlgorithm::COUPLED})
        require(global_outer_progress_v85_algorithm_compatible(algorithm,false),
                "SIMPLEC/COUPLED may use V85 with committed SST");
    require(!global_outer_progress_v85_algorithm_compatible(
                PressureVelocityAlgorithm::COUPLED,true) &&
            !global_outer_progress_v85_algorithm_compatible(
                PressureVelocityAlgorithm::SIMPLEC,true),
            "V83 provisional SST remains restricted to SIMPLE");

    // Each COUPLED Schur configuration owns an independent ledger. No history
    // may be mixed across algorithm or Schur operator changes.
    const std::array<CoupledSchurModel,6> schur_models{{
        CoupledSchurModel::BlockLocal,CoupledSchurModel::PCD,
        CoupledSchurModel::LSC,CoupledSchurModel::BFBT,
        CoupledSchurModel::SIMPLE,CoupledSchurModel::SIMPLEC}};
    for (const auto schur : schur_models) {
        GlobalOuterProgressLedgerV85 coupledLedger;
        auto coupledOne=state(.9951482239106,.580463,.281451,.003770,1e-16,
                              PressureVelocityAlgorithm::COUPLED,schur,true);
        auto coupledTwo=state(.9950336766582,.446394,.862147,.010256,1e-16,
                              PressureVelocityAlgorithm::COUPLED,schur,true);
        require(coupledLedger.evaluate(coupledOne,gates).allowed,
                "COUPLED Schur gets its own committed-SST baseline");
        coupledLedger.commit(coupledOne,gates);
        require(coupledLedger.evaluate(coupledTwo,gates).allowed,
                "COUPLED Schur progression uses the same physical guards");
    }
    GlobalOuterProgressLedgerV85 coupledLedger;
    auto coupledOne=state(.9951482239106,.580463,.281451,.003770,1e-16,
                          PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::PCD,true);
    coupledLedger.commit(coupledOne,gates);
    auto mixed=state(.9950336766582,.446394,.862147,.010256,1e-16,
                     PressureVelocityAlgorithm::COUPLED,CoupledSchurModel::LSC,true);
    must_throw([&]{(void)coupledLedger.evaluate(mixed,gates);},
               "Schur change inside one ledger rejected");
    mixed=three;
    must_throw([&]{(void)coupledLedger.evaluate(mixed,gates);},
               "algorithm change inside one ledger rejected");

    GlobalOuterProgressLedgerV85 simplecLedger;
    auto simplecOne=state(.9951482239106,.580463,.281451,.003770,1e-16,
                          PressureVelocityAlgorithm::SIMPLEC,CoupledSchurModel::BlockLocal,true);
    require(simplecLedger.evaluate(simplecOne,gates).allowed,
            "SIMPLEC gets an independent committed-SST baseline");
    simplecLedger.commit(simplecOne,gates);
    must_throw([&]{(void)ledger.evaluate(three,gates,false);},"changed physical operator rejected");
    auto noProgress=state(.99512,.5804,.282,.00377);
    GlobalOuterProgressLedgerV85 fresh;fresh.commit(one,gates);
    require(!fresh.evaluate(noProgress,gates).allowed,"no physical progress rejected");
    auto direct=state(.5,.45,.15,.002);
    GlobalOuterProgressLedgerV85 safe;safe.commit(one,gates);
    require(safe.evaluate(direct,gates).allowed,"contracting physical state accepted");
    require(ledger.committed_steps()==6,"rejection never changes ledger");
    std::cout<<"V85 PROGRESSION_LEDGER_OK measured V84 third trial accepted privately; "
                "pressure/positivity/finite/COUPLED/omega-debt guards retained\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}