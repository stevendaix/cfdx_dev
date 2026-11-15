@@ -0,0 +1,56 @@
#include "cfdx/physics/timestep_control.h"
#include "cfdx/physics/global_outer_progress_v85.h"
#include <cmath>
int main() {
    using cfdx::physics::NonlinearRetryController;
    using cfdx::physics::NonlinearRetryControls;
    NonlinearRetryControls c;
    c.max_retries=5;
    c.relaxation_shrink=0.5;
    c.minimum_alpha_u=0.015625;
    c.minimum_alpha_p=0.00625;
    NonlinearRetryController retry(c);
    const double alpha_u=0.25, alpha_p=0.1;
    retry.reject(false,true); // pressure-only failed preflight
    if(std::abs(retry.alpha_u(alpha_u)-0.25)>1e-14 ||
       std::abs(retry.alpha_p(alpha_p)-0.05)>1e-14)return 1;
    retry.reject(true,false); // transverse-only failed preflight
    if(std::abs(retry.alpha_u(alpha_u)-0.125)>1e-14 ||
       std::abs(retry.alpha_p(alpha_p)-0.05)>1e-14)return 2;
    retry.reject(); // shared global residual rejection
    if(std::abs(retry.alpha_u(alpha_u)-0.0625)>1e-14 ||
       std::abs(retry.alpha_p(alpha_p)-0.025)>1e-14)return 3;
    if(retry.retries()!=3)return 4;
    retry.reset();
    if(retry.retries()!=0 || retry.alpha_u(alpha_u)!=alpha_u ||
       retry.alpha_p(alpha_p)!=alpha_p)return 5;
    // Measured V89 state 3: momentum and k improve, but the worst
    // active component increases 1.85%. V85 default 1% rejects it.
    // The explicitly guarded V89 3% progression window can accept it,
    // WITHOUT changing final tolerances or prior V85 behavior.
    auto state=[](double worst,double global,double k,double omega){
        cfdx::physics::GlobalNonlinearResidualV8 s;
        s.algorithm=cfdx::physics::PressureVelocityAlgorithm::SIMPLE;
        s.sst=true;s.bounded=true;s.all_finite=true;
        s.momentum.global_relative=global;
        s.momentum.max_active_component_relative=worst;
        s.k=k;s.omega=omega;s.continuity=1e-16;
        s.eddy_viscosity_closure=0;
        s.merit=std::max({worst,global,k,omega});
        return s;
    };
    cfdx::physics::ConvergenceCriteria gates;
    gates.relative_tolerance=1e-8;gates.continuity_tolerance=1e-11;
    cfdx::physics::GlobalOuterProgressLedgerV85 ledger;
    const auto first=state(.995148,.580463,.281451,.003770);
    const auto second=state(.9808236,.446335,.870368,.010274);
    const auto third=state(.998931,.433873,.534251,.023165);
    ledger.commit(first,gates);
    if(!ledger.evaluate(second,gates).allowed)return 6;
    ledger.commit(second,gates);
    if(ledger.evaluate(third,gates).allowed)return 7;
    if(!ledger.evaluate(third,gates,true,1.03).allowed)return 8;
    ledger.commit(third,gates,true,1.03);
    if(ledger.committed_steps()!=3)return 9;
    return 0;
}
