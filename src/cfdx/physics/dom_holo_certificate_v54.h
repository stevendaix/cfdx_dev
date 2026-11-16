#pragma once
// V54 research-only final independent certificate for an already converged
// V53 coupled HO-SMM/AMG candidate. Never use the HO-LO identity as the
// independent angular or thermophysical validation oracle.
#include "cfdx/physics/dom_independent_angular_audit_v54.h"
#include "cfdx/physics/dom_holo_amg_global_v53.h"
#include <string>
#include <utility>
#include <vector>
namespace cfdx::physics {
struct DomV54Certificate {
    bool published=false;
    DomV53Result pilot;
    DomV54AngularAudit angular;
    std::string diagnostic;
    std::vector<double> accepted_state;
};
inline DomV54Certificate certify_dom_holo_v54(
    DomMomentSchurSystemV44& ho,DomV53Result candidate,
    double independent_tolerance,std::size_t simulated_partitions=1){
    DomV54Certificate result;
    result.pilot=std::move(candidate);
    if(!result.pilot.converged || !result.pilot.coupled_amg ||
       result.pilot.state.size()!=ho.dimension()){
        result.pilot.state.clear();
        result.diagnostic="V54: no converged HO-SMM/AMG state to certify";
        return result;
    }
    // Independent recomputation from every ordinate and every signed FV face.
    result.angular=audit_dom_angular_v54(ho,result.pilot.state,
                                        independent_tolerance,simulated_partitions);
    if(!result.angular.accepted){
        result.pilot.state.clear();
        result.pilot.converged=false;
        result.pilot.coupled_amg=false;
        result.diagnostic="V54: independent angular/gray-wall/energy gate rejected";
        return result;
    }
    result.accepted_state=result.pilot.state;
    result.published=true;
    result.diagnostic="V54: independent angular certificate accepted";
    return result;
}
inline DomV54Certificate solve_dom_holo_certified_v54(
    DomMomentSchurSystemV44& ho,const std::vector<double>& initial,
    const DomV53Options& options={},std::size_t simulated_partitions=1) {
    return certify_dom_holo_v54(ho,solve_dom_holo_amg_v53(ho,initial,options),
        options.base.base.physical_tolerance,simulated_partitions);
}
}
