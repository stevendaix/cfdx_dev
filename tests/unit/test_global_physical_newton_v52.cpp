#include "cfdx/physics/global_physical_newton_v52.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
using namespace cfdx::physics;
static void require(bool value, const char* name) {
    if(!value)throw std::runtime_error(name);
}
static PhysicalResidualV51 toy(const PhysicalCandidateV51& x) {
    const std::size_t n=x.cells;
    PhysicalResidualV51 result;
    result.signed_raw.assign(6*n,0.0);
    result.normalized.assign(6*n,0.0);
    result.eddy_viscosity.assign(n,0.0);
    result.face_flux={x.fields[0]+x.fields[1],x.fields[3*n]-x.fields[3*n+1]};
    for(std::size_t i=0;i<n;++i) {
        const auto at=[&](std::size_t block){return x.fields[block*n+i];};
        result.eddy_viscosity[i]=at(4)/at(5);
        result.signed_raw[i]=at(0)-1.0+0.03*at(4)+0.01*at(3);
        result.signed_raw[n+i]=at(1)+0.4+0.02*at(3);
        result.signed_raw[2*n+i]=at(2)-0.1;
        result.signed_raw[4*n+i]=at(4)-0.8+0.1*at(0);
        result.signed_raw[5*n+i]=at(5)-1.2+0.05*at(4);
    }
    // A redundant continuity equation, to exercise the pinning policy.
    const double div=(x.fields[3*n]-x.fields[3*n+1])+
                     0.03*(x.fields[0]-x.fields[1]);
    result.signed_raw[3*n]=div;
    result.signed_raw[3*n+1]=-div;
    result.merit=0.0;
    for(std::size_t j=0;j<result.signed_raw.size();++j) {
        result.normalized[j]=result.signed_raw[j];
        result.merit=std::max(result.merit,std::abs(result.normalized[j]));
    }
    result.flux_closed=true;
    result.flux_fixed_point_defect=0.0;
    return result;
}
int main(){
    PhysicalCandidateV51 original;
    original.cells=2;
    original.fields={0.1,0.2, 0.0,0.2, 0.0,0.0, 0.3,-0.2,
                     1.1,1.0, 1.6,1.7};
    GlobalPhysicalNewtonControlsV52 controls;
    controls.maximum_steps=3;
    controls.krylov_dimension=12;
    controls.linear_relative_tolerance=1e-4;
    controls.maximum_backtracks=9;
    const auto initial=original.fields;
    const auto report=global_physical_newton_v52(original,toy,controls);
    require(!report.steps.empty(),"no Newton step recorded");
    require(report.accepted_steps>=1,"private Newton failed to contract a coupled problem");
    require(report.last_private_residual.merit<toy(original).merit,
            "physical merit did not decrease");
    require(!report.production_committed && !report.physically_qualified,
            "private Newton was falsely promoted");
    require(initial==original.fields,"private solver modified accepted state");
    // The experimental SST 2x2 block is an explicit opt-in. Probes must
    // originate from the very same residual callback, never stale SST data.
    auto sst_block=controls;
    sst_block.sst_two_by_two_preconditioner=true;
    const auto sst_report=global_physical_newton_v52(original,toy,sst_block);
    require(!sst_report.steps.empty() &&
            sst_report.steps.front().sst_two_by_two_ready,
            "SST 2x2 diagnostic probe was not available on the toy problem");
    require(!sst_report.production_committed && !sst_report.physically_qualified,
            "SST preconditioner silently promoted a private solve");
    for(const auto& step:report.steps) {
        if(step.accepted) {
            require(step.alpha>0.0 && step.alpha<=1.0,"bad accepted alpha");
            require(step.merit_after<step.merit_before,"uncontracted step");
            require(step.linear_true_relative<=controls.linear_relative_tolerance,
                    "Krylov estimated rather than true residual certified");
        }
    }
    // V53 adds a distinct frozen Pmat builder to the SAME matrix-free Jv;
    // neither the builder nor its action owns or commits production state.
    auto newer=controls;
    newer.v53_extended_krylov=true;
    bool builder_used=false;
    auto builder=[&](const PhysicalCandidateV51& state,
                     const PhysicalResidualV51&,
                     const std::array<double,6>&,
                     const std::array<double,6>&) {
        if(!builder_used)
            require(state.fields==original.fields,"V53 builder received wrong initial state");
        builder_used=true;
        GlobalPhysicalFrozenPmatV53 frozen;
        frozen.apply=[](const std::vector<double>& v){return v;};
        frozen.probes=13;frozen.transport_nnz=28;
        return frozen;
    };
    const auto new_report=global_physical_newton_v52(
        original,toy,newer,{}, {},builder);
    require(builder_used && !new_report.steps.empty() &&
            new_report.steps[0].v53_global_block_ready &&
            new_report.steps[0].v53_probe_evaluations==13 &&
            new_report.steps[0].v53_transport_nnz==28 &&
            !new_report.production_committed,
            "V53 Pmat builder was not wired transactionally");
    bool rejected_new=false;
    try{(void)global_physical_newton_v52(original,toy,newer);}
    catch(const std::invalid_argument&){rejected_new=true;}
    require(rejected_new,"V53 extended Krylov accepted without Pmat");
    // With fixed pressure BC the auxiliary gauge is suppressed.
    auto with_boundary=controls;with_boundary.pin_pressure=false;
    const auto result2=global_physical_newton_v52(original,toy,with_boundary);
    require(!result2.physically_qualified,"fixed-pressure branch falsely qualified");
    // An inadmissible initial omega must fail before any evaluation.
    auto invalid=original;invalid.fields[10]=0.0;
    bool threw=false;
    try{(void)global_physical_newton_v52(invalid,toy,controls);}
    catch(const std::invalid_argument&){threw=true;}
    require(threw,"nonpositive omega admitted");
    // A broken preconditioner is a failure, never an identity fallback.
    threw=false;
    try{(void)global_physical_newton_v52(original,toy,controls,
              [](const std::vector<double>&){return std::vector<double>(1,0.0);});}
    catch(const std::runtime_error&){threw=true;}
    require(threw,"malformed preconditioner silently accepted");
    // Identical residual at every candidate must not produce a false step.
    const auto flat=[](const PhysicalCandidateV51& x) {
        auto residual=toy(x);
        for(double& r:residual.signed_raw)r=1.0;
        for(double& r:residual.normalized)r=1.0;
        residual.merit=1.0;
        return residual;
    };
    const auto stalled=global_physical_newton_v52(original,flat,controls);
    require(stalled.accepted_steps==0,"flat physical operator accepted as Newton");
}
