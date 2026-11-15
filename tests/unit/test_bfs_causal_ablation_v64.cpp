@@ -0,0 +1,91 @@
// V64 falsification tests: no solver mocks may be called "physical BFS".
#include "cfdx/physics/bfs_causal_ablation_v64.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::physics;
static void require(bool ok,const char* reason) {
    if(!ok)throw std::runtime_error(reason);
}
static bool close(double x,double y,double tol=1e-12){
    return std::abs(x-y)<tol;
}
int main() {
    constexpr std::size_t n=3;
    PhysicalCandidateV51 oldstate{n,std::vector<double>(6*n,0.0)};
    for(std::size_t i=0;i<n;++i) {
        oldstate.fields[i]=2.0;
        oldstate.fields[3*n+i]=1.0;
        oldstate.fields[4*n+i]=1.0;
        oldstate.fields[5*n+i]=4.0;
    }
    PhysicalCandidateV51 updated=oldstate;
    updated.fields[1]=3.0; // local Ux change
    updated.fields[4*n+1]=3.0; // local k change
    updated.fields[5*n+2]=5.0; // independent omega change
    auto physical=[](const PhysicalCandidateV51& x) {
        PhysicalResidualV51 r;
        r.signed_raw.resize(6*x.cells);
        r.normalized.resize(6*x.cells);
        r.eddy_viscosity.resize(x.cells,0.001);
        r.face_flux={x.fields[1]*x.fields[4*x.cells+1],0.2};
        r.eddy_viscosity[1]=0.001+x.fields[1]*x.fields[4*x.cells+1]*1e-5;
        r.flux_fixed_point_defect=0.0;
        r.flux_closed=true;
        r.merit=0.0;
        for(std::size_t i=0;i<x.cells;++i) {
            const double ux=x.fields[i],k=x.fields[4*x.cells+i];
            // Exact flow-SST interaction in the momentum equation.
            r.signed_raw[i]=2.0*ux*k + ux + 7*k;
            r.signed_raw[x.cells+i]=ux;
            r.signed_raw[2*x.cells+i]=1.;
            r.signed_raw[3*x.cells+i]=ux-2.;
            r.signed_raw[4*x.cells+i]=k-1.;
            r.signed_raw[5*x.cells+i]=x.fields[5*x.cells+i]-4.;
        }
        r.normalized=r.signed_raw;
        return r;
    };
    const auto after=evaluate_physical_candidate_v51(updated,physical,0.,0.);
    const auto a=audit_bfs_causal_ablation_v64(oldstate,updated,after,physical);
    require(a.flow_changed && a.sst_changed,"missing changed blocks");
    require(close(a.blocks[0].interaction_linf,4.0),"wrong 2x2 cross interaction");
    require(a.blocks[0].interaction_hotspot_cell==1,"wrong interaction cell");
    require(close(a.blocks[4].interaction_linf,0.0),"SST linear block invented an interaction");
    require(close(a.blocks[0].flow_effect_linf,3.0),"flow-only effect incorrect");
    require(close(a.blocks[0].sst_effect_linf,22.0),"SST-only effect incorrect");
    require(a.max_recomposition_error<1e-12,"2x2 recomposition failed");
    require(close(a.rc_flux_interaction,2.0),"RC coupling interaction lost");
    require(close(a.eddy_viscosity_interaction,2e-5),"nut interaction lost");
    auto same=audit_bfs_causal_ablation_v64(oldstate,oldstate,
        evaluate_physical_candidate_v51(oldstate,physical),physical);
    require(!same.flow_changed&&!same.sst_changed,"unchanged state falsely called an update");
    require(close(same.blocks[0].interaction_linf,0.),"zero delta invented interaction");
    auto bad=updated;bad.fields[5*n+1]=-1.;
    bool refused=false;
    try{(void)audit_bfs_causal_ablation_v64(oldstate,bad,after,physical);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"negative omega was not rejected");
    refused=false;
    try{auto bogus=after;bogus.flux_closed=false;
        (void)audit_bfs_causal_ablation_v64(oldstate,updated,bogus,physical);}
    catch(const std::runtime_error&){refused=true;}
    require(refused,"unclosed Rhie-Chow flux was not rejected");
    refused=false;
    try{auto badshape=oldstate;badshape.cells=2;
        (void)audit_bfs_causal_ablation_v64(badshape,updated,after,physical);}
    catch(const std::invalid_argument&){refused=true;}
    require(refused,"incompatible shapes not rejected");
    refused=false;
    auto tampered=after;
    tampered.signed_raw[0]+=1.0;
    try{(void)audit_bfs_causal_ablation_v64(oldstate,updated,tampered,physical);}
    catch(const std::runtime_error&){refused=true;}
    require(refused,"stale candidate residual was not rejected");
    std::cout<<"V64_CAUSAL_ABLATION_UNIT_PASS interaction=4 zero_interaction=0 "
                "nonsmooth_newton_claim=0 production_gate_modified=0\n";
    return 0;
}
