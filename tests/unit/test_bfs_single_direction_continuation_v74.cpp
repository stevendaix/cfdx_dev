@@ -0,0 +1,81 @@
#include "cfdx/physics/bfs_single_direction_continuation_v74.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void need(bool good,const char* msg){if(!good)throw std::runtime_error(msg);}
int main(){
    const std::array<std::size_t,3> own{0,0,1};
    const std::array<int,3> nei{-1,1,-1};
    const auto owner=[&](std::size_t f){return own.at(f);};
    const auto neighbour=[&](std::size_t f){return nei.at(f);};
    PhysicalCandidateV51 x;x.cells=2;x.fields.assign(12,0.);
    x.fields[8]=x.fields[9]=.2;x.fields[10]=x.fields[11]=2.;
    const auto mass_flux=[](const PhysicalCandidateV51& s){
        return std::vector<double>{-1.,1.8+s.fields[6]-s.fields[7],.3+s.fields[7]};};
    std::size_t calls=0;
    const auto evaluate=[&](const PhysicalCandidateV51& s){
        ++calls;
        PhysicalResidualV51 r;r.flux_closed=true;r.flux_fixed_point_defect=0.;
        r.face_flux=mass_flux(s);r.signed_raw.assign(12,0.);
        r.normalized.assign(12,0.);r.eddy_viscosity.assign(2,.1);
        r.inverse_momentum_diagonal.assign(6,.25);
        r.signed_raw[6]=r.face_flux[0]+r.face_flux[1];
        r.signed_raw[7]=-r.face_flux[1]+r.face_flux[2];
        r.normalized=r.signed_raw;r.merit=std::abs(r.signed_raw[6]);
        return r;
    };
    const auto transport=[](const PhysicalCandidateV51& s,const std::vector<double>& phi){
        need(phi.size()==3,"incomplete SST face flux");return s;};
    const auto project=[&](const PhysicalCandidateV51& s,
                           const PhysicalResidualV51& r,double alpha){
        const auto response=[&](const std::vector<double>& dp){
            BfsPressureResponseV70 q;q.state=s;
            q.state.fields[6]+=dp[0];q.state.fields[7]+=dp[1];
            q.face_flux=mass_flux(q.state);return q;};
        return bfs_project_pressure_v70(s,r,response,owner,neighbour,true,0,
                                        1e-10,alpha);
    };
    const auto baseline=evaluate(x);
    const auto saved=x.fields;const auto saved_phi=baseline.face_flux;
    PhysicalCandidateV51 last;PhysicalResidualV51 last_f;
    const auto report=bfs_continue_single_direction_v74(x,baseline,x,baseline,
        evaluate,project,transport,owner,neighbour,1,0,1,8,
        0.,0.,1e-8,1e-8,1e-10,&last,&last_f);
    need(report.rows.size()==8 && report.accepted_private_cycles==8,
         "V74 should repeat fixed direction without searching");
    need(report.final_mass<report.initial_mass &&
         report.final_mass>1e-10 && !report.production_commit_authorized &&
         !report.native_momentum_certificate && !report.native_sst_certificate,
         "V74 miscertified private pressure continuation");
    need(x.fields==saved && baseline.face_flux==saved_phi,
         "V74 mutated authoritative fields or flux");
    for(const auto& row:report.rows)
        need(row.evaluated && row.selected && row.guards_ok &&
             row.strict_mass_descent && row.merit_descends &&
             row.relevant_groups_nonincreasing && row.mass_after<row.mass_before,
             "V74 did not monitor every successive reassembly");
    need(calls>20,"V74 did not independently re-evaluate Fphys");
    bool bad_limit=false;
    try{(void)bfs_continue_single_direction_v74(x,baseline,x,baseline,
        evaluate,project,transport,owner,neighbour,1,0,1,9,0.,0.,1e-8,1e-8);}
    catch(const std::invalid_argument&){bad_limit=true;}
    need(bad_limit,"V74 exceeded bounded private continuation budget");
    bool bad_face=false;
    try{(void)bfs_continue_single_direction_v74(x,baseline,x,baseline,
        evaluate,project,transport,owner,neighbour,1,1,0,2,0.,0.,1e-8,1e-8);}
    catch(const std::invalid_argument&){bad_face=true;}
    need(bad_face,"V74 accepted reversed face-3 connectivity");
    const auto broken_transport=[](const PhysicalCandidateV51&,const std::vector<double>&)
        ->PhysicalCandidateV51{throw std::runtime_error("forced SST failure");};
    const auto failed=bfs_continue_single_direction_v74(x,baseline,x,baseline,
        evaluate,project,broken_transport,owner,neighbour,1,0,1,3,0.,0.,1e-8,1e-8);
    need(failed.rows.size()==1 && failed.accepted_private_cycles==0 &&
         failed.final_mass==failed.initial_mass &&
         !failed.production_commit_authorized,
         "V74 must stop without commit on SST failure");
    std::cout<<"V74_FIXED_DIRECTION_PASS cycles="<<report.accepted_private_cycles
             <<" mass="<<report.final_mass<<" committed=0\n";
}
