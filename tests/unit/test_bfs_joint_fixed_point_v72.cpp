// V72: true repeated private Picard / Rhie-Chow / SST sequence,
// six physical certificates, stale operator injection, no production commit.
#include "cfdx/physics/bfs_joint_fixed_point_v72.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
int main(){
    constexpr std::array<std::size_t,3> own={0,0,1};
    constexpr std::array<int,3> neigh={-1,1,-1};
    const auto owner=[&](std::size_t f){return own.at(f);};
    const auto neighbour=[&](std::size_t f){return neigh.at(f);};
    PhysicalCandidateV51 start;
    start.cells=2;start.fields.assign(12,0.);
    start.fields[0]=1.;start.fields[1]=1.;
    start.fields[8]=.3;start.fields[9]=.3;
    start.fields[10]=2.;start.fields[11]=2.;
    const auto flux=[](const PhysicalCandidateV51& s){
        return std::vector<double>{-1.,.2 + (s.fields[6]-s.fields[7]) +
                                   .2*(s.fields[0]-s.fields[1]),.3+s.fields[7]};
    };
    std::size_t evaluations=0;
    const auto evaluate=[&](const PhysicalCandidateV51& s){
        ++evaluations;
        PhysicalResidualV51 r;
        r.flux_closed=true;r.flux_fixed_point_defect=0.;
        r.face_flux=flux(s);
        r.signed_raw.assign(12,0.);r.normalized.assign(12,0.);
        r.eddy_viscosity.resize(2);r.inverse_momentum_diagonal.assign(6,.2);
        for(std::size_t i=0;i<2;++i){
            r.signed_raw[i]=s.fields[i]-(1.5-.1*s.fields[8]);
            r.signed_raw[8+i]=s.fields[8+i]-(.2+.02*s.fields[i]);
            r.signed_raw[10+i]=s.fields[10+i]-(2.+.1*s.fields[i]);
            r.eddy_viscosity[i]=s.fields[8+i]/s.fields[10+i];
        }
        r.signed_raw[6]=r.face_flux[0]+r.face_flux[1];
        r.signed_raw[7]=-r.face_flux[1]+r.face_flux[2];
        for(std::size_t i=0;i<12;++i){
            r.normalized[i]=r.signed_raw[i];
            r.merit=std::max(r.merit==std::numeric_limits<double>::infinity()?0.:r.merit,
                             std::abs(r.signed_raw[i]));
        }
        return r;
    };
    const auto start_f=evaluate(start);
    std::size_t transport_calls=0;
    const auto transport=[&](const PhysicalCandidateV51& s,
                             const std::vector<double>& candidate_flux){
        ++transport_calls;
        require(candidate_flux.size()==3,"wrong private SST flux size");
        auto x=s;
        for(std::size_t i=0;i<2;++i){
            x.fields[8+i]=.2+.02*s.fields[i];
            x.fields[10+i]=2.+.1*s.fields[i];
        }
        return x;
    };
    const auto project=[&](const PhysicalCandidateV51& s,
                           const PhysicalResidualV51& r,double fraction){
        const auto response=[&](const std::vector<double>& dp){
            BfsPressureResponseV70 out;out.state=s;
            out.state.fields[6]+=dp[0];out.state.fields[7]+=dp[1];
            out.face_flux=flux(out.state);
            return out;
        };
        return bfs_project_pressure_v70(s,r,response,owner,neighbour,true,0,1e-10,fraction);
    };
    const auto immutable_start=start;
    const auto original_f=start_f;
    PhysicalCandidateV51 last;
    PhysicalResidualV51 last_f;
    const std::vector<BfsJointChoiceV72> choices={{.25,1.},{0.,.5}};
    const auto result=bfs_joint_private_fixed_point_v72(
        start,start_f,start,start_f,evaluate,project,transport,owner,neighbour,
        choices,3,0.,0.,1e-8,1e-8,1e-10,&last,&last_f);
    require(start.fields==immutable_start.fields &&
            start_f.face_flux==original_f.face_flux,"private V72 leaked into accepted state");
    require(!result.attempts.empty() && evaluations>6 && transport_calls>0,
            "V72 did not execute physical operators");
    require(!result.commit_authorized,"V72 private continuation committed a state");
    for(const auto& c:result.certificates){
        require(c.replay_identical,"no independent residual replay");
        require(c.signed_continuity_vs_face_gap<1e-12,"wrong face continuity");
        require(!c.native_momentum_linear_recertified &&
                !c.native_sst_transport_recertified &&
                !c.eligible_for_production_commit,"fake native linear certificate");
    }
    const auto final_cert=bfs_recertify_six_field_v72(
        start,start_f,last,last_f,evaluate,owner,neighbour,1e-8,1e-8,1e-10);
    require(!final_cert.eligible_for_production_commit,"replay alone authorized commit");
    // Falsify continuity while retaining an otherwise plausible operator.
    const auto fake_continuity=[&](const PhysicalCandidateV51& s){
        auto r=evaluate(s);r.signed_raw[6]+=.1;r.normalized[6]+=.1;return r;};
    bool continuity_rejected=false;
    try {
        const auto invalid=fake_continuity(last);
        (void)bfs_recertify_six_field_v72(start,start_f,last,invalid,
            fake_continuity,owner,neighbour,1e-8,1e-8,1e-10);
    } catch(const std::exception&){continuity_rejected=true;}
    require(continuity_rejected,"fake continuity was accepted");
    // Nonreproducible candidate: repeated evaluate must be bitwise identical.
    std::size_t calls=0;
    const auto history_dependent=[&](const PhysicalCandidateV51& s){
        auto r=evaluate(s);if(++calls%2==0)r.signed_raw[0]+=.001;return r;};
    const auto first=history_dependent(last);
    bool replay_rejected=false;
    try{ (void)bfs_recertify_six_field_v72(start,start_f,last,first,
          history_dependent,owner,neighbour,1e-8,1e-8,1e-10);
    }catch(const std::exception&){replay_rejected=true;}
    require(replay_rejected,"history-dependent residual was accepted");
    bool invalid_budget=false;
    try{(void)bfs_joint_private_fixed_point_v72(start,start_f,start,start_f,
        evaluate,project,transport,owner,neighbour,choices,17,0.,0.,1e-8,1e-8);}
    catch(const std::exception&){invalid_budget=true;}
    require(invalid_budget,"unbounded fixed point permitted");
    // A partial projection must be repeated, not declared converged from
    // its frozen-coefficient pressure equation. Physical mass halves on each
    // cycle, but still fails the unchanged 1e-10 physical conservation gate.
    const std::vector<BfsJointChoiceV72> half_pressure={{0.,.5}};
    const auto half=bfs_joint_private_fixed_point_v72(start,start_f,start,start_f,
        evaluate,project,transport,owner,neighbour,half_pressure,3,
        0.,0.,1e-8,1e-8);
    require(half.accepted_private_cycles==3 &&
            std::abs(half.initial_mass-.8)<1e-12 &&
            std::abs(half.final_mass-.1)<1e-12 &&
            !half.physical_mass_closed && !half.commit_authorized &&
            half.certificates.size()==3,
            "repeated pressure/SST cycle falsely closed or changed conservation");
    std::cout<<"V72_FIXED_POINT_PASS cycles="<<result.attempted_cycles
             <<" private_accepted="<<result.accepted_private_cycles
             <<" mass="<<result.final_mass<<" commits=0\n";
}