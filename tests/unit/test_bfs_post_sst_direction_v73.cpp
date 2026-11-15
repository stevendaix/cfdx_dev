@@ -0,0 +1,140 @@
// V73: the same U,p,k,omega state and real V71 operators for every private
// search direction; verify face, rAU, momentum/SST evidence and fail closed.
#include "cfdx/physics/bfs_post_sst_direction_v73.h"
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool ok,const char* why) {if(!ok)throw std::runtime_error(why);}
int main(){
    const std::array<std::size_t,3> own{0,0,1};
    const std::array<int,3> nei{-1,1,-1};
    const auto owner=[&](std::size_t f){return own.at(f);};
    const auto neighbour=[&](std::size_t f){return nei.at(f);};
    PhysicalCandidateV51 start;start.cells=2;start.fields.assign(12,0.);
    start.fields[0]=1.;start.fields[1]=1.;
    start.fields[8]=.3;start.fields[9]=.3;
    start.fields[10]=2.;start.fields[11]=2.;
    const auto flux=[](const PhysicalCandidateV51& x){
        return std::vector<double>{-1.,.2+x.fields[6]-x.fields[7]
                            +.2*(x.fields[0]-x.fields[1]),.3+x.fields[7]};
    };
    int evaluated=0;
    const auto eval=[&](const PhysicalCandidateV51& x){
        ++evaluated;
        PhysicalResidualV51 r;r.flux_closed=true;r.flux_fixed_point_defect=0.;
        r.face_flux=flux(x);r.signed_raw.assign(12,0.);
        r.normalized.assign(12,0.);r.inverse_momentum_diagonal.assign(6,.2);
        r.eddy_viscosity.assign(2,0.);
        for(std::size_t i=0;i<2;++i) {
            r.signed_raw[i]=x.fields[i]-(1.5-.1*x.fields[8+i]);
            r.signed_raw[8+i]=x.fields[8+i]-(.2+.02*x.fields[i]);
            r.signed_raw[10+i]=x.fields[10+i]-(2.+.1*x.fields[i]);
            r.eddy_viscosity[i]=x.fields[8+i]/x.fields[10+i];
        }
        r.signed_raw[6]=r.face_flux[0]+r.face_flux[1];
        r.signed_raw[7]=-r.face_flux[1]+r.face_flux[2];
        r.merit=0.;
        for(std::size_t j=0;j<12;++j){r.normalized[j]=r.signed_raw[j];
            r.merit=std::max(r.merit,std::abs(r.signed_raw[j]));}
        return r;
    };
    const auto transport=[&](const PhysicalCandidateV51& x,
                             const std::vector<double>& f){
        check(f.size()==3,"wrong candidate flux");
        auto y=x;
        for(std::size_t i=0;i<2;++i) {
            y.fields[8+i]=.2+.02*x.fields[i];
            y.fields[10+i]=2.+.1*x.fields[i];
        }
        return y;
    };
    const auto project=[&](const PhysicalCandidateV51& x,
                           const PhysicalResidualV51& f,double alpha){
        const auto response=[&](const std::vector<double>& dp){
            BfsPressureResponseV70 p;p.state=x;
            p.state.fields[6]+=dp[0];p.state.fields[7]+=dp[1];
            p.face_flux=flux(p.state);return p;};
        return bfs_project_pressure_v70(x,f,response,owner,neighbour,true,0,
                                         1e-10,alpha);
    };
    const auto old_f=eval(start);
    const auto saved=start.fields;
    const auto flux_saved=old_f.face_flux;
    const std::vector<BfsJointChoiceV72> samples{{0.,.5},{.25,1.}};
    int observed=0;
    const auto observe=[&](const PhysicalCandidateV51& x,
                           const PhysicalResidualV51& f,
                           const BfsDirectionAttemptV73& a){
        ++observed;
        check(x.fields.size()==12 && f.face_flux.size()==3,"observed wrong FV data");
        check(a.cells[0].cell==0 && a.cells[1].cell==1,
              "face owner neighbour cell audit wrong");
        check(std::abs(a.face_flux_after-f.face_flux[1])<1e-15,
              "wrong interface flux in observation");
        check(!a.eligible_for_commit,"V73 observation forged a native certificate");
    };
    const auto result=bfs_probe_post_sst_directions_v73(
        start,old_f,start,old_f,eval,project,transport,owner,neighbour,
        1,0,1,samples,observe,0.,0.,1e-8,1e-8,1e-10);
    check(result.attempts.size()==2 && result.evaluated_count==2 &&
          observed==2,"V73 failed to evaluate real private Picard alternatives");
    check(result.commits==0 && !result.commit_authorized &&
          !result.native_momentum_linear_recertified &&
          !result.native_sst_transport_recertified,
          "V73 falsely issued native certificates");
    check(start.fields==saved && old_f.face_flux==flux_saved,
          "V73 modified the authoritative state");
    check(evaluated>10,"V73 did not re-evaluate independent physical residuals");
    check(result.attempts[0].face_flux_before==flux_saved[1],
          "V73 wrong same-state interface baseline");
    // The same interface must be real, oriented owner-to-neighbour.
    bool wrong_topology=false;
    try{(void)bfs_probe_post_sst_directions_v73(start,old_f,start,old_f,
         eval,project,transport,owner,neighbour,1,1,0,samples,observe,
         0.,0.,1e-8,1e-8);}
    catch(const std::invalid_argument&){wrong_topology=true;}
    check(wrong_topology,"V73 allowed reversing face 3 topology");
    // The F_phys replay must reject stale rAU/flux evidence even when fields
    // themselves have not changed.
    bool stale=false;auto bad_f=old_f;bad_f.face_flux[1]+=.01;
    try{(void)bfs_probe_post_sst_directions_v73(start,old_f,start,bad_f,
         eval,project,transport,owner,neighbour,1,0,1,samples,observe,
         0.,0.,1e-8,1e-8);}
    catch(const std::exception&){stale=true;}
    check(stale,"V73 accepted stale reconstructed face flux");
    bool unbounded=false;
    try{auto extra=samples;extra.resize(25,{0.,.5});
         (void)bfs_probe_post_sst_directions_v73(start,old_f,start,old_f,
         eval,project,transport,owner,neighbour,1,0,1,extra,observe,
         0.,0.,1e-8,1e-8);}
    catch(const std::invalid_argument&){unbounded=true;}
    check(unbounded,"V73 allowed unbounded expensive direction search");
    // A transport exception cannot be turned into a fake accepted state.
    const auto broken_transport=[](const PhysicalCandidateV51&,
                                   const std::vector<double>&)->PhysicalCandidateV51{
        throw std::runtime_error("forced SST failure");};
    const auto broken=bfs_probe_post_sst_directions_v73(start,old_f,start,old_f,
        eval,project,broken_transport,owner,neighbour,1,0,1,samples,observe,
        0.,0.,1e-8,1e-8);
    check(broken.evaluated_count==0 && broken.guard_accepted_count==0 &&
          broken.commits==0 && !broken.commit_authorized,
          "V73 fake accepted a failed transport");
    // The native source audit is a required part of a usable evaluation.
    // An incomplete omega decomposition cannot be counted as a success.
    const auto failed_source_audit=[](const PhysicalCandidateV51&,
                                      const PhysicalResidualV51&,
                                      const BfsDirectionAttemptV73&){
        throw std::runtime_error("forced source-budget closure failure");};
    const auto refused=bfs_probe_post_sst_directions_v73(
        start,old_f,start,old_f,eval,project,transport,owner,neighbour,
        1,0,1,samples,failed_source_audit,0.,0.,1e-8,1e-8,1e-10);
    check(refused.evaluated_count==0 && refused.guard_accepted_count==0 &&
          refused.commits==0 && !refused.commit_authorized,
          "V73 counted an incomplete SST source audit");
    std::cout<<"V73_DIRECTION_PASS evaluated="<<result.evaluated_count
             <<" native_certificates=0 committed=0\n";
}
