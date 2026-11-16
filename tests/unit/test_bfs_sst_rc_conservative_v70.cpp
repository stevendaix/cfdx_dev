#include "cfdx/physics/bfs_sst_rc_conservative_v70.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool b,const char* why){if(!b)throw std::runtime_error(why);}
template<class F>static void reject(F&& f){bool seen=false;try{f();}catch(const std::exception&){seen=true;}require(seen,"V70 falsified state must fail closed");}
int main(){
    constexpr std::array<std::size_t,3> owners={0,0,1};
    constexpr std::array<int,3> neighbours={-1,1,-1};
    const auto own=[&](std::size_t f){return owners.at(f);};
    const auto nb=[&](std::size_t f){return neighbours.at(f);};
    const auto flux_at_pressure=[](const std::vector<double>& p){
        if(p.size()!=2)throw std::runtime_error("bad pressure size");
        return std::vector<double>{-1.,.2+p[0]-p[1],.3+p[1]};
    };
    PhysicalCandidateV51 state;
    state.cells=2;state.fields.assign(12,1.);
    state.fields[6]=0.;state.fields[7]=0.; // p block
    state.fields[8]=.4;state.fields[9]=.8; // k block
    state.fields[10]=6.;state.fields[11]=8.; // omega block
    auto old=state;
    const auto response=[&](const std::vector<double>& delta){
        BfsPressureResponseV70 r;
        r.state=state;
        for(std::size_t i=0;i<2;++i)r.state.fields[6+i]+=delta.at(i);
        r.face_flux=flux_at_pressure({r.state.fields[6],r.state.fields[7]});
        return r;
    };
    PhysicalResidualV51 f;
    f.flux_closed=true;f.flux_fixed_point_defect=0.;
    f.signed_raw.assign(12,1.);f.normalized.assign(12,1.);
    f.face_flux=flux_at_pressure({0.,0.});f.eddy_viscosity={.01,.02};
    f.inverse_momentum_diagonal.assign(6,1.);f.merit=1.;
    const auto a=bfs_mass_balance_v70(2,f.face_flux,own,nb);
    require(std::abs(a.boundary_signed+.7)<1e-15,"initial boundary imbalance");
    require(std::abs(a.divergence_linf-.8)<1e-15,"initial continuity");
    const auto p=bfs_project_pressure_v70(state,f,response,own,nb,true,0,1e-12);
    require(p.pressure_unknowns==2 && p.divergence_after<1e-12,"conservative projection");
    require(std::abs(p.state.fields[6]-1.5)<1e-12 &&
            std::abs(p.state.fields[7]-.7)<1e-12,"correct pressure increments");
    require(std::abs(p.face_flux[1]-1.)<1e-12 &&
            std::abs(p.face_flux[2]-1.)<1e-12,"owner neighbour and boundary balance");
    // A coupled physical-velocity response: each pressure correction also
    // changes U with a frozen positive momentum inverse diagonal. Conserved
    // face fluxes depend on BOTH updated pressure and updated U.
    const auto coupled_response=[&](const std::vector<double>& delta){
        auto r=response(delta);
        r.state.fields[0]-=.2*delta[0];
        r.state.fields[1]-=.2*delta[1];
        const double p0=r.state.fields[6],p1=r.state.fields[7];
        const double u0=r.state.fields[0],u1=r.state.fields[1];
        r.face_flux={-1.,.2+(p0-p1)+.3*(u0-u1),.3+p1+.2*(u1-1.)};
        return r;
    };
    const auto coupled=bfs_project_pressure_v70(
        state,f,coupled_response,own,nb,true,0,1e-12);
    require(coupled.fully_projected && coupled.divergence_after<1e-12,
            "pressure/velocity combined face projection must conserve mass");
    require(coupled.state.fields[0]<state.fields[0] &&
            coupled.state.fields[1]<state.fields[1] &&
            coupled.state.fields[6]>state.fields[6],
            "the pressure correction must also update physical velocity");
    require(old.fields==state.fields,"private pressure projection changed live state");
    const auto partial=bfs_project_pressure_v70(state,f,response,own,nb,true,0,
                                                1e-12,.25);
    require(!partial.fully_projected && partial.projection_fraction==.25 &&
            partial.divergence_after<partial.divergence_before &&
            partial.divergence_after>1e-12,
            "partial pressure step must reduce but NOT claim mass closure");
    auto after=f;
    after.face_flux=p.face_flux;
    after.signed_raw.assign(12,.5);
    after.normalized.assign(12,.5);
    after.signed_raw[6]=0.;after.signed_raw[7]=0.;
    after.normalized[6]=0.;after.normalized[7]=0.;
    after.merit=.5;
    const auto gate=bfs_certify_projected_trial_v70(
        old,f,p,after,own,nb,1e-8,1e-8,1e-12);
    require(gate.accepted && gate.joint_allowed && gate.monotone_allowed,
            "independently reassembled V65 must pass");
    reject([&]{auto altered=f;altered.face_flux[1]+=.01;
        bfs_project_pressure_v70(state,altered,response,own,nb,true,0);});
    reject([&]{bfs_project_pressure_v70(state,f,response,own,nb,false,0);});
    reject([&]{bfs_project_pressure_v70(state,f,
        [&](const auto& delta){auto r=response(delta);r.face_flux={-1.,.2,.3};return r;},own,nb,true,0);});
    reject([&]{bfs_project_pressure_v70(state,f,response,own,nb,true,8);});
    reject([&]{bfs_project_pressure_v70(state,f,response,own,nb,true,0,0.);});
    reject([&]{bfs_project_pressure_v70(state,f,response,own,nb,true,0,1e-12,0.);});
    reject([&]{bfs_project_pressure_v70(state,f,response,own,nb,true,0,1e-12,1.1);});
    reject([&]{auto broken=after;broken.face_flux[2]+=.001;
        bfs_certify_projected_trial_v70(old,f,p,broken,own,nb,1e-8,1e-8,1e-12);});
    reject([&]{auto broken=after;broken.signed_raw[6]+=.01;
        bfs_certify_projected_trial_v70(old,f,p,broken,own,nb,1e-8,1e-8,1e-12);});
    reject([&]{auto broken=after;broken.signed_raw[10]=20.;
        bfs_certify_projected_trial_v70(old,f,p,broken,own,nb,1e-8,1e-8,1e-12);});
    std::cout<<"V70_CONSERVATIVE_PRESSURE_PASS p0="<<p.state.fields[6]
             <<" p1="<<p.state.fields[7]
             <<" residual="<<p.divergence_after<<" no_publication=1\n";
}