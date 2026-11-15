#include "cfdx/physics/bfs_sst_implicit_block_v76.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
using namespace cfdx::physics;
struct State {std::size_t cells=2;std::vector<double> fields=std::vector<double>(12,0.);};
struct Phys {std::vector<double> signed_raw=std::vector<double>(12,0.);};
static void need(bool b,const char* m){if(!b)throw std::runtime_error(m);}
int main(){
    // Physical analog: coupled k/omega equations with production-destruction
    // and neighbor (face) coupling. The root is at k0=2,w0=3,k1=1,w1=4.
    auto evaluate=[](const State& x){
        Phys f;
        const double k0=x.fields[8],w0=x.fields[10];
        const double k1=x.fields[9],w1=x.fields[11];
        // The off-diagonal terms mimic production and transport coupling.
        f.signed_raw[8]=k0*w0+0.15*(k0-k1)-6.15;
        f.signed_raw[10]=w0*w0+0.08*(w0-w1)-8.92;
        f.signed_raw[9]=k1*w1+0.12*(k1-k0)-3.88;
        f.signed_raw[11]=w1*w1+0.07*(w1-w0)-16.07;
        return f;
    };
    State x;x.fields[8]=2.2;x.fields[9]=1.12;x.fields[10]=3.25;x.fields[11]=4.2;
    const Phys fx=evaluate(x);
    auto guard=[](const State& x0,const Phys& f0,const State& x1,const Phys& f1){
        (void)x0;(void)x1;
        double before=0.,after=0.;
        for(const auto i:{8,9,10,11}){
            before=std::max(before,std::abs(f0.signed_raw[i]));
            after=std::max(after,std::abs(f1.signed_raw[i]));
        }
        return after<before;
    };
    State answer;
    Phys ansf;
    const auto r=bfs_private_implicit_sst_block_v76(x,fx,0,1,1e-12,1e-12,
        evaluate,guard,&answer,&ansf);
    need(r.attempted && r.selected && r.scaled_local_after<r.scaled_local_before,
         "positive 4x4 SST correction not accepted");
    need(!r.production_commit_authorized && !r.native_sst_linear_certificate &&
         !r.native_momentum_linear_certificate,"fake production certificate");
    for(const auto i:{8,9,10,11})need(answer.fields[i]>0.,"negative SST state");
    State untouched;
    auto deny=[](const State&,const Phys&,const State&,const Phys&){return false;};
    const auto rejected=bfs_private_implicit_sst_block_v76(x,fx,0,1,1e-12,1e-12,
        evaluate,deny,&untouched);
    need(!rejected.selected && untouched.fields==x.fields,"rejected step mutated state");
    bool bad=false;
    try {(void)bfs_private_implicit_sst_block_v76(x,fx,0,0,1e-12,1e-12,evaluate,guard);}
    catch(const std::invalid_argument&){bad=true;}
    need(bad,"duplicate cells accepted");
    auto badstate=x;badstate.fields[8]=0.;
    const auto invalid=bfs_private_implicit_sst_block_v76(
        badstate,evaluate(badstate),0,1,1e-12,1e-12,evaluate,guard);
    need(!invalid.selected && invalid.reason=="nonpositive_start","floor failure not detected");
    // A strongly nonsmooth residual must fail CLOSED, not be called a Jacobian.
    auto kink=[](const State& xx){
        Phys f;
        f.signed_raw[8]=std::abs(xx.fields[8]-2.2)+0.01*(xx.fields[10]-3.25)+0.4;
        f.signed_raw[10]=xx.fields[10]-3.25+0.1;
        f.signed_raw[9]=xx.fields[9]-1.0;
        f.signed_raw[11]=xx.fields[11]-4.0;
        return f;
    };
    const auto non_smooth=bfs_private_implicit_sst_block_v76(
        x,kink(x),0,1,1e-12,1e-12,kink,guard);
    need(!non_smooth.selected && !non_smooth.derivative_regular,"nonsmooth branch accepted");
    const auto shifted=bfs_private_implicit_sst_block_v76(
        x,fx,0,1,1e-12,1e-12,evaluate,guard,static_cast<State*>(nullptr),static_cast<Phys*>(nullptr),0.15);
    need(shifted.attempted && !shifted.production_commit_authorized,
         "regularized local block lacks guards");
    std::cout<<"V76_SST_BLOCK_TEST_PASS accepted_alpha="<<r.accepted_alpha
             <<" eval="<<r.residual_evaluations<<" selected="<<r.selected<<"\n";
}
