#include "cfdx/physics/physical_six_candidate_v51.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool value) {if(!value) throw std::runtime_error("V51 test failure");}
int main() {
    PhysicalCandidateV51 x;
    x.cells=2;
    x.fields={2,3, 0.2,0.4, 0,0, 4,5, 1,2, 10,11};
    auto evaluator=[](const PhysicalCandidateV51& state) {
        const std::size_t n=state.cells;
        PhysicalResidualV51 r;
        r.signed_raw.resize(6*n);
        r.normalized.resize(6*n);
        r.face_flux={state.fields[0]+state.fields[4*n],state.fields[n+1]};
        r.eddy_viscosity.resize(n);
        for(std::size_t c=0;c<n;++c) {
            // Smooth coupled proxy: momentum and SST both depend on ALL
            // candidate blocks. Catches accidentally frozen/lagged fields.
            const double u=state.fields[c], v=state.fields[n+c];
            const double p=state.fields[3*n+c], k=state.fields[4*n+c];
            const double w=state.fields[5*n+c];
            r.eddy_viscosity[c]=k/w;
            for(std::size_t j=0;j<6;++j) {
                const double val=(j+1)*u+(j+2)*v+(j+3)*p+
                                 (j+4)*k+(j+5)*w+(j+1)*u*k;
                r.signed_raw[j*n+c]=val;
                r.normalized[j*n+c]=val/10;
            }
        }
        r.flux_fixed_point_defect=0.0;
        r.flux_closed=true;
        r.merit=1;
        return r;
    };
    const auto a=evaluate_physical_candidate_v51(x,evaluator);
    const auto b=evaluate_physical_candidate_v51(x,evaluator);
    verify_candidate_replay_v51(a,b);
    auto bad=b; bad.face_flux[0]+=0.0001;
    bool threw=false;
    try{verify_candidate_replay_v51(a,bad);}catch(const std::runtime_error&){threw=true;}
    check(threw);
    std::vector<double> v(x.fields.size(),0.0);
    for(std::size_t j=0;j<6;++j) v[j*x.cells]=0.02*(j+1);
    const auto jv=physical_directional_jv_v51(x,v,evaluator,1e-5);
    const auto jv_half=physical_directional_jv_v51(x,v,evaluator,5e-6);
    for(std::size_t i=0;i<jv.size();++i)
        check(std::abs(jv_half[i]-jv[i])<1e-5);
    for(std::size_t c=0;c<x.cells;++c)
        for(std::size_t j=0;j<6;++j) {
            const double u=x.fields[c],k=x.fields[4*x.cells+c];
            const double expected=(j+1)*v[c]+(j+2)*v[x.cells+c]+
                (j+3)*v[3*x.cells+c]+(j+4)*v[4*x.cells+c]+
                (j+5)*v[5*x.cells+c]+(j+1)*(k*v[c]+u*v[4*x.cells+c]);
            check(std::abs(jv[j*x.cells+c]-expected)<1e-5);
        }
    auto invalid=x; invalid.fields[5*x.cells]=0;
    threw=false;
    try{evaluate_physical_candidate_v51(invalid,evaluator);}catch(const std::invalid_argument&){threw=true;}
    check(threw);
    auto invalid_result=[](const PhysicalCandidateV51& s){
        PhysicalResidualV51 r;
        r.signed_raw.resize(6*s.cells);
        r.normalized.resize(6*s.cells);
        r.face_flux={1};r.eddy_viscosity.resize(s.cells);
        r.flux_closed=false;r.merit=0;return r;
    };
    threw=false;
    try{evaluate_physical_candidate_v51(x,invalid_result);}catch(const std::runtime_error&){threw=true;}
    check(threw);
    auto invalid_direction=v; invalid_direction[0]=std::numeric_limits<double>::quiet_NaN();
    threw=false;
    try{physical_directional_jv_v51(x,invalid_direction,evaluator,1e-5);}catch(const std::invalid_argument&){threw=true;}
    check(threw);
    check(x.fields[0]==2 && x.fields[5*x.cells]==10); // no mutation
}
