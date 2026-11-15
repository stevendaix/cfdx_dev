#include "cfdx/physics/sparse_local_physical_jacobian_v63.h"
#include <array>
#include <stdexcept>
#include <cmath>
#include <iostream>
using namespace cfdx::physics;

#define CHECK(expr) do { if (!(expr)) { \
    std::cerr << "CHECK FAILED: " << #expr << "\n"; return 1; \
} } while(false)
int main() {
    const std::size_t n=12;
    SparseLocalPhysicalJacobianV63::Graph graph(n);
    for(std::size_t i=1;i<n;++i) {
        graph[i].push_back(i-1);graph[i-1].push_back(i);
    }
    PhysicalCandidateV51 x;
    x.cells=n;x.fields.resize(6*n);
    for(std::size_t i=0;i<n;++i) {
        x.fields[i]=2.+0.03*i;
        x.fields[n+i]=0.5+0.01*i;
        x.fields[2*n+i]=0.0;
        x.fields[3*n+i]=0.05*i;
        x.fields[4*n+i]=1.2+0.01*i;
        x.fields[5*n+i]=1.5+0.02*i;
    }
    bool inject_nonlocal=false,inject_kink=false;
    auto physical=[&](const PhysicalCandidateV51& state) {
        PhysicalResidualV51 r;
        r.signed_raw.assign(6*n,0.0);
        r.normalized.assign(6*n,0.0);
        r.face_flux.resize(n+1,1.0);
        r.eddy_viscosity.resize(n,0.);
        r.flux_closed=true;r.flux_fixed_point_defect=0.;
        for(std::size_t i=0;i<n;++i) {
            const double ux=state.fields[i],uy=state.fields[n+i],
                  uz=state.fields[2*n+i],p=state.fields[3*n+i],
                  k=state.fields[4*n+i],w=state.fields[5*n+i];
            const double left=i ? state.fields[i-1] : ux;
            const double right=i+1<n ? state.fields[i+1] : ux;
            const double mu=k/(1.+w);
            // A small manufactured FV stencil with momentum-SST cross terms,
            // viscosity, a face-flux response and a pinned pressure gauge.
            r.signed_raw[i]=ux-0.2*uy+0.07*p+0.14*mu+0.1*(2*ux-left-right);
            r.signed_raw[n+i]=uy+0.03*uz+0.08*k*ux;
            r.signed_raw[2*n+i]=uz;
            r.signed_raw[3*n+i]=0.2*ux-0.1*uy+0.03*p;
            r.signed_raw[4*n+i]=0.8*k-0.1*w+0.01*ux*k;
            r.signed_raw[5*n+i]=w+0.07*k+0.02*uy*w;
            r.eddy_viscosity[i]=std::max(0.,mu);
            r.face_flux[i]=1.+0.03*ux-0.02*p;
        }
        if(inject_nonlocal)r.signed_raw[0]+=0.5*state.fields[n-1];
        if(inject_kink)r.signed_raw[0]+=std::abs(state.fields[0]-x.fields[0]);
        for(std::size_t i=0;i<6*n;++i)r.normalized[i]=r.signed_raw[i];
        r.merit=1.0;
        return r;
    };
    const auto r0=physical(x);
    std::array<double,6> state_scales{{3,2,1,1,2,2}};
    std::array<double,6> row_scales{{1,1,1,1,1,1}};
    auto create=[&]() {
        return SparseLocalPhysicalJacobianV63(x,r0,state_scales,row_scales,
            graph,2,1e-7,true,0,0.,1e-6,physical,0.,0.);
    };
    const auto jac=create();
    const auto& stats=jac.statistics();
    CHECK(stats.trusted);
    CHECK(stats.colors<6*n && stats.nonzeros<36*n*n);
    // Momentum depends on k and omega; SST depends on momentum.
    CHECK(std::abs(jac.matrix()(0,4*n))>1e-6);
    CHECK(std::abs(jac.matrix()(0,5*n))>1e-6);
    CHECK(std::abs(jac.matrix()(4*n,0))>1e-6);
    CHECK(jac.matrix()(3*n,3*n)==1.0);
    std::vector<double> v(6*n);
    for(std::size_t i=0;i<6*n;++i)v[i]=0.01+0.0001*i;
    const auto audit=jac.audit_direction(v,physical);
    if(!(audit.fidelity<1e-6 && audit.sided<1e-6)) {
        std::cerr<<"local audit mismatch "<<audit.fidelity<<" "<<audit.sided<<"\n";
        return 1;
    }
    const auto jac_repeat=create();
    CHECK(jac_repeat.matrix().matvec(v)==jac.matrix().matvec(v));
    inject_nonlocal=true;
    const auto suspicious=create();
    const auto bad=suspicious.audit_direction(v,physical);
    CHECK(bad.fidelity>1e-5 || !suspicious.statistics().trusted);
    inject_nonlocal=false;
    inject_kink=true;
    const auto kink=create();
    CHECK(!kink.statistics().trusted);
    std::cout<<"V63_SPARSE_PHYSICAL_PASS colors="<<stats.colors
             <<" nnz="<<stats.nonzeros<<" evaluations="<<stats.physical_evaluations
             <<" audit="<<audit.fidelity<<"\n";
}