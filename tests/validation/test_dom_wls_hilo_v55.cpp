// V55 WLS arithmetic: five historical 640-cell 3D extrusions. No production
// promotion: a hi+lo residual alone cannot pass original double or HO gates.
#include "cfdx/physics/dom_wls_iterative_refinement_v55.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void ensure(bool pass,const char* text){if(!pass)throw std::runtime_error(text);}
struct Case {
    Mesh mesh;
    FvGeometry geo;
    std::vector<DomSmmSecondMomentV48::Tensor> E;
    std::vector<double> ka,kt,ib;
    DomWallBoundaryConditions bc;
    std::vector<DiscreteDirection> dirs;
    explicit Case(double depth)
        :mesh(cfdx::testing::thermal_v36::make_mesh(8,8,10,.018,0.,depth)){
        for(std::size_t i=0;i<mesh.points().size();++i){
            const double x=mesh.points().x_data()[i];
            const double y=mesh.points().y_data()[i];
            const double z=mesh.points().z_data()[i];
            const double xr=std::cos(.41)*x-std::sin(.41)*y;
            const double yr=std::sin(.41)*x+std::cos(.41)*y;
            mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,yr,
                              -std::sin(.33)*xr+std::cos(.33)*z);
        }
        geo=build_fv_geometry(mesh);
        dirs=make_dom_tensor_gauss_legendre_v39(4,16);
        for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
            const auto& name=mesh.boundary().patch(p).name;
            bc[name]={name=="bottom"?0.:name=="west"?.45:.85,
                      name=="east"?515.:name=="west"?385.:420.};
        }
        const auto n=mesh.n_cells();
        ensure(n==640,"V55 expected exactly 640 cells");
        E.resize(n);ka.resize(n);kt.resize(n);ib.resize(n);
        for(std::size_t c=0;c<n;++c){
            const double xx=geo.cell_centres[c].x;
            E[c]={1./3.+.009*xx,1./3.-.009*xx,1./3.,.014,.006,.007};
            ka[c]=.3+.02*std::abs(xx);
            kt[c]=3.+.1*std::abs(xx);
            ib[c]=blackbody_intensity(420.+10.*xx);
        }
    }
};
static void run_case(double depth,bool dense,bool iterative_enabled=false){
    Case c(depth);
    DomSmmGrayWallV49 historic(c.mesh,c.geo,c.dirs,c.bc,c.E,c.kt,c.ka);
    DomStableWlsV54 stable(c.mesh,c.geo,historic,c.E,c.kt,c.ka);
    if(depth==.06){
        auto wrong=c.E;
        wrong[0][3]+=1e-3;
        bool refused=false;
        try{DomStableWlsV54 bad(c.mesh,c.geo,historic,wrong,c.kt,c.ka);}
        catch(const std::invalid_argument&){refused=true;}
        ensure(refused,"V55 stale tensor coefficients were accepted");
    }
    const auto n=c.mesh.n_cells();
    std::vector<DomHiLoV55> state(n);
    for(std::size_t i=0;i<n;++i)
        state[i]={c.ib[i],c.ib[i]*((i%2)?1.:-1.)*1e-17};
    const auto raw=stable.evaluate_hilo_v55(state);
    const auto rounded=historic.apply(c.ib);
    std::vector<long double> x(n);
    for(std::size_t i=0;i<n;++i) x[i]=dom_to_long_double_v55(state[i]);
    const auto wide=stable.evaluate_extended(x);
    double maximum_difference=0.,original_difference=0.,max_face=0.;
    for(std::size_t i=0;i<n;++i){
        const long double d=dom_to_long_double_v55(raw.divergence_watts[i]);
        const long double e=wide.divergence_watts[i];
        maximum_difference=std::max(maximum_difference,
            static_cast<double>(std::abs(d-e)/std::max(1.L,std::abs(e))));
        original_difference=std::max(original_difference,
            std::abs(rounded[i]-static_cast<double>(d))/
            std::max({1.,std::abs(rounded[i]),std::abs(static_cast<double>(d))}));
    }
    for(std::size_t f=0;f<c.mesh.n_faces();++f){
        const long double d=dom_to_long_double_v55(raw.face_watts[f]);
        max_face=std::max(max_face,
            static_cast<double>(std::abs(d-wide.face_watts[f])));
    }
    ensure(maximum_difference<1e-12,
           "V55 split arithmetic differs from 80-bit flux evaluation");
    ensure(original_difference<1e-9,
           "V55 split action not consistent with historical WLS");
    // Independent reconstruction of face contributions with the retained hi+lo
    // flux: owner and neighbour MUST consume one identical signed value.
    std::vector<DomHiLoV55> independent(n);
    for(std::size_t f=0;f<c.mesh.n_faces();++f){
        const auto o=c.mesh.ownership().owner(f);
        const auto j=c.mesh.ownership().neighbour(f);
        const auto watts=raw.face_watts[f];
        independent[o]+=watts;
        if(j>=0)independent[static_cast<std::size_t>(j)]-=watts;
    }
    for(const auto& w:historic.walls())
        independent[w.owner]+=DomHiLoV55(w.conductance)*state[w.owner];
    for(std::size_t i=0;i<n;++i){
        independent[i]+=dom_four_pi_v55()*DomHiLoV55(c.ka[i])*
                         DomHiLoV55(c.geo.cell_volumes[i])*state[i];
        const long double err=std::abs(dom_to_long_double_v55(independent[i]-
                                                    raw.divergence_watts[i]));
        ensure(err<1e-17L,"V55 reconstructed face ledger not conservative");
    }
    std::cout<<"V55_WLS_FAMILY depth="<<depth
             <<" cells="<<n
             <<" split_vs_80="<<maximum_difference
             <<" split_vs_v48_double="<<original_difference
             <<" face_difference_watts="<<max_face<<'\n';
    if(iterative_enabled){
        const auto rhs=historic.rhs(c.ib);
        const auto iterative=solve_wls_hilo_refinement_v55(
            stable,historic,rhs,c.ib);
        ensure(iterative.amg_setup,"V55: native AMG setup failed");
        ensure(!iterative.publish,"V55: iterative research state was published");
        std::cout<<"V55_WLS_ITERATIVE depth="<<depth
                 <<" initial="<<iterative.initial_split_l2
                 <<" corrected="<<iterative.final_split_l2
                 <<" historic="<<iterative.final_original_double_l2
                 <<" updates="<<iterative.accepted_updates
                 <<" rejects="<<iterative.rejected_updates
                 <<" krylov="<<iterative.krylov_iterations
                 <<" hilo_strict="<<iterative.split_strict
                 <<" original_strict="<<iterative.original_double_strict
                 <<" publish="<<iterative.publish<<'\n';
        ensure(iterative.split_strict,
               "V55: sparse compensated refinement did not converge");
    }
    if(dense){
        std::vector<DomHiLoV55> diagnostic;
        const auto rhs=historic.rhs(c.ib);
        const auto pre=audit_dense_precision_v54(stable,historic,rhs,
                                                  1e-13,640,&diagnostic);
        ensure(pre.extended_strict,"V55: 80-bit diagnostic reference failed");
        const auto gate=audit_wls_hilo_v55(stable,historic,rhs,diagnostic);
        ensure(!gate.publish,"V55: missing physics gates published a state");
        ensure(gate.max_apply_difference<1e-12,
               "V55: split residual differs from 80-bit evaluator");
        ensure(gate.hilo_pass,"V55: compensated arithmetic not strict on extreme");
        // Falsification: a 2% perturbation carried ONLY by the low component
        // must be observed by the genuine hi+lo residual, even though V48
        // sees exactly the same rounded high-component vector.
        auto poisoned=diagnostic;
        poisoned[0].lo+=0.02*poisoned[0].hi;
        const auto poison_gate=audit_wls_hilo_v55(stable,historic,rhs,poisoned);
        ensure(!poison_gate.hilo_pass && !poison_gate.publish,
               "V55: hidden low-component perturbation escaped true residual");
        ensure(poison_gate.historic_double_l2==gate.historic_double_l2,
               "V55: original V48 changed despite identical hi-only state");
        std::cout<<"V55_WLS_FALSIFICATION poisoned_hilo="
                 <<poison_gate.hilo_l2
                 <<" same_double="<<poison_gate.historic_double_l2
                 <<" publish="<<poison_gate.publish<<'\n';
        // The historical V48 result is explicitly measured, not substituted by
        // the split norm when deciding convergence of the legacy double field.
        std::cout<<"V55_WLS_EXTREME depth="<<depth
                 <<" algebraic80="<<pre.extended_l2
                 <<" hilo="<<gate.hilo_l2
                 <<" original_v48_double="<<gate.historic_double_l2
                 <<" original_strict="<<gate.original_double_pass
                 <<" max_apply_gap="<<gate.max_apply_difference
                 <<" publish="<<gate.publish<<'\n';
    }
}
int main(int argc,char** argv){
    try{
        const std::string arg=(argc==1?"--quick":argc==2?argv[1]:"bad");
        if(arg=="--quick") run_case(.06,false);
        else if(arg=="--family")
            for(double depth:{.06,.03,.012,.006,.003})run_case(depth,false);
        else if(arg=="--extreme")run_case(.003,true,true);
        else if(arg=="--iterative-family")
            for(double depth:{.06,.03,.012,.006,.003})run_case(depth,false,true);
        else throw std::invalid_argument("V55 expects --quick, --family, --extreme, --iterative-family");
        std::cout<<"V55_WLS_PASS mode="<<arg<<'\n';
        return 0;
    }catch(const std::exception& e){
        std::cerr<<"V55_WLS_FAIL "<<e.what()<<'\n';return 1;
    }
}
