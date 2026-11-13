// V54 research-only numerical-linearity diagnostic. Never substitutes the
// original V48 WLS physical oracle or certifies a failed strict solver.
#include "cfdx/physics/dom_wls_precision_audit_v54.h"
#include "thermal_3d_mesh_v36.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static void check(bool ok,const char* what){if(!ok)throw std::runtime_error(what);}
static void run(double depth){
    auto mesh=cfdx::testing::thermal_v36::make_mesh(8,8,10,.018,0.,depth);
    for(std::size_t i=0;i<mesh.points().size();++i){
        const double x=mesh.points().x_data()[i],y=mesh.points().y_data()[i];
        const double z=mesh.points().z_data()[i];
        const double xr=std::cos(.41)*x-std::sin(.41)*y;
        const double yr=std::sin(.41)*x+std::cos(.41)*y;
        mesh.points().set(i,std::cos(.33)*xr+std::sin(.33)*z,yr,
                            -std::sin(.33)*xr+std::cos(.33)*z);
    }
    const auto geo=build_fv_geometry(mesh);
    const auto dirs=make_dom_tensor_gauss_legendre_v39(4,16);
    DomWallBoundaryConditions bc;
    for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
        const auto name=mesh.boundary().patch(p).name;
        bc[name]={name=="bottom"?0.:name=="west"?.45:.85,
                 name=="east"?515.:name=="west"?385.:420.};
    }
    const std::size_t nc=mesh.n_cells();
    std::vector<DomSmmSecondMomentV48::Tensor> E(nc);
    std::vector<double> kt(nc),ka(nc),ib(nc);
    for(std::size_t c=0;c<nc;++c){
        const double xx=geo.cell_centres[c].x;
        E[c]={1./3.+.009*xx,1./3.-.009*xx,1./3.,.014,.006,.007};
        ka[c]=.3+.02*std::abs(xx);
        kt[c]=3.+.1*std::abs(xx);
        ib[c]=blackbody_intensity(420.+10.*xx);
    }
    DomSmmGrayWallV49 old(mesh,geo,dirs,bc,E,kt,ka);
    DomStableWlsV54 stable(mesh,geo,old,E,kt,ka);
    const auto o=old.apply(ib),s=stable.apply(ib);
    double discrepancy=0.;
    for(std::size_t c=0;c<nc;++c)
        discrepancy=std::max(discrepancy,std::abs(o[c]-s[c])/
                            std::max({1.,std::abs(o[c]),std::abs(s[c])}));
    check(discrepancy<1e-10,"V54: stable WLS is inconsistent with V48 original operator");
    std::vector<double> y(nc),sum(nc);
    for(std::size_t c=0;c<nc;++c){
        y[c]=ib[c]*(.2+.02*std::sin(double(c)));
        sum[c]=ib[c]+y[c];
    }
    const auto sy=stable.apply(y),ss=stable.apply(sum);
    const auto oy=old.apply(y),os=old.apply(sum);
    double defect_stable=0.,defect_old=0.,scale=1.;
    for(std::size_t c=0;c<nc;++c){
        defect_stable=std::max(defect_stable,std::abs(ss[c]-s[c]-sy[c]));
        defect_old=std::max(defect_old,std::abs(os[c]-o[c]-oy[c]));
        scale=std::max({scale,std::abs(ss[c]),std::abs(os[c])});
    }
    defect_stable/=scale;defect_old/=scale;
    // The stencils are mathematical equivalents. A small arithmetic gain does
    // not automatically establish original-operator strict convergence.
    check(defect_stable<1e-10,"V54: stable WLS remains unacceptably nonlinear");
    const auto faces=stable.evaluate(ib);
    check(faces.face_watts.size()==mesh.n_faces() &&
          faces.divergence_watts.size()==nc,
          "V54: incomplete stable face/current ledger");
    for(std::size_t f=0;f<mesh.n_faces();++f)
        check(std::isfinite(faces.face_watts[f]),"V54: nonfinite stable face current");
    if(depth<.004 && std::numeric_limits<long double>::digits>
                     std::numeric_limits<double>::digits){
        const auto precision=audit_dense_precision_v54(stable,old,old.rhs(ib));
        check(precision.extended_strict,
              "V54: 80-bit reference did not satisfy strict algebraic threshold");
        // This test explicitly does NOT certify rounded fields using the 80-bit
        // residual. The historical operator is the only authority for that gate.
        check(precision.split_hi_lo_l2<=1e-13,
              "V54: high/low stored state failed extended arithmetic gate");
        check(precision.extended_l2<precision.authoritative_v48_l2,
              "V54: extended arithmetic did not improve the reference defect");
        std::cout<<"V54_EXTENDED_REFERENCE depth="<<depth
                 <<" algebraic80="<<precision.extended_l2
                 <<" split_hi_lo="<<precision.split_hi_lo_l2
                 <<" rounded_dense="<<precision.rounded_dense_l2
                 <<" rounded_stable="<<precision.rounded_stable_l2
                 <<" authoritative_v48="<<precision.authoritative_v48_l2
                 <<" publish="<<precision.original_strict<<'\n';
    }
    std::cout<<"V54_WLS depth="<<depth<<" cells="<<nc
             <<" difference_original="<<discrepancy
             <<" nonadditivity_original="<<defect_old
             <<" nonadditivity_stable="<<defect_stable
             <<" status=RESEARCH_ONLY_NOT_STRICT_CERTIFICATE\n";
}
int main(){try{run(.003);run(.06);std::cout<<"V54_WLS_PASS\n";return 0;}
    catch(const std::exception& e){std::cerr<<"V54_WLS_FAIL "<<e.what()<<'\n';return 1;}}
