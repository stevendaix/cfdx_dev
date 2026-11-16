#include "cfdx/physics/convective_limited_tangent_v56.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool b){if(!b) throw std::runtime_error("V56 branch tangent contract failed");}
static FaceJetV56 primal(LimitedFaceKindV56 k,double u,double d,double inc,double full,double lo,double hi) {
    return limited_high_face_v56(k,face_jet_v56(u),face_jet_v56(d),face_jet_v56(inc),
                                face_jet_v56(full),face_jet_v56(lo),face_jet_v56(hi));
}
int main() {
    const std::array<LimitedFaceKindV56,5> schemes{{
        LimitedFaceKindV56::SOU,LimitedFaceKindV56::TVD_MINMOD,
        LimitedFaceKindV56::TVD_VANLEER,LimitedFaceKindV56::TVD_VANALBADA,
        LimitedFaceKindV56::MUSCL_CELL}};
    const double h=1.e-6;
    std::size_t certified=0;
    for (auto kind:schemes) for (int i=0;i<40;++i) {
        const double up=1.0+0.09*i, down=up+0.8+0.003*i;
        const double inc=0.13+0.0003*i, full=0.7+0.0002*i;
        const double low=up-0.5, high=down+0.5;
        const double du=0.011,dd=-0.023,di=0.007,df=-0.018,dlo=0.004,dhi=-0.003;
        const auto v=limited_high_face_v56(kind,face_jet_v56(up,du),
             face_jet_v56(down,dd),face_jet_v56(inc,di),face_jet_v56(full,df),
             face_jet_v56(low,dlo),face_jet_v56(high,dhi));
        const auto plus=primal(kind,up+h*du,down+h*dd,inc+h*di,
                               full+h*df,low+h*dlo,high+h*dhi);
        const auto minus=primal(kind,up-h*du,down-h*dd,inc-h*di,
                                full-h*df,low-h*dlo,high-h*dhi);
        if(!v.nonsmooth) {
            const double fd=(plus.value-minus.value)/(2*h);
            check(std::abs(v.tangent-fd)<2.e-7);
            const auto flux=limited_convective_flux_v56(face_jet_v56(0.57,0.03),
                face_jet_v56(up,du),v,0.75);
            const double flux_plus=(0.57+h*0.03)*(up+h*du+0.75*(plus.value-up-h*du));
            const double flux_minus=(0.57-h*0.03)*(up-h*du+0.75*(minus.value-up+h*du));
            check(std::abs(flux.tangent-(flux_plus-flux_minus)/(2*h))<2.e-7);
            ++certified;
        }
    }
    check(certified>100);
    const auto kink=limited_high_face_v56(LimitedFaceKindV56::SOU,
         face_jet_v56(1.,1.),face_jet_v56(1.,-1.),face_jet_v56(0.4),
         face_jet_v56(0.),face_jet_v56(1.),face_jet_v56(1.));
    check(kink.nonsmooth);
    check(limited_convective_flux_v56(face_jet_v56(0.,1.),
        face_jet_v56(1.),face_jet_v56(2.)).nonsmooth);
    bool threw=false;
    try{(void)face_jet_v56(INFINITY);}catch(const std::invalid_argument&){threw=true;}
    check(threw);
}
