#include "cfdx/physics/global_physical_newton_v52.h"
#include <cmath>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void require(bool ok,const char* what){if(!ok)throw std::runtime_error(what);}
int main(){
    const std::vector<double> rhs{1.0,2.0,-3.0,4.0};
    const auto identity=[](const std::vector<double>& v){return v;};
    const auto exact=[](const std::vector<double>& v){
        return std::vector<double>{2*v[0]+v[1],3*v[1]-0.5*v[2],
                                   2*v[2]+0.3*v[3],v[3]};
    };
    std::vector<double> x;
    double true_r=0.,estimated=0.,linearity=0.;
    const bool accurate=global_physical_fgmres_v52(exact,identity,rhs,4,1e-9,
        x,true_r,&estimated,&linearity);
    require(accurate && true_r<1e-9 && estimated<1e-9 &&
            linearity<1e-10,"V53 linear operator must be certified");
    // A deterministic non-linear directional action: a Krylov projection is
    // not a certificate of the original equation, even with small estimate.
    const auto nonlinear=[](const std::vector<double>& v){
        const double factor=1.0+4.0*std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
        std::vector<double> result=v;
        for(double& value:result)value*=factor;
        return result;
    };
    x.clear();true_r=0.;estimated=0.;linearity=0.;
    (void)global_physical_fgmres_v52(nonlinear,identity,rhs,4,1e-6,
        x,true_r,&estimated,&linearity);
    require(linearity>0.01 && true_r>estimated,
            "V53 did not detect nonlinear directional action");
    auto state=PhysicalCandidateV51{};state.cells=1;
    state.fields={0.,0.,0.,0.,1.,2.};
    auto controls=GlobalPhysicalNewtonControlsV52{};
    controls.krylov_dimension=64;
    bool rejected=false;
    try{validate_global_physical_newton_v52(controls,state);}
    catch(const std::invalid_argument&){rejected=true;}
    require(rejected,"V52 ordinary Krylov budget relaxed by V53");
    controls.v53_extended_krylov=true;
    validate_global_physical_newton_v52(controls,state);
}
