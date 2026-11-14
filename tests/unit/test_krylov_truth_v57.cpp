#include "cfdx/physics/global_physical_newton_v52.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace cfdx::physics;
static void check(bool p,const char* msg){if(!p)throw std::runtime_error(msg);}
int main(){
    const std::vector<double> rhs{1.,2.,-3.,4.,0.5};
    const auto identity=[](const std::vector<double>& v){return v;};
    const auto A=[](const std::vector<double>& v){
        return std::vector<double>{2*v[0]+v[1],3*v[1]-0.5*v[2],
            2*v[2]+0.3*v[3],v[3]+0.2*v[4],1.7*v[4]};
    };
    std::vector<KrylovTruthSnapshotV57> history;
    std::vector<double> x;
    double t=0.,a=0.,d=0.;
    const bool ok=global_physical_fgmres_v52(A,identity,rhs,5,1e-10,
        x,t,&a,&d,{},[&](const KrylovTruthSnapshotV57& snap){history.push_back(snap);});
    check(ok && !history.empty() && history.size()>=2 &&
          history.back().true_relative<1e-10,"V57 exact Krylov did not close");
    for(std::size_t i=0;i<history.size();++i){
        check(history[i].direction==i+1 && history[i].arnoldi_gap<1e-10 &&
              history[i].action_nonadditivity_relative<1e-10,
              "V57 failed intermediate Arnoldi/true linear residual parity");
    }
    const auto nonlinear=[](const std::vector<double>& v){
        auto out=v;double sum=0.;for(double u:v)sum+=u*u;
        for(double& u:out)u*=1.+3.*std::sqrt(sum);
        return out;
    };
    history.clear();
    (void)global_physical_fgmres_v52(nonlinear,identity,rhs,5,1e-10,
        x,t,&a,&d,{},[&](const KrylovTruthSnapshotV57& snap){history.push_back(snap);});
    bool seen=false;
    for(const auto& h:history)
        seen|=h.arnoldi_gap>1e-4 || h.action_nonadditivity_relative>1e-4;
    check(seen,"V57 failed to identify false Arnoldi projection for nonlinear action");
    std::cout<<"V57_KRYLOV_PASS steps="<<history.size()<<" nonadditive_detected=1\n";
}
