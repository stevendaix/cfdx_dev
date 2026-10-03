// N16.1 — WENO5-Z finite-volume reconstruction verification.
#include "cfdx/core/numerics/high_order_reconstruction.h"
#include <array>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using cfdx::core::WENO5ZDiagnostics;
using cfdx::core::weno5z_reconstruct_right;
using cfdx::core::weno5z_reconstruct_right_bounded;
namespace {
void require(bool ok,const std::string& msg){if(!ok) throw std::runtime_error(msg);}
double order(double e0,double e1){return std::log(e0/e1)/std::log(2.0);}
double sine_average(double x,double h){
    return (std::cos(2.0*M_PI*x)-std::cos(2.0*M_PI*(x+h)))/(2.0*M_PI*h);
}
void exactness(){
    const std::array<double,5> c{{3.25,3.25,3.25,3.25,3.25}};
    WENO5ZDiagnostics d{};
    require(std::abs(weno5z_reconstruct_right(c,1e-40,2.0,&d)-3.25)<1e-13,
            "constant field is not exact");
    const std::array<double,5> linear{{-2.0,-1.0,0.0,1.0,2.0}};
    require(std::abs(weno5z_reconstruct_right(linear)-0.5)<1e-12,
            "linear cell-average field is not exact");
}
void smooth_order(){
    const std::vector<std::size_t> ns{20,40,80,160};
    std::vector<double> e;
    for(const auto n:ns){
        const double h=1.0/static_cast<double>(n);
        std::vector<double> a(n);
        for(std::size_t i=0;i<n;++i) a[i]=sine_average(i*h,h);
        double s=0.0; std::size_t count=0;
        for(std::size_t i=2;i+2<n;++i){
            const std::array<double,5> u{{a[i-2],a[i-1],a[i],a[i+1],a[i+2]}};
            const double exact=std::sin(2.0*M_PI*(static_cast<double>(i)+1.0)*h);
            const double err=weno5z_reconstruct_right(u)-exact;
            s+=err*err; ++count;
        }
        e.push_back(std::sqrt(s/static_cast<double>(count)));
        std::cout<<"N16_WENO5Z n="<<n<<" L2="<<e.back()<<"\\n";
    }
    for(std::size_t i=1;i<e.size();++i){
        const double p=order(e[i-1],e[i]);
        std::cout<<"N16_WENO5Z observed_order="<<p<<"\\n";
        require(p>4.0,"smooth WENO5-Z observed order must exceed 4");
    }
}
void bounded(){
    const std::array<double,5> u{{0.0,0.0,1.0,1.0,1.0}};
    WENO5ZDiagnostics d{};
    const double v=weno5z_reconstruct_right_bounded(u,1e-40,2.0,&d);
    require(v>=-1e-14 && v<=1.0+1e-14,"bounded WENO5-Z escaped stencil envelope");
    require(std::isfinite(v),"bounded WENO5-Z is not finite");
}
void validation(){
    const std::array<double,5> u{{1,1,1,1,1}};
    bool threw=false; try{(void)weno5z_reconstruct_right(u,0.0);}catch(const std::invalid_argument&){threw=true;}
    require(threw,"epsilon=0 must be rejected");
    threw=false; try{(void)weno5z_reconstruct_right(u,1e-40,0.0);}catch(const std::invalid_argument&){threw=true;}
    require(threw,"power=0 must be rejected");
}
}
int main(){
    try{exactness();smooth_order();bounded();validation();
        std::cout<<"N16_WENO5Z_VERIFICATION: PASS\\n"; return 0;}
    catch(const std::exception& e){std::cerr<<"N16_WENO5Z_VERIFICATION: FAIL: "<<e.what()<<"\\n"; return 1;}
}
