#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

struct Point { double x; double y; };

static std::vector<Point> naca0012(std::size_t n)
{
    if(n<20) throw std::invalid_argument("NACA0012 requires at least 20 panels");
    constexpr double t=0.12;
    std::vector<Point> p;
    p.reserve(2*n+1);
    for(std::size_t i=0;i<=n;++i) {
        const double beta=M_PI*static_cast<double>(i)/static_cast<double>(n);
        const double x=0.5*(1.0-std::cos(beta));
        const double yt=5.0*t*(0.2969*std::sqrt(x)-0.1260*x-0.3516*x*x+
                                0.2843*x*x*x-0.1015*x*x*x*x);
        p.push_back({x,yt});
    }
    for(std::size_t i=n;i>0;--i) {
        const double beta=M_PI*static_cast<double>(i)/static_cast<double>(n);
        const double x=0.5*(1.0-std::cos(beta));
        const double yt=5.0*t*(0.2969*std::sqrt(x)-0.1260*x-0.3516*x*x+
                                0.2843*x*x*x-0.1015*x*x*x*x);
        p.push_back({x,-yt});
    }
    return p;
}

int main(int argc,char** argv)
{
    const bool quick=argc==2 && std::string(argv[1])=="--quick";
    if(argc>1 && !quick) throw std::invalid_argument("usage: test_naca0012_geometry [--quick]");
    const auto p=naca0012(quick?40:160);
    if(p.front().x>1e-12 || p.back().x>1e-12)
        throw std::runtime_error("NACA0012 trailing closure contract failed");
    double max_t=0.0;
    for(const auto& q:p) max_t=std::max(max_t,std::abs(q.y));
    if(std::abs(max_t-0.0600)>0.001)
        throw std::runtime_error("NACA0012 maximum thickness contract failed");
    std::cout<<"NACA0012_GEOMETRY: PASS points="<<p.size()
             <<" max_half_thickness="<<max_t<<"\n";
    std::cout<<"NACA0012 solver qualification remains BLOCKED: external-flow force integration is not yet available.\n";
    return 0;
}
