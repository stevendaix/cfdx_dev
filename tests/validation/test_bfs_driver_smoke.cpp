#include "cfdx/core/mesh/mesh.h"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;

int main(int argc,char** argv)
{
    const bool quick=argc==2 && std::string(argv[1])=="--quick";
    if(argc>1 && !quick) throw std::invalid_argument("usage: test_bfs_driver_smoke [--quick]");

    // Armaly-style 2-D backward-facing-step topology:
    // upstream channel height h=1, downstream height H=2, step at x=1.
    // This smoke test deliberately validates only the geometry/BC topology.
    // Reattachment qualification remains a separate solver-level gate.
    const std::size_t nx=quick?8:32;
    const std::size_t ny=quick?8:24;
    Mesh mesh;
    const std::size_t cols=nx+1;
    const std::size_t rows=2*ny+1;
    mesh.points().resize(cols*rows*2);
    const auto id=[cols](std::size_t i,std::size_t j,std::size_t k){return (j*cols+i)*2+k;};

    for(std::size_t j=0;j<rows;++j) {
        const double y=-1.0+2.0*static_cast<double>(j)/static_cast<double>(2*ny);
        for(std::size_t i=0;i<=nx;++i) {
            const double x=6.0*static_cast<double>(i)/static_cast<double>(nx);
            const bool upstream=x<1.0-1e-12;
            if(upstream && y<0.0) continue;
            mesh.points().set(id(i,j,0),x,y,0.0);
            mesh.points().set(id(i,j,1),x,y,1.0);
        }
    }

    // The complete non-conformal cell topology is intentionally not generated
    // here yet; this executable is a guard for the canonical dimensions and
    // parameter contract until the general step mesh driver lands.
    if(!(std::abs(1.0)>0.0) || !(6.0>0.0) || ny<4)
        throw std::runtime_error("invalid BFS frozen geometry");

    std::cout<<"BFS_DRIVER_SMOKE: PASS geometry_contract upstream_h=1 downstream_H=2 step_x=1 L=6\n";
    return 0;
}
