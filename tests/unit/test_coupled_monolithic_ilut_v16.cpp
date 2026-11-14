#include "cfdx/physics/coupled_monolithic_ilut_v16.h"
#include "cfdx/physics/coupled_block_residual_v16.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <limits>
#define REQUIRE(p) do {if (!(p)) throw std::runtime_error("V16 ILUT failed: " #p);} while(false)
int main() {
    using namespace cfdx::core;
    SparseMatrix A(4,4);
    A.push_back(0,0,2);A.push_back(0,3,1);
    A.push_back(1,1,3);A.push_back(1,3,1);
    A.push_back(2,2,4);A.push_back(2,3,-1);
    A.push_back(3,0,1);A.push_back(3,1,1);A.push_back(3,2,-1);
    A.finalize();
    Vector b(4,0.),x(4,0.);b(0)=1;b(1)=3;b(2)=2;b(3)=0;
    cfdx::physics::CoupledMonolithicIlutV16 pre;
    REQUIRE(pre.setup(A));
    const auto r=solve_fgmres(A,b,x,4,100,1e-12,&pre);
    REQUIRE(r.status==SolverStatus::CONVERGED);
    const auto cert=cfdx::physics::independent_coupled_block_residual_v16(A,b,x,3);
    REQUIRE(cert.relative_to_global_rhs()<1e-12);
    Vector small(3,0.);
    REQUIRE(!pre.apply(small,small));
    SparseMatrix too_large(4097,4097);
    for(size_t i=0;i<4097;++i) too_large.push_back(i,i,1.);
    too_large.finalize();
    REQUIRE(!pre.setup(too_large));
    std::cout<<"V16 native ILUT converges on unmodified saddle system; >4096 rejected\n";
}
