#include "cfdx/core/linalg/communication_avoiding.h"
#include "cfdx/core/parallel/mpi_utils.h"
#include <cmath>
#include <iostream>
int main(int argc,char** argv){
 cfdx::core::parallel::mpi_init(&argc,&argv);
 const int rank=cfdx::core::parallel::mpi_rank(), size=cfdx::core::parallel::mpi_size();
 cfdx::core::ReductionPacket local{double(rank+1),double((rank+1)*(rank+1)),double(2*rank+1)};
 const auto g=cfdx::core::mpi_deterministic_fused_reduction(local);
 const bool ok=std::abs(g.dot-double(size*(size+1))/2.0)<1e-14 &&
   std::abs(g.norm2-double(size*(size+1)*(2*size+1))/6.0)<1e-14 &&
   std::abs(g.max_abs-double(2*size-1))<1e-14;
 if(rank==0&&!ok) std::cerr<<"deterministic fused reduction mismatch\n";
 cfdx::core::parallel::mpi_finalize(); return ok?0:1;
}