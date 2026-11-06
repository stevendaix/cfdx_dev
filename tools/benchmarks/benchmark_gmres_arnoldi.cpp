// Optional GMRES Arnoldi benchmark. Certifies the ORIGINAL unscaled true
// residual and holds the restart length and iteration budget fixed. Runs are
// not CI acceptance gates: compare medians on identical hardware/build flags.
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
 using namespace cfdx::core;
 std::size_t n = argc>1 ? static_cast<std::size_t>(std::stoul(argv[1])) : 2304;
 SparseMatrix A(n,n);
 for (std::size_t i=0; i<n; ++i) {
    A.push_back(i,i,2.06);
    if (i) A.push_back(i,i-1,-1.0);
    if (i+1<n) A.push_back(i,i+1,-1.0);
 }
 A.finalize();
 Vector exact(n), b(n);
 for(std::size_t i=0;i<n;++i) exact(i)=std::sin(0.171*(i+1)) + 0.21*std::cos(0.031*(i+1));
 auto bb=A.matvec(exact);
 for(std::size_t i=0;i<n;++i) b(i)=bb[i];
 KrylovControls c;
 c.adaptive_restart=false; c.restart_min=50; c.restart_max=50; c.residual_replacement=false;
 LinearOperator op;
 op.size=n;
 op.apply=[&A](const Vector& in, Vector& out) {
     const auto ax=A.matvec(in);
     if(out.size()!=ax.size())out.resize(ax.size());
     for(std::size_t i=0;i<ax.size();++i)out(i)=ax[i];
 };
 GmresWorkspace space;
 std::vector<double> ms;
 for(int k=0; k<9;++k) {
   Vector x(n,0.0);
   const auto t0=std::chrono::steady_clock::now();
   const auto r=solve_gmres(op,b,x,50,500,1e-12,nullptr,c,&space);
   const auto t1=std::chrono::steady_clock::now();
   const auto elapsed=std::chrono::duration<double,std::milli>(t1-t0).count();
   const auto ax=A.matvec(x);
   double norm=0.0, scale=0.0;
   for(std::size_t i=0;i<n;++i) {norm=std::hypot(norm,ax[i]-b(i)); scale=std::hypot(scale,b(i));}
   if (r.status != SolverStatus::CONVERGED || norm/scale > 1.01e-12) {
      std::cerr << "fail status=" << int(r.status) << " iterations=" << r.iterations << " true_res=" << norm/scale << "\n";
      return 2;
   }
   if(k>0) ms.push_back(elapsed);
   if(k==0) std::cout << "n="<<n<<" iters="<<r.iterations<<" res="<<norm/scale<<"\n";
 }
 std::sort(ms.begin(),ms.end());
 std::cout << "median_ms=" << 0.5*(ms[3]+ms[4]) << " min=" << ms[0] << " max=" << ms[7] << "\n";
}
