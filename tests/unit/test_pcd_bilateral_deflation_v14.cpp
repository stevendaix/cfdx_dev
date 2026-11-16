#include "cfdx/physics/pcd_fp_bisided_basis_v13.h"
#include "cfdx/physics/pcd_fp_deflated_quality_v14.h"
#include "cfdx/physics/pcd_fp_capture.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/advanced_preconditioners.h"
#include "cfdx/physics/pcd_fp_strong_component_anchor.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::physics;
using namespace cfdx::core;
static void require(bool v,const char* m) {if(!v) throw std::runtime_error(m);}
int main(int argc,char** argv) {
 try {
  SparseMatrix A(2,2);
  A.push_back(0,0,1);A.push_back(0,1,-1);
  A.push_back(1,0,2);A.push_back(1,1,-2);A.finalize();
  Vector w(2),z(2),b(2),x(2,0.0);
  w(0)=-2/std::sqrt(5.);w(1)=1/std::sqrt(5.);
  z(0)=1/std::sqrt(2.);z(1)=1/std::sqrt(2.);
  b(0)=1.;b(1)=0.;
  auto basis=PcdFpBiNullspaceV13::certify(A,{w},{z});
  auto rhs=basis.compatible_rhs(b);
  LinearOperator op;
  op.size=2;
  op.apply=[&](const Vector& in,Vector& out){out=basis.deflated_matvec(A,in);};
  op.apply_fp64=op.apply;
  auto result=solve_fgmres(op,rhs,x,2,30,1e-12,nullptr);
  require(result.status==SolverStatus::CONVERGED,"deflated nonsymmetric solve failed");
  const auto x_c=basis.canonical_solution(x);
  const auto independent_deflated=pcd_fp_deflated_true_residual_v14(A,basis,rhs,Vector(2,0.0),x_c);
  require(independent_deflated.relative_rhs<1e-12,"independent deflated residual not closed");
  require(basis.projected_residual_relative(A,b,x_c)<1e-12,"true projected residual not closed");
  Vector physical(2,0.0);
  A.matvec_into(x_c,physical);
  const double original=std::hypot(physical(0)-b(0),physical(1)-b(1));
  require(std::abs(original-basis.incompatibility_relative(b))<1e-11,
    "original incompatible RHS floor hidden");
  SparseMatrix changed(A);
  changed.values_data()[0]+=1e-5;
  bool stale=false;
  try {basis.deflated_matvec(changed,x);}catch(const std::invalid_argument&){stale=true;}
  require(stale,"stale Fp operator accepted");
  std::cout << "FP_V14_TOY projected=" << basis.projected_residual_relative(A,b,x_c)
            << " original=" << original << " floor=" << basis.incompatibility_relative(b)
            << " matched=" << (std::abs(original-basis.incompatibility_relative(b))<1e-11)
            << "\n";
  if (argc==3) {
    const auto captured=read_pcd_fp_capture(argv[1]);
    const auto cert=load_pcd_fp_bisided_basis_v13(argv[2],captured.matrix);
    const auto bcomp=cert.compatible_rhs(captured.rhs);
    const auto floor=cert.incompatibility_relative(captured.rhs);
    require(floor>0.01 && floor<0.2,"unexpected captured floor");
    LinearOperator fp;
    fp.size=captured.matrix.n_rows();
    fp.apply=[&](const Vector& in, Vector& out) {
      out=cert.deflated_matvec(captured.matrix,in);
    };
    fp.apply_fp64=fp.apply;
    Vector xx(captured.matrix.n_rows(),0.0);
    // This ILU0 approximates the AUXILIARY anchored operator as in the
    // real production PCD preconditioner, but Arnoldi evaluates the EXACT
    // deflated original Fp on every matvec.
    auto anchored=anchor_pcd_fp_strong_components(captured.matrix);
    ILU0Preconditioner ilu;
    require(ilu.setup(anchored.matrix),"anchored ILU0 setup failed");
    const auto run=solve_fgmres(fp,bcomp,xx,80,200,1e-8,&ilu);
    const auto independent=pcd_fp_deflated_true_residual_v14(
        captured.matrix,cert,bcomp,Vector(fp.size,0.0),xx);
    require(independent.relative_rhs<1e-8,"captured deflated true residual not closed");
    std::cout << "FP_V14_CAPTURE n=" << fp.size << " modes=" << cert.left.size()
              << " floor=" << floor
              << " deflated_true=" << independent.relative_rhs
              << " projected_true=" << cert.projected_residual_relative(
                    captured.matrix,captured.rhs,xx)
              << " iterations=" << run.iterations
              << " status=" << static_cast<int>(run.status)
              << " physical_qualified=0\n";
  } else if (argc!=1) throw std::invalid_argument("usage: [capture basis]");
  return 0;
 } catch(const std::exception& e) {std::cerr<<"FP_V14_FAIL: "<<e.what()<<"\n";return 1;}
}
