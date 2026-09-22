#include "cfdx/core/linalg/jacobi_solver.h"
#include "cfdx/core/linalg/gauss_seidel_solver.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"
#include <cmath>
using namespace cfdx::core; using namespace cfdx::testing;

SparseMatrix make_reference() {
 SparseMatrix A(3,3);
 A.push_back(0,0,4); A.push_back(0,1,1); A.push_back(0,2,1);
 A.push_back(1,0,2); A.push_back(1,1,5); A.push_back(1,2,1);
 A.push_back(2,0,1); A.push_back(2,1,1); A.push_back(2,2,3);
 A.finalize(); return A;
}
Vector make_b() { Vector b(3); b(0)=9; b(1)=15; b(2)=12; return b; }

int main() {
 run_case("jacobi_reference_solution", [] {
  auto A=make_reference(); auto b=make_b(); Vector x(3,0);
  auto r=solve_jacobi(A,b,x,1000,1e-12);
  EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
  EXPECT_TRUE(r.residual_relative<1e-11);
  EXPECT_NEAR(x(0),1,1e-10); EXPECT_NEAR(x(1),2,1e-10); EXPECT_NEAR(x(2),3,1e-10);
  EXPECT_TRUE(std::isfinite(r.residual_relative));
 });
 run_case("gauss_seidel_reference_solution", [] {
  auto A=make_reference(); auto b=make_b(); Vector x(3,0);
  auto r=solve_gauss_seidel(A,b,x,1000,1e-12);
  EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
  EXPECT_TRUE(r.residual_relative<1e-11);
  EXPECT_NEAR(x(0),1,1e-10); EXPECT_NEAR(x(1),2,1e-10); EXPECT_NEAR(x(2),3,1e-10);
  EXPECT_TRUE(std::isfinite(r.residual_relative));
 });
 run_case("gauss_seidel_uses_latest_values", [] {
  auto A=make_reference(); auto b=make_b(); Vector xj(3,0), xg(3,0);
  auto rj=solve_jacobi(A,b,xj,1000,1e-10); auto rg=solve_gauss_seidel(A,b,xg,1000,1e-10);
  EXPECT_TRUE(rj.status==SolverStatus::CONVERGED); EXPECT_TRUE(rg.status==SolverStatus::CONVERGED);
  EXPECT_TRUE(rg.iterations < rj.iterations);
 });
 run_case("stationary_solvers_reject_missing_diagonal", [] {
  SparseMatrix A(2,2); A.push_back(0,1,1); A.push_back(1,0,1); A.finalize();
  Vector b(2,1), x(2,0);
  EXPECT_TRUE(solve_jacobi(A,b,x).status==SolverStatus::NOT_APPLICABLE);
  EXPECT_TRUE(solve_gauss_seidel(A,b,x).status==SolverStatus::NOT_APPLICABLE);
 });
 run_case("stationary_solvers_report_max_iter_without_nan", [] {
  SparseMatrix A(2,2); A.push_back(0,0,1); A.push_back(0,1,2); A.push_back(1,0,2); A.push_back(1,1,1); A.finalize();
  Vector b(2,1), xj(2,0), xg(2,0);
  auto rj=solve_jacobi(A,b,xj,5,1e-14); auto rg=solve_gauss_seidel(A,b,xg,5,1e-14);
  EXPECT_TRUE(rj.status==SolverStatus::MAX_ITER_REACHED || rj.status==SolverStatus::DIVERGED);
  EXPECT_TRUE(rg.status==SolverStatus::MAX_ITER_REACHED || rg.status==SolverStatus::DIVERGED);
  for(std::size_t i=0;i<2;++i){EXPECT_TRUE(std::isfinite(xj(i))); EXPECT_TRUE(std::isfinite(xg(i)));}
 });
 run_case("stationary_solvers_validate_controls", [] {
  auto A=make_reference(); auto b=make_b(); Vector x(3,0);
  EXPECT_TRUE(solve_jacobi(A,b,x,0,1e-12).status==SolverStatus::NOT_APPLICABLE);
  EXPECT_TRUE(solve_gauss_seidel(A,b,x,10,0).status==SolverStatus::NOT_APPLICABLE);
 });
 return run_all();
}