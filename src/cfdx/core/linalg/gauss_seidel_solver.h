#pragma once
#include "sparse_matrix.h"
#include "vector.h"
#include "cg_solver.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
namespace cfdx::core {
inline SolverResult solve_gauss_seidel(const SparseMatrix& A,const Vector& b,Vector& x,std::size_t max_iter=1000,double tolerance=1e-12){
 SolverResult r;
 if(A.n_rows()!=A.n_cols()||b.size()!=A.n_rows()||x.size()!=A.n_cols()||max_iter==0||!(tolerance>0)||!std::isfinite(tolerance)) return r;
 const std::size_t n=A.n_rows(); const double* Av=A.values_data(); const auto* Ac=A.columns_data(); const auto* Ar=A.row_offsets_data();
 std::vector<double> d(n); for(std::size_t k=0;k<A.nnz();++k) if(!std::isfinite(Av[k])){r.status=SolverStatus::DIVERGED;return r;}
 for(std::size_t i=0;i<n;++i){bool f=false;for(std::size_t k=Ar[i];k<Ar[i+1];++k)if(Ac[k]==static_cast<std::uint32_t>(i)){d[i]=Av[k];f=true;break;}if(!f||d[i]==0||!std::isfinite(d[i])){r.status=SolverStatus::NOT_APPLICABLE;return r;}}
 for(std::size_t i=0;i<n;++i)if(!std::isfinite(b(i))||!std::isfinite(x(i))){r.status=SolverStatus::DIVERGED;return r;}
 const double bn=b.norm2(),tol=tolerance*std::max(bn,1e-15); Vector next(n,0.0);
 auto residual=[&](){double s=0;for(std::size_t i=0;i<n;++i){double q=0;for(std::size_t k=Ar[i];k<Ar[i+1];++k)q+=Av[k]*x(Ac[k]);double v=q-b(i);if(!std::isfinite(v))return std::numeric_limits<double>::infinity();s+=v*v;}return std::sqrt(s);};
 double res=residual(); if(res<=tol){r.status=SolverStatus::CONVERGED;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}
 for(std::size_t it=1;it<=max_iter;++it){for(std::size_t i=0;i<n;++i){double rhs=b(i);for(std::size_t k=Ar[i];k<Ar[i+1];++k)if(Ac[k]!=static_cast<std::uint32_t>(i))rhs-=Av[k]*x(Ac[k]);next(i)=rhs/d[i];if(!std::isfinite(next(i))){r.status=SolverStatus::DIVERGED;r.iterations=it;return r;}}for(std::size_t i=0;i<n;++i)x(i)=next(i);res=residual();if(!std::isfinite(res)){r.status=SolverStatus::DIVERGED;r.iterations=it;return r;}if(res<=tol){r.status=SolverStatus::CONVERGED;r.iterations=it;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}}
 r.status=SolverStatus::MAX_ITER_REACHED;r.iterations=max_iter;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;
}}
