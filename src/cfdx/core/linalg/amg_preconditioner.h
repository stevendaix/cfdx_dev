#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include "cfdx/core/linalg/preconditioner.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
namespace cfdx::core {
class MatrixFreeVcyclePreconditioner final:public Preconditioner{
public:
 MatrixFreeVcyclePreconditioner(const LinearOperatorBase& op,double omega=0.7,std::size_t pre=2,std::size_t post=2):op_(op),omega_(omega),pre_(pre),post_(post){
  if(op.rows()!=op.cols()||!(omega_>0.0&&omega_<2.0))throw std::invalid_argument("invalid V-cycle operator");setup_diagonal();
 }
 bool setup(const SparseMatrix& A)override{
  if(A.n_rows()!=op_.rows()||A.n_cols()!=op_.cols())return false;
  if(!op_.has_diagonal()){
   const auto* row=A.row_offsets_data();const auto* col=A.columns_data();const auto* val=A.values_data();inv_diag_.assign(A.n_rows(),0.0);
   for(std::size_t i=0;i<A.n_rows();++i){for(std::size_t k=row[i];k<row[i+1];++k)if(col[k]==i&&std::abs(val[k])>1e-30){inv_diag_[i]=1.0/val[k];break;}if(inv_diag_[i]==0.0)inv_diag_[i]=1.0;}
   build_coarse_diagonal();
  } return true;
 }
 bool apply(const Vector& r,Vector& z)const override{
  if(r.size()!=op_.rows())return false;if(z.size()!=r.size())z.resize(r.size());z.fill(0.0);smooth(r,z,pre_);
  op_.apply(z,ws_.fine_A);for(std::size_t i=0;i<r.size();++i)ws_.fine_r(i)=r(i)-ws_.fine_A(i);
  ws_.coarse_r.fill(0.0);for(std::size_t i=0;i<r.size();++i)ws_.coarse_r(i)+=0.5*ws_.fine_r(i);
  for(std::size_t i=0;i<ws_.coarse_x.size();++i)ws_.coarse_x(i)=ws_.coarse_r(i)*ws_.coarse_inv_diag[i];
  for(std::size_t i=0;i<r.size();++i)z(i)+=ws_.coarse_x(i/2);smooth(r,z,post_);return true;
 }
 const char* name()const override{return "matrix-free-agglomerated-vcycle";}
private:
 struct Workspace{Vector fine_r,fine_A,coarse_r,coarse_x;std::vector<double> coarse_inv_diag;};
 void setup_diagonal(){Vector d(op_.rows(),1.0);if(op_.has_diagonal())op_.diagonal(d);inv_diag_.resize(d.size());for(std::size_t i=0;i<d.size();++i)inv_diag_[i]=std::abs(d(i))>1e-30?1.0/d(i):1.0;build_coarse_diagonal();}
 void build_coarse_diagonal(){const std::size_t nc=(inv_diag_.size()+1)/2;ws_.coarse_r.resize(nc);ws_.coarse_r.fill(0.0);ws_.coarse_x.resize(nc);ws_.coarse_x.fill(0.0);ws_.coarse_inv_diag.assign(nc,1.0);for(std::size_t c=0;c<nc;++c){const std::size_t i0=2*c;double d=1.0/std::max(inv_diag_[i0],1e-30);if(i0+1<inv_diag_.size())d+=1.0/std::max(inv_diag_[i0+1],1e-30);ws_.coarse_inv_diag[c]=1.0/std::max(d,1e-30);}ws_.fine_r.resize(inv_diag_.size());ws_.fine_r.fill(0.0);ws_.fine_A.resize(inv_diag_.size());ws_.fine_A.fill(0.0);}
 void smooth(const Vector& r,Vector& x,std::size_t sweeps)const{for(std::size_t s=0;s<sweeps;++s){op_.apply(x,ws_.fine_A);for(std::size_t i=0;i<r.size();++i)x(i)+=omega_*inv_diag_[i]*(r(i)-ws_.fine_A(i));}}
 const LinearOperatorBase& op_;double omega_;std::size_t pre_,post_;std::vector<double> inv_diag_;mutable Workspace ws_;
};
using AgglomeratedAMGPreconditioner=MatrixFreeVcyclePreconditioner;
} // namespace cfdx::core
