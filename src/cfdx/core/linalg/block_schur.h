#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include <cstddef>
#include <functional>
#include <utility>
namespace cfdx::core {
class BlockSchurPreconditioner {
public:
 using Solve=std::function<bool(const Vector&,Vector&)>;
 BlockSchurPreconditioner(const LinearOperatorBase& G,const LinearOperatorBase& D,Solve solve_momentum,Solve solve_schur)
 :G_(G),D_(D),solve_momentum_(std::move(solve_momentum)),solve_schur_(std::move(solve_schur)){}
 bool apply(const Vector& rhs_u,const Vector& rhs_p,Vector& u,Vector& p)const{
  if(!solve_momentum_||!solve_schur_||rhs_u.size()!=G_.cols()||rhs_p.size()!=D_.rows())return false;
  if(u.size()!=rhs_u.size())u.resize(rhs_u.size());if(p.size()!=rhs_p.size())p.resize(rhs_p.size());
  if(!solve_momentum_(rhs_u,u))return false;
  Vector Du(rhs_p.size());D_.apply(u,Du);
  Vector schur_rhs(rhs_p.size());for(std::size_t i=0;i<p.size();++i)schur_rhs(i)=rhs_p(i)-Du(i);
  if(!solve_schur_(schur_rhs,p))return false;
  Vector Gp(rhs_u.size());G_.apply(p,Gp);Vector correction(rhs_u.size());if(!solve_momentum_(Gp,correction))return false;
  for(std::size_t i=0;i<u.size();++i)u(i)-=correction(i);
  return true;
 }
private:const LinearOperatorBase&G_;const LinearOperatorBase&D_;Solve solve_momentum_,solve_schur_;
};
} // namespace cfdx::core
