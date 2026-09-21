#pragma once
#include "cfdx/core/linalg/linear_operator.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
namespace cfdx::core {
class ChebyshevSmoother {
public:
 struct Controls{std::size_t degree=2;double lambda_min=0.05,lambda_max=1.0,damping=1.0;};
 explicit ChebyshevSmoother(Controls c={}):c_(c){if(c_.degree==0||!(c_.lambda_min>0.0)||!(c_.lambda_max>=c_.lambda_min)||!(c_.damping>0.0))throw std::invalid_argument("invalid Chebyshev controls");}
 void apply(const LinearOperatorBase& A,const Vector& rhs,Vector& x)const{
  if(rhs.size()!=A.rows()||A.rows()!=A.cols())throw std::invalid_argument("Chebyshev dimension mismatch");
  if(x.size()!=rhs.size()){x.resize(rhs.size());x.fill(0.0);}
  const double theta=0.5*(c_.lambda_max+c_.lambda_min),delta=0.5*(c_.lambda_max-c_.lambda_min);
  double alpha=1.0/theta,beta=0.0;Vector r(rhs.size()),Ar(rhs.size()),xprev=x;
  for(std::size_t k=0;k<c_.degree;++k){
   A.apply(x,Ar);for(std::size_t i=0;i<rhs.size();++i)r(i)=rhs(i)-Ar(i);
   if(k==0){alpha=1.0/theta;beta=0.0;}else{beta=std::pow(0.5*delta*alpha,2.0);alpha=1.0/(theta-beta*delta);}
   for(std::size_t i=0;i<rhs.size();++i){const double old=x(i);x(i)+=c_.damping*alpha*r(i)+beta*(x(i)-xprev(i));xprev(i)=old;}
  }
 }
private:Controls c_;
};
} // namespace cfdx::core
