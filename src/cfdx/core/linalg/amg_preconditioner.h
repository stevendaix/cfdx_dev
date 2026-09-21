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
        if(op.rows()!=op.cols()||!(omega>0.0&&omega<2.0))throw std::invalid_argument("invalid V-cycle operator");
        inv_diag_.assign(op.rows(),1.0);
    }
    bool setup(const SparseMatrix&)override{return true;}
    bool apply(const Vector& r,Vector& z)const override{
        if(r.size()!=op_.rows())return false;if(z.size()!=r.size())z.resize(r.size(),0.0);std::fill(z.data(),z.data()+z.size(),0.0);
        smooth(r,z,pre_); Vector Az(r.size(),0.0),rc(r.size(),0.0);
        op_.apply(z,Az);for(std::size_t i=0;i<r.size();++i)rc(i)=r(i)-Az(i);
        const std::size_t nc=(r.size()+1)/2;Vector ec(nc,0.0);
        for(std::size_t i=0;i<r.size();++i)ec(i/2)+=rc(i);
        for(std::size_t i=0;i<nc;++i)ec(i)/=((2*i+1<r.size())?2.0:1.0);
        for(std::size_t i=0;i<r.size();++i)z(i)+=ec(i/2);smooth(r,z,post_);return true;
    }
    const char* name()const override{return "matrix-free-vcycle";}
private:
    void smooth(const Vector& r,Vector& x,std::size_t sweeps)const{Vector Ax(r.size(),0.0);for(std::size_t s=0;s<sweeps;++s){op_.apply(x,Ax);for(std::size_t i=0;i<r.size();++i)x(i)+=omega_*(r(i)-Ax(i))*inv_diag_[i];}}
    const LinearOperatorBase& op_;double omega_;std::size_t pre_,post_;std::vector<double> inv_diag_;
};
}