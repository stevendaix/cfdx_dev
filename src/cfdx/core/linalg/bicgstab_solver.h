#pragma once
#include "sparse_matrix.h"
#include "vector.h"
#include "cg_solver.h"
#include "solver_workspace.h"
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace cfdx::core {

inline SolverResult solve_bicgstab(
    const SparseMatrix& A,const Vector& b,Vector& x,
    std::size_t max_iter=1000,double tolerance=1e-12)
{
    SolverResult result;
    if(b.size()!=A.n_rows()||x.size()!=A.n_cols()||A.n_rows()!=A.n_cols()){
        result.status=SolverStatus::NOT_APPLICABLE;return result;
    }
    const std::size_t n=A.n_rows();const double* Av=A.values_data();const auto* Ac=A.columns_data();const auto* Ar=A.row_offsets_data();
    KrylovWorkspace w;w.resize(n);
    std::vector<double> inv_diag(n,0.0);
    for(std::size_t i=0;i<n;++i){for(std::size_t k=Ar[i];k<Ar[i+1];++k)if(Ac[k]==static_cast<std::uint32_t>(i)){inv_diag[i]=1.0/Av[k];break;}if(!std::isfinite(inv_diag[i])){result.status=SolverStatus::NOT_APPLICABLE;return result;}}
    auto matvec=[&](const Vector& in,Vector& out){for(std::size_t i=0;i<n;++i){double s=0.0;for(std::size_t k=Ar[i];k<Ar[i+1];++k)s+=Av[k]*in(Ac[k]);out(i)=s;}};
    Vector Ax(n,0.0);matvec(x,Ax);
    for(std::size_t i=0;i<n;++i)w.r(i)=b(i)-Ax(i);
    for(std::size_t i=0;i<n;++i)w.r_hat(i)=w.r(i);
    const double bnorm=b.norm2(),tol=tolerance*std::max(bnorm,1e-15);double res=w.r.norm2();
    if(res<=tol){result.status=SolverStatus::CONVERGED;result.residual=res;result.residual_relative=bnorm>0?res/bnorm:0.0;return result;}
    double rho=1.0,alpha=1.0,omega=1.0;
    for(std::size_t i=0;i<n;++i){w.v(i)=0.0;w.p(i)=0.0;}
    for(std::size_t iter=1;iter<=max_iter;++iter){
        double rho_new=0.0;for(std::size_t i=0;i<n;++i)rho_new+=w.r_hat(i)*w.r(i);
        if(std::abs(rho_new)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter-1;return result;}
        const double beta=(rho!=0.0)?(rho_new/rho)*(alpha/omega):0.0;
        for(std::size_t i=0;i<n;++i)w.p(i)=w.r(i)+beta*(w.p(i)-omega*w.v(i));
        for(std::size_t i=0;i<n;++i)w.z(i)=w.p(i)*inv_diag[i];
        matvec(w.z,w.v);
        double rv=0.0;for(std::size_t i=0;i<n;++i)rv+=w.r_hat(i)*w.v(i);
        if(std::abs(rv)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
        alpha=rho_new/rv;
        for(std::size_t i=0;i<n;++i)w.s(i)=w.r(i)-alpha*w.v(i);
        const double snorm=w.s.norm2();
        if(snorm<=tol){for(std::size_t i=0;i<n;++i)x(i)+=alpha*w.z(i);result.status=SolverStatus::CONVERGED;result.iterations=iter;result.residual=snorm;result.residual_relative=bnorm>0?snorm/bnorm:0.0;return result;}
        for(std::size_t i=0;i<n;++i)w.z_s(i)=w.s(i)*inv_diag[i];
        matvec(w.z_s,w.t);
        double ts=0.0,tt=0.0;for(std::size_t i=0;i<n;++i){ts+=w.t(i)*w.s(i);tt+=w.t(i)*w.t(i);}
        if(tt<=1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
        omega=ts/tt;if(std::abs(omega)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
        for(std::size_t i=0;i<n;++i){x(i)+=alpha*w.z(i)+omega*w.z_s(i);w.r(i)=w.s(i)-omega*w.t(i);}
        res=w.r.norm2();
        if(res<=tol){result.status=SolverStatus::CONVERGED;result.iterations=iter;result.residual=res;result.residual_relative=bnorm>0?res/bnorm:0.0;return result;}
        rho=rho_new;
    }
    result.status=SolverStatus::MAX_ITER_REACHED;result.iterations=max_iter;result.residual=res;result.residual_relative=bnorm>0?res/bnorm:0.0;return result;
}
}