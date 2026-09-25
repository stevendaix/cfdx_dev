#pragma once
#include "sparse_matrix.h"
#include "vector.h"
#include "cg_solver.h"
#include "solver_workspace.h"
#include "preconditioner.h"
#include "mixed_precision.h"
#include "krylov_reductions.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace cfdx::core {

inline SolverResult solve_bicgstab(
    const SparseMatrix& A,const Vector& b,Vector& x,
    std::size_t max_iter=1000,double tolerance=1e-12,
    Preconditioner* preconditioner=nullptr, PrecisionPolicy precision={},
    KrylovReductionPolicy reduction={})
{
    SolverResult result;
    if(b.size()!=A.n_rows()||x.size()!=A.n_cols()||A.n_rows()!=A.n_cols()||
       max_iter==0||!std::isfinite(tolerance)||tolerance<=0.0){
        result.status=SolverStatus::NOT_APPLICABLE;return result;
    }
    if(preconditioner && !preconditioner->setup(A)){
        result.status=SolverStatus::NOT_APPLICABLE;return result;
    }

    const std::size_t n=A.n_rows();
    const bool mp32 = precision.enabled && precision.operator_precision == SolverPrecision::FP32;
    const SolverPrecision redp = precision.enabled ? precision.reduction_precision : SolverPrecision::FP64;
    const double* Av=A.values_data();const auto* Ac=A.columns_data();const auto* Ar=A.row_offsets_data();
    auto finite_vector = [](const Vector& v) {
        for (std::size_t i = 0; i < v.size(); ++i)
            if (!std::isfinite(v(i))) return false;
        return true;
    };
    for(std::size_t k=0;k<A.nnz();++k)
        if(!std::isfinite(Av[k])){result.status=SolverStatus::DIVERGED;return result;}
    // Without a preconditioner, retain the established applicability
    // requirement for an explicit finite diagonal. A block preconditioner
    // may provide the required scaling itself, so do not reject saddle-point
    // matrices whose continuity rows have a structural zero diagonal.
    if (!preconditioner) {
        for (std::size_t i = 0; i < n; ++i) {
            bool found_diagonal = false;
            for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k) {
                if (Ac[k] == i) {
                    if (!std::isfinite(Av[k]) || Av[k] == 0.0) {
                        result.status = SolverStatus::NOT_APPLICABLE;
                        return result;
                    }
                    found_diagonal = true;
                    break;
                }
            }
            if (!found_diagonal) {
                result.status = SolverStatus::NOT_APPLICABLE;
                return result;
            }
        }
    }
    if(!finite_vector(b)||!finite_vector(x)){result.status=SolverStatus::DIVERGED;return result;}

    auto matvec=[&](const Vector& in,Vector& out){
        if (mp32) { mixed_precision_matvec(A, in, out, SolverPrecision::FP32); return; }
        if(out.size()!=n) out.resize(n);
        for(std::size_t i=0;i<n;++i){
            double s=0.0;
            for(std::size_t k=Ar[i];k<Ar[i+1];++k)s+=Av[k]*in(Ac[k]);
            out(i)=s;
        }
    };

    KrylovWorkspace w;w.resize(n);
    Vector Ax(n,0.0);
    matvec(x,Ax);
    if(!finite_vector(Ax)){result.status=SolverStatus::DIVERGED;return result;}
    for(std::size_t i=0;i<n;++i)w.r(i)=b(i)-Ax(i);
    if(!finite_vector(w.r)){result.status=SolverStatus::DIVERGED;return result;}
    for(std::size_t i=0;i<n;++i)w.r_hat(i)=w.r(i);

    const double bnorm=krylov_norm2(b, SolverPrecision::FP64, reduction),tol=tolerance*std::max(bnorm,1e-15);
    if(!std::isfinite(bnorm)||!std::isfinite(tol)){result.status=SolverStatus::DIVERGED;return result;}
    Vector true_r(n); mixed_precision_true_residual(A,b,x,true_r); double res=krylov_norm2(true_r, redp, reduction);
    if(res<=tol){
        result.status=SolverStatus::CONVERGED;result.residual=res;
        result.residual_relative=bnorm>0?res/bnorm:0.0;return result;
    }

    double rho=1.0,alpha=1.0,omega=1.0;
    for(std::size_t i=0;i<n;++i){w.v(i)=0.0;w.p(i)=0.0;}

    auto true_residual = [&]() {
        matvec(x, Ax);
        if (!finite_vector(Ax)) return std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < n; ++i) {
            const double ri = b(i) - Ax(i);
            if (!std::isfinite(ri)) return std::numeric_limits<double>::infinity();
            w.r(i) = ri;
        }
        return krylov_norm2(w.r, redp, reduction);
    };

    for(std::size_t iter=1;iter<=max_iter;++iter){
        const double rho_new=krylov_dot(w.r_hat, w.r, redp, reduction);
        if(!std::isfinite(rho_new)||std::abs(rho_new)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter-1;return result;}
        const double beta=(rho!=0.0)?(rho_new/rho)*(alpha/omega):0.0;
        if(!std::isfinite(beta)){result.status=SolverStatus::DIVERGED;result.iterations=iter-1;return result;}
        for(std::size_t i=0;i<n;++i)w.p(i)=w.r(i)+beta*(w.p(i)-omega*w.v(i));

        if(preconditioner){
            if(!preconditioner->apply(w.p,w.z)){result.status=SolverStatus::NOT_APPLICABLE;result.iterations=iter-1;return result;}
        }else{
            for(std::size_t i=0;i<n;++i)w.z(i)=w.p(i);
        }
        matvec(w.z,w.v);
        if(!finite_vector(w.v)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}

        const double rv=krylov_dot(w.r_hat, w.v, redp, reduction);
        if(!std::isfinite(rv)||std::abs(rv)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
        alpha=rho_new/rv;
        if(!std::isfinite(alpha)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
        for(std::size_t i=0;i<n;++i)w.s(i)=w.r(i)-alpha*w.v(i);
        if(!finite_vector(w.s)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}

        const double snorm=krylov_norm2(w.s, redp, reduction);
        if(snorm<=tol){
            for(std::size_t i=0;i<n;++i)x(i)+=alpha*w.z(i);
            res=true_residual();
            if(res<=tol){
                result.status=SolverStatus::CONVERGED;result.iterations=iter;result.residual=res;
                result.residual_relative=bnorm>0?res/bnorm:0.0;return result;
            }
            if (!std::isfinite(res)) {
                result.status=SolverStatus::DIVERGED;
                result.iterations=iter;
                return result;
            }
            // The recurrence residual s can be below tolerance while the
            // recomputed true residual is not. Continue from the recomputed
            // residual; never mix the two residual states.
            for (std::size_t i = 0; i < n; ++i) {
                w.r_hat(i) = w.r(i);
                w.p(i) = 0.0;
                w.v(i) = 0.0;
            }
            rho = 1.0;
            alpha = 1.0;
            omega = 1.0;
            continue;
        }

        if(preconditioner){
            if(!preconditioner->apply(w.s,w.z_s)){result.status=SolverStatus::NOT_APPLICABLE;result.iterations=iter;return result;}
        }else{
            for(std::size_t i=0;i<n;++i)w.z_s(i)=w.s(i);
        }
        matvec(w.z_s,w.t);
        if(!finite_vector(w.t)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}

        const double ts=krylov_dot(w.t, w.s, redp, reduction);
        const double tt=krylov_dot(w.t, w.t, redp, reduction);
        if(!std::isfinite(ts)||!std::isfinite(tt)||tt<=1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
        omega=ts/tt;if(!std::isfinite(omega)||std::abs(omega)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
        for(std::size_t i=0;i<n;++i){x(i)+=alpha*w.z(i)+omega*w.z_s(i);w.r(i)=w.s(i)-omega*w.t(i);}
        if(!finite_vector(x)||!finite_vector(w.r)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}

        res=true_residual();
        if(res<=tol){
            result.status=SolverStatus::CONVERGED;result.iterations=iter;result.residual=res;
            result.residual_relative=bnorm>0?res/bnorm:0.0;return result;
        }
        rho=rho_new;
    }
    result.status=SolverStatus::MAX_ITER_REACHED;result.iterations=max_iter;result.residual=res;
    result.residual_relative=bnorm>0?res/bnorm:0.0;return result;
}
}