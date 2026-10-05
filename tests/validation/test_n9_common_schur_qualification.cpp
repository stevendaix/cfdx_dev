#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/mgr_preconditioner.h"
#include "cfdx/core/linalg/simplerc_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
using Dense = std::vector<std::vector<double>>;

SparseMatrix make_sparse(std::size_t rows, std::size_t cols,
                         const std::vector<std::tuple<std::size_t,std::size_t,double>>& e) {
    SparseMatrix A(rows, cols);
    for (const auto& [i,j,v] : e) A.push_back(i,j,v);
    A.finalize();
    return A;
}

Dense dense(const SparseMatrix& A) {
    Dense d(A.n_rows(), std::vector<double>(A.n_cols(),0.0));
    for(std::size_t i=0;i<A.n_rows();++i)
        for(std::size_t k=A.row_offsets_data()[i];k<A.row_offsets_data()[i+1];++k)
            d[i][A.columns_data()[k]] += A.values_data()[k];
    return d;
}

Dense inverse(Dense A) {
    const std::size_t n=A.size(); Dense I(n,std::vector<double>(n,0.0));
    for(std::size_t i=0;i<n;++i) I[i][i]=1.0;
    for(std::size_t k=0;k<n;++k) {
        std::size_t p=k;
        for(std::size_t i=k+1;i<n;++i) if(std::abs(A[i][k])>std::abs(A[p][k])) p=i;
        EXPECT_TRUE(std::abs(A[p][k])>1e-12);
        std::swap(A[k],A[p]); std::swap(I[k],I[p]);
        const double diag=A[k][k];
        for(std::size_t j=0;j<n;++j){A[k][j]/=diag;I[k][j]/=diag;}
        for(std::size_t i=0;i<n;++i) if(i!=k) {
            const double f=A[i][k];
            for(std::size_t j=0;j<n;++j){A[i][j]-=f*A[k][j];I[i][j]-=f*I[k][j];}
        }
    }
    return I;
}

Vector matvec(const Dense& A,const Vector& x){
    Vector y(A.size(),0.0);
    for(std::size_t i=0;i<A.size();++i)
        for(std::size_t j=0;j<A[i].size();++j) y(i)+=A[i][j]*x(j);
    return y;
}

double relerr(const Vector& a,const Vector& b){
    return (a-b).norm2()/std::max(1.0,b.norm2());
}

struct Case {
    std::size_t n;
    double convection;
};

struct Blocks { SparseMatrix A,G,D,C; };

SparseMatrix make_block_matrix(const Blocks& b) {
    const std::size_t nu = b.A.n_rows();
    const std::size_t np = b.C.n_rows();
    std::vector<std::tuple<std::size_t,std::size_t,double>> e;
    auto append = [&](const SparseMatrix& A, std::size_t row0, std::size_t col0) {
        for (std::size_t i=0;i<A.n_rows();++i)
            for (std::size_t k=A.row_offsets_data()[i];k<A.row_offsets_data()[i+1];++k)
                e.emplace_back(row0+i,col0+A.columns_data()[k],A.values_data()[k]);
    };
    append(b.A,0,0);
    append(b.G,0,nu);
    append(b.D,nu,0);
    append(b.C,nu,nu);
    return make_sparse(nu+np,nu+np,e);
}

Blocks make_case(const Case& cfg) {
    const std::size_t np=cfg.n, nv=3*np;
    std::vector<std::tuple<std::size_t,std::size_t,double>> ae,ge,de,ce;
    for(std::size_t c=0;c<np;++c) {
        for(std::size_t comp=0;comp<3;++comp) {
            const std::size_t r=comp*np+c;
            ae.emplace_back(r,r,5.0+cfg.convection);
            if(c>0) ae.emplace_back(r,comp*np+c-1,-1.0-cfg.convection);
            if(c+1<np) ae.emplace_back(r,comp*np+c+1,-1.0);
        }
        // Same non-zero-flow pressure gradient/divergence stencil for every model.
        for(std::size_t comp=0;comp<3;++comp) {
            const std::size_t r=comp*np+c;
            ge.emplace_back(r,c,0.35);
            if(c>0) ge.emplace_back(r,c-1,-0.25);
            if(c+1<np) ge.emplace_back(r,c+1,0.25);
        }
        for(std::size_t comp=0;comp<3;++comp) {
            const std::size_t r=c;
            de.emplace_back(r,comp*np+c,-0.35);
            if(c>0) de.emplace_back(r,comp*np+c-1,0.25);
            if(c+1<np) de.emplace_back(r,comp*np+c+1,-0.25);
        }
        ce.emplace_back(c,c,0.20);
    }
    return {make_sparse(nv,nv,ae),make_sparse(nv,np,ge),
            make_sparse(np,nv,de),make_sparse(np,np,ce)};
}

SparseMatrix pressure_operator(std::size_t n,double diagonal,double offdiag,double bias=0.0){
    std::vector<std::tuple<std::size_t,std::size_t,double>> e;
    for(std::size_t i=0;i<n;++i){
        e.emplace_back(i,i,diagonal+bias);
        if(i>0) e.emplace_back(i,i-1,-offdiag);
        if(i+1<n) e.emplace_back(i,i+1,-offdiag);
    }
    return make_sparse(n,n,e);
}

double true_residual(const SparseMatrix& A,const Vector& x,const Vector& b){
    const auto ax=A.matvec(x); double r2=0,b2=0;
    for(std::size_t i=0;i<b.size();++i){const double r=b(i)-ax[i];r2+=r*r;b2+=b(i)*b(i);}
    return std::sqrt(r2/std::max(b2,1e-300));
}

}

int main() {
    run_case("n9_common_oseen_schur_contract", [] {
        const std::vector<Case> cases{{4,1.0},{4,4.0},{6,1.0}};
        for(const auto& cfg:cases){
            const auto d=make_case(cfg);
            const BlockOperator blocks(d.A,d.G,d.D,d.C);
            blocks.validate();

            const Dense Ad=dense(d.A), Gd=dense(d.G), Dd=dense(d.D), Cd=dense(d.C);
            const Dense Ainv=inverse(Ad);
                Dense Soracle=Cd;
            for(std::size_t i=0;i<cfg.n;++i)
                for(std::size_t j=0;j<cfg.n;++j)
                    for(std::size_t k=0;k<3*cfg.n;++k)
                        for(std::size_t l=0;l<3*cfg.n;++l)
                            Soracle[i][j]-=Dd[i][k]*Ainv[k][l]*Gd[l][j];
            const auto P=pressure_operator(cfg.n,2.0,0.25);
            const auto F=pressure_operator(cfg.n,2.0,0.25,cfg.convection*0.15);
            const auto M=pressure_operator(cfg.n,1.0,0.0);
            const auto pinv=inverse(dense(P));
            const auto finv=inverse(dense(F));

            ExactSchurApproximation exact([Ainv](const Vector& rhs,Vector& x){x=matvec(Ainv,rhs);return true;});
            EXPECT_TRUE(exact.setup(blocks));

            SimplerSchurApproximation simple(SimplerSchurMode::SIMPLE);
            SimplerSchurApproximation simplec(SimplerSchurMode::SIMPLEC);
            EXPECT_TRUE(simple.setup(blocks));
            EXPECT_TRUE(simplec.setup(blocks));

            LscBfbtSchurApproximation lsc(
                LscBfbtSchurApproximation::Mode::LSC,
                [&](const Vector& rhs,Vector& x){x=matvec(pinv,rhs);return true;});
            LscBfbtSchurApproximation bfbt(
                LscBfbtSchurApproximation::Mode::BFBT,
                [&](const Vector& rhs,Vector& x){x=matvec(pinv,rhs);return true;});
            EXPECT_TRUE(lsc.setup(blocks));
            EXPECT_TRUE(bfbt.setup(blocks));

            PcdSchurApproximation pcd(
                M,P,F,
                [&](const Vector& rhs,Vector& x){x=matvec(inverse(dense(P)),rhs);return true;},
                [&](const Vector& rhs,Vector& x){x=matvec(finv,rhs);return true;});
            EXPECT_TRUE(pcd.setup(blocks));

            Vector q(cfg.n); for(std::size_t i=0;i<cfg.n;++i) q(i)=std::sin(0.31*(i+1));
            Vector exact_Sq(cfg.n,0.0); EXPECT_TRUE(exact.apply_schur(q,exact_Sq));
            Vector simple_Sq(cfg.n,0.0),simplec_Sq(cfg.n,0.0);
            EXPECT_TRUE(simple.apply(q,simple_Sq)); EXPECT_TRUE(simplec.apply(q,simplec_Sq));
            std::cout<<"n9_schur_contract cells="<<cfg.n
                     <<" convection="<<cfg.convection
                     <<" simple_operator_error="<<relerr(simple_Sq,exact_Sq)
                     <<" simplec_operator_error="<<relerr(simplec_Sq,exact_Sq)<<'\n';

            Vector rhs(cfg.n); for(std::size_t i=0;i<cfg.n;++i) rhs(i)=std::cos(0.17*(i+1));
            Vector exact_inv=matvec(inverse(Soracle),rhs);
            Vector z(cfg.n,0.0);
            EXPECT_TRUE(lsc.apply(rhs,z));
            std::cout<<"n9_schur_contract cells="<<cfg.n<<" convection="<<cfg.convection
                     <<" lsc_inverse_error="<<relerr(z,exact_inv)<<'\n';
            EXPECT_TRUE(bfbt.apply(rhs,z));
            EXPECT_TRUE(pcd.apply(rhs,z));
        }
    });

    run_case("n9_common_oseen_full_block_true_residual_multiplex", [] {
        const std::size_t n=6;
        const auto d=make_case({n,4.0});
        const std::size_t N=4*n;
        const SparseMatrix K = make_block_matrix(d);
        Vector exact(N); for(std::size_t i=0;i<N;++i) exact(i)=std::sin(0.071*(i+1));
        const auto kb=K.matvec(exact);
        Vector b(N,0.0);
        for(std::size_t i=0;i<N;++i) b(i)=kb[i];
        const auto P=pressure_operator(n,2.0,0.25);
        const auto F=pressure_operator(n,2.0,0.25,4.0*0.15);
        const auto M=pressure_operator(n,1.0,0.0);
        const auto P_inv=inverse(dense(P));
        const auto F_inv=inverse(dense(F));
        const std::vector<std::pair<std::string,CoupledSchurApproximationModel>> models={
            {"simple",CoupledSchurApproximationModel::SIMPLE},
            {"simplec",CoupledSchurApproximationModel::SIMPLEC},
            {"lsc",CoupledSchurApproximationModel::LSC},
            {"bfbt",CoupledSchurApproximationModel::BFBT},
            {"pcd",CoupledSchurApproximationModel::PCD},
        };
        for(const auto& [name,model]:models){
            CoupledBlockSchurOptions o; o.schur_approximation=model;
            o.diagonal_schur_sign_flip=(model==CoupledSchurApproximationModel::PCD);
            CoupledBlockSchurAMGPreconditioner pc(n,o);
            if(model==CoupledSchurApproximationModel::SIMPLE)
                pc.set_simpler_schur(std::make_unique<SimplerSchurApproximation>(SimplerSchurMode::SIMPLE));
            if(model==CoupledSchurApproximationModel::SIMPLEC)
                pc.set_simpler_schur(std::make_unique<SimplerSchurApproximation>(SimplerSchurMode::SIMPLEC));
            if(model==CoupledSchurApproximationModel::LSC || model==CoupledSchurApproximationModel::BFBT) {
                const auto mode=model==CoupledSchurApproximationModel::LSC
                    ? LscBfbtSchurApproximation::Mode::LSC : LscBfbtSchurApproximation::Mode::BFBT;
                pc.set_algebraic_schur(std::make_unique<LscBfbtSchurApproximation>(
                    mode,[P_inv](const Vector& r,Vector& x){x=matvec(P_inv,r);return true;}));
            }
            if(model==CoupledSchurApproximationModel::PCD) {
                pc.set_pcd_schur(std::make_unique<PcdSchurApproximation>(
                    M,P,F,
                    [P_inv](const Vector& r,Vector& x){x=matvec(P_inv,r);return true;},
                    [F_inv](const Vector& r,Vector& x){x=matvec(F_inv,r);return true;}));
            }
            EXPECT_TRUE(pc.setup(K));
            Vector z(N,0.0); EXPECT_TRUE(pc.apply(b,z));
            const double rr=true_residual(K,z,b);
            std::cout<<"n9_model="<<name<<" true_residual="<<rr<<'\n';
            EXPECT_TRUE(std::isfinite(rr));
        }

        std::vector<std::size_t> fine(3*n),coarse(n);
        for(std::size_t i=0;i<3*n;++i) fine[i]=i;
        for(std::size_t i=0;i<n;++i) coarse[i]=3*n+i;
        NativeMGRPreconditioner mgr(fine,coarse);
        EXPECT_TRUE(mgr.setup(K));
        Vector z(N,0.0); EXPECT_TRUE(mgr.apply(b,z));
        const double mgr_rr=true_residual(K,z,b);
        std::cout<<"n9_model=mgr true_residual="<<mgr_rr<<'\\n';
        EXPECT_TRUE(std::isfinite(mgr_rr));
    });
    return run_all();
}
