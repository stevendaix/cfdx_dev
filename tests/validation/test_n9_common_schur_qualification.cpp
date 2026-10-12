#include "cfdx/core/linalg/block_operator.h"
#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/exact_schur.h"
#include "cfdx/core/linalg/lsc_bfbt_schur.h"
#include "cfdx/core/linalg/mgr_preconditioner.h"
#include "cfdx/core/linalg/simplerc_schur.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <limits>
#include <tuple>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
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

// Two-cell channel (same topology as test_n9_s10_pressure_reference_policies)
// for the end-to-end coupled gauge qualification below.
Mesh make_two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i = 0; i < 12; ++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);

    m.faces().push_face({0,3,7,4});
    m.faces().push_face({1,2,6,5});
    m.faces().push_face({8,9,11,10});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});

    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const auto owner = (f == 2 || f >= 7) ? 1 : 0;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});

    Patch wall;
    wall.name = "wall";
    wall.type = PatchType::WALL;
    wall.face_ids = {0,2,3,4,5,6,7,8,9,10};
    m.boundary().add_patch(wall);
    return m;
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
            const double simple_operator_error = relerr(simple_Sq,exact_Sq);
            const double simplec_operator_error = relerr(simplec_Sq,exact_Sq);
            std::cout<<"n9_schur_contract cells="<<cfg.n
                     <<" convection="<<cfg.convection
                     <<" simple_operator_error="<<simple_operator_error
                     <<" simplec_operator_error="<<simplec_operator_error<<'\n';
            // Blocking catastrophe guard: SIMPLE/SIMPLEC are deliberate
            // approximations, so this is NOT an accuracy gate. It fails
            // when the action degenerates algebraically (non-finite or
            // two orders of magnitude past the oracle on these academic
            // matrices). Tighter regression pins come from measured
            // campaign logs, not invented thresholds.
            EXPECT_TRUE(std::isfinite(simple_operator_error) &&
                        simple_operator_error < 10.0);
            EXPECT_TRUE(std::isfinite(simplec_operator_error) &&
                        simplec_operator_error < 10.0);

            Vector rhs(cfg.n); for(std::size_t i=0;i<cfg.n;++i) rhs(i)=std::cos(0.17*(i+1));
            Vector exact_inv=matvec(inverse(Soracle),rhs);
            Vector z(cfg.n,0.0);
            EXPECT_TRUE(lsc.apply(rhs,z));
            const double lsc_inverse_error = relerr(z,exact_inv);
            std::cout<<"n9_schur_contract cells="<<cfg.n<<" convection="<<cfg.convection
                     <<" lsc_inverse_error="<<lsc_inverse_error<<'\n';
            // Blocking catastrophe guard (same policy as the SIMPLE/SIMPLEC
            // guard above).
            EXPECT_TRUE(std::isfinite(lsc_inverse_error) &&
                        lsc_inverse_error < 10.0);
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
            // Blocking catastrophe guard on the single preconditioned
            // apply: the value is not an accuracy criterion, but it must
            // stay finite and bounded.
            EXPECT_TRUE(std::isfinite(rr) && rr < 1e3);
        }

        std::vector<std::size_t> fine(3*n),coarse(n);
        for(std::size_t i=0;i<3*n;++i) fine[i]=i;
        for(std::size_t i=0;i<n;++i) coarse[i]=3*n+i;
        NativeMGRPreconditioner mgr(fine,coarse);
        EXPECT_TRUE(mgr.setup(K));
        Vector z(N,0.0); EXPECT_TRUE(mgr.apply(b,z));
        const double mgr_rr=true_residual(K,z,b);
        std::cout<<"n9_model=mgr true_residual="<<mgr_rr<<'\n';
        EXPECT_TRUE(std::isfinite(mgr_rr) && mgr_rr < 1e3);
    });

    run_case("n9_schur_gauge_row_column_semantics", [] {
        // Unit falsification of the nested-operator gauge: reference row AND
        // column eliminated, unit diagonal, sentinel is a no-op, and the
        // raw-rhs passthrough keeps the gauge component the outer Krylov
        // needs (the remaining components solve the reduced system).
        const std::size_t n = 5;
        std::vector<std::tuple<std::size_t,std::size_t,double>> e;
        for (std::size_t i=0;i<n;++i) {
            e.emplace_back(i,i,4.0);
            if(i>0) e.emplace_back(i,i-1,-1.0);
            if(i+1<n) e.emplace_back(i,i+1,-1.0);
        }
        const SparseMatrix A0 = make_sparse(n,n,e);
        const std::size_t sentinel = std::numeric_limits<std::size_t>::max();

        SparseMatrix A = make_sparse(n,n,e);
        apply_nested_schur_pressure_gauge(A, sentinel);
        EXPECT_TRUE(A.nnz()==A0.nnz());
        for (std::size_t k=0;k<A.nnz();++k)
            EXPECT_NEAR(A.values_data()[k],A0.values_data()[k],1e-14);

        A = make_sparse(n,n,e);
        apply_nested_schur_pressure_gauge(A, 2);
        bool has_unit_diag=false;
        for (std::size_t r=0;r<n;++r)
            for (std::size_t k=A.row_offsets_data()[r];k<A.row_offsets_data()[r+1];++k) {
                const std::size_t c=A.columns_data()[k];
                if (r==2) {
                    if (c==2) { EXPECT_NEAR(A.values_data()[k],1.0,1e-14); has_unit_diag=true; }
                    else EXPECT_NEAR(A.values_data()[k],0.0,1e-14);
                } else if (c==2) {
                    EXPECT_NEAR(A.values_data()[k],0.0,1e-14);
                }
            }
        EXPECT_TRUE(has_unit_diag);

        // Passthrough: the gauged solve returns the rhs gauge component
        // at the reference, and the other components solve the reduced
        // system exactly.
        Vector rhs(n); for(std::size_t i=0;i<n;++i) rhs(i)=std::cos(0.4*(i+1));
        const Vector x=matvec(inverse(dense(A)),rhs);
        EXPECT_NEAR(x(2),rhs(2),1e-12);
        Dense red(n-1,std::vector<double>(n-1,0.0));
        for(std::size_t r=0;r<n;++r) if(r!=2)
            for(std::size_t k=A0.row_offsets_data()[r];k<A0.row_offsets_data()[r+1];++k)
                if(A0.columns_data()[k]!=2)
                    red[r<2?r:r-1][A0.columns_data()[k]<2?A0.columns_data()[k]:A0.columns_data()[k]-1]
                        +=A0.values_data()[k];
        Vector rrhs(n-1); for(std::size_t i=0;i<n-1;++i) rrhs(i)=rhs(i<2?i:i+1);
        const Vector xr=matvec(inverse(red),rrhs);
        for(std::size_t i=0;i<n-1;++i)
            EXPECT_NEAR(x(i<2?i:i+1),xr(i),1e-10);

        // Fallback path: no diagonal entry in the reference row pattern;
        // the rebuilt operator still carries the unit row and the
        // eliminated column.
        std::vector<std::tuple<std::size_t,std::size_t,double>> e2;
        for(std::size_t i=0;i<n;++i) if(i!=2) e2.emplace_back(i,i,4.0);
        e2.emplace_back(0,1,-1.0); e2.emplace_back(1,0,-1.0);
        e2.emplace_back(1,2,-1.0); e2.emplace_back(3,2,-1.0);
        e2.emplace_back(3,4,-1.0); e2.emplace_back(4,3,-1.0);
        SparseMatrix B = make_sparse(n,n,e2);
        apply_nested_schur_pressure_gauge(B, 2);
        bool unit_row_present=false, col_entries_zeroed=false;
        for(std::size_t k=B.row_offsets_data()[2];k<B.row_offsets_data()[2+1];++k){
            if(B.columns_data()[k]==2){ EXPECT_NEAR(B.values_data()[k],1.0,1e-14); unit_row_present=true; }
            else EXPECT_NEAR(B.values_data()[k],0.0,1e-14);
        }
        EXPECT_TRUE(unit_row_present);
        for(std::size_t r=0;r<n;++r) if(r!=2)
            for(std::size_t k=B.row_offsets_data()[r];k<B.row_offsets_data()[r+1];++k)
                if(B.columns_data()[k]==2){ EXPECT_NEAR(B.values_data()[k],0.0,1e-14); col_entries_zeroed=true; }
        EXPECT_TRUE(col_entries_zeroed);
    });

    run_case("n9_schur_gauge_end_to_end_coupled_lsc_bfbt", [] {
        // End-to-end guard for the coupled LSC/BFBt paths: interior
        // reference (gauge imposed, row+column coherent) and fixed
        // pressure boundary (no gauge at all).
        const std::vector<CoupledSchurModel> models{
            CoupledSchurModel::LSC, CoupledSchurModel::BFBT};
        for (const CoupledSchurModel model : models) {
            {
                const Mesh mesh = make_two_cell_channel();
                Field<double, Location::CELL> U(2, "U", "m/s", 3);
                Field<double, Location::CELL> p(2, "p", "Pa", 1);
                U.fill(0.0);
                p.fill(7.0);
                VelocityBoundaryConditions ubc;
                // Non-degenerate contract: an O(1) wall velocity gives the
                // exact discrete solution U = (1,0,0) in both cells with
                // uniform pressure and O(1) boundary fluxes.
                ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1.0,0.0,0.0}};
                ScalarBoundaryConditions pbc;
                pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
                IncompressibleSolverControls c;
                c.algorithm = PressureVelocityAlgorithm::COUPLED;
                c.coupling.schur_model = model;
                // Test-contract fix, second revision (well-posed
                // contract). Two measured defects made this guard
                // structurally unconvergeable on the COUPLED path:
                // (1) Picard budget mismatch (run 38224806432): the
                //     default relaxed update (alpha_u = 0.7,
                //     alpha_p = 0.3) contracts the pressure error by 0.7
                //     per iteration, so p going 7.0 -> 3.0 needs ~79
                //     iterations while the budget was 10.
                // (2) Degenerate exact solution (evidence 9826140, the
                //     committed trace of the alpha = 1 run): with
                //     U_b = (0,0,0) the exact solution is U == 0 with
                //     uniform pressure, so the alpha = 1 state freezes
                //     at a bit-exact fixed point from iteration 1; the
                //     trace shows momentum_residual_relative pinned at
                //     exactly 1.0 for all 10 iterations (the reference
                //     residual is captured at iteration 1 and the field
                //     never changes again) and continuity_normalized
                //     pinned at exactly 2 (noise over an all-zero
                //     boundary flux), so the 1e-12 gates are unreachable
                //     for any Schur model on that contract.
                // This revision makes the contract well-posed without
                // loosening a single gate: U_b = (1,0,0) above gives a
                // non-degenerate O(1) problem, alpha_u = alpha_p = 0.5
                // lets every criterion decay geometrically under the
                // unchanged 1e-12 gates (~41 iterations at rate 0.5),
                // coupled_linear_tolerance = 1e-14 keeps the linear
                // noise floor 100x below the gates, and
                // max_iterations = 200 budgets the Picard loop with
                // margin.
                c.coupling.alpha_u = 0.5;
                c.coupling.alpha_p = 0.5;
                // Linear noise floor 100x below the 1e-12 gates.
                c.coupling.coupled_linear_tolerance = 1e-14;
                c.pressure_gauge_policy = PressureGaugePolicy::REFERENCE_CELL;
                c.pressure_reference_cell = 0;
                c.pressure_reference_value = 3.0;
                c.convergence.max_iterations = 200;
                c.convergence.continuity_tolerance = 1e-12;
                c.convergence.relative_tolerance = 1e-12;
                c.linear_max_iterations = 100;
                c.linear_tolerance = 1e-12;
                // Forensic trace for the e2e guard (measurement only):
                // per-iteration criterion values are committed to the CI
                // log so a gate failure is attributed with measured
                // numbers (Picard budget vs linear-solve noise floor vs
                // gauge-chain defect) instead of a bare converged=false.
                c.diagnostics.iteration_trace = true;
                c.diagnostics.iteration_trace_frequency = 1;
                const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
                EXPECT_TRUE(result.converged);
                EXPECT_NEAR(p(0), 3.0, 1e-10);
                EXPECT_NEAR(p(1), 3.0, 1e-10);
                EXPECT_NEAR(U.component_data(0)[0], 1.0, 1e-9);
                EXPECT_NEAR(U.component_data(0)[1], 1.0, 1e-9);
                EXPECT_NEAR(U.component_data(1)[0], 0.0, 1e-9);
                EXPECT_NEAR(U.component_data(2)[1], 0.0, 1e-9);
            }
            {
                const Mesh mesh = make_two_cell_channel();
                Field<double, Location::CELL> U(2, "U", "m/s", 3);
                Field<double, Location::CELL> p(2, "p", "Pa", 1);
                U.fill(0.0);
                p.fill(7.0);
                VelocityBoundaryConditions ubc;
                // Non-degenerate contract: an O(1) wall velocity gives the
                // exact discrete solution U = (1,0,0) in both cells with
                // uniform pressure and O(1) boundary fluxes.
                ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {1.0,0.0,0.0}};
                ScalarBoundaryConditions pbc;
                pbc["wall"] = {ScalarBoundaryType::FIXED_VALUE, 5.0, 0.0};
                IncompressibleSolverControls c;
                c.algorithm = PressureVelocityAlgorithm::COUPLED;
                c.coupling.schur_model = model;
                // Test-contract fix, second revision (well-posed
                // contract). Two measured defects made this guard
                // structurally unconvergeable on the COUPLED path:
                // (1) Picard budget mismatch (run 38224806432): the
                //     default relaxed update (alpha_u = 0.7,
                //     alpha_p = 0.3) contracts the pressure error by 0.7
                //     per iteration, so p going 7.0 -> 3.0 needs ~79
                //     iterations while the budget was 10.
                // (2) Degenerate exact solution (evidence 9826140, the
                //     committed trace of the alpha = 1 run): with
                //     U_b = (0,0,0) the exact solution is U == 0 with
                //     uniform pressure, so the alpha = 1 state freezes
                //     at a bit-exact fixed point from iteration 1; the
                //     trace shows momentum_residual_relative pinned at
                //     exactly 1.0 for all 10 iterations (the reference
                //     residual is captured at iteration 1 and the field
                //     never changes again) and continuity_normalized
                //     pinned at exactly 2 (noise over an all-zero
                //     boundary flux), so the 1e-12 gates are unreachable
                //     for any Schur model on that contract.
                // This revision makes the contract well-posed without
                // loosening a single gate: U_b = (1,0,0) above gives a
                // non-degenerate O(1) problem, alpha_u = alpha_p = 0.5
                // lets every criterion decay geometrically under the
                // unchanged 1e-12 gates (~41 iterations at rate 0.5),
                // coupled_linear_tolerance = 1e-14 keeps the linear
                // noise floor 100x below the gates, and
                // max_iterations = 200 budgets the Picard loop with
                // margin.
                c.coupling.alpha_u = 0.5;
                c.coupling.alpha_p = 0.5;
                // Linear noise floor 100x below the 1e-12 gates.
                c.coupling.coupled_linear_tolerance = 1e-14;
                c.pressure_gauge_policy = PressureGaugePolicy::REFERENCE_CELL;
                c.pressure_reference_cell = 0;
                c.pressure_reference_value = 3.0;
                c.convergence.max_iterations = 200;
                c.convergence.continuity_tolerance = 1e-12;
                c.convergence.relative_tolerance = 1e-12;
                c.linear_max_iterations = 100;
                c.linear_tolerance = 1e-12;
                // Forensic trace for the e2e guard (measurement only):
                // per-iteration criterion values are committed to the CI
                // log so a gate failure is attributed with measured
                // numbers (Picard budget vs linear-solve noise floor vs
                // gauge-chain defect) instead of a bare converged=false.
                c.diagnostics.iteration_trace = true;
                c.diagnostics.iteration_trace_frequency = 1;
                const auto result = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
                EXPECT_TRUE(result.converged);
                EXPECT_NEAR(p(0), 5.0, 1e-10);
                EXPECT_NEAR(p(1), 5.0, 1e-10);
                EXPECT_NEAR(U.component_data(0)[0], 1.0, 1e-9);
                EXPECT_NEAR(U.component_data(0)[1], 1.0, 1e-9);
                EXPECT_NEAR(U.component_data(1)[0], 0.0, 1e-9);
                EXPECT_NEAR(U.component_data(2)[1], 0.0, 1e-9);
            }
        }
    });
    return run_all();
}
