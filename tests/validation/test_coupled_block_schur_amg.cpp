#include "cfdx/core/linalg/coupled_amg_schur.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <algorithm>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
SparseMatrix make_coupled(std::size_t n, double velocity_diag=4.0) {
    SparseMatrix A(4*n,4*n);
    for(std::size_t c=0;c<n;++c){
        for(std::size_t d=0;d<3;++d) A.push_back(d*n+c,d*n+c,velocity_diag);
        A.push_back(3*n+c,3*n+c,1.0);
        A.push_back(3*n+c,c,1.0);
        A.push_back(c,3*n+c,1.0);
        if(c+1<n){ A.push_back(c,3*n+c+1,-1.0); A.push_back(3*n+c+1,c,-1.0); }
    }
    A.finalize(); return A;
}

SparseMatrix make_coupled_changed_pattern(std::size_t n, double velocity_diag=4.0) {
    SparseMatrix A(4*n,4*n);
    for(std::size_t c=0;c<n;++c){
        for(std::size_t d=0;d<3;++d) A.push_back(d*n+c,d*n+c,velocity_diag);
        A.push_back(3*n+c,3*n+c,1.0);
        A.push_back(3*n+c,c,1.0);
        A.push_back(c,3*n+c,1.0);
        if(c+1<n){ A.push_back(c,3*n+c+1,-1.0); A.push_back(3*n+c+1,c,-1.0); }
    }
    A.push_back(3*n,1,0.125);
    A.finalize();
    return A;
}
SparseMatrix make_gauged_coupled(std::size_t n, std::size_t reference_cell,
                                      double velocity_diag=4.0) {
    SparseMatrix A(4*n,4*n);
    for(std::size_t c=0;c<n;++c){
        for(std::size_t d=0;d<3;++d) A.push_back(d*n+c,d*n+c,velocity_diag);
        if(c == reference_cell) {
            A.push_back(3*n+c,3*n+c,1.0);
            continue;
        }
        A.push_back(3*n+c,3*n+c,1.0);
        A.push_back(3*n+c,c,1.0);
        A.push_back(c,3*n+c,1.0);
        if(c+1<n && c+1 != reference_cell){
            A.push_back(c,3*n+c+1,-1.0);
            A.push_back(3*n+c+1,c,-1.0);
        }
    }
    A.finalize();
    return A;
}

Vector matvec(const SparseMatrix&A,const Vector&x){auto v=A.matvec(x);Vector y(v.size());for(std::size_t i=0;i<v.size();++i)y(i)=v[i];return y;}
double relres(const SparseMatrix&A,const Vector&x,const Vector&b){auto ax=A.matvec(x);double r=0,q=0;for(std::size_t i=0;i<b.size();++i){double e=b(i)-ax[i];r+=e*e;q+=b(i)*b(i);}return std::sqrt(r/std::max(q,1e-300));}
}
int main(){
 run_case("coupled_4n_block_schur_preserves_pinned_pressure_gauge",[] {
   constexpr std::size_t n = 32;
   constexpr std::size_t reference_cell = 7;
   const auto A=make_gauged_coupled(n, reference_cell);
   CoupledBlockSchurAMGPreconditioner pc(n, CoupledBlockSchurOptions{}, reference_cell);
   EXPECT_TRUE(pc.setup(A));
   EXPECT_TRUE(pc.is_ready());
   EXPECT_TRUE(pc.pressure_hierarchy_builds()==1);
   // Coarsening depth is backend/graph dependent; the contract here is that
   // the reduced pressure operator is solvable and the pinned correction is zero.
   Vector rhs(4*n, 0.0), z(4*n, 0.0);
   for(std::size_t i=0;i<4*n;++i) rhs(i)=std::sin(0.037*(i+1));
   EXPECT_TRUE(pc.apply(rhs,z));
   EXPECT_TRUE(std::abs(z(3*n+reference_cell)) < 1e-14);
 });
 run_case("coupled_4n_block_schur_uses_pressure_amg",[] {
   const auto A=make_coupled(16); Vector exact(64);
   for(std::size_t i=0;i<64;++i) exact(i)=std::sin(0.07*(i+1));
   const auto b=matvec(A,exact); Vector x(64,0.0);
   CoupledBlockSchurAMGPreconditioner pc(16);
   EXPECT_TRUE(pc.setup(A)); EXPECT_TRUE(pc.is_ready());
   EXPECT_TRUE(pc.pressure_hierarchy_builds()==1);
   const auto& S=pc.schur_matrix();
   // Exact Schur for this deliberately staggered D/G stencil:
   // S = C - D M^{-1} G, M=4I. The first pressure cell has one
   // velocity-pressure coupling contribution; every later pressure row couples
   // to its own velocity and the preceding velocity, hence diag={0.75,0.5,...,0.5}.
   // Both off-diagonals are +0.25.
   for(std::size_t i=0;i<16;++i) {
     const double expected_diag = (i == 0) ? 0.75 : 0.5;
     EXPECT_TRUE(std::abs(S(i,i)-expected_diag)<1e-12);
     if(i+1<16) {
       EXPECT_TRUE(std::abs(S(i,i+1)-0.25)<1e-12);
       EXPECT_TRUE(std::abs(S(i+1,i)-0.25)<1e-12);
     }
   }
   EXPECT_TRUE(pc.apply(b,x)); EXPECT_TRUE(relres(A,x,b)<0.5);
 });
 run_case("coupled_4n_gmres_converges_with_block_schur_amg",[] {
   const auto A=make_coupled(16); Vector exact(64);
   for(std::size_t i=0;i<64;++i) exact(i)=std::cos(0.05*(i+1));
   const auto b=matvec(A,exact); Vector x(64,0.0);
   CoupledBlockSchurAMGPreconditioner pc(16);
   const auto r=solve_gmres(A,b,x,20,200,1e-10,&pc);
   EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
   EXPECT_TRUE(relres(A,x,b)<1e-9);
   EXPECT_TRUE((x-exact).norm_inf()<1e-8);
 });
 // Schur scaling ladder: verify that the 4N block construction remains
 // convergent as the pressure/velocity system grows beyond the small exact
 // algebraic contract above.
 for (const std::size_t n : {8u, 16u, 32u, 64u}) {
   const auto A=make_coupled(n);
   CoupledBlockSchurAMGPreconditioner pc(n);
   EXPECT_TRUE(pc.setup(A));
   EXPECT_TRUE(pc.is_ready());
   EXPECT_TRUE(pc.pressure_hierarchy_builds()==1);
   Vector exact(4*n);
   for(std::size_t i=0;i<4*n;++i) exact(i)=std::cos(0.031*(i+1));
   const auto b=matvec(A,exact); Vector x(4*n,0.0);
   const auto r=solve_gmres(A,b,x,20,400,1e-10,&pc);
   const double rr=relres(A,x,b);
   std::cout << "schur_scaling n_cells=" << n
             << " unknowns=" << 4*n
             << " status=" << static_cast<int>(r.status)
             << " iterations=" << r.iterations
             << " true_residual=" << rr
             << " pressure_hierarchy_builds=" << pc.pressure_hierarchy_builds()
             << " pressure_numeric_updates=" << pc.pressure_numeric_updates() << '\n';
   EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
   EXPECT_TRUE(std::isfinite(rr));
   EXPECT_TRUE(rr<1e-9);
 }

 run_case("coupled_4n_petsc_style_schur_factorizations",[] {
   const auto A=make_coupled(16); Vector exact(64);
   for(std::size_t i=0;i<64;++i) exact(i)=std::sin(0.11*(i+1));
   const auto b=matvec(A,exact);
   const CoupledSchurFactorization modes[] = {
       CoupledSchurFactorization::Diagonal,
       CoupledSchurFactorization::Lower,
       CoupledSchurFactorization::Upper,
       CoupledSchurFactorization::Full};
   Vector rp_only(64,0.0), z_diag(64,0.0), z_lower(64,0.0), z_upper(64,0.0), z_full(64,0.0);
   for(std::size_t c=0;c<16;++c) rp_only(48+c)=1.0+0.01*c;
   for(const auto mode : modes) {
     CoupledBlockSchurOptions options;
     options.factorization=mode;
     CoupledBlockSchurAMGPreconditioner pc(16, options);
     Vector x(64,0.0);
     EXPECT_TRUE(pc.setup(A));
     EXPECT_TRUE(pc.apply(b,x));
     EXPECT_TRUE(relres(A,x,b)<0.5);
     if(mode==CoupledSchurFactorization::Diagonal) { EXPECT_TRUE(pc.apply(rp_only,z_diag)); }
     if(mode==CoupledSchurFactorization::Lower) { EXPECT_TRUE(pc.apply(rp_only,z_lower)); }
     if(mode==CoupledSchurFactorization::Upper) { EXPECT_TRUE(pc.apply(rp_only,z_upper)); }
     if(mode==CoupledSchurFactorization::Full) { EXPECT_TRUE(pc.apply(rp_only,z_full)); }
   }
   // The four algebraic factorizations must not collapse to one implementation:
   // diagonal/lower leave the velocity field untouched for a pressure-only
   // RHS, whereas upper/full apply the G-correction. The pressure block is
   // expected to be non-zero for all four variants.
   auto velocity_norm = [](const Vector& v) {
     double m = 0.0;
     for (std::size_t i=0; i<48; ++i) m = std::max(m, std::abs(v(i)));
     return m;
   };
   auto pressure_norm = [](const Vector& v) {
     double m = 0.0;
     for (std::size_t i=48; i<64; ++i) m = std::max(m, std::abs(v(i)));
     return m;
   };
   EXPECT_TRUE(velocity_norm(z_diag) < 1e-14);
   EXPECT_TRUE(velocity_norm(z_lower) < 1e-14);
   EXPECT_TRUE(pressure_norm(z_diag) > 1e-12);
   EXPECT_TRUE(pressure_norm(z_lower) > 1e-12);
   EXPECT_TRUE(velocity_norm(z_upper) > 1e-12);
   EXPECT_TRUE(velocity_norm(z_full) > 1e-12);
   EXPECT_TRUE(pressure_norm(z_upper) > 1e-12);
   EXPECT_TRUE(pressure_norm(z_full) > 1e-12);
   EXPECT_TRUE((z_upper-z_full).norm_inf() < 1e-12);
 });

 run_case("coupled_4n_hypre_style_block_diagonal_schur_approximation",[] {
   SparseMatrix A(32,32);
   const std::size_t n=8;
   for(std::size_t c=0;c<n;++c) {
     const double block[3][3]={{4.0,0.5,0.0},{0.5,5.0,0.25},{0.0,0.25,3.0}};
     for(std::size_t r=0;r<3;++r)
       for(std::size_t q=0;q<3;++q)
         if(block[r][q]!=0.0) A.push_back(r*n+c,q*n+c,block[r][q]);
     A.push_back(3*n+c,3*n+c,1.0);
     A.push_back(3*n+c,c,1.0);
     A.push_back(c,3*n+c,1.0);
     if(c+1<n) {
       A.push_back(c,3*n+c+1,-1.0);
       A.push_back(3*n+c+1,c,-1.0);
     }
   }
   A.finalize();

   CoupledBlockSchurOptions block_options;
   block_options.velocity_approximation=CoupledSchurVelocityApproximation::Block;
   CoupledBlockSchurAMGPreconditioner block_pc(n,block_options);
   EXPECT_TRUE(block_pc.setup(A));

   CoupledBlockSchurOptions diag_options;
   diag_options.velocity_approximation=CoupledSchurVelocityApproximation::Diagonal;
   CoupledBlockSchurAMGPreconditioner diag_pc(n,diag_options);
   EXPECT_TRUE(diag_pc.setup(A));

   EXPECT_TRUE(block_pc.schur_matrix().nnz()==diag_pc.schur_matrix().nnz());
   bool differs=false;
   for(std::size_t r=0;r<n && !differs;++r)
     for(std::size_t c=0;c<n;++c)
       if(std::abs(block_pc.schur_matrix()(r,c)-diag_pc.schur_matrix()(r,c))>1e-12)
         { differs=true; break; }
   EXPECT_TRUE(differs);

   Vector exact(32);
   for(std::size_t i=0;i<32;++i) exact(i)=std::cos(0.09*(i+1));
   const auto b=matvec(A,exact);
   Vector x(32,0.0);
   const auto r=solve_gmres(A,b,x,20,200,1e-10,&block_pc);
   EXPECT_TRUE(r.status==SolverStatus::CONVERGED);
   EXPECT_TRUE(relres(A,x,b)<1e-9);
   EXPECT_TRUE((x-exact).norm_inf()<1e-8);
 });

 {
   // A gauged Schur operator must keep its reduced (n-1) AMG hierarchy
   // across numeric updates. This catches accidentally updating AMG with the
   // full n-cell Schur matrix after setup eliminated the reference DOF.
   constexpr std::size_t n = 8;
   constexpr std::size_t reference_cell = 3;
   const auto A = make_gauged_coupled(n, reference_cell);
   const auto A2 = make_gauged_coupled(n, reference_cell, 4.25);
   CoupledBlockSchurAMGPreconditioner pc(
       n, CoupledBlockSchurOptions{}, reference_cell);
   EXPECT_TRUE(pc.setup(A));
   EXPECT_TRUE(pc.pressure_hierarchy_builds() == 1);
   EXPECT_TRUE(pc.update_values(A2));
   EXPECT_TRUE(pc.pressure_hierarchy_builds() == 1);
   EXPECT_TRUE(pc.pressure_numeric_updates() == 1);
   Vector rhs(4*n, 1.0), z(4*n, 0.0);
   EXPECT_TRUE(pc.apply(rhs, z));
   EXPECT_TRUE(std::abs(z(3*n + reference_cell)) < 1e-14);
 }

 {
   const auto A=make_coupled(8);
   // update_values() is numeric-only: it must not implicitly create a
   // symbolic hierarchy when setup() has not been called.
   CoupledBlockSchurAMGPreconditioner uninitialized(8);
   EXPECT_TRUE(!uninitialized.update_values(A));
   EXPECT_TRUE(!uninitialized.is_ready());

   // Change only coefficients: the Schur graph is unchanged and the native
   // AMG hierarchy must be reused through the numeric-update path.
   const auto A2=make_coupled(8,4.25);
   CoupledBlockSchurAMGPreconditioner pc(8);
   EXPECT_TRUE(pc.setup(A));
   EXPECT_TRUE(pc.pressure_hierarchy_builds()==1);
   EXPECT_TRUE(pc.pressure_numeric_updates()==0);
   EXPECT_TRUE(pc.update_values(A2));
   EXPECT_TRUE(pc.pressure_hierarchy_builds()==1);
   EXPECT_TRUE(pc.pressure_numeric_updates()==1);

   // A new Schur edge invalidates the symbolic hierarchy and must be rejected
   // rather than silently rebuilding or falling back.
   const auto changed=make_coupled_changed_pattern(8,4.25);
   EXPECT_TRUE(!pc.update_values(changed));
   EXPECT_TRUE(pc.last_error().find("pattern changed") != std::string::npos);
   // Rejecting the update must preserve the old valid hierarchy; only an
   // explicit setup() is allowed to replace it.
   EXPECT_TRUE(pc.is_ready());
   EXPECT_TRUE(pc.pressure_hierarchy_builds()==1);
   EXPECT_TRUE(pc.pressure_numeric_updates()==1);
 }

 return run_all();
}
