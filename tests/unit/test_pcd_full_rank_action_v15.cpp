#include "cfdx/physics/pcd_fp_basis_registry_v15.h"
#include "cfdx/physics/pcd_fp_deflated_quality_v14.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
static void verify(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main() {
  try {
    // Rank-one nonsymmetric matrix with different left/right nullspaces.
    SparseMatrix A(2,2);
    A.push_back(0,0,1); A.push_back(0,1,-1);
    A.push_back(1,0,2); A.push_back(1,1,-2); A.finalize();
    Vector w(2), z(2);
    w(0)=-2/std::sqrt(5.);w(1)=1/std::sqrt(5.);
    z(0)=1/std::sqrt(2.);z(1)=1/std::sqrt(2.);
    auto basis=PcdFpBiNullspaceV13::certify(A,{w},{z});
    Vector incompatible=w, legacy=basis.compatible_rhs(incompatible);
    verify(basis.dot(legacy,legacy)<1e-27,"V14 projected action should remove incompatible RHS");
    LinearOperator deflated;
    deflated.size=2;
    deflated.apply=[&](const Vector& x,Vector& out) {out=basis.deflated_matvec(A,x);};
    deflated.apply_fp64=deflated.apply;
    Vector full_solution(2,0.0);
    auto solved=solve_fgmres(deflated,incompatible,full_solution,2,24,1e-12,nullptr);
    verify(solved.status==SolverStatus::CONVERGED,"V15 full-rank deflated solve failed");
    verify(std::abs(basis.dot(z,full_solution)-1.)<1e-11,
           "V15 failed to transport left-null RHS to right-null solution");
    const auto quality=pcd_fp_deflated_true_residual_v14(
        A,basis,incompatible,Vector(2,0.0),full_solution);
    verify(quality.relative_rhs<1e-11,"V15 independent auxiliary true residual not closed");
    verify(std::abs(basis.incompatibility_relative(incompatible)-1.)<1e-13,
           "V15 original Fp floor hidden");
    auto changed=A;
    changed.values_data()[0]+=1e-4;
    verify(pcd_fp_basis_registry_key_v15(changed)!=pcd_fp_basis_registry_key_v15(A),
           "V15 coefficient update did not change key");
    bool stale=false;
    try {(void)basis.deflated_matvec(changed,full_solution);}
    catch(const std::invalid_argument&) {stale=true;}
    verify(stale,"V15 stale basis used on changed Fp");
    const auto directory=std::filesystem::temp_directory_path()/"cfdx-v15-fp-basis-test";
    std::filesystem::create_directories(directory);
    const auto file=directory/(pcd_fp_basis_registry_key_v15(A)+".cfpb");
    {
        std::ofstream out(file,std::ios::binary|std::ios::trunc);
        const std::uint64_t header[]={2,1,basis.operator_hash};
        out.write("CFDXNB13",8);
        out.write(reinterpret_cast<const char*>(header),sizeof(header));
        out.write(reinterpret_cast<const char*>(w.data()),2*sizeof(double));
        out.write(reinterpret_cast<const char*>(z.data()),2*sizeof(double));
    }
    const auto matched=load_pcd_fp_basis_registry_v15(directory.string(),A);
    verify(matched.operator_hash==basis.operator_hash,"V15 registry selected wrong exact basis");
    // A new outer iteration can rebuild Fp with identical sparsity but new
    // coefficients.  Its own certificate must be loaded, never the old one.
    auto twice=A;
    for (std::size_t entry=0;entry<twice.nnz();++entry)
        twice.values_data()[entry]*=2.;
    auto new_basis=PcdFpBiNullspaceV13::certify(twice,{w},{z});
    const auto new_file=directory/(pcd_fp_basis_registry_key_v15(twice)+".cfpb");
    {
        std::ofstream out(new_file,std::ios::binary|std::ios::trunc);
        const std::uint64_t header[]={2,1,new_basis.operator_hash};
        out.write("CFDXNB13",8);
        out.write(reinterpret_cast<const char*>(header),sizeof(header));
        out.write(reinterpret_cast<const char*>(w.data()),2*sizeof(double));
        out.write(reinterpret_cast<const char*>(z.data()),2*sizeof(double));
    }
    const auto refreshed=load_pcd_fp_basis_registry_v15(directory.string(),twice);
    verify(refreshed.operator_hash==new_basis.operator_hash &&
               refreshed.operator_hash!=basis.operator_hash,
           "V15 coefficient refresh did not select the fresh certificate");
    bool missing=false;
    try {(void)load_pcd_fp_basis_registry_v15(directory.string(),changed);}
    catch(const std::runtime_error&) {missing=true;}
    verify(missing,"V15 registry silently reused previous matrix basis");
    std::filesystem::remove(file);
    std::filesystem::remove(new_file);
    std::filesystem::remove(directory);
    std::cout << "PCD_FULL_RANK_V15_PASS mapped_left_to_right=1 original_floor=1 "
                 "deflated_true=" << quality.relative_rhs << " registry_stale_rejected=1 registry_refresh=1\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr << "PCD_FULL_RANK_V15_FAIL: " << error.what() << '\n';
    return 1;
  }
}
