#include "cfdx/core/linalg/cg_solver.h"
#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/bicgstab_solver.h"
#include "cfdx/core/parallel/mpi_utils.h"
#include <cmath>
#include <iostream>

int main(int argc, char** argv) {
    cfdx::core::parallel::mpi_init(&argc, &argv);
    const int rank = cfdx::core::parallel::mpi_rank();
    const int size = cfdx::core::parallel::mpi_size();
    if (size != 2) {
        if (rank == 0) std::cerr << "requires 2 MPI ranks\n";
        cfdx::core::parallel::mpi_finalize();
        return 77;
    }
    cfdx::core::SparseMatrix A(4, 4);
    for (std::size_t i = 0; i < 4; ++i) A.add_entry(i, i, 2.0);
    A.finalize();
    cfdx::core::Vector b(4, 1.0);
    cfdx::core::KrylovReductionPolicy reduction;
    reduction.mpi_enabled = true;
    reduction.deterministic = true;
    reduction.local_begin = static_cast<std::size_t>(rank * 2);
    reduction.local_end = reduction.local_begin + 2;

    cfdx::core::Vector x1(4, 0.0), x2(4, 0.0), x3(4, 0.0);
    const auto cg = cfdx::core::solve_cg(A,b,x1,20,1e-12,{},reduction);
    const auto bi = cfdx::core::solve_bicgstab(A,b,x2,20,1e-12,nullptr,{},reduction);
    cfdx::core::KrylovControls controls;
    controls.adaptive_restart = false;
    controls.restart_min = controls.restart_max = 4;
    controls.reduction = reduction;
    cfdx::core::LinearOperator op;
    op.size = 4;
    op.apply = [&A](const cfdx::core::Vector& in, cfdx::core::Vector& out) {
        const auto y = A.matvec(in);
        out.resize(y.size());
        for (std::size_t i = 0; i < y.size(); ++i) out(i) = y[i];
    };
    const auto gm = cfdx::core::solve_gmres(op,b,x3,4,20,1e-12,nullptr,controls);

    auto good=[](const cfdx::core::SolverResult& r,const cfdx::core::Vector& x){
        if(r.status!=cfdx::core::SolverStatus::CONVERGED) return false;
        for(std::size_t i=0;i<x.size();++i) if(std::abs(x(i)-0.5)>1e-12) return false;
        return true;
    };
    const bool ok=good(cg,x1)&&good(bi,x2)&&good(gm,x3);
    if(rank==0&&!ok) std::cerr<<"deterministic Krylov reduction wiring failed\n";
    cfdx::core::parallel::mpi_finalize();
    return ok?0:1;
}
