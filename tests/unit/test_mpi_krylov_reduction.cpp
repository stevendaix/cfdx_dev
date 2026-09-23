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
        if (rank == 0) std::cerr << "test_mpi_krylov_reduction requires 2 ranks\n";
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

    cfdx::core::Vector x_cg(4, 0.0);
    const auto cg = cfdx::core::solve_cg(A, b, x_cg, 20, 1e-12, {}, reduction);

    cfdx::core::Vector x_bicg(4, 0.0);
    const auto bicg = cfdx::core::solve_bicgstab(A, b, x_bicg, 20, 1e-12, nullptr, {}, reduction);

    cfdx::core::KrylovControls controls;
    controls.adaptive_restart = false;
    controls.restart_min = 4;
    controls.restart_max = 4;
    controls.reduction = reduction;
    cfdx::core::Vector x_gmres(4, 0.0);
    const auto gmres = cfdx::core::solve_gmres(A, b, x_gmres, 4, 20, 1e-12, nullptr, controls);

    const auto good = [](const cfdx::core::SolverResult& r, const cfdx::core::Vector& x) {
        if (r.status != cfdx::core::SolverStatus::CONVERGED) return false;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (std::abs(x(i) - 0.5) > 1e-12) return false;
        return true;
    };

    const bool ok = good(cg, x_cg) && good(bicg, x_bicg) && good(gmres, x_gmres);
    if (rank == 0 && !ok) {
        std::cerr << "distributed deterministic Krylov reduction wiring failed\n";
    }

    cfdx::core::parallel::mpi_finalize();
    return ok ? 0 : 1;
}
