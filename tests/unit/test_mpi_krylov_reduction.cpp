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

    using namespace cfdx::core;
    SparseMatrix A(4, 4);
    for (std::size_t i = 0; i < 4; ++i) A.push_back(i, i, 2.0);
    A.finalize();
    Vector b(4, 1.0);
    KrylovReductionPolicy reduction;
    reduction.mpi_enabled = true;
    reduction.deterministic = true;
    reduction.local_begin = static_cast<std::size_t>(rank * 2);
    reduction.local_end = reduction.local_begin + 2;

    Vector x1(4, 0.0), x2(4, 0.0), x3(4, 0.0);
    const auto cg = solve_cg(A, b, x1, 20, 1e-12, {}, reduction);
    const auto bi = solve_bicgstab(A, b, x2, 20, 1e-12, nullptr, {}, reduction);
    KrylovControls controls;
    controls.adaptive_restart = false;
    controls.restart_min = controls.restart_max = 4;
    controls.reduction = reduction;
    LinearOperator op;
    op.size = 4;
    op.apply = [&A](const Vector& in, Vector& out) {
        const auto y = A.matvec(in);
        out.resize(y.size());
        for (std::size_t i = 0; i < y.size(); ++i) out(i) = y[i];
    };
    const auto gm = solve_gmres(op, b, x3, 4, 20, 1e-12, nullptr, controls);

    auto good = [](const SolverResult& r, const Vector& x) {
        if (r.status != SolverStatus::CONVERGED) return false;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (std::abs(x(i) - 0.5) > 1e-12) return false;
        return std::isfinite(r.residual) && r.residual_relative < 1e-12;
    };
    bool ok = good(cg, x1) && good(bi, x2) && good(gm, x3);

    // Homogeneous RHS with nonzero initial guess exercises the special GMRES
    // path. Both ranks own two entries; residual certification must use the
    // same owned-range reduction as all other Krylov norms.
    Vector zero_rhs(4, 0.0);
    Vector x_zero_rhs(4, 1.0);
    const auto zero_result = solve_gmres(
        op, zero_rhs, x_zero_rhs, 4, 20, 1e-12, nullptr, controls);
    ok = ok && zero_result.status == SolverStatus::CONVERGED &&
         std::isfinite(zero_result.residual) &&
         zero_result.residual_relative <= 1e-12;
    for (std::size_t i = 0; i < x_zero_rhs.size(); ++i)
        ok = ok && std::abs(x_zero_rhs(i)) <= 1e-12;

    // Also exercise the ordinary MPI_Allreduce path, not only deterministic
    // reductions, and compare its solution against the same exact answer.
    KrylovReductionPolicy nondeterministic = reduction;
    nondeterministic.deterministic = false;
    KrylovControls controls_fast = controls;
    controls_fast.reduction = nondeterministic;
    Vector x_fast(4, 0.0);
    const auto gm_fast = solve_gmres(
        op, b, x_fast, 4, 20, 1e-12, nullptr, controls_fast);
    ok = ok && good(gm_fast, x_fast);

    if (rank == 0 && !ok)
        std::cerr << "MPI Krylov reductions or zero-RHS certification failed\n";
    cfdx::core::parallel::mpi_finalize();
    return ok ? 0 : 1;
}
