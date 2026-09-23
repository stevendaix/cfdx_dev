#include "cfdx/core/parallel/distributed_poisson.h"
#include "cfdx/core/solvers/scalar_diffusion.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::core::parallel;

static Mesh two_cell_channel()
{
    Mesh m;
    m.points().resize(12);
    m.points().set(0,0,0,0);   m.points().set(1,0.5,0,0);
    m.points().set(2,1,0,0);   m.points().set(3,0,1,0);
    m.points().set(4,0.5,1,0); m.points().set(5,1,1,0);
    m.points().set(6,0,0,1);   m.points().set(7,0.5,0,1);
    m.points().set(8,1,0,1);   m.points().set(9,0,1,1);
    m.points().set(10,0.5,1,1); m.points().set(11,1,1,1);

    m.faces().push_face({0,3,9,6});       // left
    m.faces().push_face({1,4,10,7});      // interface
    m.faces().push_face({2,5,11,8});      // right
    m.faces().push_face({0,1,4,3});       // cell 0 y-
    m.faces().push_face({6,9,10,7});      // cell 0 y+
    m.faces().push_face({0,6,7,1});       // cell 0 z-
    m.faces().push_face({3,4,10,9});      // cell 0 z+
    m.faces().push_face({1,2,5,4});       // cell 1 y-
    m.faces().push_face({7,10,11,8});     // cell 1 y+
    m.faces().push_face({1,7,8,2});       // cell 1 z-
    m.faces().push_face({4,5,11,10});     // cell 1 z+
    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        // Faces 0..6 belong to cell 0 except the shared interface (1),
        // while faces 2 and 7..10 belong to cell 1. Every owner index must
        // therefore be in [0, n_cells).
        const std::size_t owner =
            (f == 2 || f >= 7) ? 1u : 0u;
        m.ownership().set_owner(f, owner);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1, 0);
    m.ownership().set_neighbour(1, 1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});
    return m;
}

int main(int argc, char** argv)
{
    mpi_init(&argc, &argv);
    const int rank = mpi_rank();
    const int size = mpi_size();
    if (size != 2) {
        if (rank == 0)
            std::fprintf(stderr, "test_phase8_mpi_poisson requires exactly 2 ranks\n");
        mpi_finalize();
        return 2;
    }

    const Mesh mesh = two_cell_channel();
    auto bc = PoissonBoundaryCondition::dirichlet(mesh.n_faces());
    bc.face_values[0] = 0.0;
    bc.face_values[2] = 1.0;

    Partition partition;
    partition.n_parts = 2;
    partition.cell_rank = {0, 1};
    partition.face_owner_rank.assign(mesh.n_faces(), 0);
    partition.face_ghost_rank.assign(mesh.n_faces(), -1);
    partition.face_owner_rank[1] = 0;
    partition.face_ghost_rank[1] = 1;

    DistributedPoissonConfig cfg;
    cfg.tolerance = 1e-12;
    cfg.max_iterations = 100;

    const auto parallel = solve_poisson_mpi(mesh, bc, {0.0, 0.0}, partition, cfg);
    if (parallel.linear_result.status != SolverStatus::CONVERGED) {
        std::fprintf(stderr, "rank %d: MPI Poisson did not converge\n", rank);
        mpi_finalize();
        return 1;
    }

    // Compare against the serial canonical assembly. With two cells and
    // physical Dirichlet faces at x=0 and x=1 the cell-centred solution is
    // [1/3, 2/3], i.e. the same discrete problem on every rank.
    Vector rhs;
    const auto A = assemble_cell_diffusion_matrix(mesh, bc, 1.0, {0.0,0.0}, rhs);
    Vector serial_x(2, 0.0);
    const auto serial = solve_cg(A, rhs, serial_x, 100, 1e-12);
    if (serial.status != SolverStatus::CONVERGED)
        return 1;

    const auto& ids = parallel.solution.global_ids();
    for (std::size_t i = 0; i < ids.size(); ++i) {
        const auto gid = static_cast<std::size_t>(ids[i]);
        if (std::abs(parallel.solution(i) - serial_x(gid)) > 1e-11) {
            std::fprintf(stderr, "rank %d: serial/MPI mismatch for cell %zu\n", rank, gid);
            mpi_finalize();
            return 1;
        }
    }

    const double max_err = mpi_allreduce_max(parallel.global_residual);
    if (!(max_err < 1e-11)) {
        if (rank == 0)
            std::fprintf(stderr, "MPI Poisson residual too large: %.17g\n", max_err);
        mpi_finalize();
        return 1;
    }

    if (rank == 0)
        std::printf("MPI Poisson: 2-rank serial-equivalent validation PASS\n");
    mpi_finalize();
    return 0;
}
