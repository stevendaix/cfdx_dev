#include "cfdx/core/numerics/backend_equivalence.h"
#include "cfdx/core/parallel/distributed_poisson.h"
#include "cfdx/core/parallel/mpi_utils.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::core::numerics;
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

    m.faces().push_face({0,3,9,6});
    m.faces().push_face({1,4,10,7});
    m.faces().push_face({2,5,11,8});
    m.faces().push_face({0,1,4,3});
    m.faces().push_face({6,9,10,7});
    m.faces().push_face({0,6,7,1});
    m.faces().push_face({3,4,10,9});
    m.faces().push_face({1,2,5,4});
    m.faces().push_face({7,10,11,8});
    m.faces().push_face({1,7,8,2});
    m.faces().push_face({4,5,11,10});
    m.ownership().resize(11);
    for (std::size_t f = 0; f < 11; ++f) {
        const std::size_t owner = (f == 2 || f >= 7) ? 1u : 0u;
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
            std::fprintf(stderr, "test_n13_mpi_equivalence requires exactly 2 MPI ranks\n");
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
    const auto result = solve_poisson_mpi(mesh, bc, {0.0, 0.0}, partition, cfg);

    bool local_ok = result.linear_result.status == SolverStatus::CONVERGED;
    std::vector<double> local_values;
    for (std::size_t i = 0; i < result.solution.local_size(); ++i)
        local_values.push_back(result.solution(i));

    const auto gathered = mpi_allgather(local_values);
    const std::vector<double> serial_reference{0.25, 0.75};

    if (rank == 0) {
        const auto metrics = compare_vectors(serial_reference, gathered);
        local_ok = local_ok &&
                   numerically_equivalent(metrics, 1e-12, 1e-12);
        if (!local_ok)
            std::fprintf(stderr, "N13 serial/MPI equivalence failed: L2rel=%.17g Linfrel=%.17g\n",
                         metrics.l2_relative, metrics.max_relative);
        else
            std::printf("N13 serial/MPI equivalence PASS L2rel=%.17g Linfrel=%.17g\n",
                        metrics.l2_relative, metrics.max_relative);
    }

    const int global_ok = mpi_allreduce_min(local_ok ? 1 : 0);
    if (!global_ok) {
        mpi_finalize();
        return 1;
    }

    const double local_sum = static_cast<double>(rank + 1);
    const double normal_sum = mpi_allreduce_sum(local_sum);
    const double deterministic_sum = mpi_deterministic_sum(local_sum);
    const double reduction_delta = std::abs(normal_sum - deterministic_sum);

    if (rank == 0 && reduction_delta > 1e-14)
        std::fprintf(stderr, "N13 reduction equivalence failed: %.17g\n", reduction_delta);

    const int reduction_ok = mpi_allreduce_min(
        (reduction_delta <= 1e-14) ? 1 : 0);
    mpi_finalize();
    return (reduction_ok && global_ok) ? 0 : 1;
}
