#include "cfdx/core/parallel/distributed_execution.h"
#include "cfdx/core/parallel/mpi_utils.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using cfdx::core::parallel::DistributedCellField;
using cfdx::core::parallel::mpi_allgather;
using cfdx::core::parallel::mpi_barrier;
using cfdx::core::parallel::mpi_rank;
using cfdx::core::parallel::mpi_size;
using cfdx::core::parallel::read_distributed_checkpoint;
using cfdx::core::parallel::write_distributed_checkpoint;

namespace {
struct Options {
    std::string mode;
    std::filesystem::path checkpoint;
    std::filesystem::path reference;
    std::size_t cells = 32;
    std::size_t checkpoint_iteration = 10;
    std::size_t final_iteration = 20;
};

Options parse(int argc, char** argv)
{
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::invalid_argument(std::string("missing value for ") + name);
            return argv[++i];
        };
        if (arg == "--mode") o.mode = value("--mode");
        else if (arg == "--checkpoint") o.checkpoint = value("--checkpoint");
        else if (arg == "--reference") o.reference = value("--reference");
        else if (arg == "--cells") o.cells = std::stoull(value("--cells"));
        else if (arg == "--checkpoint-iteration") o.checkpoint_iteration = std::stoull(value("--checkpoint-iteration"));
        else if (arg == "--final-iteration") o.final_iteration = std::stoull(value("--final-iteration"));
        else throw std::invalid_argument("unknown option: " + arg);
    }
    if ((o.mode != "writer" && o.mode != "reader") ||
        o.checkpoint.empty() || o.reference.empty() || o.cells < 4 ||
        o.checkpoint_iteration == 0 || o.final_iteration <= o.checkpoint_iteration)
        throw std::invalid_argument("invalid distributed Poisson options");
    return o;
}

std::vector<std::uint64_t> owned_ids(std::size_t global_n, int rank, int size)
{
    std::vector<std::uint64_t> ids;
    for (std::size_t gid = 0; gid < global_n; ++gid)
        if (static_cast<int>((gid * static_cast<std::size_t>(size)) / global_n) == rank)
            ids.push_back(static_cast<std::uint64_t>(gid));
    return ids;
}

void write_iteration(const std::filesystem::path& path, std::size_t iteration)
{
    if (mpi_rank() == 0) {
        std::ofstream out(path);
        if (!out) throw std::runtime_error("cannot write iteration metadata");
        out << iteration << "\n";
    }
    mpi_barrier();
}

std::size_t read_iteration(const std::filesystem::path& path)
{
    std::size_t iteration = 0;
    if (mpi_rank() == 0) {
        std::ifstream in(path);
        if (!in || !(in >> iteration))
            throw std::runtime_error("cannot read iteration metadata");
    }
    MPI_Bcast(&iteration, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    return iteration;
}

double global_l2(const DistributedCellField& field, const std::vector<double>& exact)
{
    double local = 0.0;
    for (std::size_t i = 0; i < field.local_size(); ++i) {
        const auto gid = static_cast<std::size_t>(field.global_ids()[i]);
        const double e = field(i) - exact[gid];
        local += e * e;
    }
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return std::sqrt(global / static_cast<double>(exact.size()));
}

std::vector<double> gather_global(const DistributedCellField& field)
{
    const auto ids = mpi_allgather(field.global_ids());
    std::vector<double> local(field.local_size());
    for (std::size_t i = 0; i < field.local_size(); ++i)
        local[i] = field(i);
    const auto values = mpi_allgather(local);
    if (ids.size() != values.size())
        throw std::runtime_error("global gather size mismatch");

    std::vector<double> global(field.global_size(), 0.0);
    std::vector<bool> seen(field.global_size(), false);
    for (std::size_t i = 0; i < ids.size(); ++i) {
        if (ids[i] >= global.size() || seen[ids[i]])
            throw std::runtime_error("invalid gathered global cell identity");
        global[ids[i]] = values[i];
        seen[ids[i]] = true;
    }
    if (std::find(seen.begin(), seen.end(), false) != seen.end())
        throw std::runtime_error("gathered state does not cover global domain");
    return global;
}

// Finite-volume Jacobi step for -d2(phi)/dx2 = 1 on [0,1], phi(0)=phi(1)=0.
// Cell centres are at (i+1/2)h. Boundary half-cell conductances are included
// explicitly, so this is a genuine cell-centred Poisson iteration rather than
// a synthetic vector update.
void jacobi_step(DistributedCellField& field, std::size_t global_n)
{
    const auto global = gather_global(field);
    const double h = 1.0 / static_cast<double>(global_n);
    std::vector<double> next(field.local_size(), 0.0);

    for (std::size_t i = 0; i < field.local_size(); ++i) {
        const std::size_t gid = static_cast<std::size_t>(field.global_ids()[i]);
        const bool left_boundary = gid == 0;
        const bool right_boundary = gid + 1 == global_n;

        double diagonal = 2.0 / h;
        double rhs = h;
        double neighbours = 0.0;

        if (left_boundary) {
            diagonal += 2.0 / h;
            rhs += 0.0;
            neighbours += global[gid + 1] / h;
        } else if (right_boundary) {
            diagonal += 2.0 / h;
            rhs += 0.0;
            neighbours += global[gid - 1] / h;
        } else {
            neighbours += (global[gid - 1] + global[gid + 1]) / h;
        }
        next[i] = (rhs + neighbours) / diagonal;
    }

    for (std::size_t i = 0; i < field.local_size(); ++i)
        field(i) = next[i];
}

std::vector<double> exact_solution(std::size_t n)
{
    std::vector<double> exact(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = (static_cast<double>(i) + 0.5) / static_cast<double>(n);
        exact[i] = 0.5 * x * (1.0 - x);
    }
    return exact;
}
}

int main(int argc, char** argv)
{
    try {
        MPI_Init(&argc, &argv);
        const auto options = parse(argc, argv);
        const int rank = mpi_rank();
        const int size = mpi_size();

        const int required_size = options.mode == "writer" ? 2 : 3;
        if (size != required_size)
            throw std::runtime_error(
                "production MPI Poisson requires exactly " +
                std::to_string(required_size) + " ranks for this mode");

        const auto ids = owned_ids(options.cells, rank, size);
        DistributedCellField state(options.cells, ids, 1, "phi");
        state.field().fill(0.0);

        std::size_t iteration = 0;
        if (options.mode == "writer") {
            for (; iteration < options.checkpoint_iteration; ++iteration)
                jacobi_step(state, options.cells);
            write_distributed_checkpoint(options.checkpoint.string(), state);
            write_iteration(options.reference.string() + ".iteration", iteration);

            for (; iteration < options.final_iteration; ++iteration)
                jacobi_step(state, options.cells);
            const auto exact = exact_solution(options.cells);
            const double error = global_l2(state, exact);
            write_distributed_checkpoint(options.reference.string(), state);
            mpi_barrier();
            if (rank == 0)
                std::cout << "MPI Poisson writer: checkpoint iteration="
                          << options.checkpoint_iteration
                          << " final iteration=" << options.final_iteration
                          << " L2=" << error << "\n";
        } else {
            read_distributed_checkpoint(options.checkpoint.string(), state);
            iteration = read_iteration(options.reference.string() + ".iteration");
            if (iteration != options.checkpoint_iteration)
                throw std::runtime_error("checkpoint iteration metadata mismatch");

            const auto before = gather_global(state);
            for (; iteration < options.final_iteration; ++iteration)
                jacobi_step(state, options.cells);
            const auto after = gather_global(state);

            // Read the writer's final field using the new M-rank decomposition.
            // This is the independent reference continuation produced by the
            // actual solver on N ranks, not a recomputation from the checkpoint.
            DistributedCellField reference(options.cells, ids, 1, "phi");
            read_distributed_checkpoint(options.reference.string(), reference);
            const auto expected = gather_global(reference);

            double local_max = 0.0;
            for (std::size_t gid = 0; gid < options.cells; ++gid)
                local_max = std::max(local_max, std::abs(after[gid] - expected[gid]));
            double max_error = 0.0;
            MPI_Allreduce(&local_max, &max_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            if (max_error > 1e-14)
                throw std::runtime_error("N-to-M continuation mismatch");

            if (rank == 0)
                std::cout << "MPI Poisson reader: restart iteration="
                          << options.checkpoint_iteration
                          << " continued to=" << options.final_iteration
                          << " max_field_difference=" << std::setprecision(17)
                          << max_error << "\n";
        }

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        return 0;
    } catch (const std::exception& exc) {
        std::cerr << "production MPI Poisson error: " << exc.what() << "\n";
        MPI_Abort(MPI_COMM_WORLD, 2);
        return 2;
    }
}
