#pragma once

#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/parallel/distributed_execution.h"
#include "cfdx/core/solvers/scalar_diffusion.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace cfdx::core::parallel {

struct DistributedPoissonConfig {
    double diffusivity = 1.0;
    std::size_t max_iterations = 2000;
    double tolerance = 1e-10;
    bool deterministic_reduction = true;
};

struct DistributedPoissonResult {
    SolverResult linear_result;
    DistributedCellField solution;
    double global_residual = std::numeric_limits<double>::infinity();
    double global_rhs_norm = 0.0;
};

namespace detail {

struct LocalPoissonRow {
    double diagonal = 0.0;
    std::vector<std::pair<std::uint64_t, double>> neighbours;
    double rhs = 0.0;
};

inline double distributed_dot(double local, bool deterministic, MPI_Comm comm) {
    return deterministic ? mpi_deterministic_sum(local, comm)
                         : mpi_allreduce_sum(local, comm);
}

inline std::vector<double> exchange_scalar_halo(
    const DistributedCellField& field,
    const DistributedHalo& halo,
    std::unordered_map<std::uint64_t, double>& values,
    MPI_Comm comm)
{
    const auto packed = exchange_distributed_cell_halo(field, halo, comm);
    std::size_t offset = 0;
    for (const auto& peer_ids : halo.recv_global_ids) {
        for (const auto gid : peer_ids) {
            if (offset >= packed.size())
                throw std::runtime_error("distributed Poisson: halo packing mismatch");
            values.emplace(gid, packed[offset++]);
        }
    }
    if (offset != packed.size())
        throw std::runtime_error("distributed Poisson: halo value count mismatch");
    return packed;
}

inline void apply_local_operator(
    const std::vector<LocalPoissonRow>& rows,
    const DistributedCellField& x,
    const std::unordered_map<std::uint64_t, double>& remote,
    std::vector<double>& y)
{
    y.assign(rows.size(), 0.0);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        double value = rows[i].diagonal * x(i);
        for (const auto& [gid, coefficient] : rows[i].neighbours) {
            double neighbour_value = 0.0;
            if (gid < x.global_size()) {
                try {
                    neighbour_value = x(x.local_index(gid));
                } catch (const std::out_of_range&) {
                    const auto it = remote.find(gid);
                    if (it == remote.end())
                        throw std::runtime_error("distributed Poisson: missing remote neighbour value");
                    neighbour_value = it->second;
                }
            }
            value += coefficient * neighbour_value;
        }
        y[i] = value;
    }
}

inline double residual_norm(
    const std::vector<double>& r,
    bool deterministic,
    MPI_Comm comm)
{
    double local = 0.0;
    for (double v : r) {
        if (!std::isfinite(v))
            throw std::runtime_error("distributed Poisson: non-finite residual");
        local += v * v;
    }
    return std::sqrt(std::max(0.0, distributed_dot(local, deterministic, comm)));
}

} // namespace detail

/// Distributed matrix-free CG for the canonical
///     -div(Gamma grad(phi)) = S
/// Poisson operator. The global mesh is available on every rank; only owned
/// unknowns are iterated. MPI interfaces are handled through explicit cell
/// halos, and all Krylov scalar products are globally reduced.
inline DistributedPoissonResult solve_poisson_mpi(
    const Mesh& mesh,
    const PoissonBoundaryCondition& boundary,
    const std::vector<double>& source,
    const Partition& partition,
    const DistributedPoissonConfig& config = {},
    MPI_Comm comm = MPI_COMM_WORLD)
{
    const int rank = mpi_rank(comm);
    const int size = mpi_size(comm);
    if (partition.n_parts != size)
        throw std::invalid_argument("solve_poisson_mpi: partition size differs from communicator");
    if (partition.cell_rank.size() != mesh.n_cells())
        throw std::invalid_argument("solve_poisson_mpi: partition cell count mismatch");
    if (source.size() != mesh.n_cells())
        throw std::invalid_argument("solve_poisson_mpi: source size mismatch");
    if (!(config.diffusivity > 0.0) || !std::isfinite(config.diffusivity))
        throw std::invalid_argument("solve_poisson_mpi: diffusivity must be positive");
    if (!(config.tolerance >= 0.0) || !std::isfinite(config.tolerance))
        throw std::invalid_argument("solve_poisson_mpi: tolerance must be finite and non-negative");
    boundary.validate(mesh.n_faces());

    const auto geometry = make_geometry_cache(mesh);
    const auto ids = owned_global_cell_ids(partition, rank);
    validate_distributed_ids(ids, mesh.n_cells(), comm);
    DistributedCellField x(mesh.n_cells(), ids, 1, "phi");
    for (std::size_t i = 0; i < x.local_size(); ++i)
        x(i) = 0.0;

    DistributedCellField work(mesh.n_cells(), ids, 1, "phi_operator");
    const auto halo = build_distributed_halo(mesh, partition, work, comm);
    std::vector<detail::LocalPoissonRow> rows(x.local_size());

    const auto& cell_faces = mesh.cells().faces();
    const auto& offsets = mesh.cells().offsets();
    const auto& own = mesh.ownership();

    bool has_dirichlet = false;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f)
        if (std::isfinite(boundary.face_values[f]) &&
            boundary.type_for_face(f) == PoissonBoundaryType::DIRICHLET)
            has_dirichlet = true;
    if (!has_dirichlet)
        throw std::invalid_argument(
            "solve_poisson_mpi: pure-Neumann problems require a dedicated gauge-aware backend");

    for (std::size_t li = 0; li < ids.size(); ++li) {
        const std::size_t c = static_cast<std::size_t>(ids[li]);
        auto& row = rows[li];
        row.rhs = source[c] * geometry.cell_volumes[c];

        for (std::size_t k = offsets[c]; k < offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const std::size_t owner = own.owner(f);
            const auto neighbour = own.neighbour(f);
            const double area = geometry.face_Sf[f].mag();
            if (!(area > 0.0) || !std::isfinite(area))
                throw std::runtime_error("solve_poisson_mpi: invalid face area");

            if (neighbour >= 0) {
                const std::size_t n = static_cast<std::size_t>(
                    owner == c ? neighbour : owner);
                const double d =
                    (geometry.cell_centres[n] - geometry.cell_centres[c]).mag();
                if (!(d > std::numeric_limits<double>::epsilon()) || !std::isfinite(d))
                    throw std::runtime_error("solve_poisson_mpi: invalid cell-centre distance");
                const double g = config.diffusivity * area / d;
                row.diagonal += g;
                row.neighbours.push_back({static_cast<std::uint64_t>(n), -g});
            } else {
                if (!std::isfinite(boundary.face_values[f]))
                    continue;
                const double d =
                    (geometry.face_centres[f] - geometry.cell_centres[c]).mag();
                if (!(d > std::numeric_limits<double>::epsilon()) || !std::isfinite(d))
                    throw std::runtime_error("solve_poisson_mpi: invalid boundary distance");
                const double g = config.diffusivity * area / d;
                if (boundary.type_for_face(f) == PoissonBoundaryType::DIRICHLET) {
                    row.diagonal += g;
                    row.rhs += g * boundary.face_values[f];
                } else if (boundary.type_for_face(f) == PoissonBoundaryType::NEUMANN) {
                    row.rhs += boundary.face_values[f] * area;
                } else {
                    throw std::invalid_argument("solve_poisson_mpi: unsupported boundary type");
                }
            }
        }
    }

    std::vector<double> rhs(rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i) rhs[i] = rows[i].rhs;

    std::vector<double> r(rows.size(), 0.0);
    std::vector<double> z(rows.size(), 0.0);
    std::vector<double> p(rows.size(), 0.0);
    std::vector<double> Ap;
    std::unordered_map<std::uint64_t, double> remote;

    auto refresh_remote = [&]() {
        remote.clear();
        detail::exchange_scalar_halo(work, halo, remote, comm);
    };
    auto apply = [&](const std::vector<double>& v, std::vector<double>& y) {
        for (std::size_t i = 0; i < work.local_size(); ++i)
            work(i) = v[i];
        refresh_remote();
        detail::apply_local_operator(rows, work, remote, y);
    };

    std::vector<double> x0(rows.size(), 0.0);
    apply(x0, Ap);
    for (std::size_t i = 0; i < rows.size(); ++i) {
        r[i] = rhs[i] - Ap[i];
        if (!(rows[i].diagonal > 0.0))
            throw std::runtime_error("solve_poisson_mpi: non-positive local diagonal");
        z[i] = r[i] / rows[i].diagonal;
        p[i] = z[i];
    }

    double local_b2 = 0.0;
    for (double v : rhs) local_b2 += v * v;
    const double b_norm = std::sqrt(std::max(0.0, detail::distributed_dot(
        local_b2, config.deterministic_reduction, comm)));
    const double tol_abs = config.tolerance * std::max(b_norm, 1e-15);

    double local_rz = 0.0;
    for (std::size_t i = 0; i < r.size(); ++i) local_rz += r[i] * z[i];
    double rz = detail::distributed_dot(local_rz, config.deterministic_reduction, comm);

    DistributedPoissonResult result{
        SolverResult{}, x, 0.0, b_norm};
    result.linear_result.residual = std::sqrt(std::max(0.0, rz));
    result.linear_result.residual_relative =
        b_norm > 0.0 ? result.linear_result.residual / b_norm : 0.0;

    if (result.linear_result.residual <= tol_abs) {
        result.linear_result.status = SolverStatus::CONVERGED;
        return result;
    }

    for (std::size_t iter = 1; iter <= config.max_iterations; ++iter) {
        std::vector<double> pv = p;
        apply(pv, Ap);

        double local_pAp = 0.0;
        for (std::size_t i = 0; i < p.size(); ++i)
            local_pAp += p[i] * Ap[i];
        const double pAp = detail::distributed_dot(
            local_pAp, config.deterministic_reduction, comm);
        if (!(pAp > 0.0) || !std::isfinite(pAp)) {
            result.linear_result.status = SolverStatus::NOT_APPLICABLE;
            result.linear_result.iterations = iter;
            return result;
        }

        const double alpha = rz / pAp;
        if (!std::isfinite(alpha)) {
            result.linear_result.status = SolverStatus::DIVERGED;
            result.linear_result.iterations = iter;
            return result;
        }

        for (std::size_t i = 0; i < x.local_size(); ++i) {
            x(i) += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        const double res = detail::residual_norm(
            r, config.deterministic_reduction, comm);
        result.linear_result.residual = res;
        result.linear_result.residual_relative =
            b_norm > 0.0 ? res / b_norm : 0.0;
        result.linear_result.iterations = iter;
        if (res <= tol_abs) {
            result.linear_result.status = SolverStatus::CONVERGED;
            break;
        }

        for (std::size_t i = 0; i < z.size(); ++i)
            z[i] = r[i] / rows[i].diagonal;
        double local_rz_new = 0.0;
        for (std::size_t i = 0; i < r.size(); ++i)
            local_rz_new += r[i] * z[i];
        const double rz_new = detail::distributed_dot(
            local_rz_new, config.deterministic_reduction, comm);
        if (!std::isfinite(rz_new) || rz == 0.0) {
            result.linear_result.status = SolverStatus::DIVERGED;
            break;
        }
        const double beta = rz_new / rz;
        for (std::size_t i = 0; i < p.size(); ++i)
            p[i] = z[i] + beta * p[i];
        rz = rz_new;
    }

    if (result.linear_result.status == SolverStatus::NOT_APPLICABLE)
        result.linear_result.status = SolverStatus::MAX_ITER_REACHED;

    // Verify the final residual with the actual distributed operator.
    std::vector<double> final_x(x.local_size());
    for (std::size_t i = 0; i < x.local_size(); ++i) final_x[i] = x(i);
    apply(final_x, Ap);
    double local_res2 = 0.0;
    for (std::size_t i = 0; i < Ap.size(); ++i) {
        const double e = rhs[i] - Ap[i];
        local_res2 += e * e;
    }
    result.global_residual = std::sqrt(std::max(0.0, detail::distributed_dot(
        local_res2, config.deterministic_reduction, comm)));
    result.linear_result.residual = result.global_residual;
    result.linear_result.residual_relative =
        b_norm > 0.0 ? result.global_residual / b_norm : 0.0;
    // The result object was initialized before the Krylov iterations to carry
    // the metadata and global RHS norm. Keep its solution synchronized with
    // the final distributed iterate; otherwise callers observe the initial
    // zero field even when the linear solve converged.
    result.solution = x;
    return result;
}

} // namespace cfdx::core::parallel
