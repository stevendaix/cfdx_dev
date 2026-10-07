#include "cfdx/core/linalg/linear_solver_models.h"
#include "cfdx/core/parallel/mpi_utils.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>

using namespace cfdx::core;

namespace {

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

} // namespace

int main(int argc, char** argv) {
    cfdx::core::parallel::mpi_init(&argc, &argv);
    const int rank = cfdx::core::parallel::mpi_rank();
    const int size = cfdx::core::parallel::mpi_size();

    int rc = 0;
    try {
        require(size == 2, "N8 MPI boundary test requires exactly 2 ranks");

        // N8's current MPI boundary is explicit: the distributed Krylov
        // reduction path is supported, while the native AMG and coupled
        // Schur production preconditioners are serial-only. This test guards
        // the capability declaration so a future change cannot silently turn
        // an unsupported production path into an implied MPI qualification.
        require(krylov_model_catalog()[static_cast<std::size_t>(KrylovModel::CG)].supports_mpi,
                 "CG must remain declared MPI-capable");
        require(krylov_model_catalog()[static_cast<std::size_t>(KrylovModel::GMRES)].supports_mpi,
                 "GMRES must remain declared MPI-capable");
        require(krylov_model_catalog()[static_cast<std::size_t>(KrylovModel::FGMRES)].supports_mpi,
                 "FGMRES must remain declared MPI-capable");

        require(!preconditioner_model_catalog()[static_cast<std::size_t>(PreconditionerModel::NativeAMG)].supports_mpi,
                 "Native AMG must remain outside the N8 MPI qualification boundary");
        require(!preconditioner_model_catalog()[static_cast<std::size_t>(PreconditionerModel::SmoothedAggregationAMG)].supports_mpi,
                 "Smoothed Aggregation AMG must remain outside the N8 MPI qualification boundary");
        require(!preconditioner_model_catalog()[static_cast<std::size_t>(PreconditionerModel::CoupledBlockSchur)].supports_mpi,
                 "Coupled Block Schur must remain outside the N8 MPI qualification boundary");
        require(!preconditioner_model_catalog()[static_cast<std::size_t>(PreconditionerModel::PCD)].supports_mpi,
                 "PCD must remain outside the N8 MPI qualification boundary");
        require(!preconditioner_model_catalog()[static_cast<std::size_t>(PreconditionerModel::MGR)].supports_mpi,
                 "MGR must remain outside the N8 MPI qualification boundary");

        // Exercise the same deterministic reduction primitive used by the
        // existing MPI numerical-equivalence campaign. This is the supported
        // parallel evidence in this boundary test, not a claim of distributed
        // Schur/AMG qualification.
        const double global =
            cfdx::core::parallel::mpi_deterministic_sum(static_cast<double>(rank + 1));
        const double expected = static_cast<double>(size * (size + 1)) / 2.0;
        require(std::abs(global - expected) < 1e-14,
                "deterministic MPI reduction mismatch");

        if (rank == 0) {
            std::printf(
                "n8_mpi_boundary ranks=%d "
                "krylov_mpi=1 "
                "native_amg_mpi=0 "
                "smoothed_aggregation_amg_mpi=0 "
                "coupled_block_schur_mpi=0 "
                "pcd_mpi=0 "
                "mgr_mpi=0 "
                "deterministic_reduction=pass\n",
                size);
        }
    } catch (const std::exception& error) {
        rc = 1;
        if (rank == 0)
            std::fprintf(stderr, "N8 MPI boundary failure: %s\n", error.what());
    }

    cfdx::core::parallel::mpi_finalize();
    return rc;
}
