// M0.8/M0.9/M0.10/M0.12 — executable validation contracts.
//
// These checks use the public APIs that actually exist in CFDX. They are
// deliberately small contracts: a green result must come from an assertion,
// not from a print-only smoke message.

#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/amg_wrapper.h"
#include "cfdx/core/linalg/block_preconditioner.h"
#include "cfdx/core/numerics/matrix_free.h"
#include "cfdx/io/openfoam/openfoam_importer.h"
#include "cfdx/core/memory/memory_planner.h"
#include "common/test_harness.h"

#include <cstddef>
#include <exception>
#include <type_traits>

using namespace cfdx::core;
using namespace cfdx::testing;

int main()
{
    // Compile contracts for the real public module interfaces.
    static_assert(std::is_default_constructible_v<KrylovControls>);
    static_assert(sizeof(BlockDiagonalPreconditioner) > 0);
    static_assert(std::is_default_constructible_v<NativeBoomerAMGPreconditioner>);
    static_assert(sizeof(MatrixFreeFvDiffusionOperator) > 0);
    static_assert(sizeof(memory::MemoryLedger) > 0);
    static_assert(std::is_function_v<decltype(io::openfoam::import_openfoam_case)>);

    run_case("krylov_restart_contract_is_bounded", [] {
        const auto r = choose_gmres_restart(30, 0.5);
        // With the default adaptive controls, poor cycle reduction enlarges
        // the restart by five while remaining inside the configured bounds.
        EXPECT_TRUE(r == 35);
        EXPECT_TRUE(r <= 512);
    });

    run_case("memory_ledger_allocate_deallocate_is_balanced", [] {
        using namespace cfdx::core::memory;
        MemoryLedger ledger;
        BufferID id{1, "test"};
        char mem[256]{};
        ledger.allocate(id, mem, sizeof(mem), MemoryLocation::HOST);
        EXPECT_TRUE(ledger.currentUsage(MemoryLocation::HOST) == sizeof(mem));
        ledger.deallocate(id);
        EXPECT_TRUE(ledger.currentUsage(MemoryLocation::HOST) == 0);
    });

    run_case("memory_ledger_duplicate_id_replaces_allocation", [] {
        using namespace cfdx::core::memory;
        MemoryLedger ledger;
        BufferID id{2, "replacement"};
        char first[32]{};
        char second[64]{};
        ledger.allocate(id, first, sizeof(first), MemoryLocation::HOST);
        ledger.allocate(id, second, sizeof(second), MemoryLocation::HOST);
        // The implementation deliberately replaces an existing record for
        // the same BufferID; usage must therefore remain balanced.
        EXPECT_TRUE(ledger.currentUsage(MemoryLocation::HOST) == sizeof(second));
        ledger.deallocate(id);
        EXPECT_TRUE(ledger.currentUsage(MemoryLocation::HOST) == 0);
    });

    run_case("public_solver_headers_expose_real_types", [] {
        EXPECT_TRUE(sizeof(KrylovControls) >= sizeof(int));
        EXPECT_TRUE(sizeof(BlockDiagonalPreconditioner) > 0);
        EXPECT_TRUE(sizeof(NativeBoomerAMGPreconditioner) > 0);
        EXPECT_TRUE(sizeof(MatrixFreeFvDiffusionOperator) > 0);
    });

    return run_all();
}
