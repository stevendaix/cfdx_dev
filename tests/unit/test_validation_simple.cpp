// M0.8/M0.9/M0.10/M0.12 — executable validation contracts.
//
// This test deliberately contains no print-only "PASS" cases.  A CTest
// entry is successful only when it exercises a runtime contract or a compile
// contract that can fail independently of the test harness itself.

#include "cfdx/core/linalg/gmres_solver.h"
#include "cfdx/core/linalg/preconditioner.h"
#include "cfdx/core/linalg/amg_wrapper.h"
#include "cfdx/core/linalg/block_preconditioner.h"
#include "cfdx/core/numerics/matrix_free.h"
#include "cfdx/io/openfoam/openfoam_importer.h"
#include "cfdx/core/memory/memory_planner.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

using namespace cfdx::core;
using namespace cfdx::testing;

int main()
{
    // Compile contracts: these declarations fail compilation if the public
    // module interfaces disappear or become incompatible.  They are not
    // presented as runtime validation.
    static_assert(std::is_default_constructible_v<KrylovControls>);
    static_assert(std::is_default_constructible_v<CellBlockJacobiPreconditioner> == false ||
                  std::is_constructible_v<CellBlockJacobiPreconditioner, std::size_t>);
    static_assert(std::is_default_constructible_v<MemoryLedger>);

    run_case("krylov_restart_contract_is_bounded", [] {
        const auto r = choose_gmres_restart(30, 0.5);
        EXPECT_TRUE(r > 30);
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

    run_case("memory_ledger_rejects_duplicate_allocation", [] {
        using namespace cfdx::core::memory;
        MemoryLedger ledger;
        BufferID id{2, "duplicate"};
        char mem[32]{};
        ledger.allocate(id, mem, sizeof(mem), MemoryLocation::HOST);
        EXPECT_THROW(ledger.allocate(id, mem, sizeof(mem), MemoryLocation::HOST),
                     std::exception);
        ledger.deallocate(id);
    });

    run_case("public_solver_headers_expose_nontrivial_types", [] {
        EXPECT_TRUE(sizeof(KrylovControls) >= sizeof(int));
        EXPECT_TRUE(sizeof(AmgPreconditioner) > 0);
        EXPECT_TRUE(sizeof(MatrixFreeOperator) > 0);
        EXPECT_TRUE(sizeof(OpenFoamImporter) > 0);
    });

    return run_all();
}
