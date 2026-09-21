#include "cfdx/core/memory/memory_planner.h"

#include <cassert>
#include <cstdio>
#include <vector>

using namespace cfdx::core::memory;

int main() {
    int passed = 0;
    int total = 0;

    // 1. Backward-compatible budget prediction.
    {
        MemoryPlanner planner;
        const auto plan = planner.plan({}, 10, 500, 1500, 3, 5);
        assert(plan.feasible);
        assert(plan.budget.bytes_per_cell > 0.0);
        assert(plan.budget.bytes_per_cell_per_iteration > 0.0);
        assert(plan.budget.peak_ram > 0);
        assert(plan.budget.peak_vram > 0);
        assert(plan.budget.topology_bytes > 0);
        ++passed;
        ++total;
    }

    // 2. Descriptor-driven lifetime reuse: two non-overlapping temporary
    // buffers must share storage while an overlapping buffer cannot.
    {
        std::vector<BufferDescriptor> buffers{
            {{1, "tmp_a"}, BufferType::TEMPORARY, 1024, MemoryLocation::DEVICE, 0, 2, {}},
            {{2, "tmp_b"}, BufferType::TEMPORARY, 1024, MemoryLocation::DEVICE, 3, 5, {}},
            {{3, "tmp_c"}, BufferType::TEMPORARY, 512, MemoryLocation::DEVICE, 1, 4, {}},
        };

        MemoryPlanner planner;
        const auto plan = planner.plan(buffers, 5, 10, 20, 1, 1);
        assert(plan.feasible);
        assert(plan.pools.size() == 1);
        const auto& pool = plan.pools.front();
        assert(pool.location == MemoryLocation::DEVICE);
        assert(pool.total_size == 1536);
        assert(pool.allocations.size() == 3);

        const auto find = [&](uint64_t id) -> const BufferReuseOptimizer::Allocation& {
            for (const auto& allocation : pool.allocations) {
                if (allocation.buffer_id.id == id) return allocation;
            }
            assert(false);
            return pool.allocations.front();
        };
        const auto& a = find(1);
        const auto& b = find(2);
        const auto& c = find(3);
        assert(a.offset == b.offset);
        assert(a.offset != c.offset);
        assert(plan.budget.temporaries_bytes == 2560);
        assert(plan.budget.peak_vram == 1536);
        ++passed;
        ++total;
    }

    // 3. Invalid lifetimes must fail instead of silently producing a plan.
    {
        MemoryPlanner planner;
        const std::vector<BufferDescriptor> invalid{
            {{7, "bad"}, BufferType::FIELD, 128, MemoryLocation::HOST, 4, 2, {}}
        };
        const auto plan = planner.plan(invalid, 5, 10, 20, 1, 1);
        assert(!plan.feasible);
        assert(!plan.error_message.empty());
        ++passed;
        ++total;
    }

    // 4. Memory ledger is instance-local; one ledger must not affect another.
    {
        MemoryLedger first;
        MemoryLedger second;
        BufferID id{42, "test_kpi"};
        char mem[1024];

        first.allocate(id, mem, sizeof(mem), MemoryLocation::HOST);
        assert(first.currentUsage(MemoryLocation::HOST) == sizeof(mem));
        assert(first.peakUsage(MemoryLocation::HOST) == sizeof(mem));
        assert(second.currentUsage(MemoryLocation::HOST) == 0);
        assert(second.peakUsage(MemoryLocation::HOST) == 0);

        first.deallocate(id);
        assert(first.currentUsage(MemoryLocation::HOST) == 0);
        assert(first.peakUsage(MemoryLocation::HOST) == sizeof(mem));
        ++passed;
        ++total;
    }

    // 5. RCM must handle empty and disconnected meshes deterministically.
    {
        RCMReorderer rcm;
        const auto empty = rcm.reorder({}, {}, 0);
        assert(empty.old_to_new.empty());

        const std::vector<uint32_t> owner{0, 2};
        const std::vector<uint32_t> neighbour{1, 3};
        const auto disconnected = rcm.reorder(owner, neighbour, 4);
        assert(disconnected.old_to_new.size() == 4);
        assert(disconnected.new_to_old.size() == 4);
        for (size_t i = 0; i < 4; ++i) {
            assert(disconnected.old_to_new[disconnected.new_to_old[i]] == i);
        }
        ++passed;
        ++total;
    }

    // 6. Applying an RCM permutation must preserve boundary sentinels.
    {
        MeshReorderer reorderer;
        std::vector<uint32_t> owner{0, 0};
        std::vector<uint32_t> neighbour{1, 0xFFFFFFFFU};
        const auto mapping = reorderer.reorder(
            owner, neighbour, 1, MeshOrdering::RCM);
        reorderer.applyReordering(owner, neighbour, mapping.old_to_new);
        assert(neighbour[1] == 0xFFFFFFFFU);
        ++passed;
        ++total;
    }

    std::printf("Memory planner validation: %d/%d PASS\n", passed, total);
    return passed == total ? 0 : 1;
}