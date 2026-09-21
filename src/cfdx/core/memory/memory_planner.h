#pragma once

#include <chrono>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace cfdx::core::memory {

enum class MemoryLocation {
    HOST,
    DEVICE,
    MANAGED
};

enum class BufferType {
    FIELD,
    GEOMETRY,
    TOPOLOGY,
    SOLVER_VECTOR,
    TEMPORARY,
    AMG_LEVEL,
    HALO
};

struct BufferID {
    uint64_t id;
    std::string name;
    bool operator==(const BufferID& other) const { return id == other.id; }
};

struct BufferIDHash {
    size_t operator()(const BufferID& b) const noexcept {
        return std::hash<uint64_t>{}(b.id);
    }
};

struct BufferDescriptor {
    BufferID id;
    BufferType type;
    size_t size_bytes;
    MemoryLocation location;
    int birth_op;
    int death_op;
    std::vector<BufferID> dependencies;
};

class ExecutionGraph;

class LifetimeAnalyzer {
public:
    struct LifetimeResult {
        std::vector<BufferDescriptor> buffers;
        int total_operations = 0;
        std::vector<std::pair<BufferID, std::pair<int, int>>> lifetimes;
    };
    LifetimeResult analyze(const ExecutionGraph& graph);
};

class BufferReuseOptimizer {
public:
    struct Allocation {
        BufferID buffer_id;
        size_t offset;
        size_t size_bytes;
        MemoryLocation location;
    };

    struct MemoryPool {
        MemoryLocation location;
        size_t total_size;
        std::vector<Allocation> allocations;
    };

    std::vector<MemoryPool> optimize(
        const LifetimeAnalyzer::LifetimeResult& lifetimes) const;
};

struct MemoryBudget {
    size_t topology_bytes = 0;
    size_t geometry_bytes = 0;
    size_t fields_bytes = 0;
    size_t solver_vectors_bytes = 0;
    size_t amg_bytes = 0;
    size_t temporaries_bytes = 0;
    size_t halo_bytes = 0;
    size_t peak_ram = 0;
    size_t peak_vram = 0;
    double bytes_per_cell = 0.0;
    double bytes_per_cell_per_iteration = 0.0;
};

class MemoryLedger {
public:
    struct AllocationRecord {
        BufferID buffer_id;
        void* ptr;
        size_t size_bytes;
        MemoryLocation location;
        std::chrono::steady_clock::time_point timestamp;
    };

    void allocate(const BufferID& id, void* ptr, size_t size, MemoryLocation loc);
    void deallocate(const BufferID& id);
    size_t currentUsage(MemoryLocation loc) const;
    size_t peakUsage(MemoryLocation loc) const;
    void printReport() const;

private:
    std::unordered_map<BufferID, AllocationRecord, BufferIDHash> allocations_;
    std::unordered_map<MemoryLocation, size_t> current_usage_;
    std::unordered_map<MemoryLocation, size_t> peak_usage_;
};

class MemoryPlanner {
public:
    struct Plan {
        std::vector<BufferReuseOptimizer::MemoryPool> pools;
        MemoryBudget budget;
        bool feasible = false;
        std::string error_message;
    };

    Plan plan(const std::vector<BufferDescriptor>& buffers,
              int total_operations,
              size_t num_cells = 1000,
              size_t num_faces = 200,
              int num_fields = 3,
              int num_solver_vectors = 5) const;
};

class RCMReorderer {
public:
    struct ReorderingResult {
        std::vector<uint32_t> old_to_new;
        std::vector<uint32_t> new_to_old;
        double bandwidth_before;
        double bandwidth_after;
        double profile_before;
        double profile_after;
    };

    ReorderingResult reorder(const std::vector<uint32_t>& owner,
                             const std::vector<uint32_t>& neighbour,
                             size_t num_cells);

private:
    double computeBandwidth(const std::vector<uint32_t>& owner,
                            const std::vector<uint32_t>& neighbour,
                            const std::vector<uint32_t>& mapping);
    double computeProfile(const std::vector<uint32_t>& owner,
                          const std::vector<uint32_t>& neighbour,
                          const std::vector<uint32_t>& mapping);
};

enum class MeshOrdering {
    Original,
    RCM,
    SFC,
    GPUOptimized,
    PartitionAware
};

class MeshReorderer {
public:
    struct ReorderingResult {
        std::vector<uint32_t> old_to_new;
        std::vector<uint32_t> new_to_old;
        MeshOrdering selected_strategy;
        double estimated_speedup;
    };

    ReorderingResult reorder(const std::vector<uint32_t>& owner,
                             const std::vector<uint32_t>& neighbour,
                             int num_partitions,
                             MeshOrdering policy);

    void applyReordering(std::vector<uint32_t>& owner,
                         std::vector<uint32_t>& neighbour,
                         const std::vector<uint32_t>& old_to_new) const;
};

} // namespace cfdx::core::memory
