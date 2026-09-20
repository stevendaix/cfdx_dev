#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <memory>

namespace cfdx::core::memory {

// ============================================
// Memory Location
// ============================================

enum class MemoryLocation {
    HOST,     // RAM
    DEVICE,   // VRAM
    MANAGED   // Unified memory
};

// ============================================
// Buffer Types (from v4 spec)
// ============================================

enum class BufferType {
    FIELD,
    GEOMETRY,
    TOPOLOGY,
    SOLVER_VECTOR,
    TEMPORARY,
    AMG_LEVEL,
    HALO
};

// ============================================
// Buffer ID & Descriptor
// ============================================

struct BufferID {
    uint64_t id;
    std::string name;
    bool operator==(const BufferID& other) const { return id == other.id; }
};

// Custom hash wrapper (avoids std namespace shadowing)
struct BufferIDHash {
    size_t operator()(const BufferID& b) const {
        return std::hash<uint64_t>{}(b.id);
    }
};

struct BufferDescriptor {
    BufferID id;
    BufferType type;
    size_t size_bytes;
    MemoryLocation location;
    // lifetime operations
    int birth_op;
    int death_op;
    std::vector<BufferID> dependencies;
};

// ============================================
// Lifetime Analyzer (A1-T01)
// ============================================

class ExecutionGraph;  // Forward declaration

class LifetimeAnalyzer {
public:
    struct LifetimeResult {
        std::vector<BufferDescriptor> buffers;
        int total_operations;
        std::map<BufferID, std::pair<int, int>> lifetimes;
    };
    LifetimeResult analyze(const ExecutionGraph& graph);
};

// ============================================
// Buffer Reuse Optimizer (A1-T03)
// ============================================

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
    std::vector<MemoryPool> optimize(const LifetimeAnalyzer::LifetimeResult& lifetimes);
};

// ============================================
// Budget Predictor (A1-T02 / A5)
// ============================================

class MemoryPlanner;  // Forward

struct MemoryBudget {
    size_t topology_bytes;
    size_t geometry_bytes;
    size_t fields_bytes;
    size_t solver_vectors_bytes;
    size_t amg_bytes;
    size_t temporaries_bytes;
    size_t halo_bytes;
    size_t peak_ram;
    size_t peak_vram;
    double bytes_per_cell;
    double bytes_per_cell_per_iteration;
};

// ============================================
// Memory Ledger (A1-T04)
// ============================================

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
};

// ============================================
// Memory Planner Facade (A1)
// ============================================

class MemoryPlanner {
public:
    struct Plan {
        std::vector<BufferReuseOptimizer::MemoryPool> pools;
        MemoryBudget budget;
        bool feasible;
        std::string error_message;
    };
    Plan plan(const std::vector<BufferDescriptor>& buffers,
               int total_operations,
               size_t num_cells = 1000,
               size_t num_faces = 200,
               int num_fields = 3,
               int num_solver_vectors = 5);
};

// ============================================
// Mesh Reordering (A3-T01 / A3-T02 / A3-T03)
// ============================================

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
