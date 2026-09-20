# CFDX — Spécification Détaillée : Memory Planner & Mesh Reordering

---

## 1. MEMORY PLANNER

### 1.1 Objectif

Le Memory Planner est le composant central de l'architecture mémoire CFDX. Il doit :
- Analyser le cycle de vie de chaque buffer (birth/death)
- Réutiliser les buffers dont les durées de vie ne se chevauchent pas
- Prédire le pic mémoire avant exécution
- Tracker la mémoire en temps réel
- Minimiser le traffic VRAM via buffer reuse

### 1.2 Architecture

```cpp
// ============================================
// CORE TYPES
// ============================================

enum class MemoryLocation {
    HOST,    // RAM
    DEVICE,  // VRAM
    MANAGED  // Unified memory
};

enum class BufferType {
    FIELD,           // Champs physiques (U, p, T)
    GEOMETRY,        // Données géométriques
    TOPOLOGY,        // Connectivité
    SOLVER_VECTOR,   // Vecteurs Krylov (r, p, Ap)
    TEMPORARY,       // Buffers temporaires d'opérateurs
    AMG_LEVEL,       // Données multigrille
    HALO             // Cellules fantômes MPI
};

struct BufferID {
    uint64_t id;
    std::string name;
    
    bool operator==(const BufferID& other) const {
        return id == other.id;
    }
};

struct BufferDescriptor {
    BufferID id;
    BufferType type;
    size_t size_bytes;
    MemoryLocation location;
    Precision precision;  // FP64, FP32, FP16
    
    // Lifetime (en numéros d'opérations)
    int birth_op;
    int death_op;
    
    // Dependencies
    std::vector<BufferID> dependencies;
};

// ============================================
// LIFETIME ANALYSIS
// ============================================

class LifetimeAnalyzer {
public:
    struct LifetimeResult {
        std::vector<BufferDescriptor> buffers;
        int total_operations;
        std::map<BufferID, std::pair<int, int>> lifetimes; // birth, death
    };
    
    // Analyse le graphe d'exécution pour déterminer birth/death
    LifetimeResult analyze(const ExecutionGraph& graph) {
        LifetimeResult result;
        result.total_operations = graph.numOperations();
        
        // Pour chaque buffer, trouver première et dernière utilisation
        for (const auto& buffer : graph.allBuffers()) {
            int birth = findFirstUse(graph, buffer);
            int death = findLastUse(graph, buffer);
            
            BufferDescriptor desc;
            desc.id = buffer.id;
            desc.type = buffer.type;
            desc.size_bytes = buffer.size_bytes;
            desc.location = buffer.location;
            desc.precision = buffer.precision;
            desc.birth_op = birth;
            desc.death_op = death;
            desc.dependencies = buffer.dependencies;
            
            result.buffers.push_back(desc);
            result.lifetimes[buffer.id] = {birth, death};
        }
        
        return result;
    }
    
private:
    int findFirstUse(const ExecutionGraph& graph, const BufferID& buffer) {
        for (int op_idx = 0; op_idx < graph.numOperations(); ++op_idx) {
            const auto& op = graph.operation(op_idx);
            if (op.usesBuffer(buffer)) {
                return op_idx;
            }
        }
        return -1;
    }
    
    int findLastUse(const ExecutionGraph& graph, const BufferID& buffer) {
        for (int op_idx = graph.numOperations() - 1; op_idx >= 0; --op_idx) {
            const auto& op = graph.operation(op_idx);
            if (op.usesBuffer(buffer)) {
                return op_idx;
            }
        }
        return -1;
    }
};

// ============================================
// BUFFER REUSE (ALIASING)
// ============================================

class BufferReuseOptimizer {
public:
    struct Allocation {
        BufferID buffer_id;
        size_t offset;      // Offset dans le pool mémoire
        size_t size_bytes;
        MemoryLocation location;
    };
    
    struct MemoryPool {
        MemoryLocation location;
        size_t total_size;
        std::vector<Allocation> allocations;
    };
    
    // Algorithme principal : détermine quelles buffers peuvent partager la mémoire
    std::vector<MemoryPool> optimize(const LifetimeAnalyzer::LifetimeResult& lifetimes) {
        std::vector<MemoryPool> pools;
        
        // Séparer par location (HOST, DEVICE)
        auto host_buffers = filterByLocation(lifetimes.buffers, MemoryLocation::HOST);
        auto device_buffers = filterByLocation(lifetimes.buffers, MemoryLocation::DEVICE);
        
        // Optimiser chaque pool séparément
        pools.push_back(allocatePool(host_buffers, MemoryLocation::HOST));
        pools.push_back(allocatePool(device_buffers, MemoryLocation::DEVICE));
        
        return pools;
    }
    
private:
    std::vector<BufferDescriptor> filterByLocation(
        const std::vector<BufferDescriptor>& buffers,
        MemoryLocation location
    ) {
        std::vector<BufferDescriptor> filtered;
        for (const auto& buf : buffers) {
            if (buf.location == location) {
                filtered.push_back(buf);
            }
        }
        return filtered;
    }
    
    // Algorithme d'aliasing : First-Fit Decreasing
    MemoryPool allocatePool(
        const std::vector<BufferDescriptor>& buffers,
        MemoryLocation location
    ) {
        MemoryPool pool;
        pool.location = location;
        pool.total_size = 0;
        
        // Trier par taille décroissante (First-Fit Decreasing)
        auto sorted = buffers;
        std::sort(sorted.begin(), sorted.end(), 
            [](const auto& a, const auto& b) {
                return a.size_bytes > b.size_bytes;
            });
        
        // Pour chaque buffer, trouver un espace libre
        for (const auto& buffer : sorted) {
            size_t offset = findFreeSlot(pool, buffer);
            
            Allocation alloc;
            alloc.buffer_id = buffer.id;
            alloc.offset = offset;
            alloc.size_bytes = buffer.size_bytes;
            alloc.location = location;
            
            pool.allocations.push_back(alloc);
            pool.total_size = std::max(pool.total_size, offset + buffer.size_bytes);
        }
        
        return pool;
    }
    
    // Trouve un offset libre pour ce buffer (pas de chevauchement temporel)
    size_t findFreeSlot(const MemoryPool& pool, const BufferDescriptor& buffer) {
        size_t offset = 0;
        
        while (true) {
            bool conflict = false;
            
            // Vérifier tous les buffers déjà alloués
            for (const auto& alloc : pool.allocations) {
                // Chevauchement spatial ?
                bool spatial_overlap = 
                    (offset < alloc.offset + alloc.size_bytes) &&
                    (offset + buffer.size_bytes > alloc.offset);
                
                if (!spatial_overlap) continue;
                
                // Chevauchement temporel ?
                auto other_buffer = findBuffer(alloc.buffer_id);
                bool temporal_overlap = 
                    (buffer.birth_op <= other_buffer.death_op) &&
                    (buffer.death_op >= other_buffer.birth_op);
                
                if (spatial_overlap && temporal_overlap) {
                    // Conflict : avancer l'offset
                    offset = alloc.offset + alloc.size_bytes;
                    conflict = true;
                    break;
                }
            }
            
            if (!conflict) {
                return offset;  // Espace libre trouvé
            }
        }
    }
    
    BufferDescriptor findBuffer(BufferID id) {
        // Recherche dans la liste globale
        // ...
    }
};

// ============================================
// BUDGET PREDICTION
// ============================================

class BudgetPredictor {
public:
    struct MemoryBudget {
        // Par catégorie
        size_t topology_bytes;
        size_t geometry_bytes;
        size_t fields_bytes;
        size_t solver_vectors_bytes;
        size_t amg_bytes;
        size_t temporaries_bytes;
        size_t halo_bytes;
        
        // Totaux
        size_t peak_ram;
        size_t peak_vram;
        
        // Métriques
        double bytes_per_cell;
        double bytes_per_cell_per_iteration;
    };
    
    // Prédit le budget avant exécution
    MemoryBudget predict(
        const Mesh& mesh,
        const PhysicsConfig& physics,
        const SolverConfig& solver,
        const std::vector<MemoryPool>& pools
    ) {
        MemoryBudget budget;
        
        size_t num_cells = mesh.numCells();
        size_t num_faces = mesh.numFaces();
        
        // Topologie (compact indices)
        budget.topology_bytes = 
            num_faces * sizeof(uint32_t) * 2 +  // owner, neighbour
            num_cells * sizeof(uint32_t) +       // cell offsets
            num_faces * sizeof(uint32_t);        // face connectivity
        
        // Géométrie (policies)
        budget.geometry_bytes = estimateGeometryMemory(mesh, physics.geometry_policies);
        
        // Champs physiques
        budget.fields_bytes = num_cells * physics.num_fields * sizeof(float);  // FP32
        
        // Vecteurs solveur
        budget.solver_vectors_bytes = num_cells * solver.num_krylov_vectors * sizeof(float);
        
        // AMG
        budget.amg_bytes = estimateAMGMemory(mesh, solver.amg_config);
        
        // Temporaires (depuis pools)
        budget.temporaries_bytes = pools[1].total_size;  // DEVICE pool
        
        // Halo MPI
        budget.halo_bytes = mesh.numHaloCells() * physics.num_fields * sizeof(float);
        
        // Totaux
        budget.peak_ram = budget.topology_bytes + budget.geometry_bytes + 
                          budget.fields_bytes + budget.halo_bytes;
        budget.peak_vram = budget.solver_vectors_bytes + budget.amg_bytes + 
                           budget.temporaries_bytes;
        
        // Métriques
        budget.bytes_per_cell = static_cast<double>(budget.peak_ram + budget.peak_vram) / num_cells;
        budget.bytes_per_cell_per_iteration = estimateTrafficPerIteration(mesh, physics);
        
        return budget;
    }
    
private:
    size_t estimateGeometryMemory(const Mesh& mesh, const GeometryPolicies& policies) {
        size_t total = 0;
        size_t num_cells = mesh.numCells();
        size_t num_faces = mesh.numFaces();
        
        if (policies.cell_volume == StoragePolicy::Stored) {
            total += num_cells * sizeof(double);
        }
        if (policies.face_area == StoragePolicy::Stored) {
            total += num_faces * sizeof(Vec3<double>);
        }
        if (policies.face_centre == StoragePolicy::Cached) {
            total += num_faces * sizeof(Vec3<double>);
        }
        // ... autres quantités
        
        return total;
    }
    
    size_t estimateAMGMemory(const Mesh& mesh, const AMGConfig& config) {
        // Estimation basée sur nombre de niveaux et coarsening ratio
        size_t num_cells = mesh.numCells();
        double coarsening_ratio = 8.0;  // typique
        int num_levels = static_cast<int>(std::log(num_cells) / std::log(coarsening_ratio));
        
        size_t total = 0;
        size_t cells_at_level = num_cells;
        
        for (int level = 0; level < num_levels; ++level) {
            // Matrice CSR à ce niveau
            size_t nnz = cells_at_level * 27;  // stencil 3D typique
            total += cells_at_level * sizeof(double) + nnz * sizeof(double) + nnz * sizeof(int);
            
            cells_at_level /= coarsening_ratio;
        }
        
        return total;
    }
    
    double estimateTrafficPerIteration(const Mesh& mesh, const PhysicsConfig& physics) {
        size_t num_faces = mesh.numFaces();
        
        // Opérations typiques par itération
        // gradient : read U, write gradU
        // interpolate : read U, gradU, write Uf
        // flux : read Uf, phi, write flux
        // divergence : read flux, write divU
        
        size_t reads_per_face = 4 * sizeof(float);   // 4 champs lus
        size_t writes_per_face = 2 * sizeof(float);  // 2 champs écrits
        
        return num_faces * (reads_per_face + writes_per_face);
    }
};

// ============================================
// MEMORY LEDGER (TRACKING TEMPS RÉEL)
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
    
    void allocate(const BufferID& id, void* ptr, size_t size, MemoryLocation loc) {
        AllocationRecord record;
        record.buffer_id = id;
        record.ptr = ptr;
        record.size_bytes = size;
        record.location = loc;
        record.timestamp = std::chrono::steady_clock::now();
        
        allocations_[id] = record;
        current_usage_[loc] += size;
        peak_usage_[loc] = std::max(peak_usage_[loc], current_usage_[loc]);
    }
    
    void deallocate(const BufferID& id) {
        auto it = allocations_.find(id);
        if (it != allocations_.end()) {
            current_usage_[it->second.location] -= it->second.size_bytes;
            allocations_.erase(it);
        }
    }
    
    size_t currentUsage(MemoryLocation loc) const {
        return current_usage_.at(loc);
    }
    
    size_t peakUsage(MemoryLocation loc) const {
        return peak_usage_.at(loc);
    }
    
    void printReport() const {
        std::cout << "=== MEMORY LEDGER REPORT ===\n";
        std::cout << "Current RAM: " << currentUsage(MemoryLocation::HOST) / 1e9 << " GB\n";
        std::cout << "Peak RAM: " << peakUsage(MemoryLocation::HOST) / 1e9 << " GB\n";
        std::cout << "Current VRAM: " << currentUsage(MemoryLocation::DEVICE) / 1e9 << " GB\n";
        std::cout << "Peak VRAM: " << peakUsage(MemoryLocation::DEVICE) / 1e9 << " GB\n";
        std::cout << "Total allocations: " << allocations_.size() << "\n";
    }
    
private:
    std::map<BufferID, AllocationRecord> allocations_;
    std::map<MemoryLocation, size_t> current_usage_;
    std::map<MemoryLocation, size_t> peak_usage_;
};

// ============================================
// MEMORY PLANNER (FACADE)
// ============================================

class MemoryPlanner {
public:
    struct Plan {
        std::vector<BufferReuseOptimizer::MemoryPool> pools;
        BudgetPredictor::MemoryBudget budget;
        bool feasible;
        std::string error_message;
    };
    
    // Planifie l'allocation mémoire pour un cas donné
    Plan plan(
        const Mesh& mesh,
        const PhysicsConfig& physics,
        const SolverConfig& solver,
        const ExecutionGraph& graph,
        const HardwareConstraints& hardware
    ) {
        Plan result;
        
        // 1. Lifetime analysis
        LifetimeAnalyzer analyzer;
        auto lifetimes = analyzer.analyze(graph);
        
        // 2. Buffer reuse optimization
        BufferReuseOptimizer optimizer;
        result.pools = optimizer.optimize(lifetimes);
        
        // 3. Budget prediction
        BudgetPredictor predictor;
        result.budget = predictor.predict(mesh, physics, solver, result.pools);
        
        // 4. Feasibility check
        result.feasible = 
            result.budget.peak_ram <= hardware.ram_capacity &&
            result.budget.peak_vram <= hardware.vram_capacity;
        
        if (!result.feasible) {
            result.error_message = "Insufficient memory: " +
                std::to_string(result.budget.peak_ram / 1e9) + " GB RAM required, " +
                std::to_string(hardware.ram_capacity / 1e9) + " GB available. " +
                std::to_string(result.budget.peak_vram / 1e9) + " GB VRAM required, " +
                std::to_string(hardware.vram_capacity / 1e9) + " GB available.";
        }
        
        return result;
    }
    
    // Alloue les buffers selon le plan
    void allocate(const Plan& plan, MemoryLedger& ledger) {
        for (const auto& pool : plan.pools) {
            void* base_ptr = allocateMemory(pool.total_size, pool.location);
            
            for (const auto& alloc : pool.allocations) {
                void* buffer_ptr = static_cast<char*>(base_ptr) + alloc.offset;
                ledger.allocate(alloc.buffer_id, buffer_ptr, alloc.size_bytes, pool.location);
            }
        }
    }
    
private:
    void* allocateMemory(size_t size, MemoryLocation location) {
        switch (location) {
            case MemoryLocation::HOST:
                return malloc(size);
            case MemoryLocation::DEVICE:
                void* ptr;
                cudaMalloc(&ptr, size);
                return ptr;
            case MemoryLocation::MANAGED:
                void* ptr;
                cudaMallocManaged(&ptr, size);
                return ptr;
        }
        return nullptr;
    }
};
```

### 1.3 Exemple d'utilisation

```cpp
// Cas d'utilisation typique

int main() {
    // Charger le cas
    auto mesh = Mesh::load("case.cfdx.h5");
    auto physics = PhysicsConfig::load("case.cfdx.h5");
    auto solver = SolverConfig::load("case.cfdx.h5");
    
    // Construire le graphe d'exécution
    ExecutionGraph graph;
    graph.addOperator(Gradient("U", "gradU"));
    graph.addOperator(Interpolate("U", "Uf"));
    graph.addOperator(Flux("Uf", "phi", "flux"));
    graph.addOperator(Divergence("flux", "divU"));
    graph.addOperator(Residual("divU", "residual"));
    
    // Contraintes hardware
    HardwareConstraints hardware;
    hardware.ram_capacity = 64ULL * 1024 * 1024 * 1024;   // 64 GB
    hardware.vram_capacity = 24ULL * 1024 * 1024 * 1024;  // 24 GB
    
    // Planifier
    MemoryPlanner planner;
    auto plan = planner.plan(mesh, physics, solver, graph, hardware);
    
    if (!plan.feasible) {
        std::cerr << "Error: " << plan.error_message << "\n";
        return 1;
    }
    
    // Afficher le budget
    std::cout << "=== MEMORY BUDGET ===\n";
    std::cout << "Peak RAM: " << plan.budget.peak_ram / 1e9 << " GB\n";
    std::cout << "Peak VRAM: " << plan.budget.peak_vram / 1e9 << " GB\n";
    std::cout << "Bytes/cell: " << plan.budget.bytes_per_cell << "\n";
    std::cout << "Bytes/cell/iteration: " << plan.budget.bytes_per_cell_per_iteration << "\n";
    
    // Allouer
    MemoryLedger ledger;
    planner.allocate(plan, ledger);
    
    // Exécuter
    Executor executor;
    executor.execute(graph, ledger);
    
    // Rapport final
    ledger.printReport();
    
    return 0;
}
```

### 1.4 Métriques et Acceptance Criteria

| Métrique | Cible | Méthode de mesure |
|----------|-------|-------------------|
| **Réduction pic mémoire** | 30-50% vs sans reuse | Comparaison peak VRAM avec/sans buffer reuse |
| **Précision prédiction** | < 5% erreur | Comparaison budget prédit vs mesuré |
| **Overhead planification** | < 1% temps total | Temps de plan() vs temps d'exécution |
| **Memory ledger accuracy** | 100% | Tracking exact vs malloc/free |
| **Bytes/cell/iteration** | Benchmarké | Mesure traffic VRAM par itération |

---

## 2. MESH REORDERING

### 2.1 Objectif

Le mesh reordering améliore la localité mémoire en rapprochant les cellules connectées dans l'ordre de stockage. Cela réduit :
- Les cache misses
- La bande passante mémoire nécessaire
- Le temps d'accès aux données

### 2.2 Algorithmes

#### 2.2.1 Reverse Cuthill-McKee (RCM)

RCM réduit la bande de la matrice (bandwidth reduction), ce qui améliore la localité.

```cpp
class RCMReorderer {
public:
    struct ReorderingResult {
        std::vector<uint32_t> old_to_new;  // Mapping ancien → nouveau
        std::vector<uint32_t> new_to_old;  // Mapping nouveau → ancien
        double bandwidth_before;
        double bandwidth_after;
        double profile_before;
        double profile_after;
    };
    
    ReorderingResult reorder(const Mesh& mesh) {
        ReorderingResult result;
        
        size_t num_cells = mesh.numCells();
        result.old_to_new.resize(num_cells);
        result.new_to_old.resize(num_cells);
        
        // 1. Calculer degrés (nombre de voisins)
        std::vector<int> degree(num_cells);
        for (size_t i = 0; i < num_cells; ++i) {
            degree[i] = mesh.cellNeighbours(i).size();
        }
        
        // 2. Trouver noeud de départ (degré minimum)
        uint32_t start = findMinimumDegreeNode(degree);
        
        // 3. Cuthill-McKee (BFS depuis start)
        std::vector<uint32_t> ordering;
        std::vector<bool> visited(num_cells, false);
        
        std::queue<uint32_t> queue;
        queue.push(start);
        visited[start] = true;
        
        while (!queue.empty()) {
            uint32_t current = queue.front();
            queue.pop();
            ordering.push_back(current);
            
            // Trier voisins par degré croissant
            auto neighbours = mesh.cellNeighbours(current);
            std::sort(neighbours.begin(), neighbours.end(),
                [&degree](uint32_t a, uint32_t b) {
                    return degree[a] < degree[b];
                });
            
            for (uint32_t neighbour : neighbours) {
                if (!visited[neighbour]) {
                    visited[neighbour] = true;
                    queue.push(neighbour);
                }
            }
        }
        
        // 4. Reverse (d'où le nom RCM)
        std::reverse(ordering.begin(), ordering.end());
        
        // 5. Construire mappings
        for (uint32_t new_idx = 0; new_idx < ordering.size(); ++new_idx) {
            uint32_t old_idx = ordering[new_idx];
            result.old_to_new[old_idx] = new_idx;
            result.new_to_old[new_idx] = old_idx;
        }
        
        // 6. Calculer métriques
        result.bandwidth_before = computeBandwidth(mesh, std::vector<uint32_t>());
        result.bandwidth_after = computeBandwidth(mesh, result.old_to_new);
        result.profile_before = computeProfile(mesh, std::vector<uint32_t>());
        result.profile_after = computeProfile(mesh, result.old_to_new);
        
        return result;
    }
    
private:
    uint32_t findMinimumDegreeNode(const std::vector<int>& degree) {
        return std::distance(degree.begin(), 
            std::min_element(degree.begin(), degree.end()));
    }
    
    double computeBandwidth(const Mesh& mesh, const std::vector<uint32_t>& mapping) {
        size_t max_bandwidth = 0;
        
        for (size_t face = 0; face < mesh.numFaces(); ++face) {
            uint32_t owner = mesh.owner(face);
            uint32_t neighbour = mesh.neighbour(face);
            
            if (neighbour == INVALID_CELL) continue;  // Boundary face
            
            uint32_t owner_new = mapping.empty() ? owner : mapping[owner];
            uint32_t neighbour_new = mapping.empty() ? neighbour : mapping[neighbour];
            
            size_t bandwidth = std::abs(static_cast<int>(owner_new) - static_cast<int>(neighbour_new));
            max_bandwidth = std::max(max_bandwidth, bandwidth);
        }
        
        return static_cast<double>(max_bandwidth);
    }
    
    double computeProfile(const Mesh& mesh, const std::vector<uint32_t>& mapping) {
        size_t profile = 0;
        
        for (size_t cell = 0; cell < mesh.numCells(); ++cell) {
            uint32_t cell_new = mapping.empty() ? cell : mapping[cell];
            
            for (uint32_t neighbour : mesh.cellNeighbours(cell)) {
                uint32_t neighbour_new = mapping.empty() ? neighbour : mapping[neighbour];
                
                if (neighbour_new < cell_new) {
                    profile += cell_new - neighbour_new;
                }
            }
        }
        
        return static_cast<double>(profile);
    }
};
```

#### 2.2.2 Space-Filling Curves (SFC)

SFC (Z-order, Hilbert) préserve la localité spatiale.

```cpp
class SFCReorderer {
public:
    enum class CurveType {
        Z_ORDER,
        HILBERT
    };
    
    struct ReorderingResult {
        std::vector<uint32_t> old_to_new;
        std::vector<uint32_t> new_to_old;
        double spatial_locality_before;
        double spatial_locality_after;
    };
    
    ReorderingResult reorder(const Mesh& mesh, CurveType curve_type) {
        ReorderingResult result;
        
        size_t num_cells = mesh.numCells();
        
        // 1. Calculer coordonnées des centres
        std::vector<Vec3<double>> centres(num_cells);
        for (size_t i = 0; i < num_cells; ++i) {
            centres[i] = mesh.cellCentre(i);
        }
        
        // 2. Normaliser coordonnées [0, 1]
        auto bounds = computeBounds(centres);
        normalizeCoordinates(centres, bounds);
        
        // 3. Calculer codes SFC
        std::vector<uint64_t> codes(num_cells);
        for (size_t i = 0; i < num_cells; ++i) {
            if (curve_type == CurveType::Z_ORDER) {
                codes[i] = zOrderCode(centres[i]);
            } else {
                codes[i] = hilbertCode(centres[i]);
            }
        }
        
        // 4. Trier par code SFC
        std::vector<uint32_t> ordering(num_cells);
        std::iota(ordering.begin(), ordering.end(), 0);
        std::sort(ordering.begin(), ordering.end(),
            [&codes](uint32_t a, uint32_t b) {
                return codes[a] < codes[b];
            });
        
        // 5. Construire mappings
        result.old_to_new.resize(num_cells);
        result.new_to_old.resize(num_cells);
        for (uint32_t new_idx = 0; new_idx < ordering.size(); ++new_idx) {
            uint32_t old_idx = ordering[new_idx];
            result.old_to_new[old_idx] = new_idx;
            result.new_to_old[new_idx] = old_idx;
        }
        
        // 6. Calculer métriques
        result.spatial_locality_before = computeSpatialLocality(mesh, std::vector<uint32_t>());
        result.spatial_locality_after = computeSpatialLocality(mesh, result.old_to_new);
        
        return result;
    }
    
private:
    struct Bounds {
        Vec3<double> min;
        Vec3<double> max;
    };
    
    Bounds computeBounds(const std::vector<Vec3<double>>& points) {
        Bounds bounds;
        bounds.min = Vec3<double>(std::numeric_limits<double>::max());
        bounds.max = Vec3<double>(std::numeric_limits<double>::lowest());
        
        for (const auto& p : points) {
            bounds.min.x = std::min(bounds.min.x, p.x);
            bounds.min.y = std::min(bounds.min.y, p.y);
            bounds.min.z = std::min(bounds.min.z, p.z);
            bounds.max.x = std::max(bounds.max.x, p.x);
            bounds.max.y = std::max(bounds.max.y, p.y);
            bounds.max.z = std::max(bounds.max.z, p.z);
        }
        
        return bounds;
    }
    
    void normalizeCoordinates(std::vector<Vec3<double>>& points, const Bounds& bounds) {
        Vec3<double> range = bounds.max - bounds.min;
        
        for (auto& p : points) {
            p.x = (p.x - bounds.min.x) / range.x;
            p.y = (p.y - bounds.min.y) / range.y;
            p.z = (p.z - bounds.min.z) / range.z;
        }
    }
    
    uint64_t zOrderCode(const Vec3<double>& p) {
        // Discrétiser en grille 21 bits par dimension (63 bits total)
        uint32_t x = static_cast<uint32_t>(p.x * 2097151);  // 2^21 - 1
        uint32_t y = static_cast<uint32_t>(p.y * 2097151);
        uint32_t z = static_cast<uint32_t>(p.z * 2097151);
        
        // Interleaver les bits
        return interleave(x, y, z);
    }
    
    uint64_t interleave(uint32_t x, uint32_t y, uint32_t z) {
        uint64_t result = 0;
        
        for (int i = 0; i < 21; ++i) {
            result |= ((x >> i) & 1ULL) << (3 * i);
            result |= ((y >> i) & 1ULL) << (3 * i + 1);
            result |= ((z >> i) & 1ULL) << (3 * i + 2);
        }
        
        return result;
    }
    
    uint64_t hilbertCode(const Vec3<double>& p) {
        // Implémentation simplifiée de courbe de Hilbert 3D
        // (voir https://en.wikipedia.org/wiki/Hilbert_curve)
        uint32_t x = static_cast<uint32_t>(p.x * 1023);  // 10 bits
        uint32_t y = static_cast<uint32_t>(p.y * 1023);
        uint32_t z = static_cast<uint32_t>(p.z * 1023);
        
        return xy2d(10, x, y, z);  // Simplifié
    }
    
    uint64_t xy2d(int n, uint32_t x, uint32_t y, uint32_t z) {
        // Algorithme de conversion coordonnées → distance Hilbert
        // ... (implémentation complète nécessaire)
        return 0;  // Placeholder
    }
    
    double computeSpatialLocality(const Mesh& mesh, const std::vector<uint32_t>& mapping) {
        double total_distance = 0;
        size_t num_faces = 0;
        
        for (size_t face = 0; face < mesh.numFaces(); ++face) {
            uint32_t owner = mesh.owner(face);
            uint32_t neighbour = mesh.neighbour(face);
            
            if (neighbour == INVALID_CELL) continue;
            
            uint32_t owner_new = mapping.empty() ? owner : mapping[owner];
            uint32_t neighbour_new = mapping.empty() ? neighbour : mapping[neighbour];
            
            total_distance += std::abs(static_cast<int>(owner_new) - static_cast<int>(neighbour_new));
            num_faces++;
        }
        
        return total_distance / num_faces;  // Distance moyenne
    }
};
```

#### 2.2.3 GPUOptimized (Auto-Selection)

```cpp
class GPUOptimizedReorderer {
public:
    enum class Strategy {
        RCM,
        SFC_Z_ORDER,
        SFC_HILBERT,
        PARTITION_AWARE
    };
    
    struct ReorderingResult {
        std::vector<uint32_t> old_to_new;
        std::vector<uint32_t> new_to_old;
        Strategy selected_strategy;
        double estimated_speedup;
    };
    
    ReorderingResult reorder(const Mesh& mesh, const HardwareInfo& hardware) {
        ReorderingResult result;
        
        // Analyser le maillage
        MeshAnalysis analysis = analyzeMesh(mesh);
        
        // Choisir stratégie
        Strategy strategy = selectStrategy(analysis, hardware);
        result.selected_strategy = strategy;
        
        // Appliquer reordering
        switch (strategy) {
            case Strategy::RCM: {
                RCMReorderer rcm;
                auto rcm_result = rcm.reorder(mesh);
                result.old_to_new = rcm_result.old_to_new;
                result.new_to_old = rcm_result.new_to_old;
                result.estimated_speedup = estimateSpeedupRCM(rcm_result);
                break;
            }
            case Strategy::SFC_Z_ORDER: {
                SFCReorderer sfc;
                auto sfc_result = sfc.reorder(mesh, SFCReorderer::CurveType::Z_ORDER);
                result.old_to_new = sfc_result.old_to_new;
                result.new_to_old = sfc_result.new_to_old;
                result.estimated_speedup = estimateSpeedupSFC(sfc_result);
                break;
            }
            case Strategy::SFC_HILBERT: {
                SFCReorderer sfc;
                auto sfc_result = sfc.reorder(mesh, SFCReorderer::CurveType::HILBERT);
                result.old_to_new = sfc_result.old_to_new;
                result.new_to_old = sfc_result.new_to_old;
                result.estimated_speedup = estimateSpeedupSFC(sfc_result);
                break;
            }
            case Strategy::PARTITION_AWARE: {
                // Pour maillages MPI partitionnés
                // ... (implémentation spécifique)
                break;
            }
        }
        
        return result;
    }
    
private:
    struct MeshAnalysis {
        size_t num_cells;
        size_t num_faces;
        double avg_neighbours_per_cell;
        double bandwidth;
        double profile;
        bool is_partitioned;
        int num_partitions;
    };
    
    MeshAnalysis analyzeMesh(const Mesh& mesh) {
        MeshAnalysis analysis;
        analysis.num_cells = mesh.numCells();
        analysis.num_faces = mesh.numFaces();
        
        // Calculer degré moyen
        size_t total_neighbours = 0;
        for (size_t i = 0; i < mesh.numCells(); ++i) {
            total_neighbours += mesh.cellNeighbours(i).size();
        }
        analysis.avg_neighbours_per_cell = static_cast<double>(total_neighbours) / mesh.numCells();
        
        // Calculer bande et profil
        RCMReorderer rcm;
        analysis.bandwidth = rcm.computeBandwidth(mesh, std::vector<uint32_t>());
        analysis.profile = rcm.computeProfile(mesh, std::vector<uint32_t>());
        
        // Détecter partitionnement
        analysis.is_partitioned = mesh.hasPartitionInfo();
        analysis.num_partitions = analysis.is_partitioned ? mesh.numPartitions() : 1;
        
        return analysis;
    }
    
    Strategy selectStrategy(const MeshAnalysis& analysis, const HardwareInfo& hardware) {
        // Règles de sélection
        if (analysis.is_partitioned && analysis.num_partitions > 1) {
            return Strategy::PARTITION_AWARE;
        }
        
        if (analysis.bandwidth > 1000) {
            // Haute bande → RCM bénéfique
            return Strategy::RCM;
        }
        
        if (hardware.is_gpu && analysis.num_cells > 1000000) {
            // Gros maillage GPU → SFC Z-order (plus rapide à calculer)
            return Strategy::SFC_Z_ORDER;
        }
        
        // Défaut : RCM
        return Strategy::RCM;
    }
    
    double estimateSpeedupRCM(const RCMReorderer::ReorderingResult& result) {
        double bandwidth_reduction = 1.0 - (result.bandwidth_after / result.bandwidth_before);
        return 1.0 + bandwidth_reduction * 0.3;  // Estimation conservative
    }
    
    double estimateSpeedupSFC(const SFCReorderer::ReorderingResult& result) {
        double locality_improvement = 1.0 - (result.spatial_locality_after / result.spatial_locality_before);
        return 1.0 + locality_improvement * 0.2;  // Estimation conservative
    }
};
```

### 2.3 Application du Reordering

```cpp
class MeshReorderer {
public:
    void applyReordering(Mesh& mesh, const std::vector<uint32_t>& old_to_new) {
        // 1. Réordonner les cellules
        reorderCells(mesh, old_to_new);
        
        // 2. Mettre à jour la connectivité faces
        updateFaceConnectivity(mesh, old_to_new);
        
        // 3. Réordonner les champs
        reorderFields(mesh, old_to_new);
        
        // 4. Mettre à jour les patches frontière
        updateBoundaryPatches(mesh, old_to_new);
    }
    
private:
    void reorderCells(Mesh& mesh, const std::vector<uint32_t>& old_to_new) {
        // Réordonner cellCentres, cellVolumes, etc.
        auto& centres = mesh.cellCentres();
        auto& volumes = mesh.cellVolumes();
        
        std::vector<Vec3<double>> new_centres(centres.size());
        std::vector<double> new_volumes(volumes.size());
        
        for (size_t old_idx = 0; old_idx < centres.size(); ++old_idx) {
            uint32_t new_idx = old_to_new[old_idx];
            new_centres[new_idx] = centres[old_idx];
            new_volumes[new_idx] = volumes[old_idx];
        }
        
        centres = std::move(new_centres);
        volumes = std::move(new_volumes);
    }
    
    void updateFaceConnectivity(Mesh& mesh, const std::vector<uint32_t>& old_to_new) {
        auto& owners = mesh.faceOwners();
        auto& neighbours = mesh.faceNeighbours();
        
        for (size_t i = 0; i < owners.size(); ++i) {
            if (owners[i] != INVALID_CELL) {
                owners[i] = old_to_new[owners[i]];
            }
            if (neighbours[i] != INVALID_CELL) {
                neighbours[i] = old_to_new[neighbours[i]];
            }
        }
    }
    
    void reorderFields(Mesh& mesh, const std::vector<uint32_t>& old_to_new) {
        for (auto& field : mesh.allFields()) {
            field.reorder(old_to_new);
        }
    }
    
    void updateBoundaryPatches(Mesh& mesh, const std::vector<uint32_t>& old_to_new) {
        for (auto& patch : mesh.boundaryPatches()) {
            auto& face_indices = patch.faceIndices();
            for (auto& idx : face_indices) {
                // Les indices de faces ne changent pas, seulement les cellules
            }
        }
    }
};
```

### 2.4 Métriques et Acceptance Criteria

| Métrique | Cible | Méthode de mesure |
|----------|-------|-------------------|
| **Bandwidth reduction** | > 50% | Comparaison bandwidth avant/après RCM |
| **Distance moyenne** | > 50% réduction | Distance moyenne entre cellules connectées |
| **Cache miss rate** | > 30% réduction | Profiling cache (perf stat) |
| **SpMV speedup** | 10-30% | Benchmark SpMV avant/après |
| **Overhead reordering** | < 5% temps total | Temps de reorder vs temps de solve |

---

## 3. INTÉGRATION

```cpp
// Intégration Memory Planner + Mesh Reordering

int main() {
    // 1. Charger le maillage
    auto mesh = Mesh::load("case.cfdx.h5");
    
    // 2. Mesh reordering (optionnel mais recommandé)
    GPUOptimizedReorderer reorderer;
    auto reorder_result = reorderer.reorder(mesh, getHardwareInfo());
    
    std::cout << "Selected strategy: " << strategyName(reorder_result.selected_strategy) << "\n";
    std::cout << "Estimated speedup: " << reorder_result.estimated_speedup << "x\n";
    
    MeshReorderer applier;
    applier.applyReordering(mesh, reorder_result.old_to_new);
    
    // 3. Planifier mémoire
    auto physics = PhysicsConfig::load("case.cfdx.h5");
    auto solver = SolverConfig::load("case.cfdx.h5");
    auto graph = buildExecutionGraph(physics, solver);
    
    MemoryPlanner planner;
    auto plan = planner.plan(mesh, physics, solver, graph, getHardwareInfo());
    
    if (!plan.feasible) {
        std::cerr << "Error: " << plan.error_message << "\n";
        return 1;
    }
    
    std::cout << "Peak RAM: " << plan.budget.peak_ram / 1e9 << " GB\n";
    std::cout << "Peak VRAM: " << plan.budget.peak_vram / 1e9 << " GB\n";
    std::cout << "Bytes/cell: " << plan.budget.bytes_per_cell << "\n";
    
    // 4. Allouer et exécuter
    MemoryLedger ledger;
    planner.allocate(plan, ledger);
    
    Executor executor;
    executor.execute(graph, ledger);
    
    ledger.printReport();
    
    return 0;
}
```

---

## 4. RÉSUMÉ

### Memory Planner
- **Lifetime analysis** : Détermine birth/death de chaque buffer
- **Buffer reuse** : Algorithme First-Fit Decreasing pour aliasing
- **Budget prediction** : Prédit pic mémoire avant exécution
- **Memory ledger** : Tracking temps réel
- **Gain attendu** : 30-50% réduction pic mémoire

### Mesh Reordering
- **RCM** : Réduit bandwidth (50%+ réduction)
- **SFC** : Préserve localité spatiale (Z-order, Hilbert)
- **GPUOptimized** : Auto-selection selon maillage/hardware
- **Gain attendu** : 10-30% speedup SpMV, 30%+ réduction cache misses

Ces deux composants sont **P0** et doivent être implémentés en priorité car ils ont un impact direct sur :
- La capacité à traiter de gros maillages
- La performance GPU (localité mémoire)
- La prédictibilité des ressources nécessaires


et ca Ton analyse est excellente et change fondamentalement la donne. Tu as identifié le vrai problème : **ce n'est pas le stockage, c'est le traffic**. Voici la roadmap v4 qui intègre tes découvertes.

---

# CFDX ROADMAP v4 — Memory-Traffic-First Architecture

## Changement de paradigme

**Ancien objectif** : "4-6× moins de RAM/VRAM"  
**Nouvel objectif** : "Minimiser le peak memory et le memory traffic par cellule, avec une architecture permettant au solveur de fonctionner sur un maillage dépassant la VRAM disponible"

---

## KPI fondamentaux (6 métriques)

| Métrique | Définition | Cible |
|----------|-----------|-------|
| **Bytes/cell stored** | Mémoire statique (topology + geometry + fields) | Benchmarké |
| **Bytes/cell/iteration** | Traffic mémoire par itération (reads + writes) | Benchmarké |
| **Peak RAM/cell** | Pic mémoire CPU pendant setup/solve | Benchmarké |
| **Peak VRAM/cell** | Pic mémoire GPU pendant solve | Benchmarké |
| **CPU-GPU bytes/iteration** | Transferts PCIe par itération | < 1% du temps |
| **MPI bytes/iteration** | Communication halo par itération | Benchmarké |

---

## Architecture mémoire CFDX

```
┌─────────────────────────────────────────────────────────────┐
│                         CASE                                │
│  (mesh + physics + numerics + solver + state + metadata)   │
└──────────────────────────┬──────────────────────────────────┘
                           │
                    ┌──────▼──────┐
                    │ Memory      │
                    │ Planner     │
                    └──────┬──────┘
                           │
              ┌────────────┼────────────┐
              │            │            │
        ┌─────▼─────┐ ┌───▼───┐ ┌─────▼─────┐
        │ Topology  │ │Fields │ │ Operators │
        │ compact   │ │ mixed │ │ matrix-   │
        │ indices   │ │ prec  │ │ free opt  │
        └─────┬─────┘ └───┬───┘ └─────┬─────┘
              │            │            │
              └────────────┼────────────┘
                           │
                    ┌──────▼──────┐
                    │ Lifetime    │
                    │ DAG         │
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │ Buffer      │
                    │ Reuse       │
                    └──────┬──────┘
                           │
              ┌────────────┴────────────┐
              │                         │
        ┌─────▼─────┐           ┌──────▼──────┐
        │    RAM    │           │    VRAM     │
        │ cold/full │           │ hot working │
        │ mesh      │           │ set         │
        └─────┬─────┘           └──────┬──────┘
              │                         │
              └────────────┬────────────┘
                           │
                    ┌──────▼──────┐
                    │ GPU         │
                    │ execution   │
                    │ (fused      │
                    │  kernels)   │
                    └─────────────┘
```

---

## PHASE A — Memory Architecture (P0, Fondations)

### A1 — Memory Planner
**Objectif** : Prédire et optimiser l'allocation mémoire

```cpp
class MemoryPlanner {
    // Lifetime analysis
    struct BufferLifetime {
        BufferID id;
        size_t size;
        int birth_op;
        int death_op;
        MemoryLocation location; // RAM/VRAM
    };
    
    // Budget prediction
    struct MemoryBudget {
        size_t topology_bytes;
        size_t geometry_bytes;
        size_t fields_bytes;
        size_t solver_bytes;
        size_t amg_bytes;
        size_t temporaries_bytes;
        size_t peak_ram;
        size_t peak_vram;
    };
    
    // Buffer reuse
    std::vector<BufferAllocation> plan(
        const ExecutionGraph& graph,
        const MemoryConstraints& constraints
    );
    
    // Aliasing: buffers with non-overlapping lifetimes share memory
    void enable_aliasing(bool enable);
};
```

**Tâches** :
- [ ] A1-T01 : Lifetime analysis (birth/death per buffer)
- [ ] A1-T02 : Budget prediction (avant exécution)
- [ ] A1-T03 : Buffer reuse (aliasing non-overlapping)
- [ ] A1-T04 : Memory ledger (tracking en temps réel)
- [ ] A1-T05 : Memory profiler per array

**Acceptance** :
- Budget mémoire prédit avant lancement
- Buffer reuse réduit pic de 30-50%
- Memory ledger trace chaque allocation

---

### A2 — Compact Topology
**Objectif** : Réduire footprint topologique

```cpp
// Avant : uint64_t partout
struct TopologyOld {
    uint64_t owner;      // 8 bytes
    uint64_t neighbour;  // 8 bytes
    uint64_t* faceCells; // 8 bytes each
};

// Après : uint32_t + local/global separation
struct TopologyCompact {
    uint32_t owner;      // 4 bytes
    uint32_t neighbour;  // 4 bytes
    uint32_t* faceCells; // 4 bytes each
    
    // Global IDs seulement si nécessaire
    uint64_t* globalCellID; // optionnel
};

// Gain : 100M faces → 400 MB économisés (800 MB → 400 MB)
```

**Tâches** :
- [ ] A2-T01 : uint32_t pour indices locaux (< 4.29B entities)
- [ ] A2-T02 : Séparation local_id / global_id
- [ ] A2-T03 : CSR offsets en uint32_t
- [ ] A2-T04 : Global IDs optionnels (I/O, debug)

**Acceptance** :
- Topologie 50% plus compacte
- Global IDs seulement aux frontières MPI
- Performance identique (pas de conversion runtime)

---

### A3 — Mesh Reordering
**Objectif** : Améliorer localité mémoire

```cpp
enum class MeshOrdering {
    Original,        // Ordre du fichier
    RCM,            // Reverse Cuthill-McKee (bandwidth reduction)
    SFC,            // Space-Filling Curve (Z-order, Hilbert)
    GPUOptimized,   // Auto-choix selon hardware
    PartitionAware  // Après partitionnement MPI
};

class MeshReorderer {
    void reorder(Mesh& mesh, MeshOrdering policy);
    
    // Metrics
    struct LocalityMetrics {
        double avg_cell_distance;
        double bandwidth;
        double cache_miss_rate;
    };
    
    LocalityMetrics measure(const Mesh& mesh);
};
```

**Tâches** :
- [ ] A3-T01 : RCM implementation
- [ ] A3-T02 : SFC (Z-order) implementation
- [ ] A3-T03 : Locality metrics
- [ ] A3-T04 : GPUOptimized auto-selection
- [ ] A3-T05 : Reordering après partitionnement MPI

**Acceptance** :
- Distance moyenne entre cellules connectées réduite de 50%+
- Cache miss rate mesuré et réduit
- Performance SpMV améliorée de 10-30%

---

### A4 — Geometry Storage Policies
**Objectif** : Contrôle fin stockage géométrique

```cpp
enum class StoragePolicy {
    Stored,      // Toujours en mémoire
    Recomputed,  // Calculé à la demande
    Cached,      // Calculé une fois, mis en cache
    Compressed   // Stocké compressé
};

class GeometryManager {
    // Per-quantity policy
    void setPolicy(GeometryQuantity qty, StoragePolicy policy);
    
    // Examples
    // cellVolume       → Stored (critique, souvent utilisé)
    // faceArea         → Stored (critique)
    // faceNormal       → Recomputed (dérivé de faceArea)
    // faceCentre       → Cached (coûteux à recalculer)
    // distance         → Cached
    // nonOrthCorrection→ Cached
    
    // Cost model
    struct CostModel {
        double memory_pressure;
        double bandwidth;
        double compute_cost;
        double reuse_frequency;
    };
    
    StoragePolicy recommendPolicy(GeometryQuantity qty, CostModel model);
};
```

**Tâches** :
- [ ] A4-T01 : StoragePolicy enum + manager
- [ ] A4-T02 : Cost model (memory vs compute tradeoff)
- [ ] A4-T03 : Auto-recommendation selon profil
- [ ] A4-T04 : Benchmark policies (stored vs recomputed vs cached)

**Acceptance** :
- Mémoire géométrique réduite de 30-40% sans perte performance
- Cost model valide sur cas réels
- Policies configurables par utilisateur

---

### A5 — Kernel/Operator Fusion
**Objectif** : Réduire traffic VRAM

```cpp
// Avant : 4 kernels séparés
gradient(U) → tmp1
interpolate(tmp1) → tmp2
flux(tmp2) → tmp3
divergence(tmp3) → result

// Après : 1 kernel fusionné
gradient_interpolate_flux_divergence(U) → result

class OperatorFusion {
    // Fusion patterns
    struct FusionPattern {
        std::vector<OperatorID> operators;
        KernelID fused_kernel;
        size_t memory_saved;
        double speedup;
    };
    
    // Auto-detection
    std::vector<FusionPattern> detectFusions(const ExecutionGraph& graph);
    
    // Manual override
    void forceFusion(std::vector<OperatorID> ops);
};
```

**Tâches** :
- [ ] A5-T01 : Fusion pattern detection
- [ ] A5-T02 : Gradient + interpolation + flux + divergence fusion
- [ ] A5-T03 : Memory saved measurement
- [ ] A5-T04 : Performance benchmark (fused vs unfused)

**Acceptance** :
- Traffic VRAM réduit de 50-70% sur séquences fusionnables
- Performance améliorée de 20-40%
- Patterns détectés automatiquement

---

### A6 — GPU-Resident Execution
**Objectif** : CPU = orchestrateur, GPU = worker

```cpp
// Avant : CPU↔GPU à chaque itération
GPU solve
CPU residual_norm
GPU
CPU convergence_check
GPU

// Après : tout sur GPU, CPU orchestre seulement
GPU {
    solve
    residual_norm
    convergence_check
    → return 4 scalaires au CPU
}

class GPUResidentExecutor {
    // Tout reste sur GPU
    void execute(ExecutionGraph& graph);
    
    // CPU reçoit seulement résultats finaux
    struct GPUResult {
        double residual_l2;
        double residual_linf;
        double mass_imbalance;
        bool converged;
        int iterations;
    };
    
    GPUResult executeAndGetResult(ExecutionGraph& graph);
};
```

**Tâches** :
- [ ] A6-T01 : GPU-residual computation (L2, Linf)
- [ ] A6-T02 : GPU-convergence check
- [ ] A6-T03 : GPU-mass imbalance
- [ ] A6-T04 : Minimal CPU-GPU transfers (4-10 scalaires)

**Acceptance** :
- CPU-GPU bytes/iteration < 1% du temps
- Résidus calculés sur GPU
- Convergence vérifiée sur GPU

---

## PHASE B — Linear Algebra Memory-Optimized (P0)

### B1 — Matrix-Free Operators
**Objectif** : Supprimer stockage matrice quand possible

```cpp
// Avant : SparseMatrix CSR
class SparseMatrix {
    std::vector<double> values;
    std::vector<int> col_indices;
    std::vector<int> row_offsets;
};

// Après : Matrix-free operator
class MatrixFreeOperator {
    const Mesh& mesh;
    const Field& coefficients;
    
    // y = A(x) sans stocker A
    void apply(const Field& x, Field& y) const {
        // Calcul direct depuis topology + geometry + coefficients
        for (auto face : mesh.faces()) {
            auto owner = mesh.owner(face);
            auto neighbour = mesh.neighbour(face);
            auto coeff = coefficients[face];
            
            y[owner] += coeff * (x[neighbour] - x[owner]);
            y[neighbour] -= coeff * (x[neighbour] - x[owner]);
        }
    }
};

// Cost model : quand utiliser matrix-free vs stored
class OperatorCostModel {
    enum class Strategy { MatrixFree, Stored };
    
    Strategy recommend(
        size_t memory_pressure,
        double bandwidth,
        int operator_reuse_count
    );
};
```

**Tâches** :
- [ ] B1-T01 : MatrixFreeOperator base class
- [ ] B1-T02 : Laplacian matrix-free
- [ ] B1-T03 : Gradient matrix-free
- [ ] B1-T04 : Cost model (quand matrix-free vs stored)
- [ ] B1-T05 : Benchmark (memory saved vs performance cost)

**Acceptance** :
- Matrice CSR supprimée pour opérateurs simples
- Mémoire réduite de 30-40%
- Performance comparable (pas de dégradation > 10%)

---

### B2 — AMG Memory Control
**Objectif** : Contrôler pic mémoire AMG

```cpp
enum class AMGMemoryPolicy {
    Low,       // Coarsening agressif, setup memory minimisé
    Balanced,  // Compromis mémoire/itérations
    Fast       // Setup memory OK, minimise itérations
};

class AMGController {
    void setPolicy(AMGMemoryPolicy policy);
    
    // Hypre configuration
    void configureHypre(AMGMemoryPolicy policy) {
        switch (policy) {
            case Low:
                // Aggressive coarsening
                HYPRE_BoomerAMGSetAggressiveCoarsening(amg, 2);
                // PMIS coarsening (moins mémoire)
                HYPRE_BoomerAMGSetCoarsenType(amg, 8);
                break;
            case Balanced:
                // Standard
                HYPRE_BoomerAMGSetCoarsenType(amg, 6); // HMIS
                break;
            case Fast:
                // Maximize convergence
                HYPRE_BoomerAMGSetCoarsenType(amg, 6);
                HYPRE_BoomerAMGSetStrongThreshold(amg, 0.5);
                break;
        }
    }
    
    // Memory tracking
    struct AMGMemory {
        size_t setup_peak;
        size_t solve_peak;
        int num_levels;
        int num_iterations;
    };
    
    AMGMemory measure();
};
```

**Tâches** :
- [ ] B2-T01 : AMGMemoryPolicy enum
- [ ] B2-T02 : Hypre configuration per policy
- [ ] B2-T03 : Memory tracking (setup vs solve)
- [ ] B2-T04 : Auto-selection selon VRAM disponible
- [ ] B2-T05 : Benchmark (mémoire vs itérations)

**Acceptance** :
- Peak mémoire AMG contrôlable
- Policy Low réduit setup memory de 50%+
- Auto-selection fonctionne sur cas réels

---

## PHASE C — Mixed Precision (P1)

### C1 — Precision Policy
**Objectif** : Précision raisonnée, pas FP16 automatique

```cpp
enum class Precision {
    FP64,  // Conservation-critical
    FP32,  // Bulk fields
    FP16   // Optional buffers (après caractérisation)
};

class PrecisionPolicy {
    // Per-field precision
    void setPrecision(FieldID field, Precision prec);
    
    // Default recommendations
    static Precision recommend(FieldType type) {
        switch (type) {
            case PRESSURE:      return FP32; // ou FP64 si conservation critique
            case VELOCITY:      return FP32;
            case TEMPERATURE:   return FP32;
            case TURBULENCE:    return FP32;
            case RESIDUAL:      return FP64; // accumulation
            case COEFFICIENTS:  return FP32;
            case KRYLOV_VECTOR: return FP32;
            case TEMPORARY:     return FP16; // après validation
        }
    }
    
    // Error characterization
    struct PrecisionError {
        double absolute_error;
        double relative_error;
        bool acceptable;
    };
    
    PrecisionError measureError(Precision prec, const Field& reference);
};
```

**Tâches** :
- [ ] C1-T01 : PrecisionPolicy per field
- [ ] C1-T02 : FP32 bulk fields (U, p, T)
- [ ] C1-T03 : FP64 residuals/accumulation
- [ ] C1-T04 : FP16 optional buffers (après caractérisation)
- [ ] C1-T05 : Error characterization framework

**Acceptance** :
- Mémoire réduite de 30-50% vs FP64 pur
- Erreur < tolérance utilisateur
- FP16 seulement après validation explicite

---

## PHASE D — Working Set Manager (P2)

### D1 — VRAM/RAM Hierarchy
**Objectif** : Gérer maillages > VRAM

```cpp
class WorkingSetManager {
    // Hot data → VRAM
    void pinToVRAM(BufferID buffer);
    
    // Cold data → RAM
    void evictToRAM(BufferID buffer);
    
    // Automatic management
    struct WorkingSet {
        std::vector<BufferID> vram_resident;
        std::vector<BufferID> ram_resident;
        size_t vram_used;
        size_t ram_used;
    };
    
    WorkingSet optimize(
        const ExecutionGraph& graph,
        size_t vram_capacity,
        size_t ram_capacity
    );
    
    // Prefetching
    void prefetch(BufferID buffer);
    
    // Metrics
    struct TransferMetrics {
        size_t bytes_cpu_to_gpu;
        size_t bytes_gpu_to_cpu;
        double transfer_time;
    };
    
    TransferMetrics measure();
};
```

**Tâches** :
- [ ] D1-T01 : WorkingSetManager base
- [ ] D1-T02 : Pin/evict operations
- [ ] D1-T03 : Automatic optimization
- [ ] D1-T04 : Prefetching strategy
- [ ] D1-T05 : Transfer metrics

**Acceptance** :
- Cas > VRAM fonctionnel (avec performance degradation contrôlée)
- Transfers CPU↔GPU minimisés
- Métriques de transfer mesurées

---

## PHASE E — Memory-Aware Partitioning (P1)

### E1 — Cost-Based Partitioning
**Objectif** : Équilibrer mémoire, pas juste cellules

```cpp
class MemoryAwarePartitioner {
    // Cost model
    struct PartitionCost {
        size_t cell_memory;
        size_t face_memory;
        size_t halo_memory;
        size_t solver_memory;
        size_t amg_memory;
        
        size_t total() const {
            return cell_memory + face_memory + halo_memory 
                 + solver_memory + amg_memory;
        }
    };
    
    // Objective : minimize max(PartitionCost) across ranks
    Partitioning partition(
        const Mesh& mesh,
        int num_ranks,
        const MemoryConstraints& constraints
    );
    
    // Metrics
    struct PartitionMetrics {
        double memory_imbalance;
        double cell_imbalance;
        double halo_fraction;
    };
    
    PartitionMetrics measure(const Partitioning& part);
};
```

**Tâches** :
- [ ] E1-T01 : PartitionCost model
- [ ] E1-T02 : Memory-aware partitioning algorithm
- [ ] E1-T03 : Metrics (memory imbalance, halo fraction)
- [ ] E1-T04 : Benchmark vs cell-based partitioning

**Acceptance** :
- Memory imbalance < 10% entre ranks
- Performance MPI améliorée sur cas hétérogènes

---

## Priorités révisées (v4)

### 🔴 P0 — Maintenant (Memory Architecture)
1. **Memory Planner** (lifetime + buffer reuse)
2. **Compact indices** (uint32 + local/global separation)
3. **Mesh reordering** (RCM/SFC)
4. **SoA + local numbering**
5. **Kernel/operator fusion**
6. **GPU-resident iteration**
7. **Memory profiler per array**
8. **Geometry storage policies**
9. **Budget prediction**
10. **Memory ledger**

### 🟠 P1 — Court terme (Linear Algebra Memory)
11. **Matrix-free operators** (avec cost model)
12. **AMG memory control** (3 policies)
13. **FP32/FP64 mixed precision** (raisonnée)
14. **GPU-resident reductions**
15. **Memory-aware partitioning**
16. **Working set manager** (base)

### 🟡 P2 — Moyen terme (Advanced Memory)
17. **FP16 selective** (après caractérisation)
18. **Working set / out-of-core** (avancé)
19. **GPU compressed memory** (hardware-dependent)
20. **AMR compressed** (si besoin)
21. **Advanced kernel fusion patterns**

### 🟢 P3 — Long terme (Advanced Solvers)
22. **rGCROT** (réduction itérations)
23. **CA-Krylov** (exascale)
24. **JFNK** (nonlinéaire avancé)
25. **Adjoint/optimization**

---

## Différenciateur CFDX v4

> **Un solveur CFD conçu autour d'un budget mémoire connu par cellule et d'un working set GPU dynamique, capable de prédire avant lancement : "ce cas nécessite X GB RAM, Y GB VRAM, Z GB halo, N GB solver", puis choisit automatiquement la représentation et stratégie d'exécution optimales.**

---

## Prochaines étapes immédiates

**Cette semaine** :
1. Memory Planner (lifetime analysis + buffer reuse)
2. Compact indices (uint32 migration)
3. Memory profiler per array

**Ce mois** :
4. Mesh reordering (RCM)
5. Kernel fusion (gradient+interpolation+flux+divergence)
6. GPU-resident residuals

**Ce trimestre** :
7. Matrix-free operators (Laplacian)
8. AMG memory control (Hypre policies)
9. Mixed precision (FP32 bulk, FP64 residuals)

---

