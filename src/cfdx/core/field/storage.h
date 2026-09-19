// M0.4-T03 — StorageHandle (Host / Device / WorkingSet)
//
// Spécification CFDX v0.7 §39, §40, §41 :
//   Un Field ne doit pas contenir directement :
//     device_ptr
//   ou une logique CUDA spécifique.
//
//   Architecture :
//     Field
//       │
//       ▼
//     StorageHandle
//       │
//       ├── HostStorage
//       ├── DeviceStorage
//       └── WorkingSetStorage
//
//   Storage states (§41) :
//     HOST_ONLY, DEVICE_ONLY, SYNCHRONIZED,
//     HOST_DIRTY, DEVICE_DIRTY

#pragma once

#include <vector>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace cfdx {
namespace core {

enum class StorageState : std::uint8_t {
    HOST_ONLY = 0,
    DEVICE_ONLY,
    SYNCHRONIZED,
    HOST_DIRTY,
    DEVICE_DIRTY
};

inline const char* to_string(StorageState s) {
    switch (s) {
        case StorageState::HOST_ONLY:     return "host_only";
        case StorageState::DEVICE_ONLY:  return "device_only";
        case StorageState::SYNCHRONIZED: return "synchronized";
        case StorageState::HOST_DIRTY:    return "host_dirty";
        case StorageState::DEVICE_DIRTY:  return "device_dirty";
        default:                          return "unknown";
    }
}

// --- Host storage (§40) ---
// Le stockage CPU est la représentation persistante principale.
struct HostStorage {
    std::vector<double> data;

    std::size_t size() const noexcept { return data.size(); }
    bool empty() const noexcept { return data.empty(); }
    void resize(std::size_t n) { data.resize(n, 0.0); }
    void clear() { data.clear(); }
    double* data_ptr() noexcept { return data.data(); }
    const double* data_ptr() const noexcept { return data.data(); }
};

// --- Device storage (§39) ---
// Placeholder pour le backend GPU (M0.12+).
// Aucune logique CUDA ici — le DeviceStorage est une abstraction.
struct DeviceStorage {
    std::size_t size_bytes = 0;
    // Le pointeur device est géré par le Runtime (M0.12).
    // Pour M0, DeviceStorage est vide.

    std::size_t size() const noexcept { return size_bytes / sizeof(double); }
    bool available() const noexcept { return false; }  // GPU non disponible en M0
};

// --- Working set storage (§47) ---
// Pour le mode GPU out-of-core (M0.12+).
struct WorkingSetStorage {
    std::size_t capacity = 0;
    std::size_t used = 0;

    bool available() const noexcept { return false; }
};

// --- StorageHandle (§39) ---
// Un Field utilise un StorageHandle pour gérer le stockage,
// indépendamment du backend matériel.
class StorageHandle {
public:
    StorageHandle() : state_(StorageState::HOST_ONLY) {}

    // --- Accès ---

    HostStorage& host() { return host_; }
    const HostStorage& host() const { return host_; }

    DeviceStorage& device() { return device_; }
    const DeviceStorage& device() const { return device_; }

    StorageState state() const noexcept { return state_; }

    // --- Gestion des états (§41) ---

    void set_state(StorageState s) noexcept { state_ = s; }

    // Marque le host comme dirty (le device doit être resynchronisé).
    void mark_host_dirty() {
        if (state_ == StorageState::SYNCHRONIZED) {
            state_ = StorageState::HOST_DIRTY;
        }
    }

    // Marque le device comme dirty (le host doit être resynchronisé).
    void mark_device_dirty() {
        if (state_ == StorageState::SYNCHRONIZED) {
            state_ = StorageState::DEVICE_DIRTY;
        }
    }

    // Synchronise host → device (M0.12+).
    void sync_host_to_device() {
        // M0 : GPU non disponible. Le state reste HOST_ONLY.
        // M0.12+ : implémentation CUDA.
    }

    // Synchronise device → host (M0.12+).
    void sync_device_to_host() {
        // M0 : GPU non disponible.
    }

    // --- Dimensions ---

    std::size_t size() const noexcept { return host_.size(); }
    bool empty() const noexcept { return host_.empty(); }

    void resize(std::size_t n) {
        host_.resize(n);
        device_.size_bytes = n * sizeof(double);
    }

    void clear() {
        host_.clear();
        device_.size_bytes = 0;
        state_ = StorageState::HOST_ONLY;
    }

private:
    HostStorage host_;
    DeviceStorage device_;
    // WorkingSetStorage réservé pour M0.12+.
    StorageState state_;
};

}  // namespace core
}  // namespace cfdx