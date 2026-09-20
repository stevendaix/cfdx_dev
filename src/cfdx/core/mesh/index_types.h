#pragma once

// T04 — Centralized index type definitions
// ===========================================================================
// CFDX uses a unified index type system to ensure consistency across the
// mesh topology, geometry, fields, and linear algebra modules.
// ===========================================================================

#include <cstddef>
#include <cstdint>
#include <limits>

namespace cfdx {
namespace core {

// ---------------------------------------------------------------------------
// Core index type configuration
// ---------------------------------------------------------------------------
// Global choice: 64-bit indices for future-proofing large-scale problems.
// For GPU compatibility, 32-bit is also available as LocalIndex.
// ---------------------------------------------------------------------------

using Index = std::uint64_t;        // Global index (mesh topology, global field access)
using LocalIndex = std::uint32_t;   // Local index (per-rank, GPU kernels, intra-block)
using GlobalIndex = std::uint64_t;  // Global index (MPI global numbering)

// Type aliases for specific mesh entities (semantic clarity)
using PointIndex = Index;
using FaceIndex = Index;
using CellIndex = Index;
using VertexIndex = Index;  // Vertex of a face (local to face)
using FaceOfCellIndex = Index;  // Face within a cell (local to cell)

// Offsets in CSR structures
using Offset = std::uint64_t;

// ---------------------------------------------------------------------------
// Sentinel values (invalid indices)
// ---------------------------------------------------------------------------
static constexpr Index INVALID_INDEX = std::numeric_limits<Index>::max();
static constexpr LocalIndex INVALID_LOCAL_INDEX = std::numeric_limits<LocalIndex>::max();
static constexpr FaceIndex INVALID_FACE_INDEX = INVALID_INDEX;
static constexpr CellIndex INVALID_CELL_INDEX = INVALID_INDEX;
static constexpr PointIndex INVALID_POINT_INDEX = INVALID_INDEX;

// ---------------------------------------------------------------------------
// Boundary sentinel (matches FaceOwnership::BOUNDARY)
// ---------------------------------------------------------------------------
static constexpr std::int64_t BOUNDARY_NEIGHBOUR = -1;

// ---------------------------------------------------------------------------
// Configuration flag: enable 64-bit index checks
// ---------------------------------------------------------------------------
#if defined(CFDX_INDEX_64BIT)
    // 64-bit is default
#else
    // Could add compile-time 32-bit only mode if needed
#endif

// ---------------------------------------------------------------------------
// Compile-time checks
// ---------------------------------------------------------------------------
static_assert(sizeof(Index) == 8, "Index must be 64-bit");
static_assert(sizeof(LocalIndex) == 4, "LocalIndex must be 32-bit");
static_assert(std::numeric_limits<Index>::max() > std::numeric_limits<LocalIndex>::max(),
              "Index must be larger than LocalIndex");

}  // namespace core
}  // namespace cfdx
