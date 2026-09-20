# CFDX Baseline - Repository State (After T00-T10 + GeometryCache + HDF5 Fixes)

## Commit
- Branch: master
- Status: All tasks T00-T10 completed + GeometryCache + HDF5 roundtrip fixes

## Build Command
```bash
# Release build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4

# Sanitizer build (ASan + UBSan)
cmake -S . -B build -DCMAKE_BUILD_TYPE=DebugSanitizers
cmake --build build -j4
```

## Completed Tasks

### Phase 0 — Stabilisation du socle ✅
- **T00 - Baseline**: Documented repository state in `docs/development/BASELINE.md`
- **T01 - CMake Cleanup**: 
  - Reorganized sources by module (core/mesh, core/geometry, core/fvm, core/linalg, core/field, io/hdf5)
  - Removed headers from `add_library()` - only .cpp files listed
  - Added `add_cfdx_test()` macro for automatic test discovery
  - Added C++17, strict warnings, Debug/Release/DebugSanitizers build types
  - Fixed HDF5 detection and linking (added zlib)
- **T02 - Sanitizers**: Added `DebugSanitizers` build type with ASan+UBSan, all 26 tests pass clean
- **T03 - CI**: Created `.github/workflows/ci.yml` with Release + DebugSanitizers matrix

### Phase 1 — Refaire proprement le modèle Mesh ✅
- **T04 - Define Index types** (`src/cfdx/core/mesh/index_types.h`):
  - `Index = uint64_t` (global)
  - `LocalIndex = uint32_t` (local/GPU)
  - Semantic aliases: `PointIndex`, `FaceIndex`, `CellIndex`, `VertexIndex`, `FaceOfCellIndex`
  - `Offset = uint64_t` for CSR offsets
  - Sentinel values: `INVALID_INDEX`, `BOUNDARY_NEIGHBOUR = -1`
- **T05 - Separate MeshTopology from Geometry**:
  - Created `index_types.h` as single source of truth
  - Updated all mesh files: `point.h`, `face.h`, `ownership.h`, `cell.h`, `boundary.h`, `mesh.h`
  - Updated geometry files: `face_geometry.h`, `cell_geometry.h`, `mesh_quality.h`, `mesh_validator.h`
  - Updated FVM numerics: `gradient.h`, `divergence.h`, `flux.h`, `laplacian.h`, `integrate.h`
  - Updated linear algebra: `sparse_matrix.h`, `vector.h`
  - Updated HDF5 I/O: `hdf5_writer.cpp`, `hdf5_reader.cpp` (64-bit support)
  - Updated all test files
- **T06 - Strengthen topology validator** (`mesh.h::topo_validate()`):
  - Validates point/face/cell counts match
  - Validates CSR consistency (strictly increasing offsets)
  - Validates index ranges (no out-of-bounds references)
  - Validates owner/neighbour ↔ cell-face consistency (both directions)
  - Validates boundary patches (no overlap, all boundary faces covered)
  - Validates cell face reference count (2×internal + 1×boundary)
- **T07 - Boundary patches validation** (`boundary.h::is_consistent()`):
  - No duplicate face IDs across patches
  - All face IDs in valid range
  - Union of patches = all boundary faces
  - Intersection of patches = ∅
- **T08 - Face orientation convention** (`face_geometry.h`):
  - Documented convention: `Sf` points from owner toward neighbour/outside
  - Added `ensure_face_orientation()` utility
  - Added `compute_face_geometry_oriented()` that auto-corrects orientation
- **T09 - Fix face geometry (Sf orientation)**:
  - `compute_face_geometry()` returns Sf via fan triangulation (CCW = +Z normal)
  - `ensure_face_orientation()` flips Sf if `Sf · (C_neigh - C_owner) < 0`
  - For boundary: `Sf · (C_face - C_owner) > 0` (points outward)
- **T10 - GeometryCache** (`geometry_cache.h`):
  - Proper `GeometryCache` class with all cached geometry
  - Face centres, Sf, areas, normals, quality metrics
  - Cell centres, volumes, surface closure error
  - Delta coefficients for non-orthogonal correction
  - 9 new tests in `test_geometry_cache.cpp`

### HDF5 Roundtrip Fixes (Bonus)
- Fixed HDF5 writer to properly serialize Field SoA layout (interleave components)
- Fixed HDF5 reader to de-interleave flat array back to SoA layout
- Fixed empty mesh handling (write empty datasets instead of skipping)
- Removed unused 32-bit read functions, using 64-bit throughout

## Test Results (Final)

### Release Build
```
Test Results: 26/26 tests pass
- test_point, test_face, test_ownership, test_cell, test_boundary, test_mesh: PASS
- test_face_geometry, test_cell_geometry, test_mesh_quality, test_geometry_cache: PASS
- test_field, test_field_storage, test_boundary_field: PASS
- test_interpolation, test_gradient, test_divergence: PASS
- test_flux, test_integrate: PASS
- test_laplacian: PASS
- test_sparse_matrix, test_vector, test_linear_system: PASS
- test_cg_solver, test_bicgstab_solver: PASS
- test_hdf5_writer: 4/4 PASS
- test_hdf5_roundtrip: 7/7 PASS
```

### DebugSanitizers Build (ASan + UBSan)
```
Test Results: 26/26 tests pass
No sanitizer errors detected
```

## Known Issues
1. HDF5-DIAG errors appear in test output (HDF5 internal diagnostics for missing objects during sequential test runs) - benign, tests pass
2. Conversion warnings from `uint64_t` → `double` / `uint32_t` in geometry code (narrowing, but safe for typical mesh sizes)
3. Minor unused variable warnings in `geometry_cache.h` (`n_points`, `neigh_centre`)

## CMake Options
- `CFDX_BUILD_TESTS=ON` (default)
- `CFDX_ENABLE_MPI=ON` (default)
- `CFDX_USE_HDF5_SERIAL=ON` (default, uses local third_party/hdf5)
- Build types: `Release`, `Debug`, `DebugSanitizers`

## Next Steps (per roadmap)
- **T11: MPI parallel** (partitioning, halo exchange)
- **T12: GPU backend** (Kokkos/CUDA/HIP)
- Physics modules (Euler, Navier-Stokes, energy equation)
