# Gap Analysis Report

## Summary

| Severity | Count |
|----------|-------|
| Supported / Mapped | 23 |
| Approximated | 0 |
| Unsupported (non-blocking) | 6 |
| Unsupported (blocking) | 0 |
| Unavailable (source info missing) | 0 |
| **Total** | 29 |

**Blocking incompatibilities:** NONE

## 1. Supported / Mapped Features

| Category | Feature | Detail | Suggestion |
|----------|---------|--------|------------|
| mesh | nodes | Parsed 30 vertices from .cas | — |
| mesh | face_nodes | Parsed 35 face-node connectivities from .cas | — |
| mesh | cells | Parsed 1 cells from .cas | — |
| boundary_zone | wall | Detected zone 'wall' type=wall | — |
| boundary_zone | inlet | Detected zone 'inlet' type=velocity-inlet | — |
| boundary_zone | outlet | Detected zone 'outlet' type=pressure-outlet | — |
| boundary_zone | symmetry | Detected zone 'symmetry' type=symmetry | — |
| boundary_zone | periodic | Detected zone 'periodic' type=periodic | — |
| boundary_zone | interior | Detected zone 'interior' type=interior | — |
| boundary_zone | far-field | Detected zone 'far-field' type=pressure-far-field | — |
| boundary_zone | overset | Detected zone 'overset' type=overset | — |
| boundary_zone | fluid | Detected zone 'fluid' type=cell | — |
| material | air | Parsed material 'air' (ρ=1.225, μ=1.789e-05) | — |
| material | water | Parsed material 'water' (ρ=1000.0, μ=0.001002) | — |
| material | steel | Parsed material 'steel' (ρ=7830.0, μ=6e-05) | — |
| boundary_condition | wall | Mapped Fluent zone 'wall' to wall | — |
| boundary_condition | inlet | Mapped Fluent zone 'inlet' to inlet | — |
| boundary_condition | outlet | Mapped Fluent zone 'outlet' to pressure_outlet | — |
| boundary_condition | symmetry | Mapped Fluent zone 'symmetry' to symmetry | — |
| boundary_condition | periodic | Mapped Fluent zone 'periodic' to periodic | — |
| boundary_condition | interior | Mapped Fluent zone 'interior' to unknown | — |
| boundary_condition | far-field | Mapped Fluent zone 'far-field' to unknown | — |
| boundary_condition | overset | Mapped Fluent zone 'overset' to interface | — |

## 2. Approximated / Fallback Mappings

_None._

## 3. Unavailable Source Information

_None._

## 4. Unsupported Features (non-blocking)

| Category | Feature | Detail | Suggestion |
|----------|---------|--------|------------|
| features | binary_cas | Binary .cas.h5 format not supported (only legacy ASCII .cas) | Export as legacy ASCII: File → Export → Case... |
| features | dpm | Discrete Phase Model (DPM) particles not supported | Use single-phase flow for CFDX conversion |
| features | udfs | UDFs and custom field functions not supported | Pre-compute and export field values via Fluent surface or volume monitor |
| features | multiphase | Multiphase models (VOF/Mixture/Eulerian) not supported | Use single-phase configuration for CFDX conversion |
| features | sliding_mesh | Sliding mesh and dynamic mesh not supported | Use steady-state or static mesh for CFDX conversion |
| features | results | Fluent .dat results import requires meshio or neutral export; not parsed natively | Export results as .cgns, .vtk, or .plt via Fluent File → Export |

## 5. Unsupported Features (blocking)

_None._
