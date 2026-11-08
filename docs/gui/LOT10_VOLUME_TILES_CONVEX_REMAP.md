# CFDX Workbench Lot 10 — Polyhedral volume tiles and conservative convex remap

## Scope and non-claims

This lot is a **bounded building block** towards an industrial CPU/GPU visualization pipeline. It does not certify native GPU computations, multi-GPU distribution, 50-million-cell rendering, or conservative remapping of arbitrary concave CFDX cells. No modification to CFDX's physical solver or convergence gates is made.

The design follows VTK's piece-aware pipeline principles, but the tile reader implemented here is a separate, explicit HDF5 tile API. It does **not** implement `vtkStreamingDemandDrivenPipeline` itself. VTK/PyVista can render the returned tile through the active CPU/OpenGL backend; device, driver and frame-rate verification are distinct acceptance gates.

## Volumetric partitions — source-truth cell indices

`cfdx.volume_partitions.read_volume_partition(path,start_cell,count=128,field=None)` reads **actual** CFDX cell topology from the HDF5 datasets:

- `points`, `face_vertices`, `face_offsets`, `cell_faces`, `cell_offsets`;
- optionally, a scalar under `/fields/...` whose size is exactly `n_cells`.

Only the requested contiguous range of cells, its referenced faces and vertices, and the selected field values are materialized. The result keeps `original_cell_ids` and `original_point_ids`, plus the full polyhedral VTK cell stream. This is **not** an interpolation or a simulated mesh.

Current hard limits: at most **512 cells per tile**, at most **16 MiB of measured tile arrays**, at most **64 faces per cell**, **32 vertices per face**, and bounded HDF5 chunk decompression. Unsafe soft/external/virtual datasets, invalid offsets, nonfinite coordinates, out-of-range indices, duplicate faces, and open/nonmanifold cell surfaces are rejected.

### CLI benchmark and reproducibility

```bash
PYTHONPATH=src/cfdx/python python -m cfdx.volume_partitions mesh.cfdx.h5 \
  --start 512 --count 128 --field /fields/scalar/T
PYTHONPATH=src/cfdx/python python -m scripts.benchmark_volume_tiling \
  --sizes 1000,10000,100000 --tile 128 \
  --output artifacts/volume_tile_io.json
```

The benchmark uses synthetic, independently formed hexahedra, times the **CPU/HDF5 tile read only**, measures process RSS, and never invents GPU throughput, render FPS or CFD convergence. It is unsuitable for declaring industrial GPU performance.

### GUI

Open **Analyze → Industrial data / Compare → Volume tiles (CPU/VTK)**. Choose a native CFDX HDF5 mesh, the original first cell ID, a tile count, and an optional scalar. The HDF5 work runs in a Qt worker. The 3D scene renders that tile without claiming that missing cells have been loaded. Existing clip/slice and visualization filters apply to the displayed tile only.

A `volume_partition_ready` event reports the first/last original cell IDs and the fraction of the source represented. Switching case/mesh invalidates older in-flight requests. The default dark Workbench theme is unchanged.

## First-order convex-polyhedron conservative remapping

`cfdx.convex_polyhedral_remap.conservative_convex_remap(source,target,field)` calculates geometric intersection volumes of **closed convex polyhedral cells**. It validates face closure and orientation, convexity, nonzero volumes, finite values, common unit/time, two-way domain coverage and conservation of the source field's volume integral.

The reference algorithm is inspired by Grandy's cell-volume overlap work and the supermesh approach, but implemented here as a **limited first-order convex prototype**. SciPy/Qhull is used for convex hull geometry; candidate intersection vertices are computed using mixed triples of supporting planes. This numeric implementation is **not** an arbitrary-precision geometric predicate system and cannot be promoted to general concave/warped polyhedra.

Current bounds: at most **256 cells per mesh**, **20,000 donor-target pairs**, **24 faces per cell**, **16 vertices per face**, and **48 MiB combined topology arrays**. On domain mismatch, degenerate or warped faces, nonmanifold surfaces, invalid cell volumes, or nonconservative integrals it **fails closed** without returning a fabricated remap.

```bash
PYTHONPATH=src/cfdx/python python -m cfdx.convex_polyhedral_remap \
  source.cfdx.h5 target.cfdx.h5 /fields/scalar/T
```

CLI/GUI reports are **diagnostics only** (`certificate: false`). The existing fast Cartesian-box mapping remains available and separate.

## Bibliography

- Jeffrey Grandy, *Conservative Remapping and Region Overlays by Intersecting Arbitrary Polyhedra*, Journal of Computational Physics 148 (1999), 433–466, DOI 10.1006/jcph.1998.6125.
- Sandeep Menon and David P. Schmidt, *Conservative interpolation on unstructured polyhedral meshes*, CMAME 200 (2011), 2797–2804, DOI 10.1016/j.cma.2011.04.025.
- VTK reference, `vtkStreamingDemandDrivenPipeline`, requesting pipeline pieces/sub-extents: https://vtk.org/doc/nightly/html/classvtkStreamingDemandDrivenPipeline.html.

## Tests, gates, and roadmap

`tests/python/test_workbench_lot10_polyhedral.py` covers skew affine hexahedra, tetrahedra remapped into four sub-tetrahedra, corrupt and incomplete domains, symlink rejection, exact original cell IDs and values, invalid topology, the construction of a real `vtkPolyhedron` with native VTK when available, and optional PyVista/Qt checks. The CI GUI workflow registers the test and a separate synthetic I/O benchmark artifact. A skipped GUI/VTK test is **not** a pass.

Remaining qualification gates before production claims:

1. Reliable direct rendering and picking of real 1M/10M/50M CFDX volumes with measured time-to-first-frame, navigation FPS, CPU peak RSS, GPU memory peak, tile cache/cancellation behavior.
2. Multi-core partition scheduling and volume LOD; hardware-backed VTK-m/GPU, then multi-GPU/distributed composition verified on a supported device.
3. Robust general-polyhedron intersections including concave/warped cells, near-tangencies, nonconforming multi-region interfaces, conservative transfer of all physical fields and face flux reprojection.
4. Physical validation across NACA/BFS/Ghia/SST/CHT/radiation; no relaxation of numerical convergence gates.
