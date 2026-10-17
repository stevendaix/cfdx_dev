# NACA0012 wall-distance comparison

## Purpose

Compare the nine CFDX wall-distance methods against an independent geometric
oracle built from the official SU2 NACA0012 mesh already stored in
tests/data/su2/mesh_NACA0012_inv.su2. This is a wall-distance benchmark, not
a CFD solution validation and not a claim that CFDX reproduces SU2 lift, drag,
or pressure coefficients.

Upstream references:
- Mesh: https://github.com/su2code/TestCases/blob/master/euler/naca0012/mesh_NACA0012_inv.su2
- Configuration: https://github.com/su2code/SU2/blob/master/TestCases/euler/naca0012/inv_NACA0012.cfg
- The upstream configuration describes Mach 0.8 and angle of attack 1.25 degrees. Those flow settings do not enter the distance calculation.

## SU2 wall-distance availability audit

The checked-in `mesh_NACA0012_inv.su2` is a mesh file: it contains coordinates,
connectivity, and boundary markers, not a precomputed wall-distance field. The
repository's `inv_NACA0012_basic.cfg` is an import/validation configuration;
there is no compatible SU2 solution or volume-output file beside it. The
existing Python import test explicitly treats the synthetic
`mismatched_solution.csv` as non-validation evidence.

SU2 does calculate and store wall distance internally at geometry nodes
(`GetWall_Distance`). However, the upstream `CFlowOutput.cpp` currently
registers the `WALL_DISTANCE` volume-output field inside the hybrid RANS-LES
(`DDES`) output group, so it is not a generally available field in the
Euler NACA0012 configuration. A compatible SU2 run could be configured and
checked to export it, but that would be a solver-generated nodal field—not an
exact geometric oracle—and would require matching node coordinates, output
format, SU2 revision, and wall-marker semantics before comparison.

For this benchmark, retain the independent analytic point-to-segment oracle
as the accuracy reference. It avoids dependence on SU2's wall-distance
algorithm and avoids comparing different sample locations. SU2's exported
field can be added later as a separate cross-code comparison when an
appropriate viscous/hybrid case and verified output are available; it must not
replace the geometric oracle.

## Method

The executable reads the MARKER_TAG= airfoil line elements directly from the
SU2 mesh. Each 2-D segment is extruded in the out-of-plane direction to form
a rectangular wall panel. At sample points in the z=0 plane, exact
point-to-panel distance equals exact point-to-segment distance in the original
airfoil plane. This exercises the existing CFDX triangle-based distance
implementation without introducing a second approximation to the wall.

The sample grid is Cartesian and surrounds the airfoil marker with a padding
of 0.2 chord on each side. The solid mask is computed by a point-in-polygon
test against the original boundary segments. The grid is extruded through
three z planes; the wall panels span that depth, so the target distance field
is invariant in z. Grid resolution is adjustable to support a quick CI smoke
test and a more detailed local comparison.

The reference is computed independently from the original 2-D marker edges
using the analytic point-to-segment projection formula. It does not call
CFDX's point-to-triangle distance function or its BVH. Because the extruded
wall panels extend beyond all three sample planes, this 2-D oracle represents
the exact distance to the discrete airfoil boundary at every sample. The
exact_geometric and search_based methods must reproduce it to relative L2 and
L-infinity error 1e-12 or better. Other methods are measured, not pre-declared
PASS: their error, convergence status, residual, iteration count, and elapsed
time are recorded independently.

## Build and run

From the repository root after configuring and building:

~~~bash
cmake -S . -B build -DCFDX_BUILD_TESTS=ON
cmake --build build --target cfdx_naca0012_wall_distance_comparison -j
./build/cfdx_naca0012_wall_distance_comparison \
  tests/data/su2/mesh_NACA0012_inv.su2 \
  build/naca0012_wall_distance_comparison.csv 61
~~~

Arguments are: input SU2 mesh, output CSV, and number of x-grid points
(minimum 9; default 61). The executable prints geometry metadata and a compact
table; the CSV contains one row per method.

The CTest registration runs a smaller resolution as an integration smoke test.
It checks that the official marker is parsed and that both geometric methods
agree with the independent segment-distance oracle. It does not impose an accuracy threshold on the
PDE methods.

## Interpretation and limitations

- l2_relative and linf_relative compare computed distance to the geometric
  oracle over fluid samples; near_wall_l2_relative restricts to two minimum
  grid spacings from the wall.
- residual_inf, iterations, converged, and stopping_reason describe the
  method's own iterative solve. A small algebraic residual is not proof of a
  correct reconstructed distance.
- elapsed_ms is only comparable across methods in the same run and machine.
- The benchmark uses a Cartesian sampling grid, not the original unstructured
  SU2 cell centres. It qualifies wall-distance algorithms on the NACA0012
  geometry, not unstructured-grid interpolation or the SU2 flow solution.
- The finite Cartesian box and its outer-boundary treatment can influence PDE
  methods far from the wall. Near-wall error is reported separately.
- No SU2 solution file or aerodynamic coefficient is used as evidence here.
  The existing reader test validates mesh/config import only.
