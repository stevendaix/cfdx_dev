# NACA0012: importer repair and verification boundary, 2026-10-25

Baseline: `71858a20fbfc67006824936d0646b3dd1cb79394`. This record distinguishes
verified input/topology repairs from physical flow qualification. No NACA force
or turbulence status is promoted by this change.

## The different NACA cases

| Path | Declared problem | Evidence supplied by that path |
|---|---|---|
| `tests/data/su2/inv_NACA0012_basic.cfg` | Euler, Mach 0.8, incidence 1.25 degrees | SU2 format-reader fixture; no CFDX compressible production solver |
| `tests/data/su2/incomp_NACA0012.cfg` | Incompressible RANS, SA, incidence 10 degrees; source values imply Re approximately 6 million | Imported case metadata; production executable rejects RANS |
| `test_naca0012_qualification` | Incompressible laminar Re=1000, incidence zero, generated O-grid | SIMPLE + upwind solve in the long test; the quick test does not solve the flow |
| N9-S8 | Laminar Re=1000, incidence zero on imported NASA/TMR meshes | Contract and planned six-algorithm/three-mesh campaign; not a turbulent experimental comparison |
| TUR-NACA0012 | Turbulent NACA validation against independent data | Still requires a production RANS path and retained physical evidence |

The NASA turbulent validation problem uses Re=6 million and a modified sharp
trailing edge. Its grid may be reused in a declared laminar study, but the
experimental turbulent drag is not the oracle for a Re=1000 laminar solve.
The generated O-grid uses the finite-thickness original trailing edge; those
geometry definitions must also be recorded before interpreting a drag difference.

Primary references:

- https://tmbwg.github.io/turbmodels/naca0012_grids.html
- https://tmbwg.github.io/turbmodels/naca0012_val.html
- https://tmbwg.github.io/turbmodels/NACA0012_grids/n0012_225-65.nmf

## Reproduced defects in the original NASA/TMR importer

The official `n0012_225-65.p2dfmt` starts with three integers: `1`, `225`, `65`.
The first integer is the number of blocks. The old parser instead interpreted
`1` and `225` as the grid dimensions and rejected the actual file with:

```text
ValueError: expected 452 numeric tokens for 1x225, got 29253
```

The advertised `.gz` files were also read as plain ASCII.

After bypassing only that header defect for diagnosis, the original `convert()`
assigned every j=0 edge to `airfoil`. On the official 225x65 grid this produced
224 wall faces instead of 128. Ninety-six of these faces are the two wake cuts;
their centres extend to x=453.135486 chords. The neutral map explicitly assigns
wall points 49..177 and a reversed one-to-one connection 1..49 to 225..177.
These are **48 pairs of internal wake faces**, not solid walls. The original
converter kept their 49 coincident vertex pairs topologically disconnected.

The i=0 and i=ni-1 sides were named `wake_lower` and `wake_upper`, although the
neutral map identifies them as downstream farfield. Extrusion patches were
named `_front_cells` and `_back_cells`, receiving generic patch type 7 from the
writer instead of EMPTY type 6.

The previous `test_n9_s8_naca0012_contract.py` checked filename dimensions and
words in a chapter. It did not open a PLOT3D file or inspect its converted
connectivity, so none of these defects were exercised.

## Repair

`import_tmr_naca0012.py` now:

1. reads the single-block formatted PLOT3D header, validates token count and
   finite coordinates, and supports gzip;
2. identifies the reversed coincident j=0 wake endpoint interval rather than
   inventing a wall-index ratio from the filename;
3. validates that the detected interval terminates at the unit-chord trailing
   edge, so a partly corrupted wake cannot silently become extra solid wall;
4. joins the paired vertices in both extrusion planes, creating ordinary
   internal owner/neighbour faces along the wake;
5. exports only the airfoil as WALL, combines all outer sides as farfield, and
   exports `front` and `back` as EMPTY;
6. retains the source SHA256 and topology counts in the HDF5 provenance. The
   status remains `imported-provenance-only`.

The implementation is deliberately scoped to the advertised single-block,
unit-chord NASA NACA C-grid family. Other topology is rejected explicitly.
It does not claim to be a generic PLOT3D multiblock connector.

## Verification actually executed

Nine added importer regression cases fail **8/9 on the baseline** and pass
**9/9 after repair**. They cover the real header form, compressed input,
malformed input, nonfinite coordinates, an independent miniature C-grid with
exact face-count oracles, missing wake matching, partial wake corruption and
provenance. They are collected by the existing full `tests/python` CI job.

The official 225x65 grid was then imported and read through the actual native
CFDX HDF5 reader, `validate_mesh` and `build_fv_geometry`, linked against the
baseline production core. The following are measured values, not expected
values copied into a numerical solution:

| Native check | Measured result |
|---|---:|
| Cells | 14,336 |
| Points after wake joining | 29,152 |
| Faces | 57,584 |
| Internal faces | 28,432 |
| Internal wake faces | 48 |
| Airfoil wall faces | 128 |
| Farfield faces | 352 |
| Front / back EMPTY faces | 14,336 each |
| Unassigned physical boundary faces | 0 |
| Minimum cell volume | 3.3070392178011831e-9 |
| Minimum oriented face dot product Sf·d | 1.6535196089005912e-9 |
| Geometry errors | 0 |
| Geometry warnings | 1, large aspect ratio |
| Maximum nonorthogonality | 12.06419507085366 degrees |
| Maximum skewness | 0.15920642368411458 |
| Maximum absolute surface closure | 3.5527136788005009e-15 |

Every cell volume and every face-normal orientation was checked for finiteness
and strict positivity. All boundary faces belonged to exactly one patch.

Official source hashes:

- `n0012_225-65.p2dfmt`: `66817d6bbc9353cb0ca9867f958da1bcce6485f7d68ec7038362a341bdb0fa8c`
- `n0012_225-65.nmf`: `c6c32668f19bbd1a89510dfeb3a43c924137ae90afa7772eeb83875f9d9b34b0`

Reproduction, with assets downloaded deliberately before the command:

```bash
python scripts/import_tmr_naca0012.py \
  --input /path/to/n0012_225-65.p2dfmt \
  --output /path/to/n0012_225-65.h5
PYTHONPATH=src/cfdx/python python -m pytest -q \
  tests/python/test_tmr_naca0012_import.py \
  tests/python/test_su2_naca0012_tutorial.py \
  tests/python/test_su2_naca0012_validation.py
```

The three selected Python files contain **12 passing tests** after the complete
change. A passed mesh import is not a passed flow simulation.

## Honest SU2 tutorial scope

The original Python command stopped at canonical conversion: JST has no faithful
registered CFDX convection counterpart. Its error printer also classified
`unsupported_nonblocking` entries as blocking because it searched for the
substring `block`. Its unused "theory" comparison treated a freestream static
pressure as stagnation pressure and did not provide an aerodynamic oracle.
The shell wrapper used unsupported options `--case`, `--output`, `--max-iter`.
The production executable accepts `--mesh`, `--output-dir`, `--iterations` and
explicitly rejects the compressible energy model.

The Python command now retains distinct reader, boundary/initial-condition,
canonical-conversion and physical-run verdicts in `reader_report.json`; the
shell command delegates to it. Its default exit code remains nonzero (2) when
full aerodynamic validation is unavailable. An explicit `--reader-only` permits
success for the restricted input-reader scope:

```bash
python validation/run_su2_naca0012.py --reader-only --output-dir /path/to/evidence
```

No JST substitution is introduced and no solver result is manufactured.

The current SU2 adapter still has these independently observed metadata defects,
which this importer-focused patch does not repair:

- `MARKER_EULER=(airfoil)` becomes `wall_no_slip` through a name-based fallback;
- `farfield` becomes `outlet`/`fixed` with zero velocity;
- initial scalar speed is 272 m/s, but the vector contains the unit incidence
  direction `[0.9997620270799091, 0.02181488503456112, 0]`.

The updated report exposes all three mismatches. Existing format-reader tests
check patch names and cardinalities, not faithful boundary/initial-condition
semantics. These defects must be fixed before that importer can be described as
physically equivalent to the original SU2 case.

## Remaining numerical limitations relevant to the laminar NACA solve

The existing generated O-grid n64 is topologically valid, but native geometry
inspection gives nonorthogonality 80.550987 degrees, maximum skewness 0.986432
and 82 warnings. Its wall-adjacent centre distance has median 0.0268116 chord,
minimum 0.00511106 and maximum 0.0326396. It is a coarse, highly skewed starting
mesh. The indexing of its circular outer contour also breaks the reflection
pairing of the airfoil surface: the maximum corresponding reflection difference
in the first layer is 0.00319299 chord. These are measured quality limitations;
their contribution to Cd/Cl must be quantified by controlled flow comparisons.

On both actual meshes, the exact manufactured affine field phi=x+2y was supplied
to the production `gauss_gradient_with_boundary`. Only cells with no physical
boundary face were included, so the result is independent of boundary-value
approximations. WLS received explicit exact Dirichlet values at physical faces
and the correct zero spanwise derivative.

| Mesh | Interior cells | Production Gauss gradient max absolute error | WLS max absolute error |
|---|---:|---:|---:|
| Generated O-grid n64 | 2,838 | 0.10987995423605813 | 1.1082784313214814e-14 |
| Repaired TMR 225x65 | 13,858 | 0.17283073646457187 | 5.4393680889484996e-9 |
| Repaired TMR 449x129 | 56,386 | 0.10357401848014633 | 1.3019123605904252e-8 |
| Repaired TMR 897x257 | 227,458 | 0.07351694706579924 | 2.9871163766996353e-8 |

This demonstrates a reconstruction error before any nonlinear flow solve. The
production helper still uses Euclidean interpolation distances, separately from
`core/numerics/gradient.h`, and has no face-skew correction. Merely correcting a
standalone gradient kernel cannot certify this production path. These numbers
are **not** drag errors and do not by themselves establish the whole cause of
an aerodynamic discrepancy.

Further requirements remain unchanged: fix and independently verify the
pressure/velocity face-flux contract; validate the actual production gradient
and nonorthogonal diffusion; freeze consistent farfield conditions; compare
pressure and viscous force parts and surface Cp/Cf; execute the declared three
meshes and six algorithms; retain true residuals, independent balances and an
appropriate physical reference. A unit test, quick contract, converged linear
solve or this mesh repair cannot replace those results.

## Refined official grids and boundary-membership cost

The same native reader/geometry checks were also completed on the two finer
official files. Their gzip payloads were retrieved from the NASA archive,
`https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip`, and checked
against its ZIP member CRC and length before decompression.

| Check | 449x129 | 897x257 |
|---|---:|---:|
| Cells | 57,344 | 229,376 |
| Points after joining | 115,648 | 460,672 |
| Faces | 229,856 | 918,464 |
| Internal faces | 114,208 | 457,792 |
| Internal wake faces | 96 | 192 |
| Airfoil wall faces | 256 | 512 |
| Farfield faces | 704 | 1,408 |
| Front/back EMPTY faces each | 57,344 | 229,376 |
| Unassigned boundary faces | 0 | 0 |
| Minimum cell volume | 5.1820719983179846e-10 | 9.579601570460766e-11 |
| Minimum oriented face dot product | 2.5910359991589918e-10 | 4.7897995986945982e-11 |
| Geometry errors | 0 | 0 |
| Geometry warnings | 1, aspect ratio | 1, aspect ratio |
| Maximum nonorthogonality, degrees | 9.2209192083970688 | 6.1731507752808854 |
| Maximum skewness | 0.10395424486342719 | 0.055277414197076311 |
| Maximum absolute surface closure | 1.7763568394002505e-15 | 4.4408920985006262e-16 |

Source SHA256:

- `n0012_449-129.p2dfmt.gz`: `6d60dfa6ef9c7544ee156276fcc5487d7eacddadbdf6e81cc6f10e78fb03eff8`
- `n0012_897-257.p2dfmt.gz`: `c603ea77017ae8b97af581237b81a8251fd2d3af1c49c8c6f2cd2687aafcb28e`

During these checks, `BoundaryPatches::is_consistent_with_mesh` was found to
rescan every patch face for each physical boundary face. With the two large
EMPTY patches, this is quadratic in the number of cells. The native HDF5
reader invokes this validation, and subsequent `validate_mesh` invokes it again.
The repair retains the duplicate/range and internal-face checks, records listed
faces in a bitmap, and then scans the boundary once. It is linear in the number
of faces plus patch entries.

A benchmark on the three actual imported grids called the exact original
function body and its replacement on the same loaded mesh. Both returned true.
These are single-process wall-clock observations for that function alone,
not an end-to-end CFD speedup guarantee:

| Grid | Boundary faces | Original, seconds | Replacement, seconds |
|---|---:|---:|---:|
| 225x65 | 29,152 | 0.121096512 | 0.000140853 |
| 449x129 | 115,648 | 2.34260247 | 0.000537623 |
| 897x257 | 460,672 | 31.2930124 | 0.002588804 |

Five membership regression cases explicitly check a complete patch cover,
an absent face, a repeated face, an internal face in a patch, and an out-of-range
face. The expanded `test_mesh` passes 16/16 both with the original and repaired
function, demonstrating preserved acceptance/rejection decisions. The existing
`test_boundary` passes 11/11 with the repair. The complete native geometry and
manufactured-gradient diagnostics on all three actual grids are byte-for-byte
identical before and after this membership-only change. No timing threshold is
used as a unit-test oracle.
