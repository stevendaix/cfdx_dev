# NACA0012 XY planar and SST case separation — patch against CFDX V5

This patch is scoped to the CFDX V5 archive (not a verified GitHub PR head).
It neither changes pressure-velocity solver convergence tolerances nor claims
that NACA N256 or SST aerodynamics have passed.

## Existing campaigns: **do not merge them**

- `test_naca0012_qualification`: generated O-grid, **laminar** Re=1000,
  alpha=0, zero-thickness two-dimensional physics using a single z extrusion.
  This is the N64/N128/N256 numerical-convergence diagnostic.
- `TUR-NACA0012`: a **separate** RANS validation case in
  `CFDX_QUALIFICATION_REGISTRY.json`. Re=6e6, Mach 0.15 is a reference
  condition for a compressible code; for a strictly incompressible CFDX
  solver this is metadata, NOT a simulated compressible Mach number. The SST
  variant and turbulence boundary values need freezing before an actual run.
- `scripts/import_tmr_naca0012.py`: existing NASA/TMR C-grid reader, with
  stitched wake and EMPTY caps. The importer explicitly marks its output
  `validation_status=imported-provenance-only` and its historic
  `validation_physics=N9-S8 laminar incompressible campaign`. Neither label
  indicates an SST solution. The imported grid may be used for a **new** SST
  campaign when the model, geometry, inflow state and physics are specified.
- `tests/data/su2/incomp_NACA0012.cfg`: a **Spalart-Allmaras** metadata fixture,
  not the CFDX SST NACA benchmark.
- `tests/python/test_sst_persisted_cli_restart_e2e.py`: SST restart test using
  a tiny synthetic mesh instead of NACA; not an aerodynamic qualification.

NASA/TMR reference: https://tmbwg.github.io/turbmodels/naca0012_val.html
NASA/TMR C-grid: https://tmbwg.github.io/turbmodels/naca0012_grids.html
The NASA geometry has a modified closed trailing edge and a farfield near 500
chords, unlike the generated, blunter trailing edge and radius 20c in the
laminar CFDX O-grid. NASA recommends Re=6 million, essentially incompressible
flow and Mach=0.15 for compressible models, and expects fully turbulent flow.
The SST / SSTm variants and boundary conditions require an explicit decision:
https://tmbwg.github.io/turbmodels/sst.html

## Applied fixes

1. `require_naca_planar_geometry` refuses non-EMPTY caps, intermediate z
   planes, missing cap faces and stacked 3-D verification meshes. This is a
   **necessary** prerequisite, not proof that every discrete operator is 2-D.
2. `require_naca_planar_state` checks final `Uz` and authoritative EMPTY mass
   flux from the solver. It must not infer the mass flux from interpolated U.
   No repair/clipping of a failing state is performed.
3. The aerodynamic force coefficients in the laminar qualification are now
   nondimensionalized with the actual reference area A=c*span, moments with
   A*c. `Cp` and local `Cf` still use dynamic pressure q. The independent
   `forces::integrate` implementation checks the historical face-force sum.
4. `validate_naca_case_separation.py` checks the existing registry and,
   optionally, an imported TMR mesh. It **never** promotes an import to a
   turbulent physical PASS. Its validation is metadata only.

## Remaining blocking work

- N64/N128/N256: real physical convergence, independent residuals and wall
  forces, proper matched mesh family, pressure/shear uncertainty and GCI.
- Final-state Uz/EMPTY is only one gate: add independent face-MMS/Jacobian
  testing of any out-of-plane pressure/velocity terms and time checkpoints.
- SST: choose **SSTm** or another named variant with exact constants;
  use consistent wall k/omega values and `sst_gradient_dimension=2` only
  after verifying the full reduced-dimensional gradient semantics. Build and
  solve a genuine TMR three-grid campaign with measured y+, transition
  assumptions, and reference Cp/Cd/Cl/Cm agreement. Do not reuse the laminar
  Re=1000 drag oracle as a turbulent Re=6e6 acceptance gate.
- Thermal 2-D EMPTY support remains out of scope; do not silently relabel
  EMPTY caps as WALL when testing energy conservation.

## Commands (from the repository root)

```bash
cmake -S . -B build -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF
cmake --build build --target test_naca_planar_contract test_naca0012_qualification -j2
./build/test_naca_planar_contract
python3 scripts/generate_2d_validation_meshes.py naca0012 --quick --output-dir build/naca_validation_meshes
./build/test_naca0012_qualification --quick build/naca_validation_meshes/naca0012_laminar_n64.h5
python3 scripts/validate_naca_case_separation.py
python3 -m pytest -q tests/python/test_naca_case_separation.py tests/python/test_tmr_naca0012_import.py
```

The new quick gates must **not** be counted as N256 or SST physical PASS.
