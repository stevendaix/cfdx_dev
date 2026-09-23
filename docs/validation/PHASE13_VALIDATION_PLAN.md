# Phase 13 — Complete numerical-model validation

**Status:** validation framework implemented; numerical closure remains evidence-driven.

## Scope
Phase 13 is the numerical verification layer for CFDX. It establishes that discrete equations converge to the mathematical model, conservation is independently satisfied, and alternative execution paths produce equivalent physical results.

## Evidence hierarchy
1. Analytical solution — preferred.
2. Manufactured solution (MMS) — preferred when no closed form exists.
3. Conservation/invariant oracle — mandatory where an exact field oracle is unavailable.
4. Reference benchmark — published/reference data.
5. Cross-backend equivalence — assembled/matrix-free, serial/MPI, CPU/GPU.
6. Residuals — solver health only; never a substitute for the above.

## 13.1 Manufactured solutions
Existing evidence includes scalar diffusion and convection-diffusion MMS foundations.
Required production closure: divergence-free Navier–Stokes MMS with analytic forcing; transient PDE MMS; energy/CHT MMS; isolated turbulence-transport MMS for k-epsilon, SST and SA; radiation MMS where applicable; L1, volume-weighted L2 and Linf errors; observed order from at least three refined meshes.

## 13.2 Spatial convergence
Each production discretisation requires a controlled mesh family. Record h, cell count, scheme/order, error norms, observed order, geometry quality, solver tolerance and iterations. Non-orthogonal/skewed cases must be separate from orthogonal verification.

## 13.3 Temporal convergence
Existing RK checks are necessary but insufficient. Add transient diffusion MMS, transient energy MMS, every production time integrator, at least three dt levels, spatial error below temporal error, and observed temporal order.

## 13.4 Independent conservation
Recompute conservation from final fields without reusing the algebraic residual. For each equation: boundary_flux + volume_source - accumulation = conservation_defect. Report absolute and normalised defects for mass, momentum and energy where implemented.

## 13.5 Assembled / matrix-free equivalence
For identical mesh, coefficients and fields compare operator action, diagonal where defined, residual, final solution and independently recomputed true residual. Equivalence tolerance must be tighter than the physical discretisation error for deterministic tests.

## 13.6 Solver and preconditioner robustness
Cover SPD/nonsymmetric systems, increasing mesh size, varying condition number, zero/near-zero diagonal detection, breakdown detection, finite-value checks, true residual recomputation, iteration tracking, and all available Krylov/preconditioners. CONVERGED alone is insufficient.

## 13.7 Pressure-velocity consistency
Add divergence-free pressure/velocity MMS, independent continuity defect, pressure-correction consistency, Rhie-Chow collocated-grid behaviour, and Couette/Poiseuille regression. Pressure null-space treatment must be explicit.

## 13.8 Thermal / CHT verification
Add conduction MMS, transient thermal MMS, conjugate interface flux balance, temperature/heat-flux norms, spatial/temporal refinement, and coupling convergence independent of the energy linear residual. Radiation/thermal coupling must conserve energy.

## 13.9 Turbulence transport verification
For each enabled model, verify transport with analytic/MMS source, production/destruction/source balance, positivity, dimensional consistency, mesh convergence, and no artificial PASS from clipping. SA, SST and k-epsilon transport must be distinguished from algebraic closure tests.

## 13.10 Radiation verification
Cover closed-form/equilibrium regression, participating-medium source balance where applicable, DOM angular refinement, radiative heat flux/source convergence, energy coupling and limiting cases. CPU emulation is not GPU evidence.

## 13.11 Parallel equivalence
Run identical cases in serial and 2+ MPI partitions; compare mapped fields, global conservation and reproducibility. MPI Poisson is foundation evidence, not full M1–M4 closure.

## 13.12 GPU equivalence
On real CUDA hardware compare CPU/GPU fields, conservation, true residual and metadata. A skipped CUDA test is BLOCKED, never PASS.

## 13.13 Out-of-core equivalence
On real CUDA hardware run an oversized working set and compare against in-core execution where feasible. Verify field/conservation equivalence, checkpoint/restart, explicit OOC mode, and absence of silent GPU-to-CPU fallback.

## Automated campaign
`scripts/run_validation.py` is the software validation driver. It runs the CTest `phase13` label, emits JUnit and JSON evidence, records hardware-only campaigns as BLOCKED_WITHOUT_HARDWARE, and supports `--strict` to fail until all closure items are evidenced.

Relevant existing validation tests are labelled `phase13;validation` in CMake.

## Evidence ledger
| Area | Current evidence | Closure |
|---|---|---|
| MMS | scalar diffusion / convection-diffusion | pending full solver MMS |
| Spatial order | refinement foundations | pending complete production matrix |
| Temporal order | RK/temporal foundations | pending transient PDE |
| Conservation | M1–M4 / CHT checks | pending independent integral campaign |
| Matrix-free | operator tests | pending production equivalence |
| Linear solvers | analytical + CFD-like systems | pending scaling/conditioning campaign |
| Pressure/velocity | solver/Couette foundations | pending pressure MMS |
| Thermal/CHT | CHT validation | pending thermal MMS/convergence |
| Turbulence | hardening + thermo validation | pending transport MMS |
| Radiation | analytical/coupling foundations | pending full transport/angular campaign |
| MPI | Poisson + MPI foundations | pending production physics equivalence |
| GPU | CUDA test | blocked without real CUDA evidence |
| OOC | CPU/emulated path | blocked without real CUDA OOC evidence |

Phase 13 closes only when every item has executable quantitative evidence. No result is accepted solely from console text or a solver residual.