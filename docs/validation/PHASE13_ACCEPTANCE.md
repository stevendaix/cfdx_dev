# Phase 13 — Numerical-model verification acceptance

## Purpose

Phase 13 is the verification layer between module-level tests (Phases 0–12) and the external VMFL/reference campaign (Phase 14). It must demonstrate that implemented numerical formulations converge to the correct mathematical solution, preserve independent conservation invariants, and remain equivalent across execution backends.

A passing unit test or solver residual is not sufficient evidence for Phase 13 closure.

## Audit baseline — 2026-09-23

Current master already contains useful verification assets: scalar-diffusion MMS and mesh refinement; convection-diffusion MMS/refinement; RK2 ODE temporal refinement; Ghia cavity reference validation; assembled/matrix-free numerical primitives; small-matrix Krylov/preconditioner robustness tests; M1–M4 closure and Level-C coupled tests; MPI Poisson equivalence; CUDA Poisson verification which skips without CUDA hardware; and CPU/emulation OOC regression.

These are valuable foundations, but they do not yet constitute a complete Phase 13 campaign.

## Acceptance matrix

| Item | Current evidence | Gap | Required closure |
|---|---|---|---|
| 13.1 Manufactured solutions | Scalar diffusion + convection-diffusion MMS | No full NS, transient PDE, turbulence, energy/CHT or radiation MMS | Add solver-level MMS with analytic source terms and independent L1/L2/Linf oracles |
| 13.2 Spatial convergence | Scalar diffusion / Poiseuille and convection-diffusion refinement; Ghia grids | Coverage is not systematic across active schemes/models | 3+ similar meshes per formulation, observed order and explicit expected order |
| 13.3 Temporal convergence | RK2 ODE refinement | No PDE transient verification; RK3 only algebraic | Transient diffusion/energy MMS with dt refinement and order gates for every production integrator |
| 13.4 Conservation independent of residuals | Continuity and energy-balance diagnostics exist | No repository-wide independent integral oracle across M1–M4 | Compute integrated boundary flux + volume source + accumulation independently of algebraic residuals |
| 13.5 Assembled vs matrix-free | Matrix-free operator exists and is tested | No complete production FVM operator equivalence matrix | Compare operator action, diagonal where applicable, residual and final solution on identical cases |
| 13.6 Solver/preconditioner robustness | Small analytical matrices and CFD-like matrix tests | Limited conditioning/mesh scaling and production-system evidence | Conditioning/scale campaign, breakdown detection, true residual, finite checks, iteration robustness |
| 13.7 Pressure-correction consistency | Rhie-Chow and pressure/velocity primitives; Couette | No independent manufactured pressure-velocity oracle | Divergence-free MMS with known p/U and pressure-correction invariant checks |
| 13.8 Coupled energy | Single-cell and simple radiation/CHT cases | No systematic manufactured energy/CHT convergence | Thermal MMS, interface flux balance, temporal/spatial refinement |
| 13.9 Turbulence transport | Positivity and closure checks | No manufactured transport verification for k-epsilon/SST/SA | Isolated transport MMS with known source, units and convergence |
| 13.10 Radiation | Closed-form radiation and equilibrium coupling | No systematic transport/angular convergence campaign | DOM/P1/Rosseland verification with independent oracles and angular refinement |
| 13.11 Parallel equivalence | MPI infrastructure and Poisson equivalence | Physics solvers not systematically checked serial vs MPI | Same physical case, same discretization, mapped fields and quantitative difference gates |
| 13.12 GPU equivalence | CUDA Poisson test exists but skips without device | No CUDA-capable CI evidence; no M1–M4 GPU campaign | Hardware evidence; compare CPU/GPU fields and residuals |
| 13.13 OOC equivalence | CPU/emulated OOC tests | No real CUDA beyond-VRAM CFD evidence | Real GPU-OOC case exceeding VRAM, compare against in-core execution |

## Required verification protocol

Every Phase-13 campaign shall record: case identifier; formulation; mesh family and characteristic h; dt when transient; numerical/model controls; independent oracle; L1, volume-weighted L2 and Linf errors where applicable; observed order; integrated conservation defect independently recomputed from final fields; algebraic residuals separately; solver status and true residual; backend/execution mode; and reproducibility metadata where available.

## Evidence rules

- A residual below tolerance never substitutes for an error-to-oracle gate.
- A converged linear solve never substitutes for physical conservation.
- A single mesh never establishes convergence order.
- CPU/emulation never proves real CUDA equivalence.
- A local test never proves MPI partition invariance.
- A skipped hardware test is BLOCKED/PARTIAL, not PASS.
- OpenFOAM/reference comparison is secondary to analytic or manufactured verification when an oracle exists.

## Immediate implementation in this PR

This PR adds the Phase-13 acceptance contract and executable verification coverage for claims already supported by the current implementation. It deliberately does not convert missing campaigns into false passes.

## Pending campaign identifiers

- P13-MMS-NS: divergence-free Navier–Stokes MMS.
- P13-MMS-TD: transient diffusion/energy MMS.
- P13-MMS-TURB: k-epsilon/SST/SA transport MMS.
- P13-MMS-RAD: participating-media radiation verification.
- P13-CONS: independent M1–M4 integral conservation.
- P13-LINALG: production-system solver/preconditioner robustness.
- P13-PV: pressure-velocity manufactured consistency.
- P13-EQUIV-MF: assembled/matrix-free production equivalence.
- P13-MPI-PHYS: serial/MPI physics equivalence.
- P13-GPU: hardware CPU/GPU equivalence.
- P13-OOC: real CUDA beyond-VRAM equivalence.

## Closure gate

Phase 13 can be marked green only when every 13.1–13.13 item has an executable test/case, independent oracle or invariant, quantitative tolerance/order gate, negative/degenerate-input coverage where applicable, CI evidence for software-only requirements, explicit hardware evidence for GPU/OOC requirements, and a documented result artifact.

Until then the phase remains IMPLEMENTED / VALIDATION IN PROGRESS.
