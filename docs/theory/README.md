# CFDX Theory — Fluid Mechanics and Numerical CFD Course

The Theory domain is the scientific course and mathematical reference for CFDX. Quantitative scientific content is authored as executable Jupytext percent-format Python (chapter.py); README files provide navigation and scope.

## Master learning path

\[
\boxed{
\text{fluid mechanics}
\rightarrow
\text{conservation laws}
\rightarrow
\text{finite volume}
\rightarrow
\text{reconstruction/fluxes}
\rightarrow
\text{time integration}
\rightarrow
\text{coupled solution}
\rightarrow
\text{linear algebra}
\rightarrow
\text{physical closures}
\rightarrow
\text{V\&V}
}
\]

### A — Fluid mechanics

00 Foundations → 01 Conservation laws

### B — Spatial discretisation

02 Finite-volume method → 03 Meshes → 04 Gradients & reconstruction → 05 Fluxes

\[
\text{PDE}
\rightarrow
\int_V\text{PDE}\,dV
\rightarrow
\text{face/volume contributions}
\rightarrow
\text{algebraic coefficients}.
\]

### C — Time and flow solution

06 Time integration → 07 Pressure–velocity coupling → 08 Linear algebra → 09 AMG/MGR/Schur

### D — Physical closures

10 Turbulence → 11 Heat transfer → 12 Radiation → 13 Multiphysics

### E — Numerical evidence

14 Numerical analysis → 15 Verification & validation

### F — CFDX execution

16 Data model & file formats → 17 Computational chain → 18 Source-tree physics audit

## Standard scientific subsection contract

Every important equation or numerical method should contain, where applicable:

1. physical motivation;
2. variable definitions, dimensions and SI units;
3. assumptions and sign conventions;
4. continuous equation;
5. derivation;
6. integral/control-volume form;
7. discrete formulation;
8. algebraic contribution;
9. numerical properties;
10. failure modes and sensitivity;
11. exact CFDX source mapping;
12. executable verification;
13. analytical/benchmark reference;
14. limitations;
15. improvement paths;
16. bibliography.

A summary table may follow the derivation; it must never replace it.

## Scientific integrity

- No invented CFDX results.
- No reference value presented as a CFDX result.
- No residual-only physical validation.
- No implementation-to-qualification shortcut.
- No tolerance relaxation to manufacture a green result.
- No failing benchmark silently converted into a pass.
- Quantitative claims record mesh, physics, numerical settings and software revision where relevant.

## Status vocabulary

- **Implemented:** code path exists.
- **Verified:** a defined mathematical/software property has executable evidence.
- **Validated:** independent physical/reference evidence supports the declared case.
- **Qualified:** the declared population and acceptance gates are complete.

Documentation completeness never changes numerical maturity.

## Chapter navigation

| # | Chapter | Role |
|---|---|---|
| 00 | Foundations | Fluid-mechanics course |
| 01 | Conservation laws | Conservative governing equations |
| 02 | Finite-volume method | Cell-integrated discretisation |
| 03 | Meshes | Topology and geometry |
| 04 | Gradients & reconstruction | Derivatives and face reconstruction |
| 05 | Fluxes | Conservative interface transport |
| 06 | Time integration | Physical and pseudo-time schemes |
| 07 | Pressure–velocity coupling | Incompressible flow solution |
| 08 | Linear algebra | Sparse systems and Krylov solvers |
| 09 | AMG/MGR/Schur | Scalable coupled linear solution |
| 10 | Turbulence | RANS/LES/hybrid closures |
| 11 | Heat transfer | Energy, conduction and CHT |
| 12 | Radiation | S2S/P1/DOM and radiative coupling |
| 13 | Multiphysics | Coupled residual/Jacobian systems |
| 14 | Numerical analysis | Accuracy, stability and error |
| 15 | V&V | Verification, validation and qualification |
| 16 | Data model | Case/checkpoint/output semantics |
| 17 | Computational chain | PDE-to-solver workflow |
| 18 | Source-tree audit | Equation-to-code-to-test traceability |

See docs/references/bibliography.bib for the maintained bibliography.


## Second-pass closure audit — preparation for the next numerical-method PRs

The chapters are now a coherent mathematical baseline, but the documentation must stay ahead of implementation without claiming qualification. The next pass should close the following explicit gaps, using the open engineering issues as the planning authority.

### 1. Mathematical depth to add before claiming a method is fully documented

| Area | Documentation additions needed | Main engineering dependency |
|---|---|---|
| Mesh geometry | orientation contract, signed/positive volume, non-planar faces, quality metrics, invalid-cell policy | #59, #423 |
| Gradients/reconstruction | WLS weighting, rank/conditioning criteria, boundary stencil policy, face reconstruction independence, limiter interaction | #461 N2 |
| Convection | multidimensional boundedness, limiter assumptions, deferred correction, consistency with mass flux, high-resolution failure modes | #461 N4 |
| Diffusion | explicit non-orthogonal/skew decomposition, correction limiting, high-quality vs difficult-mesh regimes | #461 N3, #423 |
| Time integration | BDF2 history/bootstrap, variable-step consistency, restart history and CFL semantics | #60, #461 N5/N6 |
| Pressure–velocity | exact algorithmic distinction between SIMPLE/SIMPLEC/PISO/PIMPLE/projection/coupled, pressure gauge and collocated flux consistency | #390, #461 N9 |
| Linear algebra | true residual vs preconditioned residual, null spaces, conditioning, block systems, AMG diagnostics | #492, #461 N8 |
| Turbulence | equation-level contracts, constants, wall-distance dependency and model-specific admissibility | #473, #461 N17 |
| Thermal/CHT | energy-form consistency, property derivatives, interface conservation, nonlinear source linearisation | #60, #382, #388 |
| Radiation | angular quadrature accuracy, optical-thickness regimes, S2S/P1/DOM model boundaries and coupling error | #461, #118 |
| Data/restart | schema compatibility, provenance, case/mesh/physics/numerics revision coupling and restart determinism | #60 |
| V&V | claim → oracle → metric → criterion → artifact traceability and qualification population boundaries | #118, #531 |

### 2. Source-to-document traceability

Chapter source mappings must be treated as auditable links, not illustrative filenames. Before the next documentation release, every mapped path should be checked against the repository and each important numerical claim should identify its executable verification or validation artifact.

### 3. Executable theory experiments

The current chapters contain small executable invariants. The next step is to add quantitative experiments where they materially improve understanding:

- observed-order regression for gradients, diffusion and convection;
- amplification-factor plots/tables for temporal schemes;
- Schur/AMG energy-contraction experiments;
- limiter response on smooth and steep fields;
- conservation/antisymmetry checks;
- radiation view-factor closure/reciprocity;
- CHT interface balance;
- dimensionless-number sanity checks.

These experiments are educational/verification evidence only; they must not be presented as solver qualification unless they execute the real CFDX implementation.

### 4. Scientific reference closure

Each chapter should distinguish:
- analytical oracle;
- manufactured oracle;
- discrete-exact oracle;
- published benchmark;
- external-code comparison.

The reference value itself is never CFDX evidence.

### 5. Planned next theory pass

The recommended order is:

**geometry/orientation → N2 gradients → N4 convection → N5/N6 time integration → N7–N9 solver/coupling → N10 difficult meshes → N11–N12 conservation/MMS → N17 turbulence → thermal/radiation → data/restart → final source-tree audit.**

This order follows the numerical dependency graph rather than the visual chapter numbering.
