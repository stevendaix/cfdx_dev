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
