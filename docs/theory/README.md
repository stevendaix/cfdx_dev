# CFDX Theory — Numerical Methods Course

CFDX Theory is the scientific course and mathematical reference for CFDX. It is deliberately separate from the Developer implementation contract and from V&V acceptance status.

The Theory site is organised as a **course**, not as a feature list. Each chapter moves from physics to mathematics, then to finite-volume discretisation, numerical analysis, implementation, executable evidence and limitations.

## Master learning path

```
00 Foundations
   ↓
01 Conservation laws
   ↓
02 Finite-volume method
   ↓
03 Meshes
   ↓
04 Gradients & reconstruction
   ↓
05 Fluxes
   ↓
06 Time integration
   ↓
07 Pressure–velocity coupling
   ↓
08 Linear algebra
   ↓
09 AMG / MGR / Schur
   ↓
10 Turbulence
   ↓
11 Heat transfer
   ↓
12 Radiation
   ↓
13 Multiphysics
   ↓
14 Numerical analysis
   ↓
15 Verification & validation
   ↓
16 Data model & file formats
   ↓
17 Computational chain
   ↓
18 Source-tree physics audit
```

The chapters are intentionally sequential: later numerical methods rely on definitions established earlier.

## Detailed chapter navigation

| # | Chapter | Main purpose |
|---|---|---|
| 00 | Foundations | Build the continuum-mechanics, dimensional and physical basis |
| 01 | Conservation laws | Derive mass, momentum, energy and species equations |
| 02 | Finite-volume method | Transform conservation laws into locally conservative algebra |
| 03 | Meshes | Define topology, geometry, orientation and mesh-quality metrics |
| 04 | Gradients & reconstruction | Derive Green–Gauss, LS/WLS, reconstruction and limiting |
| 05 | Fluxes | Build mass, convection, diffusion, pressure and energy fluxes |
| 06 | Time integration | Derive Euler, CN, BDF2, CFL and timestep control |
| 07 | Pressure–velocity coupling | Derive segregated and coupled incompressible algorithms |
| 08 | Linear algebra | Explain sparse systems, Krylov methods, residuals and preconditioning |
| 09 | AMG / MGR / Schur | Explain block systems, Schur complements and multigrid |
| 10 | Turbulence | Derive RANS, LES and hybrid turbulence modelling |
| 11 | Heat transfer | Derive thermal transport, conduction and CHT |
| 12 | Radiation | Derive S2S, P1 and DOM radiation models |
| 13 | Multiphysics | Explain coupled nonlinear systems and interfaces |
| 14 | Numerical analysis | Establish consistency, stability, convergence and error theory |
| 15 | Verification & validation | Establish the evidence chain from equations to qualification |
| 16 | Data model | Explain the authoritative case/checkpoint/output representation |
| 17 | Computational chain | Trace a case through the complete CFDX execution pipeline |
| 18 | Source-tree audit | Connect physical equations to actual implementation files and tests |

Each chapter README contains its **full subsection register**. The executable `chapter.py` is the quantitative source for designated theory experiments.

## Standard structure of every subsection

Every scientific subsection must contain, in this order where applicable:

1. **Physical motivation**
2. **Problem definition**
3. **Variables and fields**
4. **Dimensions and SI units**
5. **Coordinate/sign conventions**
6. **Assumptions**
7. **Continuous governing equation**
8. **Integral/control-volume form**
9. **Finite-volume discretisation**
10. **Algebraic formulation**
11. **Numerical properties**
12. **Implementation consequences**
13. **CFDX source mapping**
14. **Executable verification**
15. **Benchmark/reference comparison**
16. **Limitations and failure modes**
17. **Improvement paths**
18. **Bibliography**

A table may summarise a result, but it must never replace the derivation.

## Mandatory equation contract

For every important equation, the documentation must answer:

- Where does the equation come from physically?
- What assumptions are required?
- What does every symbol mean?
- What are the units?
- What is the sign/orientation convention?
- How does the continuous equation become a control-volume balance?
- How is each face/volume term approximated?
- What algebraic coefficient or residual contribution results?
- Is the discretisation conservative?
- What consistency/order is expected?
- What stability/boundedness restrictions apply?
- Which CFDX source files implement it?
- Which executable test verifies it?
- Which benchmark or independent reference validates it, if applicable?
- What remains unverified or unqualified?

## Numerical-method traceability

The intended chain is:

```
physical law
    ↓
continuous equation
    ↓
integral balance
    ↓
geometric/discrete operator
    ↓
algebraic system
    ↓
nonlinear/linear solver
    ↓
diagnostics
    ↓
verification
    ↓
validation
    ↓
qualification
```

This is the central organising principle of the Theory documentation.

## Status vocabulary

- **Implemented** — an actual code path exists.
- **Verified** — a defined mathematical property has executable evidence.
- **Validated** — an appropriate physical/reference comparison has been completed.
- **Qualified** — the declared population and acceptance gates support the intended engineering scope.

These terms must never be collapsed into a single “supported” label.

## Executable documentation

Every chapter has an executable Python source where appropriate. The source uses Jupytext percent cells and is rendered by the Sphinx/MyST-NB documentation pipeline.

The executable source is the source of truth for numerical examples. Generated figures, tables and HTML are publication outputs.

## Scientific integrity rules

1. No invented CFDX results.
2. No reference-oracle value presented as a CFDX result.
3. No residual-only physical validation.
4. No implementation-to-qualification shortcut.
5. No tolerance relaxation to manufacture a green result.
6. No disabled failing benchmark to manufacture maturity.
7. Every quantitative result records its mesh/resolution, physics, numerical settings and software revision where relevant.
8. External code comparisons are documented as methodological comparisons unless an independent validation oracle exists.

## Relationship with the other documentation domains

```
THEORY
  mathematical meaning
       ↓
DEVELOPER
  implementation contract
       ↓
V&V
  executable evidence / status
       ↓
USER
  supported workflow
```

Information should be linked across these domains rather than copied.

## Current scope

The repository now contains a detailed subsection structure for all 19 Theory chapters and executable chapter sources. The next refinement of each subsection is expected to replace generic scaffolding with repository-grounded derivations, figures, quantitative experiments and exact code/test references.

**Important:** documentation maturity does not change the N1–N17 numerical maturity recorded by the engineering roadmap, in particular N2 remains subject to its existing V&V gates.
