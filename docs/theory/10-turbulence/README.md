# 10 — Turbulence

**Status: REPOSITORY-GROUNDED; qualification remains a V&V question.**

Reynolds decomposition is

$$
u_i=\overline u_i+u_i'.
$$

Averaging introduces the Reynolds stress

$$
-\rho\overline{u_i'u_j'},
$$

which requires closure.

## Eddy-viscosity RANS

A representative closure is

$$
-\rho\overline{u_i'u_j'}
=
2\mu_tS_{ij}-\frac23\rho k\delta_{ij}.
$$

The definition of turbulent viscosity and transported variables is model-specific.

## Spalart–Allmaras

SA transports a modified viscosity variable with production, destruction, diffusion and wall-dependent terms. Current source family: physics/spalart_allmaras.h and turbulence transport/model files.

## k–ω SST

SST transports k and ω and blends near-wall and free-shear behaviour. Current source family: physics/sst_solver.h, turbulence_models.h and turbulence_transport.h.

## Wall distance

Wall distance enters near-wall closure and wall functions. Relevant files are wall_distance.h, wall_distance_turbulence.h and wall_functions.h. The wall-distance algorithm must be verified separately from turbulence-model validation.

## LES and hybrid methods

LES resolves part of the turbulent spectrum and models subgrid stresses. Hybrid methods combine resolved and modelled regions and are strongly mesh/time-resolution dependent.

## Numerical coupling

Turbulence equations use the common reconstruction, convection, diffusion, temporal and linear-algebra contracts. Source stiffness and positivity make their convergence requirements more restrictive in some regimes.

## V&V

Separate equation verification, transport MMS where meaningful, positivity/boundedness, canonical turbulent benchmarks and independent physical validation. Every benchmark records model constants, wall treatment, mesh resolution and inlet turbulence specification.


## Scientific explanation standard

This chapter is part of the CFDX Theory course. A short descriptive statement or comparison table is not sufficient for a major numerical or physical method.

For every important equation or method, the final documentation must follow this chain:

$$
\boxed{
\text{motivation}
\rightarrow
\text{definitions}
\rightarrow
\text{derivation}
\rightarrow
\text{FVM/discretisation}
\rightarrow
\text{CFDX algorithm}
\rightarrow
\text{example}
\rightarrow
\text{error/limitations}
\rightarrow
\text{tests}
\rightarrow
\text{benchmark/V\&V}
\rightarrow
\text{bibliography}
}
$$

### Required explanation

1. Explain why the equation or method is needed and what physical/mathematical problem it solves.
2. Define every symbol, tensor/vector/scalar, unit and sign convention.
3. Derive the formula sufficiently for a reader to reproduce the result.
4. Show the finite-volume or discrete transformation where applicable.
5. Explain the actual CFDX computational sequence, not only the textbook algorithm.
6. Give a small analytical, manufactured-solution or numerical example whenever meaningful.
7. Explain truncation error, consistency, stability, conditioning, boundedness, conservation and sensitivity as applicable.
8. Identify the exact implementation files and the data passed between stages.
9. Link the mathematical property to executable verification tests.
10. Identify the benchmark and V&V evidence, including scope and limitations.
11. Cite the scientific literature and record stable bibliographic identifiers in the project bibliography.

Comparison tables remain useful, but they are summaries **after** the mathematical explanation and never substitutes for it.

### Evidence vocabulary

- **Implemented** — an executable code path exists.
- **Verified** — a defined mathematical/software property has executable evidence.
- **Validated** — comparison exists against an independent physical or trusted reference.
- **Qualified** — the declared capability is demonstrated over an explicit scope.

A chapter may be scientifically complete while a capability remains unqualified. Documentation must never promote numerical maturity.
