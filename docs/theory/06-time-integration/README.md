> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 06 — Time Integration

**Status: REPOSITORY-GROUNDED.**

After spatial discretisation,

$$
M\frac{d\mathbf u}{dt}=\mathbf R(\mathbf u,t).
$$

The temporal method acts on this semi-discrete system, so spatial and temporal errors must be separated.

## Backward Euler

$$
\frac{u^{n+1}-u^n}{\Delta t}=R(u^{n+1},t^{n+1}).
$$

It is first order in time.

## Crank–Nicolson

$$
\frac{u^{n+1}-u^n}{\Delta t}
=\frac12[R^{n+1}+R^n].
$$

For smooth problems it is second order under the usual regularity assumptions.

## BDF2

$$
\frac{3u^{n+1}-4u^n+u^{n-1}}{2\Delta t}=R^{n+1}.
$$

Variable-step coefficients must depend on the actual step history.

## Runge–Kutta

$$
u^{(i)}=u^n+\Delta t\sum_{j<i}a_{ij}R(u^{(j)}),
\qquad
u^{n+1}=u^n+\Delta t\sum_i b_iR(u^{(i)}).
$$

Stability depends on the amplification polynomial and the spatial eigenvalues.

Current numerical contracts are in core/numerics/temporal.h; higher-level temporal and adaptive controls are in physics/temporal.h, adaptive_cfl.h, local_time_stepping.h and low_storage_time_integration.h.

## Restart

Multi-step methods require their history state. A restart campaign must compare a continuous run with a run split at a checkpoint and demonstrate equivalent physical and temporal state.

## V&V

Use analytical decay, transient diffusion MMS, smooth advection, temporal refinement at fixed spatial resolution, mixed space/time refinement and restart equivalence. Report observed temporal order separately from nonlinear convergence.


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
