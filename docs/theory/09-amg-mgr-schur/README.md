> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 09 — AMG, MGR and Schur Complement Methods

**Status: REPOSITORY-GROUNDED.**

Multigrid targets different error frequencies. A relaxation step damps high-frequency error while coarse-grid correction removes smooth error.

A generic cycle is

$$
x\leftarrow S_{pre}(x,b)
\rightarrow Rr
\rightarrow A_c^{-1}
\rightarrow Px_c
\rightarrow S_{post}.
$$

## Galerkin construction

$$
A_c=RAP.
$$

The exact restriction, prolongation, coarsening and smoothing define the method.

Current CFDX files include amg_preconditioner.h, amg_wrapper.h and chebyshev_smoother.h, plus coupled and block infrastructure.

## Energy contraction

For an SPD operator,

$$
\rho_E=\frac{\lVert e_{k+1}\rVert_A}{\lVert e_k\rVert_A}.
$$

An AMG gate should use an independently recomputed measure such as energy contraction, not only an internal residual.

## MGR

MGR reduces selected variables while retaining block structure. It is conceptually different from scalar AMG coarsening.

Current implementation: mgr_preconditioner.h and the coupled/block operators.

## Schur complement

For

$$
A=\begin{bmatrix}A_{uu}&A_{up}\\A_{pu}&A_{pp}\end{bmatrix},
$$

$$
S=A_{pp}-A_{pu}A_{uu}^{-1}A_{up}.
$$

CFDX provides exact_schur.h, block_schur.h, schur_approximation.h, simplerc_schur.h and coupled_amg_schur.h.

The exact Schur path is an oracle for approximation quality, not a production performance claim.

## V&V

Use mesh-size, anisotropy and coefficient-scaling sweeps. Record hierarchy diagnostics, true residual, energy contraction where applicable, iterations, setup/solve cost and failure cases. Do not infer PETSc/Hypre parity from option names.


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
