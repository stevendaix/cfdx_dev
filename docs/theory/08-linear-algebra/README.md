> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 08 — Linear Algebra

**Status: REPOSITORY-GROUNDED.**

Finite-volume assembly produces

$$
A\mathbf x=\mathbf b.
$$

For an approximate solution,

$$
r=b-Ax,
$$

and the algebraic error satisfies

$$
Ae=r.
$$

A small residual does not by itself imply a small solution error when A is poorly conditioned.

## Krylov methods

CG is intended for symmetric positive-definite systems. BiCGStab and GMRES address broader nonsymmetric systems. Current implementations include cg_solver.h, bicgstab_solver.h, gmres_solver.h, krylov_controls.h and krylov_reductions.h.

## Preconditioning

With M,

$$
M^{-1}Ax=M^{-1}b.
$$

The preconditioner changes the numerical path while preserving the target solution when used consistently.

CFDX separates solver and preconditioner contracts and contains AMG, block, Schur, matrix-free and mixed-precision infrastructure.

## Conditioning

$$
\kappa(A)=\lVert A\rVert\lVert A^{-1}\rVert.
$$

Conditioning affects the relation between residual, perturbation and solution error.

For SPD operators, the energy norm

$$
\lVert e\rVert_A=\sqrt{e^TAe}
$$

is especially relevant to multigrid contraction.

## Required solver contract

Every production linear solve should identify operator class, scaling, Krylov method, preconditioner, stopping criteria, true residual recomputation and failure/stagnation reason.

V&V must report true residual reduction, iterations, setup cost, solve cost, memory and reproducibility. Functional CPU/MPI/GPU execution must not be confused with numerical equivalence or performance.


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
