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
