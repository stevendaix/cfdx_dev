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
