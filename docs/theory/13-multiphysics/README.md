# 13 — Multiphysics

**Status: REPOSITORY-GROUNDED.**

A coupled problem can be written

$$
F(y)=0,
\qquad y=(u,p,T,\alpha,k,\ldots).
$$

Its Jacobian has block structure:

$$
J=
\begin{bmatrix}
J_{uu}&J_{up}&J_{uT}\\
J_{pu}&J_{pp}&J_{pT}\\
J_{Tu}&J_{Tp}&J_{TT}
\end{bmatrix}.
$$

## Segregated coupling

A fixed-point iteration can be written

$$
y^{m+1}=G(y^m).
$$

Under-relaxation gives

$$
y^{m+1}\leftarrow(1-\alpha)y^m+\alpha G(y^m).
$$

Relaxation affects convergence; it does not correct an inconsistent discretisation.

## Strong coupling

A monolithic or block method targets the coupled Jacobian directly. Schur and field-split preconditioners exploit this block structure.

## Interfaces

For conjugate heat transfer, for example,

$$
T_1=T_2,
\qquad q_{1,n}+q_{2,n}=0.
$$

Other physics require different interface conditions and must not be reduced to this generic example.

## CFDX families

The current source tree contains CHT, FSI, radiation, VOF, multicomponent and turbulence coupling paths. These are distinct mathematical models.

## V&V

Each coupled model requires isolated component verification, interface conservation, coupling convergence, relaxation sensitivity, analytical/MMS subproblems where possible and independent benchmark data. Acceptance must identify independently recomputed coupling fluxes or residuals.
