# 07 — Pressure–Velocity Coupling

**Status: REPOSITORY-GROUNDED.**

For incompressible flow,

$$
\nabla\cdot\mathbf u=0.
$$

A generic linearised coupled system is

$$
\begin{bmatrix}A_u&G\\D&0\end{bmatrix}
\begin{bmatrix}u\\p\end{bmatrix}
=
\begin{bmatrix}b_u\\b_p\end{bmatrix}.
$$

Pressure acts as the constraint variable.

## Segregated algorithms

SIMPLE and SIMPLEC derive a pressure-correction system from an approximation of the momentum inverse. PISO performs multiple correction stages. PIMPLE combines outer nonlinear iterations with pressure corrections.

The exact CFDX update sequence, relaxation and face-pressure treatment must be documented from the implementation, not inferred from the algorithm name.

## Fractional step

A simplified projection method is

$$
u^*=u^n+\Delta t R(u^n,p^n),
$$

followed by

$$
\nabla^2p^{n+1}=\frac{\rho}{\Delta t}\nabla\cdot u^*,
$$

and velocity correction. The actual discrete operator must be traced to CFDX code.

## Schur complement

For

$$
A_uu+Gp=b_u,\qquad Du=b_p,
$$

eliminating u gives a pressure Schur system

$$
Sp=b_p-DA_u^{-1}b_u,
\qquad S=-DA_u^{-1}G
$$

up to the chosen sign convention.

CFDX has exact and approximate Schur infrastructure in core/linalg.

## Pressure null space

For a pure Neumann pressure problem,

$$
p'=p+C.
$$

A reference, mean constraint or explicit null-space treatment is required.

Relevant implementation families include physics/pressure_velocity.h, pressure_velocity_algorithms.h, steady_incompressible_solver.h and core/linalg/null_space.h.

## V&V

Report continuity imbalance, momentum residual, nonlinear and linear iterations, pressure reference policy, relaxation/CFL sensitivity, conservation and canonical QoIs. Couette and Poiseuille provide analytical baselines; cavity and external-flow cases test multidimensional coupling.
