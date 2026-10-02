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

eliminating u gives

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
