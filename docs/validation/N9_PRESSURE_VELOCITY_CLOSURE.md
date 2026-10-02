# N9 — Pressure–velocity coupling closure

## Purpose

N9 closes the numerical-method chain between the finite-volume momentum/continuity operators and the pressure–velocity algorithms that consume them.

The target is a common pressure–velocity architecture in which SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step projection and the fully coupled (U-p) formulation share the same momentum assembly, pressure operator, conservative face-flux representation, convergence diagnostics and pressure-gauge policy.

This document is a **qualification contract**. It does not promote an algorithm to qualified merely because an implementation exists.

## 1. Discrete coupled problem

For an incompressible collocated formulation the assembled linearized system is represented as

[
\begin{bmatrix}
A_u & G \\
D & C
\end{bmatrix}
\begin{bmatrix}
u \\
p
\end{bmatrix}
=
\begin{bmatrix}
b_u \\
b_p
\end{bmatrix}.
]

For segregated pressure correction,

[
A_u u = b_u-Gp,
]

with

[
u=A_u^{-1}(b_u-Gp),
]

and the pressure Schur complement

[
S p =
b_p-D A_u^{-1}b_u.
]

CFDX must keep the following quantities distinct:

- (rAU=1/A_P): inverse integrated momentum diagonal;
- (dAU=V/A_P): physical velocity response coefficient;
- pressure-gradient contribution;
- conservative face mass flux;
- reconstructed velocity face flux.

The pressure correction and Rhie–Chow flux must use the same discrete response represented by (dAU), or an algebraically equivalent formulation.

## 2. Algorithm contracts

### SIMPLE

One momentum predictor followed by pressure correction and velocity/flux correction:

[
p^{k+1}=p^k+alpha_p\delta p,
]

[
u^{k+1}=u^k+alpha_u\delta u.
]

The pressure equation must be assembled from the same continuity flux used for the corrected solution.

### SIMPLEC

SIMPLEC must use its own consistent velocity-correction coefficient. It must not be accepted solely because the algorithm enum dispatches to the SIMPLE path.

Required evidence:

- same physical problem as SIMPLE;
- continuity and momentum gates;
- solution/QoI comparison;
- measured iteration history.

### PISO

PISO performs multiple pressure corrections within an outer iteration:

[
p^{(0)}\rightarrow p^{(1)}\rightarrow u^{(1)}
\rightarrow p^{(2)}\rightarrow u^{(2)}\rightarrow\cdots.
]

The number of pressure correctors is an explicit numerical parameter. The implementation must not silently reduce PISO to a single SIMPLE correction.

### PIMPLE

PIMPLE combines outer nonlinear iterations and PISO-like pressure corrections:

[
\text{time/outer iteration}
\rightarrow
\{\text{momentum + pressure corrections}\}
\rightarrow
\text{convergence test}.
]

The implementation must reuse the same momentum and pressure components rather than maintain a second solver stack.

### Fractional step / projection

For an incompressible projection,

[
u^*=u^n+\Delta t\,R_u,
]

[
\nabla^2\phi =
\frac{1}{\Delta t}\nabla\cdot u^*,
]

[
u^{n+1}=u^*-\Delta t\nabla\phi.
]

The exact pressure convention used by CFDX must be explicit. Projection is a separate algorithmic contract and must not be treated as an alias for PISO.

### Fully coupled (U-p)

The coupled path solves

[
\begin{bmatrix}
A_u&G\\D&C
\end{bmatrix}x=b
]

with the native Krylov/FieldSplit/MGR/Schur infrastructure from N8.

The coupled path must expose:

- velocity block;
- pressure block;
- Schur approximation;
- pressure gauge treatment;
- true residual (b-Ax);
- iteration/setup diagnostics.

## 3. Pressure–velocity consistency

The collocated-grid mass flux must be pressure-consistent.

For an internal face,

[
\phi_f =
\rho\,U_f\cdot S_f
-
\rho\,d_f\,(\nabla p)_f\cdot S_f
+
\text{non-orthogonal/skew correction}.
]

CFDX's Rhie–Chow implementation must be treated as a distinct numerical operator and verified independently from arithmetic interpolation of cell velocity.

Required invariants:

1. constant pressure does not create a pressure-driven flux;
2. adding a constant pressure offset does not change velocity or mass flux;
3. internal-face flux is antisymmetric between owner and neighbour;
4. pressure-corrected continuity is independently recomputed;
5. reconstructed cell-velocity continuity is reported separately;
6. the pressure-gauge operation does not introduce a physical pressure gradient.

## 4. Pressure null space and reference policies

For a pure-Neumann pressure problem,

[
p\rightarrow p+C
]

is a null-space transformation.

CFDX must explicitly select and report one gauge policy:

- reference cell/value;
- mean-pressure constraint;
- null-space projection.

A test must start from a non-zero arbitrary pressure field and demonstrate that the selected gauge removes only the null-space component.

A fixed-pressure boundary is not interchangeable with a pure-Neumann gauge constraint.

## 5. Common acceptance matrix

Every production pressure–velocity algorithm must be exercised through the same diagnostic contract.

| Case | SIMPLE | SIMPLEC | PISO | PIMPLE | Fractional step | Coupled |
|---|---:|---:|---:|---:|---:|---:|
| Couette | required | required | required | required | required | required |
| Poiseuille | required | required | required | required | required | required |
| Lid-driven cavity | required | required | required | required | required | required |
| External flow | required | required | required | required | N/A where inappropriate | required |
| Skew/non-orthogonal mesh | required | required | required | required | required where applicable | required |

A method that is physically inapplicable to a benchmark is marked **N/A with justification**, never converted into PASS.

## 6. Required evidence per run

Each run records at least:

### Algebraic convergence

[
\|r_{linear}\|_2,qquad
\|b-Ax\|_2,qquad
N_{linear}.
]

### Nonlinear convergence

[
\|R_u\|,qquad
\|R_c\|,qquad
\|\Delta u\|_\infty,qquad
\|\Delta p\|_\infty,qquad
N_{nonlinear}.
]

### Physical conservation

[
\epsilon_m =
\frac{|\sum_f\phi_f + S_m|}
{\phi_{ref}},
]

with cell-local and global forms reported independently.

For momentum,

[
R_m=\sum_f\mathbf F_f+\mathbf F_{source}.
]

### Solution/QoI evidence

Where an oracle exists:

[
L_2,quad L_\infty,quad Q,quad
\text{observed order}.
]

Residual convergence alone is never sufficient.

## 7. Sensitivity campaign

The qualification matrix must vary:

- mesh resolution;
- skewness;
- non-orthogonality;
- aspect ratio;
- Reynolds number;
- pressure-correction count;
- outer corrector count;
- under-relaxation;
- CFL/time step for transient algorithms;
- linear solver and preconditioner.

The result is a domain-of-validity statement, not an iteration-count ranking.

## 8. Reproducibility

A pressure–velocity result is reproducible only when the following are recorded:

- CFDX revision;
- mesh and mesh revision;
- physical properties;
- boundary conditions;
- algorithm;
- corrector counts;
- relaxation/CFL;
- linear solver/preconditioner;
- tolerances;
- initial state;
- execution backend.

Parallel reproducibility must be classified explicitly as:

1. bitwise;
2. numerical-equivalence within a documented tolerance;
3. physical-equivalence only.

## 9. Performance

Performance is diagnostic and is evaluated only after numerical correctness.

Report:

[
T=T_{assembly}+T_{linear}+T_{nonlinear}+T_{IO},
]

plus:

- nonlinear iterations;
- linear iterations;
- pressure/Schur setup cost;
- AMG setup and solve cost;
- memory;
- communication cost where available.

No performance result may weaken numerical gates.

## 10. Current CFDX status

The repository already contains implementations for:

- SIMPLE;
- SIMPLEC;
- PISO;
- PIMPLE;
- fractional-step;
- fully coupled (U-p);
- Schur/block preconditioning;
- pressure reference/gauge handling;
- Rhie–Chow flux construction;
- independent continuity and momentum diagnostics;
- Couette and cavity acceptance paths.

The executable closure now contains a common physical matrix in `tests/validation/test_n9_physical_matrix.cpp`. It runs SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step and coupled U-p on:

- Couette flow with an analytical profile gate;
- pressure-driven Poiseuille flow with an analytical profile gate;
- lid-driven cavity with physical convergence and cross-algorithm equivalence gates;
- a controlled skew/non-orthogonal Couette mesh.

The remaining boundary of N9 is the external-flow member of the matrix: the repository does not yet expose a pressure–velocity-coupled external-flow benchmark with a production-quality oracle for all six algorithms. That item therefore remains explicitly PARTIAL rather than being converted into a synthetic PASS.

Therefore N9 remains **PARTIAL** until the external-flow gate is executable and green on the exact final HEAD.

## 11. Definition of Done

N9 can move from PARTIAL to QUALIFIED only when:

- every applicable algorithm passes the common benchmark contract;
- pressure gauge/null-space behaviour is independently verified;
- Rhie–Chow pressure consistency has an independent executable gate;
- at least one controlled skew/non-orthogonal mesh is included;
- linear, nonlinear and physical convergence are reported separately;
- conservation is independently recomputed;
- analytical/reference QoIs are checked where applicable;
- refinement evidence exists where an order claim is made;
- no tolerance is relaxed to hide a defect;
- the exact final HEAD passes the required CI and validation-total gates;
- `NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json` agrees with the executable evidence.

## Related requirements

- #461 — Numerical methods maturity
- #390 — pressure–velocity coupling / Couette
- #405 — native AMG / FieldSplit
- #118 — validation campaign
- #530 — V&V governance

