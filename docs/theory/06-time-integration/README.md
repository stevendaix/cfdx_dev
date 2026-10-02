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
