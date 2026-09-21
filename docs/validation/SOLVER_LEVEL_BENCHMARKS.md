# CFDX solver-level benchmark suite

This document separates executable solver comparisons from component-level reference checks.

## Executable solver comparisons

| Case | CFDX path | Reference | Current CI role |
|---|---|---|---|
| Couette | finite-volume scalar diffusion | exact linear profile | PASS |
| Plane Poiseuille | finite-volume scalar diffusion/source | exact quadratic profile | PASS, observed order 2.0 |
| 1-D conduction | finite-volume scalar diffusion | exact linear profile | PASS within numerical tolerance |
| Ghia cavity Re=100/400 | steady incompressible SIMPLE | Ghia et al. centreline tables | active solver benchmark; tolerance is mesh-dependent |\n\nFor Ghia, the executable gate now includes Re=100 at 32/64/128 cells per side plus Re=400 at 64 cells per side. The Re=100 run reports observed mesh-convergence orders for the centreline velocity errors and explicitly checks zero normal velocity on the moving lid. The current momentum convection operator is first-order upwind; its `bounded_convection` option is a conservative continuity correction, not a second-order LUD/TVD scheme. Therefore the validation gate requires positive mesh convergence and does not claim second-order accuracy. A future LUD/TVD implementation should add a separate order-of-accuracy gate.
| two-region CHT | coupled energy solver/interface matching | resistance/flux continuity | PASS |
| radiation-energy equilibrium | DOM + energy coupling | uniform blackbody equilibrium | PASS |

The Couette/Poiseuille cases are solver-level verification of the finite-volume transport path used by momentum diffusion. They are not a complete pressure-velocity Navier-Stokes validation.

## Reference datasets and planned full CFD cases

### Taylor-Green vortex

For the 2-D laminar Taylor-Green vortex,

u = -U0 cos(x) sin(y) exp(-2 nu t)
v =  U0 sin(x) cos(y) exp(-2 nu t)

and the pressure field has the corresponding exact viscous-decay solution. A full CFDX transient Navier-Stokes comparison requires a transient pressure-velocity driver and periodic boundaries. Until those capabilities are present, Taylor-Green remains an operator/transient-verification target rather than a solver PASS.

### Blasius flat plate

The reference quantities are the classical laminar boundary-layer similarity solution, including skin friction and displacement/thickness quantities. A valid CFDX comparison must extract the wall shear and velocity profile from an actual external-flow mesh. The existing correlation checks are reference-oracle checks only.

### Turbulent channel

The target datasets are the Moser-Kim-Mansour DNS cases around Re_tau=180, 395 and 590. A valid comparison must compute statistically converged mean velocity and Reynolds stresses from a CFDX turbulent simulation. Reading the reference table alone is not a solver validation.

### Backward-facing step

The benchmark requires separation and reattachment location, wall pressure and wall shear. It should be added only when CFDX has a sufficiently general inlet/outlet boundary implementation and a converged multi-cell Navier-Stokes driver.

### Differentially heated cavity

The de Vahl Davis family is the reference for natural convection. A true CFDX comparison requires buoyancy in the momentum equation coupled to energy, plus Nusselt-number extraction. The current M1-M4 implementation does not yet expose that complete coupling.

### Graetz thermal entrance

The benchmark requires a developing thermal field in a channel/tube and comparison of local/mean Nusselt number against the analytical entrance-region solution. The current energy solver can verify the diffusion operator, but a developing-flow driver is still required.

### Cylinder

The low-Re cylinder cases require an actual 2-D external-flow mesh and pressure/velocity solver. Recommended QoIs are drag, lift and separation/recirculation length.

### NACA 0012 / RAE 2822

These should be enabled only after the compressible/external-flow infrastructure is available. Required comparisons include pressure coefficient, lift/drag and shock location where applicable.

### Compressible cases

Sod, isentropic nozzle and oblique shock become solver-level benchmarks once the compressible equations, EOS coupling, shock capturing and appropriate boundary conditions are executable.

## Acceptance rules

For every executable benchmark:

1. record mesh, physical parameters and solver controls;
2. report the raw QoI and reference value;
3. report absolute and relative error;
4. perform at least three grid levels where a spatial convergence study is meaningful;
5. report observed order;
6. check mass/energy conservation independently;
7. never convert an unavailable solver capability into a PASS.

A benchmark may therefore have three states:

- **PASS**: CFDX produced the requested physical solution and met the quantitative acceptance criterion;
- **FAIL**: CFDX executed but disagreed with the reference or violated conservation/convergence requirements;
- **BLOCKED**: the required physical/numerical capability is not implemented yet.

This status model is intentional: a component oracle is never promoted to a solver-validation PASS.


### Solver diagnostics

The incompressible solver now distinguishes the linear-system residuals of the intermediate momentum/pressure solves from the **final physical momentum equation residual** after pressure correction. The latter is reassembled from the corrected U/p fields and is part of the nonlinear convergence gate. The iteration history also exposes a dimensionless continuity metric normalized by a characteristic mass flux, so benchmark cases can compare conservation quality without depending only on raw kg/s values.
