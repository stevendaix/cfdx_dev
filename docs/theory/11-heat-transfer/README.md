> **Executable source:** [chapter.py](chapter.py)  
> This README is navigation and chapter contract. Quantitative theory belongs to the Python percent source.

# 11 — Heat Transfer

**Status: REPOSITORY-GROUNDED.**

A representative thermal transport equation is

$$
\rho c_p\left(
\frac{\partial T}{\partial t}+\mathbf u\cdot\nabla T
\right)
=
\nabla\cdot(k\nabla T)+S_T.
$$

The conserved variable and equation change for compressible formulations; Theory follows the actual CFDX equation path.

## Conduction

Fourier's law is

$$
\mathbf q=-k\nabla T.
$$

The outward face heat rate is

$$
Q_f=-k_f(\nabla T)_f\cdot\mathbf S_f.
$$

The sign convention must remain consistent at interfaces.

## Convection

Thermal convection has the same finite-volume structure as scalar transport:

$$
\Phi_{h,f}=\dot m_fh_f.
$$

Therefore thermal accuracy depends on the common flux and reconstruction contracts.

## Properties

Current thermophysical families include thermophysical_models.h, transport_models.h and equation_of_state.h. Temperature-dependent properties introduce nonlinear coupling.

## Conjugate heat transfer

At an ideal fluid-solid interface,

$$
T_f=T_s,
\qquad
q_{f,n}+q_{s,n}=0
$$

with the sign defined by a common interface orientation.

Current solver family: cht_solver.h.

## Buoyancy

The Boussinesq approximation can represent density variation only where it materially affects buoyancy. The reference temperature and approximation assumptions must be explicit. Current source: boussinesq.h.

## V&V

Use analytical conduction, transient diffusion MMS, convection-diffusion, interface flux balance, property-regression tests and global energy conservation. Validation results must identify material data and boundary conditions.


## Scientific explanation standard

Every major equation or method in this chapter must explain: physical/mathematical motivation; variables, units and sign conventions; derivation; FVM/discrete formulation; actual CFDX algorithm/code path; analytical or numerical example; errors, limitations and sensitivities; exact source files; executable verification tests; benchmark/V&V evidence; and bibliography with stable identifiers.

Comparison tables are summaries after the mathematics, never substitutes for it. Status must remain **Implemented / Verified / Validated / Qualified** and documentation must never promote numerical maturity.