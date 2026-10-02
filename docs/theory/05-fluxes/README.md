# 05 — Fluxes

**Status: REPOSITORY-GROUNDED.**

For a transported scalar,

$$
\frac{\partial(\rho\phi)}{\partial t}
+\nabla\cdot(\rho\mathbf u\phi)
=\nabla\cdot(\Gamma\nabla\phi)+S.
$$

The finite-volume equation is

$$
\frac{d}{dt}\int_{V_P}\rho\phi\,dV+
\sum_f\dot m_f\phi_f
=
\sum_f\Gamma_f(\nabla\phi)_f\cdot\mathbf S_f+S_PV_P.
$$

## Convective flux

$$
\dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f,
\qquad
F_{\phi,f}=\dot m_f\phi_f.
$$

The current implementation separates flux, interpolation, convection and convection assembly in core/numerics/flux.h, interpolation.h, convection.h and convection_assembly.*.

A convection scheme therefore includes mass-flux definition, face-state reconstruction, limiter policy and boundary treatment.

## Diffusive flux

$$
F_{d,f}=\Gamma_f\nabla\phi_f\cdot\mathbf S_f.
$$

On an orthogonal mesh this reduces to the centre-to-centre normal difference. On non-orthogonal meshes the correction depends on geometry and face-gradient reconstruction. The current Laplacian path is core/numerics/laplacian.h.

## Conservation

For an internal face the two cell contributions must cancel:

$$
F_{f,P}+F_{f,N}=0.
$$

This property is independent of nonlinear convergence and must be independently tested.

## Boundedness

Accuracy and boundedness are separate properties:

$$
\text{boundedness}\not\Rightarrow\text{second-order accuracy}.
$$

Limiter tests must therefore include both admissibility and smooth-field accuracy.

## V&V

Each flux family requires constant/linear consistency where applicable, face antisymmetry, conservation, MMS accuracy, boundary consistency, skew/non-orthogonal behaviour and boundedness/positivity where physically required.


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
