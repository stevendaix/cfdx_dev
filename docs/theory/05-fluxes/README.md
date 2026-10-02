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
