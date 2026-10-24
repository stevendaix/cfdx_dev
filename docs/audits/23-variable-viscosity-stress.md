# F11 — Variable-viscosity momentum stress

Status: implementation and unit qualification added; exact-HEAD CI pending.

## Governing term

For a Newtonian incompressible momentum equation with effective dynamic viscosity \(\mu_{eff}(x)\), the viscous stress divergence is

\[
\nabla\cdot\left[\mu_{eff}\left(\nabla\mathbf U+(\nabla\mathbf U)^T\right)\right].
\]

The scalar diffusion operator \(\nabla\cdot(\mu_{eff}\nabla U_i)\) does not contain the transpose contribution

\[
\left[\nabla\cdot(\mu_{eff}(\nabla\mathbf U)^T)\right]_i
=\sum_j \frac{\partial}{\partial x_j}
\left(\mu_{eff}\frac{\partial U_j}{\partial x_i}\right).
\]

Even when \(\nabla\cdot\mathbf U=0\), a spatially variable viscosity leaves the term
\(\sum_j (\partial_j\mu_{eff})(\partial_i U_j)\), which is generally nonzero.

## Discrete contract

- Reconstruct cell gradients of each velocity component with the same vector velocity boundary conditions converted component-wise to scalar boundary conditions.
- Interpolate \(\mu_{eff}\) and velocity gradients to each face using the owner/neighbour normal distances.
- Compute the integrated face flux \(F_i^{T}=\mu_f\sum_j(\partial_i U_j)_f S_{f,j}\).
- Assemble the divergence conservatively: add \(F_i^T\) to the owner RHS and subtract it from the neighbour RHS. Internal-face contributions must cancel globally.
- Include the same term in the final independently assembled momentum residual; otherwise convergence acceptance would evaluate a different equation from the production solve.
- Keep the scalar diffusion contribution implicit and the transpose term explicit. Recompute it from the current velocity each nonlinear iteration; do not silently omit it in production.

## Independent manufactured field

Use \(\mathbf U=(a y,0,0)\), so \(\nabla\cdot\mathbf U=0\), and \(\mu_{eff}(x)=1+x\). The y-momentum transpose contribution is

\[
\left[\nabla\cdot(\mu_{eff}(\nabla\mathbf U)^T)\right]_y
=\frac{\partial}{\partial x}(\mu_{eff}a)=a.
\]

On interior cells of an orthogonal unit-width channel, the independent face oracle is
\(F_y^T(x_f)=\mu_{eff}(x_f)a S_x\), and the cell divergence must equal \(aV_c\). The registered test checks both face fluxes and the RHS assembled into the momentum equation.

## Limits

This test qualifies the variable-viscosity transpose term on an orthogonal affine manufactured field. It does not yet establish second-order convergence on skewed/polyhedral meshes, coupled turbulence-model convergence, or a full variable-viscosity Navier–Stokes MMS campaign. Those remain separate qualification work.
