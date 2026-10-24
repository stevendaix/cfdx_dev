# F06 — Non-orthogonal diffusion in scalar transport

Status: shared transport assembly/reconstruction correction added; exact-HEAD CI pending.

## Discrete flux

For an internal face with owner-to-neighbour vector \(\mathbf d=\mathbf C_N-\mathbf C_P\) and area vector \(\mathbf S_f\), the transport diffusion flux is decomposed into the same implicit normal conductance used by the scalar matrix and an explicit tangential correction:

\[
F_{diff,f}=D_f(\psi_N-\psi_P)
+\Gamma_f\nabla\psi_f\cdot
\left(\mathbf S_f-\alpha_f\mathbf d\right),
\qquad
\alpha_f=\frac{|\mathbf S_f|^2}{\mathbf S_f\cdot\mathbf d}.
\]

The implementation uses \(D_f=|S_f|/(d_P/\Gamma_P+d_N/\Gamma_N)\) for variable diffusion and the corresponding constant-coefficient form. The tangential term uses the same face conductance and owner/neighbour distance weighting as the implicit matrix. The explicit correction is assembled conservatively: add the oriented correction to the owner RHS and subtract it from the neighbour RHS.

## Shared implementation

- `resolve_scalar_boundary_condition` centralizes face/patch boundary-value resolution.
- `compute_scalar_gradient_with_boundary` reconstructs a scalar gradient using the same face-specific values and conditions as the equation assembly.
- `compute_nonorthogonal_diffusion_face_fluxes` returns one oriented tangential diffusion flux per face.
- `assemble_scalar_equation` uses the correction as deferred RHS, and fails closed on a non-orthogonal diffusive mesh when the current scalar field is not supplied.
- `reconstruct_scalar_transport_flux` uses the same correction helper and the same normal conductance, so the reconstructed flux and assembled operator represent the same discrete diffusion term.

Zero-gradient boundary faces do not receive an independently imposed tangential flux; fixed-value and fixed-gradient boundaries use the resolved boundary state. Empty patches are excluded.

## Verification

The N3 transport campaign adds an affine scalar field on a skewed, nonuniform hexahedral mesh. It checks:
- the corrected assembled diffusion residual is near zero for an affine field with exact Laplacian zero;
- the reconstructed transport flux has the same cellwise divergence as the corrected operator;
- the orthogonal-only residual is nonzero on the same mesh, demonstrating that the tangential correction is active.

The test uses exact operator residuals, not relaxed solver tolerances.

## Limits

The first regression covers an inclined, nonuniform hexahedral mesh. A separate polyhedral MMS, observed-order study, and qualification of every Poisson/reconstruction/deferred-correction caller remain required before the method can be promoted beyond DIAGNOSTIC/PARTIAL.
