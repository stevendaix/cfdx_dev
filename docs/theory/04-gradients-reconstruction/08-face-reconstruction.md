# 8. Face-Value Reconstruction

Gradient computation and face-value reconstruction are related but distinct operations.

A typical linear reconstruction is

\[
\phi_f=
\phi_P+

\nabla\phi_P\cdot(\mathbf x_f-\mathbf x_P).
\]

The gradient may come from Green–Gauss, least squares, weighted least squares or another method.

This separation is important: changing the gradient algorithm should not silently change the definition of the face reconstruction API.

Face reconstruction affects convective and diffusive fluxes, boundedness, conservation and nonlinear convergence. It therefore requires verification evidence of its own.
