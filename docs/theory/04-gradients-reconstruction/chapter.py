# %% [markdown]
# Gradients and Reconstruction
#
# ## 4.1 Taylor foundation
#
# For a smooth scalar,
# \[
# \phi(\mathbf x_P+\mathbf d)=
# \phi_P+\nabla\phi_P\cdot\mathbf d+
# \frac12\mathbf d^T\mathbf H_P\mathbf d+O(h^3).
# \]
# Therefore a neighbour relation is
# \[
# \Delta\phi_N=\nabla\phi_P\cdot\mathbf d_N+O(h^2).
# \]
#
# ## 4.2 Green--Gauss
#
# Gauss' theorem gives
# \[
# \int_{V_P}\nabla\phi\,dV
# =\oint_{\partial V_P}\phi\mathbf n\,dA.
# \]
# Hence
# \[
# \boxed{\nabla\phi_P\approx
# \frac1{V_P}\sum_f\phi_f\mathbf S_f}.
# \]
# With linear face interpolation,
# \[
# \phi_f=(1-w_f)\phi_P+w_f\phi_N.
# \]
# Constant-field exactness follows from \(\sum_f\mathbf S_f=0\).
#
# ## 4.3 Least squares
#
# Define
# \[
# \mathbf d_N=\mathbf x_N-\mathbf x_P,\qquad
# \Delta\phi_N=\phi_N-\phi_P.
# \]
# Minimise
# \[
# J(\mathbf g)=\sum_Nw_N
# (\Delta\phi_N-\mathbf g\cdot\mathbf d_N)^2.
# \]
# The normal equations are
# \[
# A\mathbf g=\mathbf b,
# \]
# \[
# A=\sum_Nw_N\mathbf d_N\mathbf d_N^T,\qquad
# \mathbf b=\sum_Nw_N\mathbf d_N\Delta\phi_N.
# \]
# In d dimensions,
# \[
# \operatorname{rank}(A)=d
# \]
# is required for an unconstrained unique solution.
#
# ## 4.4 Weighted LS and conditioning
#
# A common family is
# \[
# w_N=|\mathbf d_N|^{-p}.
# \]
# For SPD A,
# \[
# \kappa(A)=\frac{\lambda_{\max}(A)}{\lambda_{\min}(A)}.
# \]
# Perturbation sensitivity obeys the classical estimate
# \[
# \frac{\|\delta\mathbf g\|}{\|\mathbf g\|}
# \lesssim\kappa(A)\frac{\|\delta\mathbf b\|}{\|\mathbf b\|}.
# \]
# A robust implementation must diagnose rank deficiency and severe conditioning.
#
# ## 4.5 Vertex reconstruction
#
# A generic weighted vertex value is
# \[
# \phi_v=
# \frac{\sum_{c\in C(v)}w_{cv}\phi_c}
# {\sum_{c\in C(v)}w_{cv}},
# \qquad
# w_{cv}=|\mathbf x_v-\mathbf C_c|^{-p}.
# \]
# Weight choice and stencil definition are part of the numerical method.
#
# ## 4.6 Boundary reconstruction
#
# A physical boundary does not supply an ordinary neighbouring cell. Possible policies are constrained reconstruction, one-sided stencils or ghost/equivalent values. The policy must be compatible with the mathematical boundary condition and verified separately from the interior stencil.
#
# ## 4.7 Face reconstruction
#
# Gradient calculation and face reconstruction are separate contracts:
# \[
# \phi_f=\mathcal R_f(\phi_P,\phi_N,\nabla\phi_P,\nabla\phi_N,\mathbf x_f,\ldots).
# \]
# This separation is essential for MUSCL/TVD extensions.
#
# ## 4.8 Limiting
#
# A generic limited reconstruction is
# \[
# \phi_f=\phi_P+
# \alpha_P\nabla\phi_P\cdot(\mathbf x_f-\mathbf x_P),
# \qquad0\le\alpha_P\le1.
# \]
# A limiter must preserve the declared smooth-region order while enforcing its declared local bounds.
#
# ## 4.9 Verification
#
# For a linear manufactured field
# \[
# \phi=a+b x+c y+d z,
# \qquad
# \nabla\phi=(b,c,d)^T.
# \]
# For a smooth nonlinear field:
# \[
# E_h=
# \left[
# \frac{\sum_PV_P|\nabla\phi_h-\nabla\phi|^2}
# {\sum_PV_P|\nabla\phi|^2}
# \right]^{1/2},
# \qquad
# p_{obs}=\frac{\ln(E_h/E_{h/2})}{\ln2}.
# \]
# Controlled refinement is mandatory for an order claim.
#
# ## 4.10 CFDX traceability
#
# Gradient: src/cfdx/core/numerics/gradient.h  
# Stencil: src/cfdx/core/numerics/gradient_stencil.h  
# Least squares: src/cfdx/core/fvm/least_squares_gradient.h  
# Interpolation: src/cfdx/core/numerics/interpolation.h  
# Verification: tests/validation/test_gradient_verification.cpp and tests/validation/test_polyhedral_gradient_campaign.cpp
#
# N2 remains PARTIAL / not qualified until the declared polyhedral accuracy, conditioning, boundary-reconstruction and face-reconstruction gates are demonstrated.
#
# %%
from __future__ import annotations
import numpy as np
d=np.array([[1.,0.],[0.,1.],[-1.,0.],[0.,-1.]])
g=np.array([2.,-3.])
A=d.T@d
b=d.T@(d@g)
assert np.linalg.matrix_rank(A)==2
assert np.allclose(np.linalg.solve(A,b),g)
