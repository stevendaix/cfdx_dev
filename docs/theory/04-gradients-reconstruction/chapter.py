# %% [markdown]
# Gradients and Reconstruction
#
# ## 4.1 Taylor foundation
#
# For a smooth scalar,
# [
# phi(mathbf x_P+mathbf d)=
# phi_P+
ablaphi_Pcdotmathbf d+
# rac12mathbf d^Tmathbf H_Pmathbf d+O(h^3).
# ]
# Therefore a neighbour relation is
# [
# Deltaphi_N=
ablaphi_Pcdotmathbf d_N+O(h^2).
# ]
#
# ## 4.2 Green--Gauss cell gradient
#
# Gauss' theorem gives
# [
# int_{V_P}
ablaphi,dV
# =oint_{partial V_P}phimathbf n,dA.
# ]
# Hence
# [
# oxed{
ablaphi_Papprox
# rac1{V_P}sum_fphi_fmathbf S_f}.
# ]
# With linear face interpolation,
# [
# phi_f=(1-w_f)phi_P+w_fphi_N.
# ]
# Constant-field exactness follows from (sum_fmathbf S_f=0).
#
# The stronger linear-consistency identity is
# [
# sum_fmathbf C_fotimesmathbf S_f=V_P I
# ]
# for a closed polyhedron under the corresponding geometric convention. The identity does not by itself make a two-point face interpolation exact: on skewed polyhedra (mathbf C_f) need not lie on the owner-neighbour line.
#
# ## 4.3 Least squares
#
# Define
# [
# mathbf d_N=mathbf x_N-mathbf x_P,qquad
# Deltaphi_N=phi_N-phi_P.
# ]
# Minimise
# [
# J(mathbf g)=sum_Nw_N
# (Deltaphi_N-mathbf gcdotmathbf d_N)^2.
# ]
# The normal equations are
# [
# Amathbf g=mathbf b,
# ]
# [
# A=sum_Nw_Nmathbf d_Nmathbf d_N^T,qquad
# mathbf b=sum_Nw_Nmathbf d_NDeltaphi_N.
# ]
# In (d) dimensions,
# [
# operatorname{rank}(A)=d
# ]
# is required for an unconstrained unique solution.
#
# ## 4.4 Weighted LS, conditioning and fallback
#
# A common family is
# [
# w_N=|mathbf d_N|^{-p}.
# ]
# For SPD (A),
# [
# kappa_2(A)=rac{lambda_{max}(A)}{lambda_{min}(A)}.
# ]
# A robust implementation should report rank and conditioning rather than silently applying an unstable inverse. In a degenerate stencil the fallback policy must be explicit: enlarge the stencil, use a boundary-aware stencil, or reject the gradient.
#
# ## 4.5 Vertex Green--Gauss
#
# A generic weighted vertex value is
# [
# phi_v=
# rac{sum_{cin C(v)}w_{cv}phi_c}
# {sum_{cin C(v)}w_{cv}},
# qquad
# w_{cv}=|mathbf x_v-mathbf C_c|^{-p}.
# ]
# Face values can then be reconstructed from vertex values and integrated in the Green--Gauss sum. This changes the reconstruction contract; it must therefore be qualified independently from cell-based Green--Gauss.
#
# ## 4.6 Boundary reconstruction
#
# A physical boundary does not supply an ordinary neighbouring cell. Possible policies are constrained reconstruction, one-sided stencils or ghost/equivalent values. For a Dirichlet condition,
# [
# phi|_Gamma=phi_Gamma,
# ]
# the boundary value must enter the reconstruction consistently with the distance from (mathbf C_P) to the face. For a Neumann condition,
# [
# 
ablaphicdotmathbf n=q_Gamma,
# ]
# the normal derivative constraint must not be confused with an arbitrary zero-gradient fallback.
#
# ## 4.7 Face reconstruction independent of gradient
#
# [
# phi_f=mathcal R_f(phi_P,phi_N,
ablaphi_P,
ablaphi_N,mathbf x_f,ldots).
# ]
# The face reconstruction can therefore be first-order, linear, MUSCL or limited without changing the gradient API. This separation is important for N2/N3 qualification.
#
# ## 4.8 Limiting
#
# A generic limited reconstruction is
# [
# phi_f=phi_P+
# alpha_P
ablaphi_Pcdot(mathbf x_f-mathbf x_P),
# qquad0lealpha_Ple1.
# ]
# For a smooth solution the limiter should approach unity at the expected asymptotic rate; near extrema or discontinuities it may intentionally reduce the formal order to maintain the declared boundedness property.
#
# ## 4.9 Observed order
#
# For a smooth nonlinear field:
# [
# E_h=
# left[
# rac{sum_PV_P|
ablaphi_h-
ablaphi|^2}
# {sum_PV_P|
ablaphi|^2}
# ight]^{1/2},
# qquad
# p_{obs}=rac{ln(E_h/E_{h/2})}{ln2}.
# ]
# With arbitrary refinement ratio (r=h_c/h_f),
# [
# p_{obs}=rac{ln(E_c/E_f)}{ln r}.
# ]
# Three or more levels are preferable because a two-level slope can be pre-asymptotic.
#
# ## 4.10 CFDX traceability and status
#
# Gradient: [gradient.h](../../../src/cfdx/core/numerics/gradient.h)  
# Stencil: [gradient_stencil.h](../../../src/cfdx/core/numerics/gradient_stencil.h)  
# Least squares: [least_squares_gradient.h](../../../src/cfdx/core/fvm/least_squares_gradient.h)  
# Interpolation: [interpolation.h](../../../src/cfdx/core/numerics/interpolation.h)  
# Registry: [numerical_method_registry.h](../../../src/cfdx/core/numerics/numerical_method_registry.h)
#
# N2 remains PARTIAL / not qualified. In particular, a successful constant/linear test on an affine mesh does not qualify second-order polyhedral accuracy. The qualification must cover skewed/polyhedral meshes, rank/conditioning diagnostics, boundary stencils, face reconstruction and controlled refinement.
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
