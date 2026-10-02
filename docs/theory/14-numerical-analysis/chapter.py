# %% [markdown]
# Numerical Analysis
#
# ## 14.1 Consistency
# \[
# L_hu=Lu+\tau_h,\qquad \tau_h\to0.
# \]
# Truncation error is a local operator property.
#
# ## 14.2 Accuracy and observed order
# \[
# E_h\approx Ch^p,\qquad
# p_{obs}=\frac{\ln(E_h/E_{h/r})}{\ln r}.
# \]
# A controlled mesh sequence and negligible iterative error are required.
#
# ## 14.3 Stability
# For \(U_t=AU\), the discrete method has an amplification operator G. Stability requires bounded powers over the declared regime.
#
# ## 14.4 Boundedness and monotonicity
# A bounded reconstruction may require
# \[
# \min_N\phi_N\le\phi_f\le\max_N\phi_N.
# \]
# Monotonicity concerns the full discrete operator.
#
# ## 14.5 Conditioning
# \[
# \kappa(A)=\|A\|\|A^{-1}\|.
# \]
# It measures perturbation sensitivity and is not itself a convergence prediction for non-normal systems.
#
# ## 14.6 Error budget
# \[
# E_{total}\lesssim E_{model}+E_{space}+E_{time}+E_{iteration}+E_{roundoff}+E_{input}.
# \]
#
# ## 14.7 CFDX traceability
# src/cfdx/core/numerics/  
# src/cfdx/core/linalg/  
# tests/validation/test_mesh_refinement_order.cpp  
# tests/validation/test_mms_scalar_diffusion.cpp  
# tests/unit/test_conservation_boundedness.cpp
#
# %%
from __future__ import annotations
import numpy as np
E=np.array([1e-2,2.5e-3,6.25e-4])
assert np.allclose(np.log(E[:-1]/E[1:])/np.log(2),2.0)
