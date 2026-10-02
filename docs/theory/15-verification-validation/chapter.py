# %% [markdown]
# Verification and Validation
#
# ## 15.1 Definitions
# Verification asks whether the numerical implementation solves the specified mathematical model correctly. Validation asks whether the mathematical/physical model represents the target physical system adequately. Qualification is a declared engineering scope with explicit acceptance criteria and retained evidence.
#
# ## 15.2 Code verification
# Use unit tests, algebraic identities, dimensional checks, exact polynomial fields and manufactured solutions.
#
# ## 15.3 Solution verification
# Separate spatial, temporal and iterative errors. If
# \[
# E_h\approx Ch^p,
# \]
# then
# \[
# p_{obs}=\frac{\ln(E_h/E_{h/r})}{\ln r}.
# \]
# A GCI-style estimate is
# \[
# GCI_{12}=F_s\frac{|u_1-u_2|}{|u_1|}\frac1{r^p-1}.
# \]
# The safety factor, refinement ratio, p and assumptions must be retained.
#
# ## 15.4 MMS
# Select \(u_{exact}\) and derive
# \[
# S_{MMS}=L(u_{exact}).
# \]
# The numerical solution is then compared directly with the exact field.
#
# ## 15.5 Conservation verification
# \[
# R_C=\text{inflow}-\text{outflow}+\text{source}-dQ/dt.
# \]
# Report a normalized defect and its reference scale.
#
# ## 15.6 Validation
# Independent analytical or experimental references and uncertainty are required. Agreement with another CFD code alone is not physical validation.
#
# ## 15.7 Evidence
# Retain mesh, physical parameters, numerical settings, stopping criteria, software revision, raw diagnostics and post-processing definitions.
#
# ## 15.8 CFDX paths
# tests/unit/  
# tests/numerical/  
# tests/validation/  
# docs/validation/
#
# Qualification must be inferred only from retained evidence and declared acceptance gates.
#
# %%
from __future__ import annotations
import numpy as np
E=np.array([1e-2,2.5e-3,6.25e-4])
assert np.allclose(np.log(E[:-1]/E[1:])/np.log(2),2.0)
