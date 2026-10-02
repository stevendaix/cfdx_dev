# %% [markdown]
# Radiation
#
# ## 12.1 Emission
# \[
# E_b=\sigma T^4,\qquad E=\varepsilon\sigma T^4.
# \]
#
# ## 12.2 Radiosity
# \[
# J=E+(1-\varepsilon)G,\qquad q''=J-G.
# \]
#
# ## 12.3 View factors
# \[
# A_iF_{ij}=A_jF_{ji},\qquad \sum_jF_{ij}=1
# \]
# for a closed enclosure. These identities are independent verification targets.
#
# ## 12.4 S2S
# Surface-to-surface radiation forms a coupled radiosity system. Geometry, areas and view factors must satisfy reciprocity and closure.
#
# ## 12.5 P1
# Under a common convention:
# \[
# -\nabla\cdot\left(\frac{1}{3\beta}\nabla G\right)+\beta G=4\beta\sigma T^4.
# \]
# The exact source/scattering terms depend on the optical model.
#
# ## 12.6 DOM
# Discrete ordinates replace angular transport by a finite set of directions and quadrature weights. Angular moments are reconstructed from these directional solutions.
#
# ## 12.7 CFDX traceability
# src/cfdx/physics/radiation.h  
# src/cfdx/physics/radiation_s2s.h  
# src/cfdx/physics/radiation_solver.h  
# src/cfdx/physics/radiation_models.h  
# Tests: tests/validation/test_radiation_vv.cpp, tests/validation/test_s2s_radiation_vv.cpp, tests/validation/test_phase12_radiation_vv.cpp
#
# %%
from __future__ import annotations
import numpy as np
sigma=5.670374419e-8
assert np.isclose(sigma*300.0**4,459.300326939)
