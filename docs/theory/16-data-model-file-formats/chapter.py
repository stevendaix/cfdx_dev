# %% [markdown]
# Data Model and File Formats
#
# ## 16.1 Semantic separation
# CFDX separates authoritative setup, runtime state and visualization:
# \[
# \text{case.cfdx.h5}=\text{setup},\quad
# \text{case.dat.h5}=\text{checkpoint},\quad
# \text{case}_{t}.vtu=\text{visualization}.
# \]
#
# ## 16.2 Case model
# \[
# \mathcal C=\{\mathcal M,\mathcal B,\mathcal P,\mathcal N,\mathcal F,\mathcal V\},
# \]
# representing mesh, boundaries, physical models, numerical configuration, fields and metadata/provenance.
#
# ## 16.3 HDF5 contract
# A dataset is conceptually
# \[
# D=(shape,dtype,data,attributes).
# \]
# Schema must define units, component ordering, topology, field location, boundary identifiers, version and compatibility.
#
# ## 16.4 Restart
# A checkpoint must contain every field/history quantity needed by the selected temporal and nonlinear algorithms. Restart correctness is a numerical property, not merely successful file parsing.
#
# ## 16.5 Round trip
# For supported data:
# \[
# R(W(x))\equiv x
# \]
# within explicitly documented representation limits.
#
# ## 16.6 Provenance
# A result should identify case revision, mesh, physics, numerics, solver configuration and software revision.
#
# ## 16.7 CFDX paths
# src/cfdx/io/hdf5/  
# src/cfdx/io/restart/  
# src/cfdx/io/vtu/  
# Tests: tests/unit/test_case_hdf5_io.cpp, tests/unit/test_hdf5_roundtrip.cpp, tests/validation/test_solver_vtu_provenance.cpp
#
# %%
from __future__ import annotations
import numpy as np
x=np.arange(12,dtype=float).reshape(3,4)
assert np.array_equal(x,x.copy())
