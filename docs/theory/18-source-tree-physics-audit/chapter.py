# %% [markdown]
# Source-Tree Physics Audit
#
# ## 18.1 Traceability chain
# \[
# equation\rightarrow discrete\ operator\rightarrow implementation
# \rightarrow unit\ test\rightarrow verification\rightarrow validation
# \rightarrow qualification.
# \]
# A source file proves existence only. A test proves only the property it measures.
#
# ## 18.2 Audit questions
# What continuous equation is implemented? What assumptions are made? What discrete equation is assembled? What geometry/interpolation enters? Which file owns the implementation? Which test exercises it? Which independent reference supports it? What is the acceptance criterion? What limitations remain?
#
# ## 18.3 Source-tree map
# Mesh: src/cfdx/core/mesh/  
# Geometry: src/cfdx/core/geometry/  
# Fields: src/cfdx/core/field/  
# Boundary: src/cfdx/core/boundary/  
# Numerics: src/cfdx/core/numerics/  
# Linear algebra: src/cfdx/core/linalg/  
# Physics: src/cfdx/physics/  
# HDF5: src/cfdx/io/hdf5/  
# Restart: src/cfdx/io/restart/  
# VTU: src/cfdx/io/vtu/  
# Unit tests: tests/unit/  
# Numerical tests: tests/numerical/  
# Validation: tests/validation/
#
# ## 18.4 Status rules
# Implemented means code exists. Verified means an executable test establishes a defined property. Validated means independent physical/reference evidence supports the declared case. Qualified means the declared population and acceptance package are complete.
#
# ## 18.5 N2 traceability
# \[
# \nabla\phi\rightarrow\{Green\!-\!Gauss,LS,WLS,vertex\}
# \rightarrow face\ reconstruction\rightarrow diffusion/viscous\ terms
# \rightarrow V\&V.
# \]
# Primary paths:
# src/cfdx/core/numerics/gradient.h  
# src/cfdx/core/numerics/gradient_stencil.h  
# src/cfdx/core/fvm/least_squares_gradient.h  
# src/cfdx/core/numerics/interpolation.h
#
# Campaigns:
# tests/validation/test_gradient_verification.cpp  
# tests/validation/test_polyhedral_gradient_campaign.cpp
#
# N2 remains PARTIAL / not qualified until its declared acceptance matrix is demonstrated.
#
# %%
from __future__ import annotations
statuses=("Implemented","Verified","Validated","Qualified")
assert len(set(statuses))==4
