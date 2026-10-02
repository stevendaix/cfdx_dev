# %% [markdown]
# Computational Chain and Code Map
#
# ## 17.1 End-to-end chain
# \[
# case\rightarrow schema\rightarrow mesh\rightarrow geometry
# \rightarrow fields/BC\rightarrow numerics\rightarrow physics
# \rightarrow assembly\rightarrow solve\rightarrow diagnostics
# \rightarrow output/restart.
# \]
# Every arrow is a contract boundary.
#
# ## 17.2 Case ingestion
# Validate schema, units, references and required data before solver construction.
#
# ## 17.3 Mesh and geometry
# Build topology and geometry and verify
# \[
# \sum_f\mathbf S_f=0
# \]
# for closed cells and positive volumes under the declared convention.
#
# ## 17.4 Fields and boundaries
# Assign dimensions, locations and mathematical boundary conditions. Boundary conditions enter the algebraic operator explicitly.
#
# ## 17.5 Numerical selection
# A numerical-method identifier must resolve to a concrete implementation and contract. Silent substitution by a fallback method breaks provenance.
#
# Registry paths: src/cfdx/core/numerics/numerical_method_registry.h and src/cfdx/core/numerics/numerical_method_contract.h
#
# ## 17.6 Assembly
# Discretisation creates
# \[
# Ax=b
# \]
# or
# \[
# F(x)=0.
# \]
#
# ## 17.7 Solve and diagnostics
# Linear/nonlinear convergence, conservation, boundedness and provenance are separate diagnostics.
#
# ## 17.8 Output
# Checkpoints serve restart; VTU serves visualization. They must retain their distinct semantic roles.
#
# Application paths: src/cfdx/application/simulation_controller.cpp, apps/production_solver/main.cpp, src/cfdx/physics/finite_volume_transport.h
#
# %%
from __future__ import annotations
stages=["case","schema","mesh","geometry","fields/BC","numerics","physics","assembly","solve","diagnostics","output"]
assert stages[0]=="case" and stages[-1]=="output"
