# %% [markdown]
# 17 — Computational Chain: From PDE to CFDX Code
#
# ## 17.1 Complete chain
#
# $
# physical\ problem
# \rightarrow mathematical\ model
# \rightarrow PDE
# \rightarrow integral\ balance
# \rightarrow control\ volumes
# \rightarrow reconstruction
# \rightarrow face\ fluxes
# \rightarrow algebraic\ system
# \rightarrow nonlinear\ iterations
# \rightarrow linear\ solves
# \rightarrow diagnostics
# \rightarrow output/restart
# \rightarrow V\&V.
# $
#
# ## 17.2 Case ingestion
#
# The setup defines geometry, mesh, materials, boundary conditions, physics and numerical methods. Validation must occur before expensive solver construction:
# $
# schema\rightarrow units\rightarrow references\rightarrow required\ data.
# $
#
# ## 17.3 Mesh construction
#
# Build topology, then geometry:
# $
# topology\rightarrow \{\mathbf S_f,\mathbf C_P,V_P,\mathbf C_f\}.
# $
# Verify closure
# $
# \sum_f\mathbf S_f=0
# $
# and positive volumes before assembling any PDE.
#
# ## 17.4 Field and boundary construction
#
# Each field has location, dimensions and semantic meaning. Boundary conditions are mathematical operators, not merely labels:
# $
# Dirichlet,\quad Neumann,\quad Robin,\quad coupled/interface.
# $
#
# ## 17.5 Numerical selection
#
# A selected method resolves to an implementation contract:
# $
# method\ id\rightarrow scheme\ object\rightarrow discrete\ operator.
# $
# Silent fallback breaks reproducibility and provenance.
#
# ## 17.6 Spatial discretisation
#
# For a generic scalar:
# $
# \frac{d}{dt}(\rho\phi V_P)
# +\sum_f\dot m_f\phi_f
# =
# \sum_f\Gamma_f\nabla\phi_f\cdot\mathbf S_f
# +S_PV_P.
# $
# The gradient/reconstruction layer supplies $\nabla\phi$ and/or $\phi_f$; flux layers assemble conservative interface contributions.
#
# ## 17.7 Algebraic assembly
#
# Linearised equations become
# $
# A(U^k)\delta U=b(U^k)
# $
# or, for a direct fixed-point row,
# $
# AU=b.
# $
# Coupled pressure/velocity systems have block form
# $
# \begin{bmatrix}A_u&G\\D&C\end{bmatrix}
# \begin{bmatrix}u\\p\end{bmatrix}
# =
# \begin{bmatrix}b_u\\b_p\end{bmatrix}.
# $
#
# ## 17.8 Nonlinear solution
#
# Newton:
# $
# J(U^k)\delta U=-F(U^k).
# $
# Picard/fixed-point:
# $
# U^{k+1}=S(U^k).
# $
# Under-relaxation:
# $
# U^{k+1}\leftarrow(1-\omega)U^k+\omega S(U^k).
# $
#
# ## 17.9 Linear solution
#
# $
# Ax=b,\qquad r=b-Ax.
# $
# Krylov solvers operate on the algebraic system; preconditioners change the convergence properties without changing the intended target when applied consistently.
#
# ## 17.10 Diagnostics
#
# Track independently:
# $
# \|r_{linear}\|,\quad
# \|F_{nonlinear}\|,\quad
# \epsilon_{mass},\quad
# \epsilon_{momentum},\quad
# \epsilon_{energy},
# $
# boundedness, physical admissibility and step/time history. A single residual is never a complete convergence certificate.
#
# ## 17.11 Output and restart
#
# Checkpoint state must contain enough history to resume the selected algorithm. VTU is for visualisation; authoritative case and state remain HDF5.
#
# ## 17.12 Code map
#
# Case/application: [src/cfdx/application](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/application/). Mesh/geometry: [src/cfdx/core/mesh](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/mesh/) and [src/cfdx/core/geometry](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/geometry/). Numerics: [src/cfdx/core/numerics](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/numerics/). Linear algebra: [src/cfdx/core/linalg](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/linalg/). Physics: [src/cfdx/physics](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/physics/). I/O: [src/cfdx/io](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/io/).
#
# ## 17.13 V&V insertion points
#
# Each arrow has an evidence target:
# $
# mesh\rightarrow geometry\ verification,
# \quad reconstruction\rightarrow MMS,
# \quad flux\rightarrow conservation,
# \quad solver\rightarrow residual/gauge,
# \quad physics\rightarrow analytical/experimental\ benchmark.
# $
#
# ## 17.14 Executable chain invariant
# %%
stages=["case","model","PDE","integral balance","mesh","reconstruction","flux","assembly","nonlinear","linear","diagnostics","output","V&V"]
assert len(stages)==13
assert stages[0]=="case" and stages[-1]=="V&V"
