# %% [markdown]
# 18 — Source-Tree Physics Audit
#
# This chapter is the living register connecting equations to implementation and evidence. It deliberately does not convert source existence into qualification.
#
# ## 18.1 Audit relation
#
# $
# equation
# \rightarrow discrete\ equation
# \rightarrow implementation
# \rightarrow unit/numerical\ test
# \rightarrow verification
# \rightarrow validation
# \rightarrow qualification.
# $
#
# ## 18.2 What must be recorded
#
# For every important numerical/physical feature record:
# 1. continuous equation;
# 2. assumptions and units;
# 3. discrete equation;
# 4. geometric/interpolation inputs;
# 5. exact implementation path;
# 6. test path;
# 7. independent oracle/benchmark;
# 8. acceptance criterion;
# 9. current maturity;
# 10. limitations and next improvement.
#
# ## 18.3 Geometry family
#
# $
# \{\mathbf x_i\}
# \rightarrow
# \{\mathbf S_f,\mathbf C_f,\mathbf C_P,V_P\}
# \rightarrow
# discrete\ operators.
# $
# Source family: [core/geometry](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/geometry/) and [core/mesh](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/mesh/).
#
# ## 18.4 Gradient/reconstruction family
#
# $
# \nabla\phi
# \in\{Green\!-\!Gauss,LS,WLS,vertex\}
# $
# feeds face reconstruction and non-orthogonal diffusion. Primary sources: [gradient.h](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/numerics/gradient.h), [gradient_stencil.h](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/numerics/gradient_stencil.h), [least_squares_gradient.h](../../../src/cfdx/core/fvm/least_squares_gradient.h), [interpolation.h](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/numerics/interpolation.h).
#
# ## 18.5 Flux/FVM family
#
# $
# \partial_t(\rho\phi)+\nabla\cdot F=S
# \rightarrow
# \frac{d}{dt}(\rho\phi V)+\sum_fF_f=S_PV.
# $
# Sources include divergence, convection, laplacian, flux and finite-volume transport.
#
# ## 18.6 Pressure/linear algebra family
#
# $
# \begin{bmatrix}A_u&G\\D&C\end{bmatrix}
# \rightarrow
# Schur
# \rightarrow
# AMG/MGR/Krylov.
# $
# Source families: [core/linalg](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/core/linalg/) and [physics/pressure_velocity.h](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/physics/pressure_velocity.h).
#
# ## 18.7 Turbulence family
#
# RANS, LES and hybrid models introduce closure equations and/or constitutive relations. Wall distance is a model dependency for SA/SST/hybrid methods. Source family: [physics](https://github.com/stevendaix/cfdx_dev/tree/master/src/cfdx/physics/), especially turbulence and wall-distance modules.
#
# ## 18.8 Thermal/radiation family
#
# Thermal:
# $
# \rho c_pDT/Dt=\nabla\cdot(k\nabla T)+S_T.
# $
# Radiation:
# $
# E_b=\sigma T^4,\qquad
# J=E+(1-\varepsilon)G.
# $
# Coupled interfaces require both field and flux conservation.
#
# ## 18.9 Data/I/O family
#
# $
# setup\rightarrow HDF5\ case,\qquad
# state\rightarrow checkpoint,\qquad
# fields\rightarrow VTU.
# $
# Provenance is part of reproducibility.
#
# ## 18.10 Evidence family
#
# Unit tests establish local properties. Numerical tests establish controlled numerical properties. Validation tests establish benchmark evidence. Documentation must link the smallest relevant test rather than a generic “CI passed”.
#
# ## 18.11 Maturity vocabulary
#
# **Implemented:** code exists.
#
# **Verified:** a defined executable property is demonstrated.
#
# **Validated:** independent physical/reference comparison supports the declared case.
#
# **Qualified:** the declared population and acceptance package are complete.
#
# Documentation maturity never promotes numerical maturity.
#
# ## 18.12 N2 status
#
# N2 is currently **PARTIAL / not qualified** until the declared polyhedral gradient accuracy, conditioning/rank-deficiency, boundary reconstruction and face reconstruction gates are demonstrated on controlled campaigns.
#
# ## 18.13 Audit rule
#
# If an exact path, equation, benchmark or acceptance criterion has not been checked against the current repository, it must be marked unresolved rather than guessed.
#
# ## 18.14 Executable vocabulary check
# %%
statuses=("Implemented","Verified","Validated","Qualified")
assert len(statuses)==4 and len(set(statuses))==4
