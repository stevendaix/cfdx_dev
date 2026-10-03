# %% [markdown]
# # VMFL036 axisymmetric sphere
#
# VMFL036 must be treated as a frozen external-validation definition, not as a generic sphere example.
#
# ## Physical definition
#
# The canonical case is 2-D axisymmetric flow around a sphere of diameter (D) in an external circular domain. The exact CFDX qualification record must freeze domain extent, inlet velocity, density, viscosity, outlet/far-field condition and axis treatment.
#
# A documented reference inconsistency must remain explicit: the Fluent properties printed for the VMFL036 reproduction imply (Re=50) for (ho=1), (U=1), (D=1), (mu=0.02), while the commonly cited literature datum (C_Dapprox1.0895) is associated with (Re=100). These must be separate populations rather than silently mixed.
#
# ## Governing quantity
#
# The drag is obtained from pressure and viscous traction:
#
# [
# mathbf F_D=int_S(-pmathbf n+oldsymbol	aumathbf n)_x,dA,
# qquad
# C_D=rac{F_D}{	frac12ho U_infty^2(pi D^2/4)}.
# ]
#
# Retain pressure drag and viscous drag separately before forming total (C_D).
#
# ## Required numerical evidence
#
# Before physical qualification:
#
# 1. axisymmetric geometry/metric verification;
# 2. axis regularity;
# 3. axisymmetric diffusion/operator verification;
# 4. mass-flux conservation;
# 5. pressure-velocity verification;
# 6. wall traction/force integration verification;
# 7. coarse/medium/fine mesh campaign;
# 8. independent comparison with the frozen reference.
#
# Mesh records must include near-wall spacing, skewness, non-orthogonality and refinement ratio. A single mesh agreement is diagnostic, not a convergence study.
#
# ## Status boundary
#
# No generic 3-D sphere, injected drag value, or reference-code output substituted for CFDX execution may be used as VMFL036 evidence.
