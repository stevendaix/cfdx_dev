# %% [markdown]
# # 03 — Solution Verification
#
# Solution verification estimates whether the numerical solution is sufficiently close to the mathematical solution for the declared calculation.
#
# ## Error budget
#
# \[
# e_{total}\approx e_{space}+e_{time}+e_{iter}+e_{model}+e_{other}.
# \]
#
# Model-form error belongs to validation/model assessment and should not be hidden inside discretisation error.
#
# ## Iterative error
#
# For \(Ax=b\),
#
# \[
# r=b-Ax=-A(x-x^*).
# \]
#
# A small residual is not a universal error bound; conditioning and the chosen norm matter.
#
# ## Spatial and temporal isolation
#
# Spatial studies hold temporal and iterative errors below the spatial signal. Temporal studies do the converse. Nonlinear convergence must be tighter than the discretisation signal being measured.
#
# ## Reporting
#
# Give mesh/time-step sequence, norms, stopping criteria, observed order, asymptotic evidence and sensitivity to solver tolerances.
