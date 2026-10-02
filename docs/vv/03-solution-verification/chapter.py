# %% [markdown]
# # 03 — Solution Verification
#
# Solution verification separates discretisation error from iterative error.
#
# \[
# e_{total}\approx e_{space}+e_{time}+e_{iter}+e_{other}.
# \]
#
# The decomposition is conceptual; the terms should be isolated experimentally where possible.
#
# ## Iterative error
#
# For a linear system,
#
# \[
# r=b-Ax=-A(x-x^*).
# \]
#
# Residual magnitude alone is therefore not a universal error bound. Conditioning and the norm must be considered.
#
# ## Spatial verification
#
# Refine the mesh while controlling other errors, calculate observed order and demonstrate an asymptotic trend before claiming the formal order.
#
# ## Temporal verification
#
# Refine delta-t while keeping spatial error below the temporal signal. Nonlinear and linear solve errors must remain below that signal.
#
# ## Reporting
#
# Report mesh/time-step sequence, solver tolerances, measured quantity, errors, observed order and pre-asymptotic behaviour.
