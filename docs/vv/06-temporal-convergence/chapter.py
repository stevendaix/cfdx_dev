# %% [markdown]
# # 06 — Temporal Convergence
#
# For a method of formal order p,
#
# \[
# E_{\Delta t}=C\Delta t^p.
# \]
#
# ## Campaign
#
# Keep the spatial error below the temporal signal, refine \(\Delta t\) by a known ratio, compare at the same physical time and keep nonlinear/linear errors controlled.
#
# ## Schemes
#
# Verify explicit Euler, implicit Euler, Crank–Nicolson, BDF2 and any Runge–Kutta implementation against their expected stability and order properties.
#
# Variable-step methods require the actual time-step ratio in the analysis; constant-step formulae must not be reused blindly.
#
# ## Restart
#
# A restart is part of verification. Multistep schemes must restore the required history and reproduce the expected trajectory within the declared numerical criterion.
