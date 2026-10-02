# %% [markdown]
# # 06 — Temporal Convergence
#
# For a temporal method of order p,
#
# \[
# E_{\Delta t}=C\Delta t^p.
# \]
#
# ## Campaign
#
# Keep the spatial discretisation fixed and sufficiently fine. Reduce delta-t by a known ratio and measure a quantity at a consistent physical time.
#
# Nonlinear and linear solve errors must remain below the temporal discretisation signal.
#
# ## Restart
#
# A restarted run must preserve the mathematical state required by the time integrator. Multistep methods such as BDF2 also require the history needed by the method.
