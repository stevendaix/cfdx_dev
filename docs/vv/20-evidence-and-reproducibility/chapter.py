# %% [markdown]
# # 20 — Evidence and Reproducibility
#
# A V&V result should be replayable by another developer from retained information.
#
# ## Minimum manifest
#
# \[
# E=(revision,setup,mesh,physics,numerics,solver,environment,raw,postprocess,criterion).
# \]
#
# ## Artifact classes
#
# Retain machine-readable raw results, derived metrics, logs, mesh identifiers, configuration and the exact post-processing definition. Human-readable reports are secondary to raw evidence.
#
# ## Reproducibility
#
# Record compiler/build configuration, relevant dependency versions, hardware and parallel configuration when they can affect results.
#
# ## Replay
#
# A replay should regenerate the metric from the retained raw result without manual intervention beyond documented setup.
#
# ## Integrity
#
# Never overwrite a failed campaign with a later pass. Preserve the historical result and link the new campaign as a new evidence record.
