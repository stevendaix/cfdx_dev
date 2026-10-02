# %% [markdown]
"""# Turbulence

Reynolds decomposition creates the Reynolds-stress closure problem. RANS, LES and hybrid methods make different modelling assumptions and require model-specific verification and validation.

The chapter must derive each supported closure, define wall treatment and length scales, map the equations to CFDX source files, and link every quantitative claim to executable evidence.
"""

# %%
MODEL_FAMILIES = ('SA', 'SST', 'LES', 'DES')
assert len(MODEL_FAMILIES) == 4
