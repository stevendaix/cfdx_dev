# %% [markdown]
"""# Radiation

Blackbody emission is E_b = sigma*T**4. Surface-to-surface radiation uses view factors with closure and reciprocity under the corresponding closed-enclosure assumptions. Radiation is nonlinear in temperature and must be coupled consistently with the thermal energy balance.

Verification of geometric identities and radiative exchange precedes coupled thermal-radiation validation.
"""

# %%
sigma = 5.670374419e-8
T = 300.0
assert sigma * T**4 > 0.0
