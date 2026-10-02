# %% [markdown]
# # 14 — Turbulence Verification
#
# Turbulence verification begins with the closure equations and their admissibility constraints.
#
# ## Models
#
# Verify each implemented closure separately: k–epsilon variants, SST k–omega, Spalart–Allmaras and LES/DES-family models where implemented.
#
# ## Checks
#
# Verify positivity of turbulent kinetic energy and dissipation/frequency variables, wall-distance dependence, production/destruction signs, clipping policies and limiting.
#
# ## Wall distance
#
# Since wall distance enters several closures, verify it independently on analytical geometries before coupling it to turbulence.
#
# ## Qualification
#
# Model qualification is model- and regime-specific. A passing case does not qualify all Reynolds numbers, geometries or adverse-pressure-gradient regimes.
