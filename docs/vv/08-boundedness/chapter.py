# %% [markdown]
# # 08 — Boundedness
#
# Boundedness prevents numerical overshoots and undershoots from violating declared admissible ranges.
#
# ## Discrete requirement
#
# For a scalar with bounds \(\phi_{min}\le\phi\le\phi_{max}\),
#
# \[
# \phi_{min}\le\phi_P\le\phi_{max}.
# \]
#
# Positivity is the special case \(\phi_{min}=0\).
#
# ## Verification
#
# Test constant fields, monotone data, discontinuities, skewed meshes and limiting activation. Report minimum, maximum, number of violations and their magnitude.
#
# ## Important distinction
#
# Boundedness is not the same as accuracy. Excessive diffusion can be bounded but inaccurate. A limiter must therefore be verified for both its admissibility property and its expected formal behaviour in smooth regions.
#
# ## Failure analysis
#
# Identify whether violations originate in interpolation, flux linearisation, source treatment, boundary conditions or solver convergence.
