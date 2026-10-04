# %% [markdown]
# # 13 — Convection Verification
#
# Convective fluxes have both conservation and accuracy requirements.
#
# ## Generic flux
#
# \[
# F_f=(\rho\mathbf u\cdot\mathbf S)_f\,\phi_f.
# \]
#
# ## Schemes
#
# Verify upwind, central and higher-order/reconstructed schemes separately. The scheme identifier must correspond to the algorithm actually executed.
#
# ## Campaign
#
# Use constant-field preservation, linear-field tests, smooth MMS, transport benchmarks and boundedness tests. Include mesh distortion when reconstruction is sensitive to it.
#
# ## Stability versus accuracy
#
# A bounded scheme may be strongly dissipative. A high-order scheme may require limiting. Report both accuracy and admissibility rather than collapsing them into one score.
