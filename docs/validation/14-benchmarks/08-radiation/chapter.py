# %% [markdown]
# # Radiation benchmark family
#
# S2S, P1 and DOM/RTE are separate populations because they introduce different physical approximations and numerical errors.
#
# ## Fundamental oracles
#
# Blackbody emission:
# [
# E_b=sigma T^4.
# ]
#
# For a closed enclosure:
# [
# sum_jF_{ij}=1,qquad A_iF_{ij}=A_jF_{ji},qquad 0le F_{ij}le1.
# ]
#
# Closure and reciprocity are geometry/model checks independent of the final heat-transfer result.
#
# ## S2S
#
# For radiosities,
# [
# G_i=sum_jF_{ij}J_j,
# ]
# and the linear system is
# [
# [I-operatorname{diag}(1-arepsilon)F]J=arepsilonsigma T^4.
# ]
# Verify matrix residual, view-factor closure, reciprocity and net energy exchange separately.
#
# ## P1/DOM
#
# P1 campaigns should identify absorption, scattering and optical thickness. DOM campaigns must declare angular quadrature and test angular refinement in addition to spatial refinement.
#
# [
# 	au=eta L,qquad eta=kappa+sigma_s.
# ]
#
# Thin, intermediate and optically thick regimes should not be treated as one validation population.
#
# ## Coupled thermal-radiation
#
# Report radiation-model residual, thermal residual, interface/global energy balance and radiation iteration error independently. A converged coupled calculation requires all relevant balances to close.
