# %% [markdown]
# # 07 — Conservation
#
# A finite-volume conservation claim must be demonstrated at operator and global levels.
#
# For an internal face:
#
# \[
# F_{P,f}+F_{N,f}=0.
# \]
#
# For a global balance:
#
# \[
# R_C=
# \frac{d}{dt}\sum_Pq_PV_P+
# \sum_{f\in\partial\Omega}F_f-
# \sum_Ps_PV_P.
# \]
#
# Report signed and normalised defects:
#
# \[
# \epsilon_C=|R_C|/Q_{ref}.
# \]
#
# Required checks include face-flux antisymmetry, closed-domain cancellation, boundary flux accounting, source accounting, transient storage and local cell residuals where applicable.
#
# Conservation is necessary but does not prove accuracy or boundedness.
