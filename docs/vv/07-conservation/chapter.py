# %% [markdown]
# # 07 — Conservation
#
# Finite-volume conservation must be checked locally and globally.
#
# ## Internal faces
#
# \[
# F_{P,f}+F_{N,f}=0.
# \]
#
# ## Global balance
#
# \[
# R_C=
# \frac{d}{dt}\sum_P q_PV_P+
# \sum_{f\in\partial\Omega}F_f-
# \sum_Ps_PV_P.
# \]
#
# Normalise only with a declared reference scale:
#
# \[
# \epsilon_C=|R_C|/Q_{ref}.
# \]
#
# ## Campaign
#
# Check internal-face antisymmetry, closed-domain cancellation, boundary flux accounting, source accounting, transient storage and local residuals. Repeat on serial and partitioned calculations where parallelism is supported.
#
# Conservation is necessary but does not prove accuracy, stability or boundedness.
