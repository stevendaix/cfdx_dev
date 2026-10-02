# %% [markdown]
# # 17 — Multiphysics Verification
#
# A coupled problem can be represented as
#
# \[
# F(y)=0,
# \]
#
# with block Jacobian
#
# \[
# J=
# \begin{bmatrix}
# J_{11}&J_{12}\\
# J_{21}&J_{22}
# \end{bmatrix}.
# \]
#
# ## Coupling ladder
#
# Verify in this order:
#
# 1. each physics independently;
# 2. one-way coupling;
# 3. two-way coupling;
# 4. fully coupled nonlinear solve;
# 5. coupled MMS or manufactured block problem where practical.
#
# ## Interface checks
#
# Verify conservation of exchanged quantities, sign conventions, units, lagging/linearisation choices and convergence of the coupled residual.
#
# ## Error localisation
#
# A coupled failure should be reduced to the smallest subsystem reproducing it before changing global tolerances.
