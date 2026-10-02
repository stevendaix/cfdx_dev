# %% [markdown]
"""# Heat Transfer

A generic thermal equation is rho*cp*(dT/dt + u.grad(T)) = div(k*grad(T)) + S_T. Fourier conduction gives q = -k*grad(T). The finite-volume heat flux is -k_f*grad(T)_f.S_f. Boundary conditions must preserve the selected outward-normal sign convention.

Verification starts with one-dimensional conduction and manufactured sources before coupled flow benchmarks.
"""

# %%
k = 10.0
area = 2.0
length = 0.5
delta_t = 20.0
assert abs(k * area * delta_t / length - 800.0) < 1e-14
