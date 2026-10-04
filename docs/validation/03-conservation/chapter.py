# %% [markdown]
# # Conservation and boundedness diagnostics
#
# This chapter is the executable mathematical contract for the N11 diagnostic
# layer. It verifies invariant definitions independently of solver residuals.
#
# It does not claim solver-level qualification. Production qualification requires
# CFDX execution with authoritative face fluxes and retained evidence.

# %%
from math import isfinite

def internal_face_balance(owner_flux: float, neighbour_flux: float) -> float:
    return owner_flux + neighbour_flux

def bounded(values, lower: float, upper: float, tolerance: float = 0.0):
    finite = all(isfinite(v) for v in values)
    lower_violation = max(0.0, lower - min(values))
    upper_violation = max(0.0, max(values) - upper)
    return {
        "finite": finite,
        "min": min(values),
        "max": max(values),
        "lower_violation": lower_violation,
        "upper_violation": upper_violation,
        "bounded": finite and lower_violation <= tolerance and upper_violation <= tolerance,
    }

# %%
for phi in (0.0, 1.0, -3.5, 1.0e12):
    assert internal_face_balance(phi, -phi) == 0.0

# %%
report = bounded([0.0, 0.2, 1.0], 0.0, 1.0)
assert report["bounded"]
report = bounded([0.0, -1.0e-8, 1.0], 0.0, 1.0)
assert not report["bounded"]
assert report["lower_violation"] == 1.0e-8
report = bounded([0.0, float("nan"), 1.0], 0.0, 1.0)
assert not report["finite"]
assert not report["bounded"]

# %% [markdown]
# ## Qualification boundary
#
# The repository implementation evidence remains in
# test_conservation_boundedness, test_transport_conservation and
# test_conservation_assembly. These checks do not qualify production CFD cases.
