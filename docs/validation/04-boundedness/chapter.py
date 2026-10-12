# %% [markdown]
# # Boundedness and positivity diagnostics
#
# Executable N11 contract for transported scalar bounds and positive fields.
# Bounds are diagnostics: values are never clipped to make a test pass.

# %%
from math import isfinite

def check_bounds(values, lower, upper, tolerance=0.0):
    finite = all(isfinite(v) for v in values)
    lo = min(values)
    hi = max(values)
    lower_violation = max(0.0, lower - lo)
    upper_violation = max(0.0, hi - upper)
    return finite and lower_violation <= tolerance and upper_violation <= tolerance

# %%
assert check_bounds([0.0, 0.25, 1.0], 0.0, 1.0)
assert not check_bounds([-1.0e-8, 0.25, 1.0], 0.0, 1.0)
assert not check_bounds([0.0, 0.25, 1.0 + 1.0e-8], 0.0, 1.0)
assert not check_bounds([0.0, float("nan"), 1.0], 0.0, 1.0)

# %% [markdown]
# A positive turbulence-like variable uses an explicit admissibility floor.
# A zero floor is only the mathematical positivity specialization; model-specific
# strictly-positive requirements must declare their own floor.

# %%
assert check_bounds([0.0, 1.0e-12, 2.0], 0.0, float("inf"))
assert not check_bounds([-1.0e-12, 1.0, 2.0], 0.0, float("inf"))
