# %% [markdown]
"""# Numerical Analysis

## Scientific objective

Consistency means L_h(u)-L(u) tends to zero as h tends to zero. If E_h is approximately C*h^p, observed order is log(E_h/E_hr)/log(r). Accuracy, stability, conservation, boundedness and conditioning are distinct. Error interpretation must separate discretisation, iterative, roundoff, model and input effects.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

import numpy as np; e=np.array([4e-2,1e-2,2.5e-3]); p=np.log(e[:-1]/e[1:])/np.log(2.); assert np.allclose(p,2.)
