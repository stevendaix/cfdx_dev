# %% [markdown]
"""# Verification and Validation

## Scientific objective

Code verification checks implementation against intended mathematics. Solution verification controls numerical error for a stated question. Validation compares with independent physical or trusted reference evidence. Qualification demonstrates a declared capability over a defined scope. MMS derives forcing from a chosen exact field and tests the discrete solver independently.

## Status

Implementation, verification, validation and qualification remain separate claims.
"""

# %%
from __future__ import annotations

import numpy as np; u=np.array([1.,2.]); assert np.linalg.norm(u-u)==0.
