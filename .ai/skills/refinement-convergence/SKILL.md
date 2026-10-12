---
name: cfdx-refinement-convergence
description: Design mesh and time refinement studies with explicit convergence evidence.
version: 1
---

# Procedure

1. Define the QoI/error norm.
2. Define the refinement parameter and ratio.
3. Keep physical/model parameters controlled.
4. Ensure iterative solver error is sufficiently smaller than discretisation error.
5. Run multiple refinement levels.
6. Check monotonicity and asymptotic behavior.
7. Compute observed order only when the refinement data support it.
8. Report raw values as well as derived slopes.

# Required checks

- solver convergence;
- discretisation error;
- refinement consistency;
- observed order;
- sensitivity to boundary/corner treatment where relevant.

A single fine-grid result is not an accuracy study.
