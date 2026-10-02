# 12. Comparison of Methods

| Method | Main information | Main sensitivity |
|---|---|---|
| Green–Gauss | face geometry and face values | face reconstruction and mesh geometry |
| Least Squares | cell-centre geometry | stencil conditioning |
| Weighted LS | cell-centre geometry plus weights | conditioning and weight policy |
| Vertex-based | vertex connectivity and reconstruction | vertex treatment and topology |

This is a descriptive comparison, not a ranking. Appropriate selection depends on the discretisation, mesh quality, boundary treatment and solver requirements.

A controlled comparison should use the same fields, mesh families, error definitions and convergence criteria wherever methods are intended to be compared.
