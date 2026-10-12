# N14 — Benchmark cross-reference

N14 consumes evidence from other numerical-method packages; it does not replace their acceptance criteria.

| N14 evidence type | Existing CFDX source | Role |
|---|---|---|
| Gradient/reconstruction | `tests/validation/test_gradient_verification.cpp`, `tests/validation/test_polyhedral_gradient_campaign.cpp` | formulation and observed-order evidence for N2 |
| Convection | `tests/validation/test_convection_scheme_verification.cpp`, `test_convection_3d_verification.cpp`, `test_convection_polyhedral_campaign.cpp` | scheme-level evidence for N4 |
| Diffusion | `tests/validation/test_nonorthogonal_laplacian_campaign.cpp`, `test_nonorthogonal_skew_campaign.cpp` | N3 formulation/robustness evidence |
| Temporal | `tests/validation/test_temporal_order_matrix.cpp` | N5 order evidence |
| Linear/AMG/Schur | `tests/unit/test_simplerc_schur.cpp`, `docs/development/HYPRE_PETSC_CFDX_AUDIT.md` | N8 algebraic comparison |
| Pressure-velocity | `tests/validation/test_phase9_acceptance.cpp` | N9 algorithm evidence |
| Cross-backend | `tests/unit/test_mpi_deterministic.cpp`, `test_matrix_free_fv_operator.cpp`, restart tests | N13 equivalence evidence |
| External physical reference | `docs/validation/FLUENT_VMFL_MATRIX.md` | physical/reference cases, not method equivalence |

## Benchmark rule

A reference-code result can be recorded only when:
1. governing equations and physical data are comparable;
2. geometry and boundary conditions are comparable;
3. numerical method differences are documented;
4. the CFDX result has independent convergence/conservation evidence;
5. the comparison quantity and tolerance are explicitly stated.

A close result is evidence for that benchmark point only. It is not a global parity claim.
