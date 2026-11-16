# CFDX V68 - SST at identical state with two mass fluxes

Scope: an opt-in, strictly read-only forensic ablation of V65 iteration 23.
No new Schur, Jacobian, line search, turbulence physics or boundary closure.

## Exact comparison

Let X=(U,p,k,omega) be the *same* six-field candidate from the existing
COUPLED/V47/SST transaction. In the existing V51 candidate evaluator, compute

    R_carrier(X) = SST_V22(X, phi_carrier)
    R_physical(X) = SST_V22(X, phi_RhieChow(X))

using the same values of k, omega, strain, distance, F1/F2, wall omega,
production G, gamma_k/w, cross diffusion and boundary values. The physical
six-field operator remains R_physical; the counterfactual R_carrier is NEVER
used to accept or solve an iterate.

The reader-only `CFDX_FPHYS_V68_FLUX_SWAP=1` requires both
`CFDX_FPHYS_CANDIDATE_V51=1` and `CFDX_FPHYS_V65_SAME_OPERATOR=1`. It samples
V65 iteration 23 retries 0, 8 and 10. Other iterations and all default
paths are unchanged.

The independent face oracle evaluates, at each face and for each SST variable,

    delta R_i = [phi_P psi_up(P) - phi_C psi_up(C)]
                - [phi_P - phi_C] psi_i

with owner and neighbour orientations, the real boundary condition and both
upwind signs evaluated independently. Wall owner omega constraints are
recognized explicitly. Every cell's predicted delta must equal the **actual**
change in the two freshly assembled FV residuals to a strict relative 1e-9
identity. A second V16 signed-row oracle at cells 0 and 12 checks unchanged
diffusion and physical source budgets and parity with V22.

The reported `omega_effect_fraction=(||Romega_phys||inf-||Romega_carrier||inf)
/||Romega_phys||inf` measures ONLY the relative raw SST defect change under a
flux substitution. It is not a merit change, does not solve the carrier SST
again, and cannot certify convergence. A large/small effect alone does not
establish whether the historical flux is right: it identifies the impact
of the lagged flux at X.

## Execution

```bash
cmake -S . -B build-v68 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v68 --target cfdx_production_solver \
  test_bfs_sst_flux_swap_v68 -j2
ctest --test-dir build-v68 --output-on-failure \
  -R 'test_bfs_sst_flux_swap_v68|test_bfs_flux_swap_v68_evidence'
CFDX_FPHYS_CANDIDATE_V51=1 CFDX_FPHYS_V65_SAME_OPERATOR=1 \
CFDX_FPHYS_V68_FLUX_SWAP=1 CFDX_BFS_STAGE_TRACE=1 \
  build-v68/cfdx_production_solver --mesh bfs152.cfdx.h5 \
  --experimental-sst --iterations 23 --output-dir bfs-v68 \
  2> bfs-v68.log
python3 scripts/audit_bfs_v68_flux_swap.py --log bfs-v68.log \
  --baseline bfs-v67.log --output bfs-v68-audit.json
```

A failed 23rd nonlinear iteration is EXPECTED and must not be treated as
successful physics. The diagnostic passes only if all three probes are
reproduced, the independent FV face/budget identities close, and the previous
149 V65 decisions / 22 accepted states are unchanged. As before this is the
TH03_debug_only BFS152 pilot with historical SST walls, **not** the full
Driver-Seegmiller experiment or the BFS712/BFS20540 campaign.