# CFDX GitHub CI architecture

The repository deliberately separates four CI responsibilities.

## 1. CFDX CI

`.github/workflows/cfdx-ci.yml`

Mandatory software build/regression layer:
- Release and Debug/UBSan builds;
- fast CTest regression;
- Python regression;
- diagnostics and artifacts before the final gate.

It does not contain N-level qualification logic.

## 2. Numerical maturity

`.github/workflows/cfdx-numerical.yml`

One workflow exposes independent jobs for N1-N16 from issue #461. The jobs are intentionally independent: a failure in N7 does not skip N8.

Each N-level job:
1. configures the same Release numerical build;
2. runs only its declared numerical test set;
3. prints diagnostics before the final gate;
4. uploads evidence.

N16 has no executable acceptance gate because advanced methods are explicitly deferred in the current maturity contract; its job validates the machine-readable audit instead of pretending that deferred work is qualified.

The test mapping is kept in this workflow so GitHub remains the orchestration layer while CTest remains the source of executable test definitions.

For active N8/N9 work, `.github/workflows/cfdx-n8-n9.yml` provides label-driven dedicated gates:
- GitHub label `n8` runs the exact N8 numerical selection;
- GitHub label `n9` runs the exact N9 numerical selection;
- the labels persist across `synchronize`, so pushing a correction reruns the dedicated gate.

## 3. Total validation

`.github/workflows/cfdx-total-validation.yml`

Long/full qualification campaigns only:
- all non-long CTest tests;
- all long CTest tests;
- wall-distance campaigns;
- validation report;
- thermal/radiation regression;
- application E2E;
- N16 maturity-scope audit;
- retained evidence.

The non-long and long CTest partitions are complementary: together they execute the complete CTest registration, so N1-N15 executable tests are not silently omitted from total validation.

It is intentionally not a second N1 gate.

## 4. Maturity records

`.github/workflows/cfdx-maturity.yml`

Validates the machine-readable #461 maturity and requirements audits on every PR/push.

## Required-check policy

Required checks should point only to workflows that always create a deterministic check:
- CFDX CI;
- the N-level checks that the project chooses to make blocking;
- CFDX numerical maturity audit.

Do not use path-filtered workflows for required checks. A required workflow must exist on every relevant PR; otherwise GitHub can report a missing/skipped check even when no numerical test failed.

## Removed duplication

The old standalone N1 workflow, the overlapping MPI/HDF5 build workflow, and the path-filtered numerical maturity workflow are replaced by the four-layer structure above. The former `cfdx-validation.yml` trigger is retained only as a manual legacy workflow while the fresh `cfdx-total-validation.yml` owns the active label/schedule triggers. This separation also avoids reusing the anomalous historical workflow identity.
