# CFDX CI policy

The workflow at `.github/workflows/cfdx-ci.yml` is a required repository
invariant. It must remain enabled and must build and run the complete CTest
suite on pull requests targeting `master`.

Required triggers:
- `pull_request` targeting `master`
- `push` to `master`
- `workflow_dispatch` for manual recovery/verification

Required checks:
- CMake/Ninja configuration in both `Release` and `DebugSanitizers`
- MPI enabled (`CFDX_ENABLE_MPI` defaults to `ON`)
- HDF5 enabled (detected; `libhdf5-dev` is installed by the workflow)
- Python mesh I/O dependencies installed in the repository virtualenv
- Full `ctest --test-dir build --output-on-failure`, long tests excluded in
  `Release` and long/benchmark/campaign excluded in `DebugSanitizers`
- The `tests/python` regression suite in `Release`

Repository protection should require this workflow's check before merging to
`master`, and CODEOWNERS review should be required for changes under
`.github/workflows/` and `.github/CODEOWNERS`.

Do not remove or weaken the workflow to make a failing test pass. Fix the
underlying implementation or validation instead.

## Enforcement status

Measured on 2026-10-04. The requirements above are not currently enforced:

- `master` has **no branch protection**
  (`GET /branches/master/protection` returns 404 "Branch not protected"), so no
  check is required before merging and no review is required.
- `continue-on-error: true` on the build and CTest steps is intentional: the
  `Final CI gate` step reads `steps.<id>.outcome`, which stays `failure` under
  `continue-on-error`, and fails the job. Removing the flag makes each step's
  own conclusion honest; it does not change the verdict.

CODEOWNERS assigns `@stevendaix` to `.github/workflows/`, `.github/CODEOWNERS`
and this file. Without branch protection that assignment is advisory.