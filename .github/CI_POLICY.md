# CFDX CI policy

The workflow at `.github/workflows/cfdx-build.yml` is a required repository
invariant. It must remain enabled and must build and run the complete CTest
suite on pull requests targeting `master`.

Required triggers:
- `pull_request` targeting `master`
- `push` to `master`
- `workflow_dispatch` for manual recovery/verification

Required checks:
- CMake/Ninja Release configuration
- MPI enabled
- HDF5 enabled
- Python mesh I/O dependencies installed in the repository virtualenv
- Full `ctest --test-dir build --output-on-failure`

Repository protection should require this workflow's check before merging to
`master`, and CODEOWNERS review should be required for changes under
`.github/workflows/` and `.github/CODEOWNERS`.

Do not remove or weaken the workflow to make a failing test pass. Fix the
underlying implementation or validation instead.
