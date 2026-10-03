---
name: cfdx-cmake-build
description: Configure, build and diagnose CFDX CMake targets without masking dependency or compiler failures.
version: 1
---

# Procedure

1. Inspect CMake options, presets/toolchain files and target dependencies.
2. Configure with the repository-supported compiler/toolchain.
3. Build the smallest affected target first.
4. Inspect the first meaningful diagnostic rather than downstream cascades.
5. Run the relevant tests.
6. Re-run the required CI configuration when the change affects build infrastructure.

# Rules

Do not disable targets, tests or validation gates merely to obtain a successful build.
