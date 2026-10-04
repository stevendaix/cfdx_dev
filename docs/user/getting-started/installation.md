# Installation and environment

CFDX is built from the repository source tree. The exact supported dependency set is defined by the repository CMake configuration and CI; this page describes the reproducible workflow rather than pinning undocumented local packages.

## 1. Configure

Create an out-of-source build directory and configure CMake with the options required by the target platform. Record the compiler, CMake version, enabled MPI/GPU/HDF5 features, and the Git revision.

## 2. Build

Build the configured targets before running any validation campaign. A dependency or compiler failure is an environment/build failure, not a numerical result.

## 3. Test

Run the repository's CTest suite. Start with fast tests, then run the relevant numerical verification and benchmark campaigns.

## 4. Reproducibility

For a reportable calculation retain:

- Git revision;
- build configuration;
- compiler and dependency versions;
- operating system and hardware;
- case/setup files;
- mesh and geometry;
- numerical schemes;
- solver and stopping criteria;
- raw logs and machine-readable results.

## Supported-status rule

A documented capability is not necessarily a qualified production capability. Check the V&V and validation evidence before relying on a model, solver option, or benchmark for a new engineering calculation.
