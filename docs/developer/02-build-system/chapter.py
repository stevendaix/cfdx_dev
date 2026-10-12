# %% [markdown]
# # 02 — Build System
#
# CMake is the build-system boundary. Configuration, compilation and verification are distinct stages.
#
# \[
# \text{configure}\rightarrow\text{build}\rightarrow\text{unit tests}\rightarrow\text{numerical tests}\rightarrow\text{validation campaigns}.
# \]
#
# Compilation proves that the selected target can be built. It does not prove numerical correctness.
#
# ## Configuration
#
# Record material choices affecting numerical results: compiler, build type, optional accelerators, linear-algebra backends and feature flags.
#
# ## Dependency failures
#
# A missing dependency is an infrastructure result. It must not be converted into a numerical pass by silently skipping the associated evidence.
#
# ## HDF5 without root
#
# HDF5 is optional in CMake, and when it is absent the configuration reports
# `HDF5 not found. HDF5-dependent features disabled.` and continues. That
# continuation hides more than it appears to: the `.cfdx.h5` interchange reader
# is excluded from the build, so `cfdx_convert` and `cfdx_production_solver`
# fail to *link* on an undefined `read_case_cfdx_h5`, and every qualification
# test that reads a `.h5` mesh aborts with `cannot read ... mesh`. The suite
# still reports a result, so an unverified change can look validated.
#
# Verify that the HDF5-dependent targets actually built before treating a local
# run as complete:
#
# ```bash
# test -x build/cfdx_convert && test -x build/cfdx_production_solver
# ```
#
# On a machine without root, the headers can be unpacked into a writable prefix
# and combined with the shared libraries the runtime package already installs.
# The `libhdf5-dev` package ships only static archives, so the versioned shared
# objects are taken from the installed `libhdf5-*-1` packages:
#
# ```bash
# apt-get download libhdf5-dev          # no root required
# dpkg-deb -x libhdf5-dev_*.deb $PREFIX/hdf5raw
# mkdir -p $PREFIX/hdf5/include $PREFIX/hdf5/lib
# cp -a $PREFIX/hdf5raw/usr/include/hdf5/serial/* $PREFIX/hdf5/include/
# cd $PREFIX/hdf5/lib
# ln -sf /usr/lib/x86_64-linux-gnu/libhdf5_serial.so.103     libhdf5.so
# ln -sf /usr/lib/x86_64-linux-gnu/libhdf5_serial_hl.so.100 libhdf5_hl.so
# cmake -S . -B build -DHDF5_ROOT=$PREFIX/hdf5 \
#       -DHDF5_INCLUDE_DIR=$PREFIX/hdf5/include -DCMAKE_PREFIX_PATH=$PREFIX/hdf5
# ```
#
# A successful configuration prints `Using system HDF5 <version>`.
#
# ## CI principle
#
# Diagnostics should be available before the final gate. Failure output should identify the violated contract rather than only returning a generic non-zero status.
