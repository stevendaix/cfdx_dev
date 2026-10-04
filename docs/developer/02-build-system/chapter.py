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
# ## CI principle
#
# Diagnostics should be available before the final gate. Failure output should identify the violated contract rather than only returning a generic non-zero status.
