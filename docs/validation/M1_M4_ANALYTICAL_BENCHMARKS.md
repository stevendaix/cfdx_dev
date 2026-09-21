# CFDX M1-M4 analytical verification

This directory is the Level-A verification layer for M1-M4. These tests use closed-form solutions independently of the solver implementation wherever possible.

## Cases

| Module | Case | Reference | Quantities |
|---|---|---|---|
| M1 | Couette | (u=Uy/H) | L2/Linf velocity |
| M1 | Plane Poiseuille | (u=G y(H-y)/(2\mu)) | L2/Linf velocity, observed order |
| M3 | 1-D conduction | (T=T_0+(T_1-T_0)y/H) | L2/Linf temperature |
| M3 | CHT resistance series | (q=\Delta T/\sum L_i/k_i) | heat flux, interface temperature |
| M4 | black/gray two-surface exchange | Stefan-Boltzmann + radiosity resistance | net exchange |

## Convergence policy

Poiseuille is the primary spatial-order gate because the exact solution is quadratic and the source term is constant.

Meshes:
- 8
- 16
- 32
- 64 cells across the channel

Observed order:

[
p = \frac{\log(E_h/E_{h/2})}{\log 2}.
]

The acceptance threshold is (p \ge 1.90) for every refinement step, with an expected second-order scheme.

Couette and 1-D conduction are linear exact solutions for the orthogonal diffusion operator, so their acceptance criterion is an absolute error threshold rather than an artificial convergence-order requirement.

## Running

Configure and build normally, then:

```bash
python scripts/run_validation.py --build-dir build
```

or directly through CTest:

```bash
ctest --test-dir build -R test_analytical_benchmarks --output-on-failure
```

## Interpretation

A PASS here means the current discretization reproduces the Level-A analytical benchmark within the explicit numerical tolerances. It does **not** by itself establish correctness of the complete coupled M1-M4 stack.

Level-B work should add independent numerical references (e.g. lid-driven cavity, turbulent channel and flat plate). Level-C work should then exercise coupled turbulence/energy/radiation/CHT systems.
