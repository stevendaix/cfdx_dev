# Recent bibliography — dual-time / pseudo-time integration

| Year | Reference | Relevance to CFDX |
|---|---|---|
| 2023 | Nguyen & Park, *A Review of Preconditioning and Artificial Compressibility Dual-Time Navier–Stokes Solvers for Multiphase Flows*, Fluids 8(3), 100. | Review of dual-time, artificial compressibility and preconditioning. |
| 2024 | Zhang & Barakos, *Assessment of implicit adaptive mesh-free CFD modelling*, IJNMF 96, 670–700, DOI 10.1002/fld.5266. | Implicit dual-time with adaptive compressible RANS discretisation. |
| 2024 | *A general positivity-preserving algorithm for implicit high-order finite volume schemes solving the Euler and Navier–Stokes equations*, JCP 508, 112999, DOI 10.1016/j.jcp.2024.112999. | BDF2 dual-time plus positivity-preserving residual/flux corrections. |
| 2025 | Sayyari & Yamaleev, *Implicit dual time-stepping positivity-preserving entropy-stable schemes for the compressible Navier-Stokes equations*, arXiv:2504.11333. | BDF1/BDF2 + dual time + positivity/entropy stability. |
| 2025 | Zandbergen, Van Noorden & Heinlein, *Improving pseudo-time stepping convergence for CFD simulations with neural networks*, Computers & Mathematics with Applications 196, 64–83, DOI 10.1016/j.camwa.2025.07.006. | Recent local pseudo-time prediction; research direction rather than core dependency. |
| 2025 | Wang & Liu, *Development and application of fourth-order-accurate semi-implicit scheme for Navier-Stokes equations*, PCFD 25(3), 123–141, DOI 10.1504/PCFD.2025.146000. | Higher-order semi-implicit physical time + dual time for compressible flows. |
| 2026 | Ghidoni, Massa & Noventa, *Coupling between the time-step size adaptation and pseudo-transient continuation algorithms for an efficient time integration of steady solutions*. | Recent coupling of time-step adaptation and pseudo-transient continuation using residual/unsteadiness information. |

## CFDX consequence

The recent literature does not establish one universally superior replacement for dual-time stepping. The recurring architecture is a validated physical-time scheme, a nonlinear solve at every physical step, pseudo-time/Newton globalization, deterministic adaptation and independent convergence evidence. Compressible solvers additionally need physics-specific safeguards such as positivity and entropy stability.

For CFDX, the next research step worth evaluating is a **residual + temporal-error-aware controller** rather than a neural-network controller. The 2026 work is particularly relevant to that direction, while the 2025 ML work should remain a research benchmark rather than a core dependency.
