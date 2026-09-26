# CFDX — Ansys Fluent VMFL Verification Campaign

Reference: **Ansys Fluid Dynamics Verification Manual, Release 2026 R1**.

Official manual:
https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/pdf/Ansys_Fluid_Dynamics_Verification_Manual.pdf

The manual defines 78 Fluent/CFX verification cases (VMFL001–VMFL078). CFDX will reproduce the same configurations whenever its physics capabilities permit it. An analytical formula alone is not a Fluent validation PASS.

## Status definitions

- **PASS** — CFDX executes the relevant physical case and passes quantitative criteria.
- **READY** — required CFDX physics/reference infrastructure exists; implementation is next.
- **PARTIAL** — some required physics exists but the complete Fluent configuration is not represented.
- **BLOCKED** — a required physical/model capability is absent.

## VMFL matrix

| ID | Fluent case | CFDX status | Main comparison |
|---|---|---|---|
| VMFL001 | Rotating/stationary concentric cylinders | BLOCKED | velocity profile |
| VMFL002 | Laminar pipe, uniform heat flux | PARTIAL | pressure drop, outlet temperature |
| VMFL003 | Turbulent pipe pressure drop | PARTIAL | Δp, friction factor |
| VMFL004 | Plain Couette flow with pressure gradient | READY | velocity profile |
| VMFL005 | Poiseuille flow in pipe | READY | velocity profile, flow rate |
| VMFL006 | Multicomponent species pipe flow | BLOCKED | species profiles |
| VMFL007 | Non-Newtonian pipe flow | BLOCKED | velocity, Δp |
| VMFL008 | Rotating cavity | BLOCKED | radial/swirl velocity |
| VMFL009 | Natural convection concentric annulus | PARTIAL | temperature, Nu |
| VMFL010 | Laminar 90° tee | PARTIAL | pressure/velocity |
| VMFL011 | Laminar triangular cavity | PARTIAL | centreline velocity |
| VMFL012 | Turbulent wavy channel | PARTIAL | velocity profile |
| VMFL013 | Turbulent heat-transfer backward-facing step | PARTIAL | local Nu |
| VMFL014 | Species mixing coaxial turbulent jets | BLOCKED | species/velocity |
| VMFL015 | Engine inlet valve | BLOCKED | 3-D velocity |
| VMFL016 | Turbulent transition duct | PARTIAL | pressure coefficient |
| VMFL017 | Transonic RAE 2822 | BLOCKED | Cp/shock |
| VMFL018 | Supersonic shock reflection | BLOCKED | shock location |
| VMFL019 | Transient moving-wall flow | PARTIAL | transient velocity |
| VMFL020 | Adiabatic piston compression | BLOCKED | p-T history |
| VMFL021 | Cavitation, high inlet pressure | BLOCKED | cavity extent |
| VMFL022 | Cavitation, low inlet pressure | BLOCKED | cavity extent |
| VMFL023 | Oscillating cylinder | BLOCKED | lift/drag/St |
| VMFL024 | Immiscible liquids rotating cylinder | BLOCKED | interface |
| VMFL025 | Swirling methane combustion | BLOCKED | temperature/species |
| VMFL026 | Real-gas shock tube | BLOCKED | shock/p-T |
| VMFL027 | Turbulent backward-facing step | PARTIAL | reattachment |
| VMFL028 | Turbulent heat-transfer pipe expansion | PARTIAL | Nu, Δp |
| VMFL029 | Anisotropic conduction | PARTIAL | heat flux |
| VMFL030 | Turbulent 90° pipe bend | PARTIAL | Δp/secondary flow |
| VMFL031 | Turbulent V-gutter | PARTIAL | separated flow |
| VMFL032 | Turbulent axisymmetric afterbody | PARTIAL | separation |
| VMFL033 | Viscous heating annulus | PARTIAL | temperature |
| VMFL034 | Particle aggregation stirred tank | BLOCKED | particle distribution |
| VMFL035 | 3-D axial compressor | BLOCKED | pressure ratio |
| VMFL036 | Laminar flow past sphere | READY | drag coefficient at literature Re=100 |
| VMFL037 | Turbulent forward-facing step | PARTIAL | pressure/forces |
| VMFL038 | Falling film inclined plane | BLOCKED | film thickness |
| VMFL039 | Boiling pipe / CHF | BLOCKED | CHF |
| VMFL040 | Separated turbulent diffuser | PARTIAL | separation/pressure |
| VMFL041 | Transonic airfoil | BLOCKED | Cp |
| VMFL042 | Turbulent mixing of different densities | BLOCKED | mixing layer |
| VMFL043 | Laminar-turbulent flat-plate transition | PARTIAL | transition location |
| VMFL044 | Supersonic nozzle | BLOCKED | Mach/thrust |
| VMFL045 | Oblique shock | BLOCKED | shock angle |
| VMFL046 | Normal shock CD nozzle | BLOCKED | shock location |
| VMFL047 | Asymmetric diffuser separation | PARTIAL | separation |
| VMFL048 | Turbulent 180° bend | PARTIAL | Δp |
| VMFL049 | Axisymmetric natural-gas furnace | BLOCKED | temperature/species |
| VMFL050 | Transient conduction semi-infinite slab | READY | T(x,t) |
| VMFL051 | Isentropic expansion convex corner | BLOCKED | Mach/Cp |
| VMFL052 | Turbulent natural convection tall cavity | PARTIAL | Nu/velocity |
| VMFL053 | Compressible turbulent mixing layer | BLOCKED | mixing profile |
| VMFL054 | Trapezoidal cavity | PARTIAL | velocity |
| VMFL055 | Transitional ventilation enclosure | BLOCKED | recirculation |
| VMFL056 | Conduction + radiation square cavity | PARTIAL | T/heat flux |
| VMFL057 | Radiation + conduction composite layers | READY | heat flux/T |
| VMFL058 | Turbulent axisymmetric diffuser | PARTIAL | separation |
| VMFL059 | Composite solid block conduction | READY | heat flux/T |
| VMFL060 | Transitional supersonic rearward step | BLOCKED | shock/separation |
| VMFL061 | S2S radiation concentric cylinders | PARTIAL | net radiation |
| VMFL062 | Fully developed turbulent flow over hill | PARTIAL | velocity/shear |
| VMFL063 | Separated laminar flow over blunt plate | READY | separation/drag |
| VMFL064 | Low-Re asymmetric expansion | READY | pressure/recirculation |
| VMFL065 | Swirling turbulent diffuser | BLOCKED | swirl/pressure |
| VMFL066 | Participating-medium enclosure radiation | PARTIAL | radiative flux |
| VMFL067 | Boiling pipe CHF | BLOCKED | CHF |
| VMFL068 | Eccentric annulus | PARTIAL | pressure/velocity |
| VMFL069 | Two-phase Poiseuille | BLOCKED | phase profiles |
| VMFL070 | Radiation between parallel surfaces | READY | net q'' |
| VMFL071 | Goldman stator blade | BLOCKED | Cp/forces |
| VMFL072 | Water film over flat plate | BLOCKED | film/heat transfer |
| VMFL073 | Turbulent separated axisymmetric diffuser | PARTIAL | separation |
| VMFL074 | Plug-flow atomizer | BLOCKED | spray |
| VMFL075 | Supersonic circular-arc bump | BLOCKED | Cp/shock |
| VMFL076 | Forced convection flat plate | PARTIAL | Nu/Cf |
| VMFL077 | Free-surface ship flow | BLOCKED | wave resistance |
| VMFL078 | Polyhedral mesh accuracy | PARTIAL | profile/error |

## First executable wave

The first implementation wave uses the same Fluent case definitions and reported quantities for:

- VMFL004 — Couette + pressure gradient
- VMFL005 — Poiseuille
- VMFL036 — laminar sphere
- VMFL050 — transient semi-infinite conduction
- VMFL057 — radiation/conduction composite layers
- VMFL059 — composite solid conduction
- VMFL061 — concentric-cylinder S2S radiation
- VMFL063 — separated laminar blunt plate
- VMFL064 — low-Re asymmetric expansion
- VMFL070 — parallel-surface radiation
- VMFL078 — mesh accuracy/grid convergence

Existing CFDX analytical tests may be reused when they represent the same governing equations, but a test is not promoted to VMFL PASS until its Fluent geometry, material properties, boundary conditions and reported quantity are reproduced.

## Required result record

Every implemented VMFL case must report geometry/mesh, physical properties, boundary conditions, solver/discretisation controls, linear and nonlinear residuals, conservation error, quantity of interest, reference value/curve, absolute and relative error, and mesh/time refinement where applicable.

This keeps analytical component tests separate from complete CFD validation.