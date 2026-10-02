# CFDX Theory — Fluid Mechanics and Numerical CFD Course

CFDX Theory is the scientific course and mathematical reference for CFDX. It is deliberately separate from the Developer implementation contract and from V&V acceptance status.

The course has two complementary threads:

1. **Fluid mechanics:** continuum mechanics → kinematics → conservation → constitutive laws → Navier–Stokes → energy → compressibility → boundary layers → dimensionless analysis.
2. **Numerical CFD:** PDE → integral conservation law → mesh → reconstruction → fluxes → discretisation → algebraic system → pressure/velocity coupling → linear/nonlinear solution → convergence → verification → validation.

The executable Python files in this directory are the scientific source of truth. They use Jupytext percent format and are rendered by Sphinx/MyST-NB.

## Master learning path

```
A — FLUID MECHANICS
00 Foundations / Fluid mechanics
   ↓
01 Conservation laws

B — SPATIAL DISCRETISATION
02 Finite-volume method
   ↓
03 Meshes
   ↓
04 Gradients & reconstruction
   ↓
05 Fluxes

C — TIME + FLOW SOLUTION
06 Time integration
   ↓
07 Pressure–velocity coupling
   ↓
08 Linear algebra
   ↓
09 AMG / MGR / Schur

D — PHYSICAL CLOSURES
10 Turbulence
   ↓
11 Heat transfer
   ↓
12 Radiation
   ↓
13 Multiphysics

E — NUMERICAL EVIDENCE
14 Numerical analysis
   ↓
15 Verification & validation

F — CFDX EXECUTION
16 Data model & file formats
   ↓
17 How CFDX solves the equations
   ↓
18 Source-tree physics audit
```

## The course students should be able to follow

### Part I — Understand the fluid

Before touching a numerical scheme, the reader should be able to derive:

[
rac{Dphi}{Dt}
=
rac{partialphi}{partial t}
+mathbf ucdot
ablaphi,
]

[
rac{partialho}{partial t}
+
ablacdot(homathbf u)=0,
]

[
horac{Dmathbf u}{Dt}
=
-
abla p+
ablacdotoldsymbol	au+homathbf f,
]

and, for a Newtonian incompressible fluid,

[
ho
left(
rac{partialmathbf u}{partial t}
+mathbf ucdot
ablamathbf u
ight)
=
-
abla p+mu
abla^2mathbf u+homathbf f,
qquad

ablacdotmathbf u=0.
]

The meaning, units, assumptions and limits of every term are introduced before discretisation.

### Part II — Understand how the PDE becomes a CFD problem

The central transformation is

[
	ext{PDE}
ightarrow
int_V	ext{PDE},dV
ightarrow
	ext{face/volume contributions}
ightarrow
A(U)U=b(U).
]

For a generic transported quantity:

[
rac{partial(hophi)}{partial t}
+
ablacdot(homathbf uphi)
=

ablacdot(Gamma
ablaphi)+S.
]

Finite volume produces:

[
rac{d}{dt}int_Vhophi,dV
+
sum_f F_{c,f}
=
sum_fF_{d,f}
+
int_VS,dV.
]

The subsequent chapters explain how every term is approximated.

### Part III — Understand how the algebraic problem is solved

The course explicitly follows:

[
oxed{
A(U)U=b(U)
ightarrow
	ext{nonlinear iteration}
ightarrow
	ext{linear system}
ightarrow
	ext{preconditioned Krylov/AMG solve}
}
]

For incompressible flow:

[
egin{bmatrix}
A_u&G\
D&C
end{bmatrix}
egin{bmatrix}
u\p
end{bmatrix}
=
egin{bmatrix}
b_u\b_p
end{bmatrix},
]

with Schur complement

[
S=C-DA_u^{-1}G.
]

This explains why pressure–velocity coupling, Schur methods and AMG/MGR are not independent features: they arise from the mathematical structure of Navier–Stokes.

## Numerical-method traceability

The intended chain is:

```
physical law
    ↓
continuous equation
    ↓
integral balance
    ↓
geometric/discrete operator
    ↓
algebraic system
    ↓
nonlinear/linear solver
    ↓
diagnostics
    ↓
verification
    ↓
validation
    ↓
qualification
```

## Standard subsection contract

Where applicable, every important scientific subsection should contain:

1. physical motivation;
2. variables, dimensions and SI units;
3. assumptions;
4. continuous equation;
5. derivation;
6. integral/control-volume form;
7. discrete formulation;
8. algebraic contribution;
9. numerical properties;
10. failure modes;
11. exact CFDX source mapping;
12. executable verification;
13. benchmark/reference;
14. limitations;
15. improvement paths;
16. bibliography.

## Scientific integrity

- No invented CFDX results.
- No reference value presented as a CFDX result.
- No residual-only physical validation.
- No implementation-to-qualification shortcut.
- No tolerance relaxation to manufacture a green result.
- No failing benchmark is silently converted into a pass.
- Quantitative claims record the mesh, physics, numerical settings and software revision where relevant.

## Status vocabulary

- **Implemented:** code path exists.
- **Verified:** a defined mathematical/software property has executable evidence.
- **Validated:** independent physical/reference evidence supports the declared case.
- **Qualified:** declared population and acceptance gates are complete.

These states remain independent from documentation completeness.

## Main references

The course is aligned with the classical finite-volume presentation of Versteeg & Malalasekera and with modern finite-volume/numerical practice. OpenFOAM's current documentation likewise describes the general scalar transport equation, control-volume integration, face interpolation, gradient/divergence/laplacian schemes and algebraic solution as successive numerical layers. citeturn0search1turn0search10

See `docs/references/bibliography.bib` for the maintained bibliography.

## Chapter navigation

| # | Chapter | Role |
|---|---|---|
| 00 | Foundations | Complete fluid-mechanics introduction |
| 01 | Conservation laws | Derive conservative governing equations |
| 02 | Finite-volume method | Turn conservation laws into cell equations |
| 03 | Meshes | Geometry and topology |
| 04 | Gradients & reconstruction | Derivative and face reconstruction |
| 05 | Fluxes | Conservative transport fluxes |
| 06 | Time integration | Physical and pseudo-time discretisation |
| 07 | Pressure–velocity coupling | SIMPLE/PISO/PIMPLE/fractional-step/coupled methods |
| 08 | Linear algebra | Sparse systems and Krylov methods |
| 09 | AMG/MGR/Schur | Scalable coupled linear solution |
| 10 | Turbulence | RANS/LES/hybrid closure |
| 11 | Heat transfer | Thermal transport and CHT |
| 12 | Radiation | S2S/P1/DOM |
| 13 | Multiphysics | Coupled residual/Jacobian systems |
| 14 | Numerical analysis | Accuracy, stability and error |
| 15 | V&V | Evidence and qualification |
| 16 | Data model | Case/checkpoint/output semantics |
| 17 | Computational chain | Complete numerical solution workflow |
| 18 | Source-tree audit | Equation-to-code-to-test traceability |
