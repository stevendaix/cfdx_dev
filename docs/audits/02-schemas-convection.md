# FICHE MÉTHODE 02 — Schémas de convection (UPWIND, SOU, MUSCL/TVD, QUICK, CENTRAL, BLENDED)

**En-tête d'évidence**
- Scope : interpolation convective aux faces pour tous les scalaires et la QDM. Non couvert : flux de masse (fiche 03/04) et limiteurs d'ordre élevé WENO (fiche séparée N16).
- Statut : APPROVED (2026-10-12)
- Issues liées : #461, #663 (audit schémas convection), #34 ; PR #853 (opérateurs PCD indépendants du schéma)
- Dernier commit vérifié : 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : audit code ↔ fiche (section 7) à exécuter en Phase 1.

## 1. Références sources
- B. P. Leonard, « A stable and accurate convective modelling procedure based on quadratic upstream interpolation », *Comput. Methods Appl. Mech. Eng.* 19(1):59–98, 1979 — DOI : 10.1016/0045-7825(79)90034-3. (QUICK)
- B. van Leer, « Towards the ultimate conservative difference scheme V. A second-order sequel to Godunov's method », *J. Comput. Phys.* 32(1):101–136, 1979 — DOI : 10.1016/0021-9991(79)90145-1. (MUSCL)
- P. K. Sweby, « High resolution schemes using flux limiters for hyperbolic conservation laws », *SIAM J. Sci. Stat. Comput.* 5(2):305–321, 1984. (TVD)
- Moukalled, Mangani, Darwish (2016), ISBN 978-3-319-16873-9, ch. convection ; Ferziger & Perić (2002), DOI 10.1007/978-3-642-56026-2, § 4.2–4.5.

## 2. Formulation continue
$$ \\nabla\\cdot(\\rho\\,\\mathbf{u}\\,\\phi) \\quad\\text{sur une face} : \\quad \\dot m_f\\,\\phi_f, \\qquad \\dot m_f = \\rho_f\\,\\mathbf{u}_f\\cdot\\mathbf{S}_f $$
Le schéma définit $\\phi_f$ à partir des valeurs $\\phi_P, \\phi_N$ et des gradients.

## 3. Discrétisation Volume Fini (contrat de correction différée)
L'implémentation CFDX suit le paradigme **upwind implicite + correction explicite** :
1. La **matrice** est toujours l'opérateur upwind 1er ordre conservatif : $a_{on} = D_f + \\max(\\dot m_f, 0)$, $a_{no} = D_f + \\max(-\\dot m_f, 0)$.
2. Les schémes d'ordre ≥ 2 (SECOND_ORDER_UPWIND, TVD/MUSCL minmod, CENTRAL, BLENDED, QUICK, QUICK_BOUNDED) ajoutent **uniquement un terme explicite au second membre** (`deferred_rhs’) évalué sur le champ convecté.
3. Le limiteur TVD/minmod garantit la bornitude : $\\phi_f \\in [\\min(\\phi_P,\\phi_N),\\, \\max(\\phi_P,\\phi_N)]$.
Conséquence contractuelle : **tout opérateur ne dépendant que de la matrice (ex. PCD) est rigoureusement indépendant du schéma** — statut confirmé par la PR #853 : l'UPWIND « en dur » dans `assemble_pcd_pressure_operators’ est déjà l'opérateur exact, et propager le knob serait une erreur.

## 4. Traitement des frontières
- Advection entrante (Dirichlet) : $\\phi_f = \\phi_b$ implicite dans les coefficients upwind.
- Advection sortante : $\\phi_f = \\phi_P$ (extrapolation d'ordre 0 — déclarer tout schéma d'extrapolation s'il est utilisé).
- La correction différée doit être **nulle aux faces frontières où $\\phi_f$ est imposé** (sinon dédoublement du flux).

## 5. Pièges connus
- QUICK non borné → oscillations : seul QUICK_BOUNDED (limiteur) peut être utilisé en production ; QUICK nu est un schéma de vérification d'ordre.
- Le nombre de Péclet de maille $Pe = \\dot m_f/D_f$ : la correction différée doit être ré-évaluée sur le champ courant à chaque itération externe ; la stocker figée casse la conservation.
- BLENDED : la fraction de mélange doit être déclarée dans le rapport de cas (traçabilité).
- Convection + correction différée et relaxation : appliquer la relaxation à la solution, jamais à la correction figée.

## 6. Ordre de précision attendu et validation exigée
- UPWIND : $O(h)$. SOU/MUSCL/QUICK : $O(h^2)$ sur champs lisses ; TVD : ordre 1 aux extrema (théorème de Godunov — documenter, ne pas « corriger »).
- Validation exigée : MMS advectif ≥ 3 niveaux avec pente log–log par schéma ; test de bornitude (chapeau/carré advecté, clipping interdit → limiteurs doivent contenir sans couper) ; bilan $\\sum \\dot m_{in}-\\sum \\dot m_{out} < 10^{-10}$ ; Ghia Re=100 par schéma admis en production.

## 7. Statut d'implémentation CFDX
- Sources : `src/cfdx/physics/finite_volume_transport.h’ (`assemble_scalar_equation’, `ConvectionScheme’), `src/cfdx/core/numerics/…’.
- Tests : `tests/validation/…’ (N12 MMS QUICK, N16 WENO5Z pour l'ordre élevé) ; audit #663 → #853.
- Écarts code ↔ fiche : **à auditer en Phase 1** — notamment : correction nulle aux faces Dirichlet, ré-évaluation de la correction à chaque itération, bornitude QUICK_BOUNDED.
- Gaps : benchmark TVD quantitatif manquant (héritage #419 : « a quantitative TVD benchmark is still required »).

## 8. Journal
- 2026-10-12 — création ; incorpore la conclusion de la PR #853 (opérateurs de pression indépendants du schéma).
- 2026-10-24 — audit des équations F10 : accès au gradient restreint aux schémas qui construisent ce champ (PR #1070, validation CI en cours). Cette correction de sûreté mémoire ne vaut pas qualification numérique des schémas ni clôture des écarts code↔fiche listés ci-dessus.
