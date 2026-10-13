# FICHE MÉTHODE 05 — Complément de Schur et préconditionnement : PCD, LSC, BFBT, MGR, AMG

**En-tête d'évidence**
- Scope : solveurs de pression et préconditionneurs des systèmes couplés (ségrégé et COUPLED). Couvre les 6 modèles Schur sélectionnables, les exigences de qualification par modèle et la parité MPI.
- Statut : APPROVED (2026-10-12)
- Issues liées : #481 (N8), #601 (N9), #663 ; PRs #804, #813, #853 ; `docs/validation/PCD_SCHUR_QUALIFICATION.md’, `N8_QUALIFICATION_TASKS.md’, `N8_MPI_QUALIFICATION_BOUNDARY.md’
- Dernier commit vérifié : 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : politique d'acceptation quantitative (section 6) à définitivement figer en Phase 1 (N8).

## 1. Références sources
- D. Kay, D. Loghin, A. Wathen, « A preconditioner for the steady-state Navier–Stokes equations », *SIAM J. Sci. Comput.* 24(1):239–253, 2002. (PCD)
- D. Loghin, A. J. Wathen, « Schur complement preconditioners for the Navier–Stokes equations », *Int. J. Numer. Meth. Fluids* 40:403–412, 2002 — DOI : 10.1002/fld.296. (LSC)
- H. C. Elman, D. J. Silvester, A. J. Wathen, *Finite Elements and Fast Iterative Solvers: with applications in incompressible fluid dynamics*, 2e éd., Oxford University Press, 2014. (BFBT, analyse coût)
- U. Trottenberg, C. W. Oosterlee, A. Schüller, *Multigrid*, Academic Press, 2001 — ISBN 978-0-12-701070-0. (AMG)

## 2. Formulation continue
Pression de correction pour le système de Stokes/Oseen ; le complément de Schur exact $S = C - D M^{-1} G$ est inaccessible ; on approxime $S^{-1}$.

## 3. Discrétisation Volume Fini (contrats par modèle)
- **PCD** : $S^{-1} \\approx (M_p^{-1} F_p)^{-1} F_p \\, L_p^{-1}$ — les opérateurs de pression ($M_p$, $L_p$, $F_p$) sont **rigoureusement indépendants du schéma de convection** car assemblés sur la matrice upwind conservatrice commune (matrice seule, correction différée au second membre) — décision #663/#853 : ne **pas** propager le knob `convection_scheme’ dans le PCD.
- **LSC** : $S^{-1} \\approx (D M_u^{-1} G)^{-1} (D M_u^{-1} F_u M_u^{-1} G)(D M_u^{-1} G)^{-1}$ — exige les moyennes de faces $A_f$.
- **BFBT** : factorisation $L D^T$ de masse — invariance par jauge à démontrer (leçon fix4 : jauge inconditionnelle corrigée par sentinelle).
- **MGR / AMG**échelonnage/hierarchies : sweep anisotrope documenté ; le passage de MGR est une correction finie → l'efficacité doit être **mesurée**, pas supposée.
- Chaque modèle doit être sélectionnable **dans le chemin de production complet** (`.cfdx.h5 → numerical_selection_report → cfdx_production_solver’, #813), pas seulement en test.

## 4. Traitement des frontières
Null-space : le préconditionneur doit préserver l'espace nul de pression (conditions tout-Neumann) ; pressions de référence et pressions imposées testées explicitement (N9 S10).

## 5. Pièges connus
- Jauge LSC/BFBT inconditionnelle → sentinelle (fix4) ; faux null-space → rejet (N9.7/#851).
- Un préconditionneur « fini » n'est pas un préconditionneur efficace : la réduction de résidu doit être **mesurée** sur des systèmes représentatifs de la population CFD (pas 4–6 cellules synthétiques seulement).
- MPI : AMG/Schur/PCD/MGR **non qualifiés** en décomposition de domaine (limite #481) ; aucune promotion MPI sans campagne halos/rangs.
- Restart FGMRES : jamais écrêté ; tolérances validées par modèle et par échelle (leçon #419 C4/C5).

## 6. Politique d'acceptation quantitative (exigée pour fermer #481)
Chaque modèle Schur doit fournir, sur des systèmes représentatifs (maillage réel, Re varié, n64/n128) :
1. erreur d'action vs complément de Schur exact quand la comparaison est pertinente ;
2. vrai résidu après application (recalculé indépendamment) ;
3. taux de réduction de résidu mesuré (pas seulement « fini ») ;
4. critères de convergence explicites sur cas physiques (Ghia, Poiseuille) ;
5. tests pression de référence / pression imposée ;
6. itérations et temps mesurés sur ≥ 2 tailles de maillage ;
7. parité MPI (rangs ≥ 2, même maillage, 1e-12) — **avant** toute qualification MPI.
Les seuils sont définis à partir des propriétés mathématiques et des cas de référence, **pas ajustés après coup**.

## 7. Statut d'implémentation CFDX
- Sources : `src/cfdx/core/numerics/…schur…’, `assemble_pcd_pressure_operators’ (#853) ; `src/cfdx/core/linalg/…’ (AMG, FGMRES).
- Tests : `test_n9_common_schur_qualification.cpp’ (à renforcer selon §6), `docs/validation/AMG_PRECONDITIONER_QUALIFICATION.md’, `MGR_PRECONDITIONER_QUALIFICATION.md’, matrice #481.
- Écarts code ↔ fiche : à auditer en Phase 1 (N8/N9) ; notamment qualification sur population CFD réelle.
- Gaps : matrice de qualification #481 incomplète ; MPI non qualifié.

## 8. Journal
- 2026-10-12 — création ; intègre #804, #813, #853, #851 et la revue projet (qualification N9 à renforcer).
- 2026-10-13 — durcissement N11 B1/B2/B3 (#480) : pivot relatif à l'échelle du bloc (BlockJacobi), garde d'entrée fail-closed sur `diagonal_schur_scale_` (Schur Diagonal), gardes GS/ILU(0) relatives à la rangée (pattern `factor_schur_ilu0`). Rejets fermés uniquement — aucune tolérance, gate ou politique d'acceptation modifiée.
- 2026-10-13 — correctif B3 (#480) : l'échelle de référence des gardes ILU(0) provient de la partie U de la rangée du pivot, jamais des multiplicateurs L (sans dimension), et chaque pivot est jugé dans sa propre rangée factorisée. Tests ajoutés : rescaling uniforme multi-échelles, rangées d'échelles contrastées, coefficients non finis (BlockJacobi fail-closed), état sûr après échec de setup. Rejets fermés uniquement — aucune tolérance, gate ou politique modifiée.
