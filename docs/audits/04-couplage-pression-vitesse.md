# FICHE MÉTHODE 04 — Couplage pression-vitesse : SIMPLE, SIMPLEC, PISO, PIMPLE, COUPLED

**En-tête d'évidence**
- Scope : algorithmes de couplage pression-vitesse pour incompressible. Couvre la boucle externe, la relaxation, les correcteurs, l'équation de pression et la jauge.
- Statut : APPROVED (2026-10-12)
- Issues liées : #601 (N9), #481 (N8) ; PRs #804, #813, #821, #851
- Dernier commit vérifiéà: 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : audit code ↔ fiche (section 7) à exécuter en Phase 1 (N9).

## 1. Références sources
- S. V. Patankar, D. B. Spalding, « A calculation procedure for heat, mass and momentum transfer in three-dimensional parabolic flows », *Int. J. Heat Mass Transfer* 15:1787–1806, 1972. (SIMPLE)
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI (rééd.) : 10.1201/9781482234213.
- J. P. van Doormaal, G. D. Raithby, « Enhancements of the SIMPLE method for predicting incompressible fluid flows », *Numerical Heat Transfer* 7:147–163, 1984 — DOI : 10.1080/01495728408961817. (SIMPLEC)
- R. I. Issa, « Solution of the implicitly discretised fluid flow equations by operator-splitting », *J. Comput. Phys.* 62(1):40–60, 1986 — DOI : 10.1016/0021-9991(86)90099-9. (PISO)
- Moukalled, Mangani, Darwish (2016), ISBN 978-3-319-16873-9, ch. 11–13 ; Ferziger & Perić (2002), DOI 10.1007/978-3-642-56026-2.

## 2. Formulation continue
Navier–Stokes incompressible :
$$ \\nabla\\cdot\\mathbf{u}=0, \\qquad \\partial\\mathbf{u}/\\partial t + \\nabla\\cdot(\\mathbf{u}\\otimes\\mathbf{u}) = -\\nabla p/\\rho + \\nabla\\cdot(\\nu\\nabla\\mathbf{u}) + \\mathbf{f} $$
Système discret cible (contrat N9.0) : $[M\\;G \\;;\\;D\\;C]$ — les opérateurs doivent être **testables indépendamment du choix SIMPLE/PISO/COUPLED**.

## 3. Discrétisation Volume Fini (contrats par algorithme)
**SIMPLE** : (1) prédiction $A\\,\\mathbf{u}^* = b - G\\,p^{n}$ (relaxation **implicite** sur la matrice : $a_P/\\alpha$) ; (2) équation de pression $D\\,(G\\,p') = \\nabla\\cdot\\mathbf{u}^*$ avec $rAU = V/a_P$ calculé sur le $a_P$ **non relaxé** ; (3) correction $\\mathbf{u} = \\mathbf{u}^* - rAU\\,G\\,p'$ ; (4) mise à jour des flux.
**SIMPLEC** : $rAU$ avec dénominateur consistant $(a_P - \\Sigma a_{nb})$ — convention signée par la PR #821 ; diagnostics exigés (dénominateur effondré/négatif = rejet, #851).
**PISO** : prédiction + 2 correcteurs (gain $1-\\Sigma$ des correcteurs ; politique du nombre de correcteurs déclarée, rejet si gain invalide #851).
**PIMPLE** : PISO à l'intérieur d'une boucle externe SIMPLE — même contrats.
**COUPLED**équivalence algébrique avec le ségrégé doit être **démontrée sur le même problème discret** avant toute promotion (N9.0) ; ce n'est pas encore le cas.

## 4. Traitement des frontières
- Sortie de pression : p imposé, correction $p'=0$ sur la face (masse conservée par reconstruction du flux).
- WALL : pas de correction de flux normal.
- EMPTY : ligne/vecteur vide du Laplacien de pression.
- Périodique : couplage dans $D\\,G$.

## 5. Pièges connus
- **Jauge de pression / null-space** : opérateur tout-Neumann singulier → résolution sur l'espace projeté + contrôle du vrai résidu ; refuser les RHS incompatibles (moyenne non nulle) — N9.7/#851 verrouillent ces rejets. Aucun ancrage arbitraire d'une cellule sans déclaration.
- Relaxation explicite de la QDM + $rAU$ non relaxé → divergence (leçon A-B2 #419).
- Convergence « MAX_ITER » ≠ convergence : la quantité de continuité doit décroître (leçon N9 cavity #538).
- SIMPLEC : dénominateur négatif sur maillage grossier/très convectif → rejet explicite (pas de fallback silencieux).

## 6. Ordre de précision attendu et validation exigée
- Le couplage n'introduit **aucune erreur de troncature propre** : convergence = exactitude (l'erreur stationnaire doit → erreur de discrétisation, pas erreur d'algorithme).
- Validation exigée (critère de fermeture #601) : campagne complète Couette, Poiseuille, **Ghia Re=100** (DOI : 10.1016/0021-9991(82)90058-4, voir fiche 07), Couette skewed, pour les 6 modèles de Schur, sur le HEAD exact, vrais résidus, flux reconstruits indépendamment, chemin complet `.cfdx.h5 → numerical_selection_report → cfdx_production_solver’, bilan de masse $<10^{-10}$, GCI sur Poiseuille.

## 7. Statut d'implémentation CFDX
- Sources : `src/cfdx/physics/steady_incompressible_solver.h’, `src/cfdx/core/numerics/schur*.h’, chemin `.cfdx.h5’.
- Tests : `tests/validation/test_n9_common_schur_qualification.cpp’ (à renforcer : cf. projet chat — exiger erreur d'action vs Schur exact, réduction de résidu mesurée), `test_n9_negative_diagnostics’ (#851), fix4b/fix4c.
- Écarts code ↔ fiche : **à auditer en Phase 1** — équivalence COUPLED/ségrégé non démontrée ; campagne physique non verte sur HEAD exact.
- Gaps : fermeture #601 conditionnée par la campagne 6 modèles × 4 cas.

## 8. Journal
- 2026-10-12 — création ; intègre #821 (SIMPLEC), #851 (N9.7), leçons #419.
