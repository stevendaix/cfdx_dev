# FICHE MÉTHODE 23 — Compatibilité AMG/CG après mise à jour numérique

**En-tête d’évidence**
- Scope : réutilisation numérique de la hiérarchie AMG native, stabilité du lisseur Jacobi pondéré et préconditions d’applicabilité de CG. Ne couvre pas la qualification physique complète du cas NACA0012.
- Statut : DRAFT — CI exact-HEAD et campagnes ciblées en cours.
- Registre : audit runtime NACA0012 du 2026-10-25 ; PR #1149.
- Sources : `src/cfdx/core/linalg/amg_preconditioner.h`, `src/cfdx/core/linalg/cg_solver.h`, `src/cfdx/core/linalg/linear_solver_context.h`.
- Contrat : aucun changement de tolérance, budget d’itérations, discrétisation CFD ou gate physique.

## 1. Références sources (obligatoire : DOI ou ISBN)

- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd edition, SIAM, 2003, ISBN 978-0-89871-534-7 — méthodes de Krylov et préconditionnement.
- R. S. Varga, *Matrix Iterative Analysis*, 2nd edition, Springer, 2000, ISBN 978-3-540-66311-8 — rayon spectral, disques de Gershgorin et convergence des méthodes stationnaires.

## 2. Formulation continue

CG résout le système linéaire discret (Ax=b) sous l’hypothèse que l’opérateur est symétrique défini positif (SPD). Le préconditionneur doit également fournir une action compatible avec cette hypothèse. En particulier, le produit initial (r_0^T M^{-1}r_0) doit être strictement positif. Une valeur finie mais non positive ne doit pas être utilisée pour former un pas CG.

## 3. Discrétisation Volume Fini

La hiérarchie AMG utilise des opérateurs de Galerkin (A_c=P^T A_f P) et un lisseur Jacobi pondéré. Pour (D=mathrm{diag}(A)), le facteur modal du lisseur est (1-omegalambda), où (lambda) est une valeur propre de (D^{-1/2}AD^{-1/2}). Une condition de stabilité suffisante est
[
0 < omega < rac{2}{lambda_{max}}.
]
Lorsque les diagonales sont positives et finies, le rayon spectral est borné par la norme infinie et par la borne conservative
[
lambda_{max}leq max_isum_j
left|rac{A_{ij}}{sqrt{A_{ii}A_{jj}}}ight|.
]
La PR #1149 recalcule cette borne pour chaque niveau après actualisation des coefficients et réduit (omega) à (min(omega_{mathrm{ancien}},1.8/lambda_{mathrm{borne}})). Si la borne n’est pas calculable, la mise à jour est refusée.

## 4. Traitement des frontières

Cette correction est algébrique : elle ne modifie aucune condition limite, aucun flux de face, aucune jauge de pression ni aucune ligne de matrice CFD. Elle ne modifie pas la discrétisation SIMPLE/SIMPLEC/COUPLED. Si une matrice de pression ou ses diagonales ne satisfont pas les hypothèses nécessaires au calcul de la borne, l’actualisation numérique échoue explicitement et le contexte propriétaire tente un setup complet ; aucune méthode de solveur de substitution n’est introduite.

## 5. Pièges connus

- L’audit NACA0012 a reproduit une hiérarchie réutilisée dont les niveaux grossiers sont devenus incompatibles avec le facteur Jacobi fixe. Le produit préconditionné initial est devenu négatif et CG a retourné `NOT_APPLICABLE`.
- Sur la matrice capturée, le résidu relatif vrai indépendant était environ (3.14\times10^{-2}), alors que les champs de résultat enregistrés étaient nuls.
- La trace montre aussi une divergence non linéaire momentum antérieure à l’échec linéaire. Le défaut AMG/CG est réel mais ne prouve pas que la qualification complète NACA0012 sera réparée.
- La borne de Gershgorin est conservative et peut réduire le facteur de relaxation plus que nécessaire, augmentant le coût des cycles.
- Le résidu récursif ou le statut seul ne remplace jamais le calcul indépendant de (b-Ax). Aucun seuil n’est relâché.

## 6. Ordre de précision attendu et validation exigée

- Tests de régression pour un préconditionneur à action négative : statut `NOT_APPLICABLE`, itération zéro, solution inchangée et résidu vrai correctement rapporté.
- Test de cycle de vie : refus d’une mise à jour numérique suivi d’un setup complet explicite par le contexte, avec compteurs distincts.
- CI Release et DebugSanitizers sur le SHA exact, tests numériques et gate d’audit.
- Rejouer NACA0012 N64/N128/N256 et contrôler statut CG, résidu vrai de pression, résidus momentum RHS-normalisés, continuité, flux-vitesse et coefficients.
- Rejouer l’étude de profils Ghia séparément ; enregistrer temps des solves de pression/AMG et décroissance momentum. Un timeout n’est pas à lui seul un échec de comparaison des profils.
- La réussite algébrique ne vaut pas preuve de convergence non linéaire, de conservation, de convergence en maillage ou d’accord physique du Cd/cisaillement.

## 7. Statut d’implémentation

- Sources modifiées : `src/cfdx/core/linalg/amg_preconditioner.h`, `src/cfdx/core/linalg/cg_solver.h`, `src/cfdx/core/linalg/linear_solver_context.h`.
- Tests ajoutés : `tests/validation/test_amg_preconditioner_qualification.cpp`.
- PR candidate : #1149.
- État des preuves : CI en cours ; la première exécution de la règle de gouvernance a identifié un défaut de structure de cette fiche, corrigé dans cette version. Les résultats numériques de la campagne exacte restent en attente.
- Décision : **PENDING**. Aucun merge avant CI vert et revue des campagnes ciblées.

## 8. Journal des modifications

- 2026-10-29 — création dans #1149 pour documenter la stabilité du lisseur sur hiérarchie AMG rafraîchie, le diagnostic de positivité CG et les critères de requalification.
