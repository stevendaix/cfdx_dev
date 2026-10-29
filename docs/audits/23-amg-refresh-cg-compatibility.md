# FICHE MÉTHODE 23 — Compatibilité AMG/CG après mise à jour numérique

**En-tête d’évidence**
- Scope : réutilisation numérique de la hiérarchie AMG native, stabilité du lisseur Jacobi pondéré et préconditions d’applicabilité de CG. Ne couvre pas la qualification physique complète du cas NACA0012.
- Statut : DRAFT — CI exact-HEAD et campagnes ciblées en cours.
- Registre : audit runtime NACA0012 du 2026-10-25 ; PR #1149.
- Sources : `src/cfdx/core/linalg/amg_preconditioner.h`, `src/cfdx/core/linalg/cg_solver.h`, `src/cfdx/core/linalg/linear_solver_context.h`.
- Contrat : aucun changement de tolérance, budget d’itérations, discrétisation CFD ou gate physique.

## 1. Références sources
- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd edition, SIAM, 2003, ISBN 978-0-89871-534-7 — méthodes de Krylov et préconditionnement.
- R. S. Varga, *Matrix Iterative Analysis*, 2nd edition, Springer, 2000, ISBN 978-3-540-66311-8 — rayon spectral, disques de Gershgorin et convergence des méthodes stationnaires.

## 2. Formulation numérique

Pour CG appliqué à (Ax=b), l’opérateur doit être SPD et le préconditionneur doit avoir une action compatible avec une matrice SPD. En particulier, le produit initial
[
r_0^T M^{-1}r_0
]
doit être strictement positif. Une valeur finie mais non positive ne doit pas être utilisée pour former un pas CG.

Pour le lisseur Jacobi pondéré, le facteur modal est (1-omegalambda), où (lambda) est une valeur propre de (D^{-1/2} A D^{-1/2}). Une condition de stabilité suffisante est
[
0 < omega < rac{2}{lambda_{max}}.
]
Le rayon spectral est borné par la norme infinie, donc par la borne conservative
[
lambda_{max}leq max_isum_j
left|rac{A_{ij}}{sqrt{A_{ii}A_{jj}}}ight|
]
lorsque les diagonales sont positives et finies.

## 3. Défaut observé et correction

L’audit NACA0012 a reproduit un cas où la hiérarchie d’interpolation conservée, suivie d’une mise à jour des opérateurs de Galerkin, produit des niveaux grossiers incompatibles avec le facteur Jacobi fixe. Le premier produit préconditionné devient négatif et CG retourne `NOT_APPLICABLE`. Sur la matrice capturée, le résidu relatif vrai indépendant était environ (3.14\times10^{-2}), alors que les champs de résultat enregistrés étaient nuls. La trace montre également une divergence non linéaire momentum antérieure : le défaut AMG/CG est réel, mais n’explique pas à lui seul la qualification NACA complète.

La PR #1149 :
- recalcule la borne normalisée sur chaque niveau rafraîchi ;
- réduit (omega) à (min(omega_{mathrm{ancien}},1.8/lambda_{mathrm{borne}})) si les diagonales permettent ce calcul ;
- refuse la mise à jour si un niveau ne peut pas être borné, afin que le contexte propriétaire tente un `setup(A)` complet ;
- fait rejeter par CG un (r^T M^{-1}rleq0) avant la première mise à jour de la solution, en conservant le résidu (b-Ax) indépendamment recalculé.

Le contexte incrémente les compteurs de setup complet ou de mise à jour numérique selon le chemin réellement exécuté. Aucun solveur alternatif silencieux n’est substitué.

## 4. Tests et preuves exigés

- Régression : préconditionneur à action négative, statut `NOT_APPLICABLE`, itération zéro et résidu vrai non nul correctement rapporté.
- Régression : refus d’une mise à jour numérique suivi d’un setup complet explicite par le contexte.
- Qualification AMG native et tests Release / DebugSanitizers sur le SHA exact.
- Rejouer la campagne NACA0012 complète N64/N128/N256, contrôler résidu vrai de pression, statut CG, résidus momentum RHS-normalisés, continuité, flux-vitesse et coefficients.
- Rejouer l’étude de profils Ghia séparément et mesurer les coûts AMG/solveur et la décroissance du résidu momentum ; un timeout ne constitue pas à lui seul un échec de comparaison de profil.

## 5. Limites et décision

La borne de Gershgorin est conservative ; la baisse de (omega) peut augmenter le coût des cycles. Le succès du test algébrique ne prouve ni convergence non linéaire, ni conservation, ni convergence en maillage, ni accord physique du Cd ou du cisaillement. Si la campagne NACA reste divergente avant l’échec linéaire, la correction de discrétisation / reconstruction près du profil doit être investiguée séparément.

Statut final : **PENDING** jusqu’à CI vert, true-residual evidence et revue des campagnes ciblées. Aucun merge avant ces preuves.
