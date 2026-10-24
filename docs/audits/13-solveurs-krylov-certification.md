# FICHE MÉTHODE — Solveurs Krylov, certification du résidu et précision mixte

**En-tête d'évidence (gouvernance §16.5)**
- Scope : CG préconditionné, GMRES/FGMRES, BiCGStab, réductions Krylov/MPI, précision mixte et oracle ExactSchur. Ne couvre pas la qualification complète des solveurs de production.
- Statut : DRAFT (2026-10-24)
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
- Issue / PR liées : audit G01–G09, PR #1092
- Dernier commit vérifié : correctifs post-audit en cours ; aucune qualification complète n'est revendiquée
- Limitations connues : la validation multi-rang MPI et la campagne QUICK quadratique multi-échelle doivent être exécutées sur le HEAD exact ; les opérateurs matrix-free sans callback FP64 ne peuvent être certifiés que dans la précision de leur callback.

## 1. Références sources (obligatoire : DOI ou ISBN)

- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, 2003. DOI: [10.1137/1.9780898718003](https://doi.org/10.1137/1.9780898718003).
- H. A. van der Vorst, *Iterative Krylov Methods for Large Linear Systems*, Cambridge University Press, 2003. ISBN: 978-0-521-81828-5.
- R. Barrett et al., *Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods*, SIAM, 1994. ISBN: 978-0-89871-328-2.

## 2. Formulation continue

Le problème algébrique est (Ax=b), obtenu après discrétisation des équations physiques. Pour un itéré (x_k), le résidu vrai est
\[
r_k=b-Ax_k.
\]
Une norme préconditionnée ou une estimation de Hessenberg est un outil de récurrence, pas une preuve de convergence du système original.

## 3. Discrétisation Volume Fini

Après assemblage FV, le système creux (A x=b) est résolu par CG (matrice SPD), GMRES/FGMRES (opérateur général) ou BiCGStab (opérateur non symétrique). Le critère d'acceptation doit porter sur
\[
\|b-Ax_k\|_2 \le \max(\mathrm{atol},\mathrm{rtol}\,\|b\|_2).
\]
Le résidu de certification doit être recalculé avec l'opérateur d'origine en FP64 lorsque celui-ci est disponible. Pour rtol seul, la référence est ||b||₂ si b est non nul ; uniquement pour b=0, elle est ||r₀||₂ afin qu'un démarrage non nul puisse converger vers la solution nulle. Le champ `residual_relative` utilise la même référence que le critère d'arrêt. Pour les vecteurs distribués, les normes et produits scalaires utilisent exclusivement le domaine owned, puis une réduction globale ; les ghosts sont actualisés mais ne sont pas recomptés.

## 4. Traitement des frontières

Les conditions aux limites physiques sont incorporées lors de l'assemblage FV dans (A) et (b), avant l'appel au solveur. Cette fiche traite le système assemblé et ne redéfinit pas les contributions Dirichlet/Neumann de chaque équation. Les lignes de jauge et espaces nuls doivent respecter le contrat propre à chaque système ; un second membre nul ne justifie pas de remplacer arbitrairement une solution appartenant au noyau par zéro.

## 5. Pièges connus (Gotchas)

- CG : (sqrt{r^T M^{-1}r}) n'est pas (|r|_2) ; un seuil de convergence doit rester invariant sous la multiplication commune de (A) et (b).
- BiCGStab : les seuils absolus appliqués aux produits scalaires dépendent des unités et de l'échelle ; les ruptures doivent être évaluées relativement aux normes des vecteurs concernés.
- GMRES/FGMRES : si (b=0) et (x_0\ne0), le solveur doit essayer de réduire (Ax_0), pas déclarer immédiatement une divergence. Si (b\ne0), la référence reste (||b||₂), même lorsque (||r₀||₂) est très supérieur ; employer le résidu initial dans ce cas relâcherait silencieusement la tolérance.
- Précision mixte : calculer une norme FP64 d'un résidu obtenu avec un matvec FP32 ne restaure pas les bits perdus ; le résidu vrai doit utiliser l'opérateur FP64 d'origine.
- Arnoldi : projections et réorthogonalisation doivent partager la même politique owned/MPI que les normes.
- NaN/Inf : les valeurs non finies sont rejetées avant les chemins de second membre nul ; le résidu de sortie ne doit pas rester à zéro par défaut après échec.
- Matrice-free : sans application FP64 séparée, la certification de l'opérateur original n'est pas disponible et doit être signalée comme telle.

## 6. Ordre de précision attendu et validation exigée

- Ordre attendu : les solveurs linéaires ne sont pas une discrétisation spatiale ; ils doivent satisfaire le contrat de résidu pré-déclaré et résoudre la même solution sous changement d'échelle, dans les limites du conditionnement et de FP64.
- Oracles quantitatifs : solutions exactes de petits systèmes SPD et non symétriques ; comparaison de (\|b-Ax\|_2/\|b\|_2) avec la tolérance déclarée.
- GCI : non applicable au noyau linéaire isolé ; les études de convergence FV doivent être conduites au niveau de l'équation discrétisée et de ses maillages.
- Conservation : le résidu linéaire ne remplace pas les bilans indépendants de masse, quantité de mouvement, énergie ou scalaire.
- Parallèle : Release et DebugSanitizers ; pour les chemins MPI, partitions différentes, rafraîchissement des halos, réductions collectives et comparaison CPU/MPI à tolérance déclarée. La qualification multi-rang n'est pas encore attestée par cette PR.
- Régressions requises : (sA x=s b) pour plusieurs (s), (b=0,x_0\ne0), coefficient FP64 non représentable en FP32, plage owned excluant un ghost, et RHS NaN/Inf pour l'oracle.

## 7. Statut d'implémentation CFDX

- Fichiers source concernés : src/cfdx/core/linalg/cg_solver.h, bicgstab_solver.h, gmres_solver.h, mixed_precision.h, exact_schur.h, krylov_reductions.h, src/cfdx/core/numerics/gradient.h, interpolation.h.
- Tests exécutables : tests/unit/test_cg_solver.cpp, test_bicgstab_solver.cpp, test_gmres_solver.cpp, test_exact_schur.cpp, test_gradient.cpp, test_interpolation.cpp.
- Écarts connus code ↔ fiche : l'API GMRES matrix-free sans callback FP64 ne peut certifier que son opérateur fourni ; le test MPI multi-rang et la campagne quadratique QUICK sur plusieurs échelles doivent encore être exécutés sur le HEAD exact.
- Gaps et prochaines actions : attendre les résultats CI Release/DebugSanitizers ; corriger tout échec ; lancer une qualification multi-rang pour les réductions distribuées ; exécuter le test QUICK quadratique sur plusieurs tailles géométriques et documenter l'erreur L2/L∞.

## 8. Journal des modifications de la fiche

- 2026-10-24 — PR #1092 — création de la fiche pour l'audit G01–G09 ; statut DRAFT jusqu'à la qualification sur HEAD exact.
