# FICHE MÉTHODE — Résidus non linéaires du solveur incompressible COUPLED

**En-tête d'évidence**
- Scope : normalisation des résidus de quantité de mouvement et critère de convergence non linéaire dans `solve_steady_incompressible`, en particulier l'algorithme COUPLED.
- Statut : CORRECTION CANDIDATE — tests locaux Release passés ; validation GitHub Actions de la branche requise avant fusion.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json`.
- Issue / PR liée : PR #1144.
- Dernier commit vérifié avant cette fiche : `940635ee1ab7853002030680548525be692d2b97` (master, 2026-10-28).
- Limitation : la fiche traite la normalisation de convergence, pas la qualification physique NACA0012 ni le couplage RANS/énergie/CHT/rayonnement.

## 1. Références sources (DOI ou ISBN)
- Ferziger, J. H., Perić, M. et Street, R. L., *Computational Methods for Fluid Dynamics*, 4e éd., ISBN 978-3-319-99692-4. Référence de méthodes numériques CFD et de résidus de volume fini.
- Nocedal, J. et Wright, S. J., *Numerical Optimization*, 2e éd., ISBN 978-0-387-30303-1. Référence pour les critères de convergence et la mise à l'échelle des résidus.

## 2. Formulation continue
La convergence non linéaire doit être évaluée sur les équations discrètes de quantité de mouvement et de continuité, avec une échelle physique stable. Une normalisation par le résidu de la première itération n'est pas une mesure fiable si ce résidu initial est déjà inférieur à la tolérance physique : le ratio peut rester proche de 1 même lorsque l'équation discrète est satisfaite à la précision demandée.

## 3. Discrétisation Volume Fini
Pour chaque composante de vitesse, CFDX conserve :
- le résidu absolu de l'équation discrète ;
- le résidu relatif à l'échelle du second membre (`momentum_equation_residual_rhs_relative`) ;
- le résidu relatif au résidu de référence de la première itération (`momentum_equation_residual_relative`).

Le résidu RHS-normalisé est le critère physique principal. Le ratio relatif à la première itération reste disponible pour l'historique et le diagnostic.

## 4. Traitement des frontières
La normalisation ne change aucune condition aux limites ni aucun flux de face. Elle ne modifie ni l'assemblage du système couplé ni le préconditionneur de Schur. Les bilans de continuité, le flux de masse corrigé et les bilans indépendants de quantité de mouvement restent des portes séparées.

## 5. Pièges connus
- Si le résidu de référence est très petit mais non nul, diviser par ce résidu crée un ratio proche de 1 qui peut dominer le moniteur.
- Un résidu relatif de référence égal à 1 ne signifie pas, à lui seul, que l'équation physique est mal satisfaite.
- À l'inverse, un petit résidu RHS-normalisé ne dispense pas des portes de continuité, de bornitude, de bilan indépendant et de comparaison avec la solution analytique.
- La correction ne doit ni augmenter la tolérance utilisateur, ni supprimer un critère physique, ni accepter un résidu non fini.

## 6. Ordre de précision attendu et validation exigée
Le seuil de résolution du résidu de référence est
`max(100 epsilon_machine, tolerance_relative) * momentum_rhs_scale`.
Si le résidu de référence est sous ce seuil, son ratio n'est pas utilisé comme critère indépendant du moniteur ; le résidu RHS-normalisé reste obligatoire. Si le résidu de référence est résolu au-dessus de ce seuil, le ratio de référence continue à contribuer au moniteur. Les tests doivent couvrir les cas COUPLED et les méthodes ségrégées sans relâcher les tolérances.

## 7. Statut d'implémentation CFDX
- Implémenté dans la PR #1144 : seuil de résolution du résidu de référence fondé sur la tolérance physique déjà déclarée et plancher machine.
- Contrat N9 mis à jour : les méthodes COUPLED sont évaluées par le résidu RHS-normalisé ; SIMPLE/PISO/PIMPLE et autres voies ségrégées conservent le ratio de référence.
- Sortie d'évidence : les deux normalisations sont affichées pour permettre une revue indépendante.
- Résultats locaux Release : `test_phase9_acceptance --quick` PASS (5/5 modèles), `test_phase9_acceptance` PASS (12/12 modèles), `test_dual_time_navier_stokes` PASS (1/1).
- Statut final : en attente des contrôles CI de la tête exacte de la PR ; aucun résultat CI distant n'est déclaré réussi dans cette fiche.

## 8. Journal des modifications de la fiche
- 2026-10-28 — PR #1144 — documente le défaut de normalisation, le correctif sans changement de tolérance et les résultats de reproduction locaux.
