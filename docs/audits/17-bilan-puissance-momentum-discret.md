# FICHE MÉTHODE 17 — Validation indépendante de puissance momentum

**Statut :** DRAFT — campagne de validation logicielle, pas une preuve de stabilité ni une qualification physique.

## Objectif

Comparer deux évaluations algébriquement équivalentes mais construites par des chemins distincts :

1. **Oracle face brut** : forces vectorielles intégrées [N] stockées une fois par face, orientées owner → neighbour, et sources cellulaires intégrées [N]. Pour une face interne, le travail est \((\mathbf U_P-\mathbf U_N)\cdot\mathbf F_f\). Pour une face frontière, il est \(\mathbf U_P\cdot\mathbf F_f\).
2. **Scatter cellulaire indépendant** : assemblage signé des forces de face vers les cellules, puis produit scalaire par la vitesse cellulaire.

Les deux routes doivent donner la même puissance totale, à l'arrondi flottant près. L'oracle expose séparément travail des faces internes et frontières pour éviter de confondre échange interne conservatif et travail imposé par les frontières.

## Contrat et unités

- Forces de face et sources cellulaires déjà intégrées : [N].
- Vitesse : [m/s].
- Puissance : [W].
- Chaque face interne apparaît une seule fois.
- Le voisin d'une face frontière vaut exclusivement `FaceOwnership::BOUNDARY` (-1).
- La traction de Reynolds, la diffusion moléculaire et la traction visqueuse transposée restent des canaux distincts.
- Les sources pression, force extérieure et accumulation sont des forces intégrées cellulaires ; aucune multiplication supplémentaire par le volume n'est permise.

## Couverture de tests ajoutée

- Maillage deux cellules avec faces internes et frontières.
- Comparaison oracle saut de vitesse / scatter cellulaire indépendant.
- Décomposition explicite interne / frontière.
- Invariance au renversement de l'orientation d'une face interne si la force est également inversée.
- Rejet des tailles incompatibles, propriétaire/voisin invalides, sentinelle négative inconnue, NaN et débordement.
- Pas de seuil physique de dissipation imposé : pression, convection, frontières et sources peuvent légitimement injecter ou extraire de l'énergie.

## Limites scientifiques

Cette campagne vérifie une identité de comptabilité discrète sur des forces intégrées fournies au contrat. Elle ne prouve pas encore que les forces passées à l'oracle sont exactement celles utilisées par les opérateurs de production à l'état accepté. Le raccordement aux sorties brutes des opérateurs reste requis.

Elle ne prouve pas :
- la positivité de la dissipation visqueuse sur maillage non orthogonal ;
- la stabilité du schéma temporel ou du couplage pression-vitesse ;
- la réalisabilité d'une contrainte de Reynolds ;
- la précision d'une fermeture turbulente.

La suite devra capturer les flux bruts réels dans les chemins de production, comparer les puissances opérateur par opérateur, puis réaliser des tests de convergence maillage/temps et des comparaisons a-posteriori sur cas de référence tenus à l'écart de la calibration.

## Lien avec Lean

Les contrats Lean owner/neighbour et de puissance source constituent l'oracle mathématique abstrait. Les tests C++ vérifient des cas runtime concrets. La correspondance générale entre topologie Lean, `Mesh::ownership()` et flux de production n'est pas encore formalisée.
