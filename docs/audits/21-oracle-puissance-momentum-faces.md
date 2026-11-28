# FICHE MÉTHODE — Oracle indépendant de puissance momentum par faces

**En-tête d’évidence (gouvernance §16.5)**

- Scope : oracle d'audit algébrique du travail mécanique des forces momentum intégrées aux faces et aux cellules. Cette fiche couvre la reconstruction indépendante par scatter owner/neighbour et la validation de la topologie utilisée par cet oracle. Elle ne définit ni un schéma de convection/diffusion, ni un critère de convergence CFD.
- Statut : DRAFT — contrat d'audit documenté ; validation exacte-head et revue de la preuve encore requises.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` et exigences numériques du dépôt.
- Issue / PR liées : #1220 / #1230.
- Dernier commit vérifié : `f4ad5cb3a6851f1a2913a2174091360ca44cba10` — correction du sentinelle négatif ; CI de ce SHA en cours au moment de la création.
- Limitations connues : l'oracle vérifie une identité de travail discrète à partir de forces déjà intégrées. Il ne valide pas à lui seul la qualité physique des flux, l'ordre de précision, la stabilité, la positivité ou la convergence du solveur.

## 1. Références sources

- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI: 10.1201/9781482234213. Référence générale pour la formulation volumes finis et les bilans conservatifs.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN: 978-3-319-16874-6. Référence pour l'assemblage face-cellule et les flux conservatifs.

## 2. Formulation continue

Pour une quantité de mouvement \\(\rho\mathbf{u}\\), on considère le résidu physique \\(\mathbf{R}\\) de l'équation de quantité de mouvement et son appariement avec la vitesse :
\\[
\mathcal{P}=\int_\Omega \mathbf{u}\cdot\mathbf{R}\,d\Omega.
\\]
Cette quantité est un bilan de travail mécanique. La fiche ne suppose pas que \\(\mathcal{P}\\) soit positif ou nul : sa valeur dépend de l'état, des flux, des sources et des conventions de signe. L'objectif est de vérifier que la reconstruction discrète des contributions est cohérente avec les forces intégrées effectivement fournies à l'oracle.

## 3. Discrétisation Volume Fini

Pour chaque face interne orientée owner \\(P\\) vers neighbour \\(N\\), une force intégrée unique \\(\mathbf{F}_f\\) est associée à la face. Le scatter conservatif ajoute \\(+\mathbf{F}_f\\) à la cellule owner et \\(-\mathbf{F}_f\\) à la cellule neighbour. Son travail discret vaut donc exactement :
\\[
\mathbf{u}_P\cdot\mathbf{F}_f-\mathbf{u}_N\cdot\mathbf{F}_f
=(\mathbf{u}_P-\mathbf{u}_N)\cdot\mathbf{F}_f.
\\]
Pour une face frontière, seule la contribution owner est présente :
\\[
\mathcal{P}_f=\mathbf{u}_P\cdot\mathbf{F}_f.
\\]
Les sources volumiques déjà intégrées \\(\mathbf{S}_c\\) contribuent par \\(\sum_c\mathbf{u}_c\cdot\mathbf{S}_c\\). L'oracle calcule séparément convection, viscosité, viscosité transpose, pression, force volumique et terme temporel, puis leur somme physique. Il ne reconstitue pas les flux à partir du résidu cellulaire : les forces de face doivent provenir indépendamment des opérateurs acceptés.

## 4. Traitement des frontières

La convention topologique de `FaceOwnership` est stricte :
- face interne : \\(0\le N<n_{cells}\\), avec \\(N\ne P\\) ;
- face frontière : `FaceOwnership::BOUNDARY == -1` ;
- toute autre valeur négative est invalide et doit être rejetée ;
- tout owner hors domaine, neighbour hors domaine ou auto-connexion est rejeté.

Une face frontière utilise la vitesse de l'owner uniquement. Aucune valeur de vitesse fantôme ni condition limite supplémentaire n'est reconstruite par cet oracle. Le contrôle porte sur l'intégrité de la topologie et la cohérence du travail associé aux forces fournies.

## 5. Pièges connus

- Traiter arbitrairement tout \\(N<0\\) comme une frontière masque les sentinelles invalides ; seul `-1` est autorisé.
- Inverser l'orientation d'une face interne impose d'inverser simultanément la force orientée ; le travail doit rester invariant.
- Réutiliser le scatter du budget momentum comme oracle rendrait la comparaison circulaire. La reconstruction doit rester indépendante.
- Des forces non finies, une forme de tableau incorrecte, un owner/neighbour invalide ou un travail accumulé non fini doivent échouer explicitement.
- Un bilan de travail discret cohérent ne prouve ni l'exactitude des flux constitutifs, ni la conservation locale de chaque équation, ni la convergence physique d'une simulation.

## 6. Ordre de précision attendu et validation exigée

Aucun ordre spatial ou temporel n'est revendiqué : l'oracle teste une identité algébrique discrète, pas une approximation de l'EDP.

Validation exigée :
- tests unitaires des contributions internes, frontières et sources, avec reconstruction indépendante par scatter cellule ;
- contrôle d'invariance par inversion simultanée de l'orientation et de la force d'une face interne ;
- contrôles négatifs pour owner hors domaine, neighbour hors domaine, auto-connexion, sentinelle négatif autre que `-1`, dimensions incohérentes et données non finies ;
- exécution Release et DebugSanitizers sur le SHA exact ;
- aucune tolérance relâchée et aucun cas négatif supprimé pour obtenir le succès.

La qualification N9/S9 et les campagnes physiques restent des preuves distinctes ; elles ne sont pas remplacées par cet oracle.

## 7. Statut d'implémentation

- Fichier source : `src/cfdx/physics/momentum_vector_power_budget.h`.
- Tests exécutables : `tests/unit/test_conservation_boundedness.cpp`.
- Le test `raw_face_power_oracle_matches_independent_cell_scatter` vérifie les puissances agrégées par reconstruction indépendante.
- Le test `raw_face_power_oracle_fails_closed_on_invalid_topology_and_data` vérifie notamment le sentinelle `-2` et doit rester actif.
- Correctif courant : l'oracle rejette toute valeur négative différente de `FaceOwnership::BOUNDARY`.
- Écarts connus : le contrat est testé unitairement ; aucune qualification de la campagne longue N9/S9 ni de l'oracle Schur n'est déduite de ces tests.
- Gaps : terminer la CI exacte-head, puis exécuter explicitement le test long N9/S9 et examiner les diagnostics C-through-G indépendamment.

