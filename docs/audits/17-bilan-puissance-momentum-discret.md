# FICHE MÉTHODE 17 — Validation indépendante de puissance momentum

**Statut :** DRAFT — identité de comptabilité discrète ; aucune preuve de stabilité ni qualification physique.

## Objectif et identité

Pour une face orientée owner vers neighbour, une force intégrée unique \(\mathbf F_f\) [N] donne :
- face interne : \((\mathbf U_P-\mathbf U_N)\cdot\mathbf F_f\) ;
- face frontière : \(\mathbf U_P\cdot\mathbf F_f\).

Les sources cellulaires intégrées \(\mathbf S_c\) contribuent par \(\mathbf U_c\cdot\mathbf S_c\). Avec \(\mathbf U\) en [m/s], la puissance est en [W].

La validation compare deux chemins :
1. oracle face brut utilisant le saut de vitesse ;
2. scatter signé des forces de face vers les cellules, puis produit scalaire cellulaire.

Les puissances internes et frontières sont rapportées séparément. La comparaison est effectuée opérateur par opérateur afin qu'une erreur de signe dans un canal ne soit pas masquée par compensation entre canaux.

## Contrat et unités

- Forces de face et sources cellulaires déjà intégrées [N] ; aucun volume supplémentaire.
- Vitesse cellulaire [m/s] ; puissance [W].
- Chaque face interne apparaît une fois ; le voisin frontière vaut exclusivement `FaceOwnership::BOUNDARY` (-1).
- Convection, diffusion moléculaire, traction visqueuse transposée, pression, force extérieure et accumulation sont des canaux distincts.
- Aucune hypothèse de signe dissipatif n'est imposée à la pression, convection, travail frontière ou sources.

## Tests ajoutés

- Maillage deux cellules, faces internes et frontières.
- Comparaison oracle/scatter indépendant pour convection, viscous, transpose-viscous, pression, force et accumulation séparément.
- Décomposition interne/frontière pour chaque opérateur face.
- Invariance au renversement de l'orientation d'une face interne si la force est inversée simultanément.
- Propriété de signe sur un graphe à conductances positives : le travail du Laplacien discret avec frontière de Dirichlet homogène est égal à une somme de termes \(g_f\lVert\Delta U_f\rVert^2\) et doit être non négatif.
- Rejet des tailles incohérentes, propriétaire/voisin invalide, sentinelle négative inconnue, NaN et overflow.

## Limites scientifiques

L'oracle est un vérificateur de contrat pour des forces intégrées fournies par son appelant. La capture des forces brutes directement depuis les opérateurs de production à l'état accepté n'est pas encore raccordée. Il ne prouve donc pas encore que les opérateurs réellement assemblés satisfont l'identité sur les cas de production.

Le test de positivité ajouté ne couvre que le Laplacien de graphe à conductances positives et une frontière de Dirichlet homogène ; il ne prouve pas la positivité du flux diffusif CFDX complet sur maillage non orthogonal, avec viscosité variable et correction différée. Il ne prouve ni stabilité temporelle/couplage pression-vitesse, ni réalisabilité d'une contrainte de Reynolds, ni précision d'une fermeture turbulente.

Étape suivante : brancher les forces réelles de production, réconcilier chaque canal avec les budgets cellulaires et les lignes CSR, puis mener des études de convergence maillage/temps et des validations a-posteriori sur cas de référence indépendants de la calibration.

## Relation à Lean

Les contrats Lean owner/neighbour et source fournissent l'identité abstraite. Les tests C++ vérifient des maillages runtime déterministes. La correspondance générale entre topologie Lean, `Mesh::ownership()` et sorties des opérateurs de production reste à formaliser.
