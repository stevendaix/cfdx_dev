# FICHE MÉTHODE 11 — Solveur d’énergie et convergence linéaire

**En-tête d’évidence**
- Scope : équation scalaire d’énergie en volumes finis, tolérance du solveur Krylov, acceptation de l’état thermique et bilans énergétiques indépendants.
- Statut : PROPOSED — validation T2 en cours ; cette fiche décrit le contrat numérique et ne vaut pas qualification du cas chauffé.
- Issues liées : #880, #885 ; PR #1059.
- Dernier commit vérifié : à renseigner après validation CI et campagne T2.
- Limitations : solveur d’énergie de fluide à flux massique gelé ; la stagnation du solveur d’écoulement est hors de ce contrat et doit être qualifiée séparément.

## 1. Références sources
- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2e éd., SIAM, 2003 — ISBN 978-0-89871-534-7.
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — ISBN 978-0-89116-522-4.

## 2. Formulation continue
Pour une température (T), une conductivité (k), une capacité thermique (c_p), une densité (ho), un flux massique (dot m_f) et une source volumique (S_T), le solveur résout le bilan stationnaire ou transitoire :
[
\nabla\cdot(\rho c_p\mathbf{u}T)=\nabla\cdot(k\nabla T)+S_T.
]
Pour un flux thermique imposé positif vers le domaine, la condition de Neumann doit contribuer au bilan global avec ce même signe.

## 3. Discrétisation Volume Fini
L’intégration par cellule produit un système creux (A T=b). Pour un flux massique, une diffusivité et une source gelés, ce système est linéaire. La solution n’est acceptée que si le statut Krylov est convergé et si le résidu relatif ainsi que l’erreur arrière de l’état accepté satisfont leurs critères configurés.

La tolérance relative passée au solveur Krylov est mise à l’échelle par la norme du second membre dans l’implémentation du solveur. Elle ne doit donc pas être multipliée une seconde fois par une échelle de tolérance locale sans contrat explicite.

## 4. Traitement des frontières
- Dirichlet : température imposée dans l’opérateur et le second membre.
- Neumann : flux thermique imposé dans le second membre, avec convention de signe explicitement vérifiée par le bilan indépendant.
- Paroi adiabatique : flux diffusif nul.
- Les bilans global et local sont reconstruits indépendamment de l’équation assemblée et restent des critères bloquants ; le résidu linéaire ne remplace jamais ces bilans.

## 5. Pièges connus
- Confondre le statut du solveur linéaire avec le verdict complet du calcul d’énergie.
- Appliquer deux fois une tolérance relative, ce qui peut imposer un seuil artificiellement strict et rendre le solveur incapable d’accepter un prédicteur.
- Répéter un système linéaire inchangé après un statut Krylov d’échec alors que le prédicteur rejeté laisse l’état accepté inchangé.
- Déclarer la convergence sur le seul résidu algébrique sans vérifier le bilan thermique indépendant.
- Confondre une erreur de bilan thermique avec la stagnation du solveur d’écoulement.

## 6. Ordre de précision attendu et validation exigée
- Le système linéaire doit satisfaire le statut Krylov, le résidu relatif et l’erreur arrière de l’état accepté.
- Le bilan énergétique global et les bilans locaux indépendants doivent satisfaire leurs tolérances existantes ; aucune tolérance ne peut être augmentée pour faire passer la qualification.
- Pour T2, exiger sur chaque maillage une variation de température cohérente avec la puissance de paroi imposée, des statuts/résidus linéaires interprétables, la fermeture des bilans indépendants et la qualification de la grandeur thermique de sortie.
- Un échec linéaire doit conserver son diagnostic et terminer le solveur d’énergie sans produire des centaines d’itérations identiques.

## 7. Statut d’implémentation
- Sources : `src/cfdx/physics/energy_solver.h`, `src/cfdx/physics/production_energy_stage.h`, `apps/production_solver/main.cpp`.
- Tests et preuves : campagne `T2 heated channel`, historique thermique JSON, CI Release/Debug et tests de validation de l’étape énergie.
- Statut actuel : correction de tolérance et observabilité poussées dans PR #1059 ; la campagne T2 et la CI doivent encore confirmer le résultat.
- Écarts code ↔ fiche : qualification exacte-HEAD et revue du cas T2 encore en attente.

## Journal
- 2026-10-23 — première version pour documenter le contrat de tolérance Krylov, l’arrêt sur échec linéaire et les gates de conservation ; qualification en attente.
