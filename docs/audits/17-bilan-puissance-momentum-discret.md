# FICHE MÉTHODE 17 — Bilan de puissance du momentum en volumes finis

**En-tête d'évidence**
- Scope : identité discrète de puissance mécanique du résidu de quantité de mouvement incompressible CFDX, formulée à partir des flux owner/neighbour et des sources cellulaires intégrées. L'implémentation C++ de #1224 n'en réalise pour l'instant que l'agrégation vectorielle des budgets cellulaires déjà scatterés ; elle ne reconstruit pas encore indépendamment les flux bruts. Couvre la distinction résidu physique / résidu de l'opérateur borné et la puissance temporelle. Ne couvre pas une preuve de stabilité du solveur couplé, une fermeture de Reynolds, ni la qualification physique d'un modèle turbulent.
- Statut : DRAFT (2026-11-28) — contrat d'audit en construction ; PR #1224 non fusionnée (remplacement de #1223).
- Registre autoritaire : issue #1216, fondations Lean #1197, contrat formel #1222, implémentation C++ #1224.
- Sources de production : `src/cfdx/physics/momentum_face_budget.h`, `src/cfdx/physics/nonlinear_face_residual_certificate.h`, `src/cfdx/physics/steady_incompressible_solver.h`.
- Tests concernés : `tests/unit/test_conservation_boundedness.cpp`, `tests/audit/test_sst_integration_gaps.cpp`.
- Limites de promotion : la correspondance C++/Lean reste bornée aux contrats abstraits et aux fixtures déterministes ; aucune inégalité d'énergie ni stabilité du solveur complet n'est établie.

## 1. Références sources

- H. K. Versteeg, W. Malalasekera, *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*, 2e éd., Pearson, 2007 — ISBN 978-0-13-127498-3 (flux intégrés, volumes de contrôle, conservation).
- J. H. Ferziger, M. Perić, *Computational Methods for Fluid Dynamics*, 3e éd., Springer, 2002 — DOI : 10.1007/978-3-642-56026-2 (équations discrètes, volumes finis et bilans).
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN 978-3-319-16873-9 (assemblage conservatif et flux de face).
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI : 10.1201/9781482234213 (formulation intégrée et traitement des sources).

## 2. Formulation continue

Pour un fluide incompressible, le résidu de quantité de mouvement peut être écrit, avec une convention de signes fixée par l'assemblage CFDX, sous la forme :

\[
\frac{\partial(\rho\mathbf U)}{\partial t}
+ \nabla\cdot(\rho\mathbf U\otimes\mathbf U)
- \nabla\cdot\boldsymbol\tau_v
+ \nabla p - \rho\mathbf f = \mathbf 0.
\]

Le travail mécanique s'obtient par produit scalaire avec la vitesse. Les termes de pression, de frontière et de force extérieure ne sont pas nécessairement dissipatifs. Une contrainte symétrique ne suffit pas à garantir une puissance non positive ; le signe dépend de la loi constitutive, de la convention de traction et de l'état cinématique.

La fiche distingue trois niveaux :
- identité de comptabilité (réarrangement exact de sommes finies) ;
- propriété de signe ou inégalité énergétique (hypothèses supplémentaires) ;
- stabilité/convergence du solveur non linéaire et du couplage pression-vitesse (problème distinct).

## 3. Discrétisation Volume Fini

Pour une face orientée owner vers neighbour, soit \(\mathbf F_f\) la force intégrée unique partagée par les cellules incidentes. La contribution intégrée de face à la cellule \(c\) est signée par l'incidence owner/neighbour. Pour une source volumique \(\mathbf s_c\), le résidu intégré s'écrit :

\[
\mathbf R_c =
\sum_{f\in\partial c} s_{cf}\,\mathbf F_f
+ V_c\,\mathbf s_c,
\qquad
s_{cf}=+1\;\text{(owner)},\quad -1\;\text{(neighbour)}.
\]

Le produit global par la vitesse donne l'identité :

\[
\sum_c \mathbf U_c\cdot\mathbf R_c
=
\sum_{f\in\mathcal F_i}
(\mathbf U_P-\mathbf U_N)\cdot\mathbf F_f
+
\sum_{f\in\mathcal F_b}
\mathbf U_P\cdot\mathbf F_f
+
\sum_c V_c\,\mathbf U_c\cdot\mathbf s_c.
\]

Les flux internes apparaissent une fois et se réduisent au saut de vitesse ; les faces frontières restent explicites. Le contrat Lean #1219 formalise le travail vectoriel owner/neighbour, et #1222 étend l'identité au terme source cellulaire.

Dans le code, `diagnose_momentum_face_budget` fournit déjà les contributions intégrées par cellule et par composante : convection reconstruite à partir du flux massique conservatif autoritaire, transport visqueux (y compris la correction non orthogonale reconstruite), traction visqueuse transposée, gradient de pression, force volumique et accumulation temporelle. La puissance s'obtient en appariant ces forces intégrées [N] à la vitesse [m/s], sans remultiplier par le volume : le résultat est en watts [W].

La PR #1224 ajoute uniquement l'agrégation des trois budgets scalaires acceptés. Cette étape ne reconstitue pas les flux bruts de faces : elle agrège les contributions déjà scatterées aux cellules. Une comparaison indépendante face/source nécessitera de conserver ou de reconstruire les flux de chaque opérateur.

## 4. Traitement des frontières

L'orientation est celle de `Mesh::ownership()` : la face est stockée une fois, le propriétaire reçoit la contribution positive et le voisin interne la contribution négative. Une face frontière ne possède pas de voisin ; son travail est \(\mathbf U_P\cdot\mathbf F_f\) et ne doit pas être annulé artificiellement.

Les conditions de vitesse influencent les flux visqueux et convectifs reconstruits. Les sources de pression et de force volumique sont cellulaires dans le budget existant. Toute extension à une traction de Reynolds devra définir explicitement :
- traction intégrée et convention de signe ;
- orientation owner/neighbour et contribution frontière ;
- convention de contrainte dynamique [Pa] ou cinématique [m²/s²] ;
- conversion éventuelle par la densité ;
- cohérence entre puissance de contrainte, production de turbulence et transport de \(k\).

Aucune de ces obligations ne peut être remplacée par une simple symétrie tensorielle.

## 5. Pièges connus

- **Ne pas appliquer le volume deux fois.** Les champs `MomentumFaceCellBudget` contiennent déjà des contributions intégrées [N] ; seules les sources volumétriques ont été multipliées par \(V_c\) lors de leur construction.
- **Ne pas confondre résidu physique et résidu borné.** `physical_residual` exclut `bounded_divergence`; `bounded_operator_residual` l'inclut. Le travail de cette correction reste un diagnostic algorithmique distinct.
- **Ne pas interpréter `collocated_difference` comme une puissance Rhie–Chow pure.** Le diagnostic compare le flux conservatif autoritaire à un flux collocatif alternatif et inclut aussi le splitting de relaxation de pression.
- **Ne pas remplacer la réconciliation CSR.** `reconcile_momentum_budget_with_matrix` compare les valeurs signées cellule par cellule entre \(A\mathbf U-\mathbf b\) et le résidu borné reconstruit. L'agrégation vectorielle conserve son statut mais ne la refait pas.
- **Ne pas imposer un signe dissipatif global.** Pression, travail des frontières, force extérieure et backscatter physique peuvent donner une puissance positive.
- **Ne pas appeler l'identité Lean une preuve de stabilité.** La topologie abstraite ne prouve pas que tous les maillages runtime satisfont les hypothèses, ni que l'algorithme couplé converge.
- **Ne pas assimiler la viscosité transposée à la contrainte de Reynolds.** La fixture existante valide un opérateur visqueux ; elle n'est pas une fermeture turbulente.

## 6. Ordre de précision attendu et validation exigée

Cette fiche spécifie une identité de comptabilité, pas un ordre de précision spatial. La validation logicielle doit couvrir :
- les trois composantes de vitesse et les trois budgets issus du même état accepté ;
- contributions non nulles de convection, diffusion visqueuse, traction transposée, pression, force volumique et accumulation ;
- contribution bornée séparée du résidu physique ;
- données non finies et tailles incohérentes rejetées ;
- réconciliation CSR signée conservée indépendamment pour chaque composante ;
- tests Release et DebugSanitizers sur le HEAD exact.

La validation physique d'une fermeture exige ensuite des campagnes a-posteriori indépendantes : DNS de canal Moser–Kim–Mansour, couche limite, NACA0012, Driver–Seegmiller BFS, NASA hump et periodic hill selon le domaine du modèle. Les profils moyens, contraintes de Reynolds, frottement pariétal, séparation/réattachement, pertes de pression, conservation, convergence et coût doivent être comparés sur des géométries/Reynolds tenus à l'écart de la calibration.

Aucun seuil ne doit être relâché et aucun test de qualification ne doit être désactivé pour faire passer l'agrégation.

## 7. Statut d'implémentation CFDX

- Fichiers source concernés : `src/cfdx/physics/momentum_face_budget.h`, `src/cfdx/physics/nonlinear_face_residual_certificate.h`, `src/cfdx/physics/steady_incompressible_solver.h`, `src/cfdx/physics/momentum_vector_power_budget.h`.
- Contrats formels : `CFDX.FiniteVolume.MeshVectorWork`, `CFDX.FiniteVolume.BackwardEulerEnergy`, `CFDX.FiniteVolume.MomentumPower` (PR #1219 et #1222).
- Tests : `tests/unit/test_conservation_boundedness.cpp` pour la comptabilité et l'agrégation ; `tests/audit/test_sst_integration_gaps.cpp` pour la fixture de traction visqueuse transposée.
- État : #1219 a passé les contrôles exact-head consignés dans l'issue #1216 ; #1222 a passé ses sept workflows exact-head au SHA `2da1c351b4c8ca50aef6fb75893561a1cacc1e95`. PR #1224 ajoute l'agrégateur C++ et sa fixture arithmétique ; Release et les autres contrôles sont verts au HEAD initial `33a459b3e50258d960ec550427fbe99fbd1452f7`, DebugSanitizers encore en cours au moment de la mise à jour. Une correction supplémentaire contre le débordement des sommes de reconstruction et son test sont en cours sur la branche de la PR.
- Écarts connus code ↔ contrat : agrégation initiale à partir des budgets cellulaires scatterés ; reconstruction indépendante de puissance brute par face/source non encore livrée ; correspondance générale entre CSR runtime et topologie Lean non démontrée.
- Gaps et prochaines actions : valider #1224 ; tester l'agrégateur sur un état réellement accepté issu du solveur ; ajouter une reconstruction indépendante face/source ; formaliser séparément le BDF2 à pas variable pour le type d'opérateur réellement applicable ; introduire ensuite seulement un contrat de traction de Reynolds distinct. Aucune fermeture turbulente n'est validée ou activée.

## 8. Journal des modifications de la fiche

- 2026-11-28 — correction du suivi : #1224 remplace #1223 ; précision que les contrôles verts du SHA initial ne couvrent pas les commits ultérieurs. Audit supplémentaire : rejet explicite des débordements pendant la reconstruction des résidus et test associé ajoutés sur la branche #1224.
- 2026-11-28 — création dans la PR #1224 (issue initialement portée par #1223) après déclenchement de la règle d'or 1 ; fiche DRAFT, à confronter au code et aux résultats exact-head avant toute approbation.
