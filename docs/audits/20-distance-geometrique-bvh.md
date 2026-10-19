# FICHE MÉTHODE 20 — Distance géométrique aux parois et accélération BVH

**En-tête d'évidence**
- Scope : distances exactes point-triangle et intersection rayon-triangle utilisées par les méthodes de distance aux parois. Cette fiche qualifie le contrat géométrique et les optimisations de recherche ; elle ne valide pas à elle seule les EDP de Poisson, Eikonal, Hamilton–Jacobi ou advection-diffusion.
- Statut : DRAFT (2026-10-19) — contrat géométrique documenté ; validation numérique de la traversée BVH en cours.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` (campagnes) et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json` (exigences).
- Issue / PR liées : PR #979 ; backlog wall-distance #34.
- Dernier commit vérifié : `6fed1d42fc69bd5236c39fc1b48dcbe5e112f767` (branche PR #979).
- Limitations connues : les résultats N=128 antérieurs sont une baseline, pas une preuve du gain de cette modification. Le benchmark complet et la comparaison à la référence exhaustive restent requis.

## 1. Références sources

- T. Akenine-Möller, E. Haines, N. Hoffman, *Real-Time Rendering*, 4e éd., CRC Press, 2018, ISBN 978-1-138-62839-4 — hiérarchies de volumes englobants et intersection rayons/primitives.
- S. Woop, C. Benthin, I. Wald, « Watertight Ray/Triangle Intersection », *Journal of Computer Graphics Techniques*, 2(1), 2013 — intersection robuste rayon-triangle, DOI: 10.1111/cgf.12462.
- C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann, 2004, ISBN 978-1-55860-732-3 — boîtes AABB et requêtes géométriques.

## 2. Formulation continue

Pour un point fluide \(\mathbf{x}\) et une surface triangulée \(\Gamma_w\), la distance géométrique non signée est
\[
d(\mathbf{x})=\min_{\mathbf{y}\in\Gamma_w}\|\mathbf{x}-\mathbf{y}\|_2.
\]
Pour une face coupée, la fermeture géométrique utilise l'intersection du rayon \(\mathbf{x}+t\hat{\mathbf{r}}\) avec les triangles de la paroi, où \(\|\hat{\mathbf{r}}\|_2=1\) et \(0\le t\le t_{max}\). Le résultat est le plus petit paramètre admissible \(t\), ou \(+\infty\) en l'absence d'intersection.

## 3. Discrétisation / algorithme

La surface est indexée par un BVH binaire d'AABB. Chaque nœud borne tous les triangles de son sous-arbre ; les feuilles référencent un intervalle de triangles dans le tableau d'indices. Une requête rayon :
1. normalise la direction, en rejetant une direction nulle ou une distance maximale non positive ;
2. teste le rayon contre les AABB, avec une intersection par intervalles sur les trois axes ;
3. visite les feuilles candidates et applique l'intersection rayon-triangle ;
4. retourne le plus petit \(t\) valide dans la limite demandée.

La traversée est une accélération de recherche uniquement. Le test rayon-triangle, les tolérances barycentriques et la définition du premier impact doivent rester équivalents à la référence exhaustive. Le BVH ne doit jamais éliminer un triangle qui peut fournir un impact admissible.

## 4. Conditions aux limites et domaine de validité

- Le rayon est limité à une distance finie positive fournie par l'appelant.
- Les impacts à l'intérieur des tolérances existantes autour de l'origine et de la limite restent traités comme dans la référence historique.
- Une direction nulle, une limite non positive ou une absence d'impact renvoie une distance infinie, jamais une distance artificielle.
- Le BVH est construit à partir de la même surface triangulée que celle interrogée ; la surface ne doit pas être mutée pendant les requêtes.

## 5. Pièges connus (Gotchas)

- Faux négatifs dus à des tests AABB trop stricts, à des directions quasi parallèles ou aux impacts sur les arêtes/sommets.
- Modification involontaire du premier impact en visitant les triangles dans un ordre différent ; seule la valeur minimale valide doit déterminer le résultat.
- Confusion entre distance au point de surface le plus proche et intersection le long d'une normale de face coupée : ces quantités ne sont pas interchangeables aux arêtes et aux coins.
- Une amélioration du temps par requête ne prouve ni une amélioration du nombre d'itérations ni une meilleure précision de l'EDP.
- Les temps mesurés doivent être comparés sur le même HEAD de base, le même runner, le même maillage et les mêmes options.

## 6. Précision attendue et validation exigée

- Exactitude géométrique : écart à la recherche exhaustive inférieur à la tolérance numérique pré-déclarée ; pour les géométries analytiques simples, comparaison aux intersections exactes.
- Couverture minimale : impact intérieur à un triangle, absence d'impact, clipping par distance maximale, direction non unitaire, impacts sur arêtes/sommets, rayons parallèles aux axes et directions quasi parallèles.
- Validation de performance : benchmark N=128 des neuf méthodes, CSV individuels et mêmes conditions d'exécution. Rapporter temps mur, nombre d'itérations, résidu final, raison d'arrêt, L1/L2/L∞ et erreur près des parois.
- Aucun seuil de convergence, tolérance ou budget ne peut être relâché pour qualifier l'accélération.
- GCI des méthodes EDP : reste soumis aux fiches et campagnes propres à chaque méthode ; cette fiche géométrique ne remplace pas ces exigences.

## 7. Statut d'implémentation CFDX

- Fichier source : `src/cfdx/physics/wall_distance.h`, classe `WallDistanceBvh`.
- Tests exécutables : `tests/unit/test_wall_distance_methods.cpp` ; tests rayons ciblés et équivalence exhaustive à compléter/vérifier dans PR #979.
- Écarts connus code ↔ fiche : la branche #979 introduit une traversée AABB/BVH dans `ray_distance()`, mais son équivalence numérique et son gain de performance ne sont pas encore qualifiés.
- Prochaines actions : ajouter les tests de régression, passer le gate documentaire, vérifier Release et DebugSanitizers, puis exécuter le benchmark de référence.

## 8. Journal des modifications

- 2026-10-19 — PR #979 — création de la fiche pour cadrer l'optimisation BVH et ses exigences de qualification.
