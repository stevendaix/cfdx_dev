# FICHE MÉTHODE 20 — Distance géométrique aux parois et accélération BVH

**En-tête d'évidence**
- Scope : distances exactes point-triangle et intersection rayon-triangle utilisées par les méthodes de distance aux parois. Cette fiche qualifie le contrat géométrique et les optimisations de recherche ; elle ne valide pas à elle seule les EDP de Poisson, Eikonal, Hamilton–Jacobi ou advection-diffusion.
- Statut : DRAFT (2026-10-19) — tests d'équivalence géométrique ajoutés ; benchmark rayons BVH/brute-force intégré au cas complexe, mesures et qualification EDP encore requises.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` (campagnes) et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json` (exigences).
- Issue / PR liées : PR #979 ; backlog wall-distance #34.
- Dernier HEAD analysé avant l'ajout du benchmark : `d7068463ae27d2c5db0fafbd9043986f0b4bcbb4`. Les contrôles CI doivent être relus sur le nouveau HEAD de la branche.
- Limitations connues : les temps N=128 disponibles sont des mesures de campagnes différentes et ne prouvent pas le gain causal du BVH. Le benchmark rayons synthétique isole le coût d'intersection, mais ne remplace pas une comparaison N=128 appariée sur le même runner.

## 1. Références sources

- T. Akenine-Möller, E. Haines, N. Hoffman, *Real-Time Rendering*, 4e éd., CRC Press, 2018, ISBN 978-1-138-62839-4 — hiérarchies de volumes englobants et intersection rayons/primitives.
- S. Woop, C. Benthin, I. Wald, « Watertight Ray/Triangle Intersection », *Journal of Computer Graphics Techniques*, 2(1), 2013 — référence pour les cas d’intersection proches des arêtes et sommets.
- C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann, 2004, ISBN 978-1-55860-732-3 — boîtes AABB et requêtes géométriques.

## 2. Formulation continue

Pour un point fluide \(\mathbf{x}\) et une surface triangulée \(\Gamma_w\), la distance géométrique non signée est
\[
d(\mathbf{x})=\min_{\mathbf{y}\in\Gamma_w}\|\mathbf{x}-\mathbf{y}\|_2.
\]
Pour une face coupée, la fermeture géométrique utilise l'intersection du rayon \(\mathbf{x}+t\hat{\mathbf{r}}\) avec les triangles de la paroi, où \(\|\hat{\mathbf{r}}\|_2=1\) et \(0\le t\le t_{max}\). Le résultat est le plus petit paramètre admissible \(t\), ou \(+\infty\) en l'absence d'intersection.

## 3. Discrétisation Volume Fini

La surface est indexée par un BVH binaire d'AABB. Chaque nœud borne tous les triangles de son sous-arbre ; les feuilles référencent un intervalle de triangles dans le tableau d'indices. Une requête rayon :
1. normalise la direction, en rejetant une direction nulle ou une distance maximale non positive ;
2. teste le rayon contre les AABB, avec une intersection par intervalles sur les trois axes ;
3. visite les feuilles candidates et applique l'intersection rayon-triangle ;
4. retourne le plus petit \(t\) valide dans la limite demandée.

La traversée est une accélération de recherche uniquement. Le test rayon-triangle, les tolérances barycentriques et la définition du premier impact doivent rester équivalents à la référence exhaustive. Le BVH ne doit jamais éliminer un triangle qui peut fournir un impact admissible.

## 4. Traitement des frontières

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
- Validation de performance : (1) exécuter `benchmark_wall_distance_complex` pour mesurer BVH vs oracle exhaustif sur les mêmes requêtes de faces coupées du cas complexe, avec résultats identiques ; (2) répéter le benchmark N=128 des neuf méthodes sur le même runner et le même HEAD de base, en conservant les CSV et les métriques temps mur, itérations, résidu, raison d'arrêt, L1/L2/L∞ et erreur près des parois.
- Aucun seuil de convergence, tolérance ou budget ne peut être relâché pour qualifier l'accélération.
- GCI des méthodes EDP : reste soumis aux fiches et campagnes propres à chaque méthode ; cette fiche géométrique ne remplace pas ces exigences.

## 7. Statut d'implémentation CFDX

- Fichier source : `src/cfdx/physics/wall_distance.h`, classe `WallDistanceBvh`.
- Tests exécutables : `tests/unit/test_wall_distance_methods.cpp` couvre les cas géométriques et l'oracle ; `apps/verification/wall_distance_complex_benchmark.cpp` compare le coût BVH et brute-force sur les requêtes réelles de faces fluides/solides du cas complexe, vérifie l'équivalence et écrit `wall_distance_ray_benchmark.csv` en plus du CSV des méthodes. Le test reste un benchmark de validation et ses mesures doivent être archivées.
- Écarts connus code ↔ fiche : les tests de régression couvrent des cas choisis, et le benchmark rayons est représentatif d'un cas complexe synthétique mais borné ; il faut encore enregistrer les mesures sur runner contrôlé et qualifier le gain end-to-end des méthodes EDP. Le benchmark n'est pas une preuve d'accélération globale.
- Prochaines actions : (1) vérifier Release et DebugSanitizers sur le nouveau HEAD ; (2) lancer `ctest --test-dir build -R '^benchmark_wall_distance_complex (3) lancer un A/B N=128 apparié, mêmes runner, maillage et options ; (4) traiter séparément les échecs NACA préexistants et le test N9 proche de 1800 s, sans relâcher leurs gates.

## 8. Journal des modifications

- 2026-10-19 — PR #979 — création de la fiche, tests rayons comparés à l'oracle exhaustif et intégration d'une mesure BVH vs brute-force au benchmark géométrique complexe ; les mesures contrôlées et le benchmark EDP apparié restent à exécuter.
 -V` et archiver `wall_distance_complex_benchmark.csv` ainsi que `wall_distance_ray_benchmark.csv` ; (3) lancer un A/B N=128 apparié, mêmes runner, maillage et options ; (4) traiter séparément les échecs NACA préexistants et le test N9 proche de 1800 s, sans relâcher leurs gates.

## 8. Journal des modifications

- 2026-10-19 — PR #979 — création de la fiche, tests rayons comparés à l'oracle exhaustif et ajout d'un benchmark reproductible BVH vs brute-force ; les mesures contrôlées et le benchmark EDP apparié restent à exécuter.
