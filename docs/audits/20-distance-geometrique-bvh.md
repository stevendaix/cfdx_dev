# FICHE MÉTHODE 20 — Distance géométrique aux parois et accélération BVH

**En-tête d'évidence**
- Scope : distances exactes point-triangle et intersection rayon-triangle utilisées par les méthodes de distance aux parois. Cette fiche qualifie le contrat géométrique et les optimisations de recherche ; elle ne valide pas à elle seule les EDP de Poisson, Eikonal, Hamilton–Jacobi ou advection-diffusion.
- Statut : APPROVED pour le contrat géométrique et algorithmique décrit aux sections 1–13 (revue PR #1172, 2026-11-04). Cette approbation ne vaut pas qualification de performance end-to-end, GCI des EDP, ni validation physique SST/NACA ; ces éléments restent explicitement ouverts.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` (campagnes) et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json` (exigences).
- Issue / PR liées : PR #979 ; backlog wall-distance #34.
- HEAD de qualification : à renseigner après la prochaine campagne CI exacte ; ne pas réutiliser les résultats d’un HEAD antérieur.
- Limitations connues : les temps N=128 disponibles sont des mesures de campagnes différentes et ne prouvent pas le gain causal du BVH. Le benchmark rayons synthétique isole le coût d'intersection, mais ne remplace pas une comparaison N=128 appariée sur le même runner.

## 1. Références sources

- T. Akenine-Möller, E. Haines, N. Hoffman, *Real-Time Rendering*, 4e éd., CRC Press, 2018, ISBN 978-1-138-62839-4 — hiérarchies de volumes englobants et intersection rayons/primitives.
- S. Woop, C. Benthin, I. Wald, « Watertight Ray/Triangle Intersection », *Journal of Computer Graphics Techniques*, 2(1), 2013 — référence pour les cas d’intersection proches des arêtes et sommets.
- C. Ericson, *Real-Time Collision Detection*, Morgan Kaufmann, 2004, ISBN 978-1-55860-732-3 — boîtes AABB et requêtes géométriques.
- Y. Zhao, « A fast sweeping method for Eikonal equations », *Mathematics of Computation* 74(250), 2005, 603–627 — DOI : 10.1090/S0025-5718-04-01678-3. — schémas fast-sweeping et résolution de l'équation Eikonal.

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
- Tests exécutables : `tests/unit/test_wall_distance_methods.cpp` couvre les cas géométriques et l'oracle ; `apps/verification/wall_distance_complex_benchmark.cpp` compare le coût BVH et brute-force sur les requêtes réelles de faces fluides/solides du cas complexe, vérifie l'équivalence et écrit `wall_distance_ray_benchmark.csv` en plus du CSV des méthodes. Le test `benchmark_wall_distance_bvh_rays` exécute uniquement la comparaison géométrique et archive `wall_distance_ray_benchmark_only.csv`. Le test `benchmark_wall_distance_complex` conserve le benchmark des neuf méthodes et archive aussi `wall_distance_ray_benchmark.csv`.
- Écarts connus code ↔ fiche : les tests de régression couvrent des cas choisis, et le benchmark rayons est représentatif d'un cas complexe synthétique mais borné ; il faut encore enregistrer les mesures sur runner contrôlé et qualifier le gain end-to-end des méthodes EDP. Le benchmark n'est pas une preuve d'accélération globale.
- Prochaines actions : (1) exécuter Release et DebugSanitizers sur le HEAD exact ; (2) lancer `ctest --test-dir build -R '^benchmark_wall_distance_bvh_rays$' -V` et archiver `wall_distance_ray_benchmark_only.csv` ; (3) lancer le benchmark N=128 apparié sur le même runner ; (4) comparer les échecs NACA à master, sans relâcher les critères.


## 8. Journal des modifications

- 2026-10-19 — PR #979 — fiche créée, tests d’équivalence géométrique et benchmark BVH/brute-force ajoutés.
- 2026-10-20 — PR #979 — benchmark rayons isolé dans un test CTest rapide ; sortie console corrigée pour produire un vrai saut de ligne ; qualification exacte du HEAD toujours requise.


## 9. Extraction des parois et propagation par point de mur (PR #1172)

La surface d'entrée est la triangulation des faces des patches 'WALL' du
maillage. L'extraction doit rejeter les indices invalides, les faces non
planaires, les polygones concaves ou auto-intersectants et les géométries
non finies. Les sommets collinéaires redondants sont autorisés ; les triangles
de fan de surface nulle sont omis. Une surface auto-intersectante ne peut pas
être « réparée » par une triangulation en éventail, car elle créerait des
triangles qui n'appartiennent pas à la frontière physique.

Pour chaque nœud fluide, la référence géométrique est
\[
d_i=\min_{T\in\mathcal T_w}\operatorname{dist}(x_i,T),
\]
où \(\mathcal T_w\) est l'ensemble des triangles de la surface d'entrée.
L'oracle exhaustif et le BVH doivent retourner la même distance à la
tolérance flottante déclarée. Les triangles dégénérés acceptés par l'API
géométrique restent traités par la primitive de distance correspondante ;
les indices et coordonnées invalides sont toujours rejetés.

'MESH_WAVE' propage désormais le point de surface le plus proche associé à
chaque seed, puis évalue la distance euclidienne entre le point fluide et ce
point de mur propagé. Il ne somme plus les longueurs des arêtes du graphe :
cette ancienne métrique Manhattan dépendait du chemin et n'était pas la
distance euclidienne. La méthode reste une approximation discrète sur le
graphe du maillage : les changements de point de mur, les composantes
déconnectées, les tie-breaks et les géométries non alignées peuvent produire
une erreur. Elle ne devient pas une méthode exacte générale parce qu'elle
est exacte sur un cube aligné.

## 10. Méthodes Eikonal et garde d'exactitude (PR #1172)

Le problème continu de distance satisfait l'équation Eikonal
\[
\|\nabla d\|_2=1,\qquad d|_{\Gamma_w}=0.
\]
La méthode 'EIKONAL' est un schéma discret de premier ordre. La méthode
'EIKONAL_SECOND_ORDER' utilise une mise à jour fast-sweeping d'ordre deux
lorsque les voisins amont requis sont disponibles et retombe sur une mise à
jour de premier ordre lorsque le stencil d'ordre deux n'est pas admissible.
Le résultat reste une approximation de l'EDP discrétisée ; l'ordre formel
local ne garantit ni l'ordre global sur géométrie polyédrique arbitraire, ni
l'exactitude vis-à-vis de la triangulation.

La garde 'compute_wall_distance_accuracy_checked' compare **chaque nœud
fluide** à la distance BVH exacte de la triangulation d'entrée. Pour chaque
échantillon, le critère est
\[
|d_h-d_{BVH}|\leq \tau_{abs}+\tau_{rel}d_{BVH}.
\]
Elle rapporte le nombre total d'échantillons, de violations, de corrections,
l'erreur relative L2 et l'erreur absolue L∞. Deux politiques seulement sont
permises : correction explicite par la valeur BVH exacte, marquée dans le nom
de méthode et le rapport ; ou rejet fail-closed si une violation subsiste.
Une convergence EDP ou un résidu faible ne constitue jamais une preuve
d'exactitude géométrique.

La garde certifie la borne sur les nœuds fluides effectivement testés, par
rapport à la surface triangulée fournie. Elle ne certifie pas la distance à
une CAO lisse non représentée, les gradients pariétaux, y+, le cisaillement,
ni l'ordre de convergence d'une solution RANS. Elle a un coût de vérification
supplémentaire proportionnel au nombre de nœuds fluides et aux requêtes BVH.

## 11. Contrat SST : normales signées et frontières

Pour les conditions SST qui nécessitent une distance normale à une face
murale, la projection signée doit être strictement positive dans la convention
géométrique du solveur. Le remplacement par 'abs(dot(n, dx))' est interdit :
il masquerait une normale inversée et convertirait une géométrie incohérente
en distance apparemment valide. La validation doit échouer explicitement
lorsque la projection est nulle, négative ou non finie. Ce contrat complète
la fiche 'docs/audits/09-turbulence-rans.md' ; il ne remplace pas la
qualification de la fermeture SST.

## 12. Validation et limites de la PR #1172

Les tests ajoutés/durcis couvrent notamment : cube analytique et cube tourné,
plan décalé, sphère, échelles géométriques de 1 à 1e-10 pour les rayons,
polygone concave et étoile auto-intersectante, point de mur propagé, erreur
de reconstruction Poisson par cancellation, Eikonal d'ordre deux, garde
exhaustive avec et sans correction BVH, absence de seeds, coordonnées NaN,
indices de triangles invalides et durée de vie du BVH.

La campagne unitaire/analytique est une vérification du code et du contrat
géométrique. Les campagnes EDP N=128 appariées, les mesures end-to-end BVH,
les études GCI et la qualification physique SST/NACA restent des tâches
séparées. Aucun statut 'PASS' de validation physique n'est attribué ici.

## 13. Journal des modifications

- 2026-11-04 — PR #1172 : contrats de propagation par point de mur, Eikonal
  d'ordre deux, garde d'exactitude et normales SST signées ajoutés.
