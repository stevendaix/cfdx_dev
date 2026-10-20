# FICHE MÉTHODE 20 — Cache des offsets de faces coupées (wall-distance PDE)

**En-tête d'évidence**
- Scope : optimisation de calcul des offsets géométriques `delta` pour les stencils wall-distance Hamilton–Jacobi et hybride Poisson–Eikonal ; ne redéfinit ni la géométrie ni les équations discrètes.
- Statut : APPROUVEE (2026-10-20) — qualification exact-HEAD terminée (A/B pinned run 39397470405, campagne validation-total 39374525195, profiling variantes run 39463278884) ; feu vert merge Steven 2026-10-20 avec la régression HJ acceptée comme compromis de production documenté (voir §8).
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json`.
- Issue / PR liées : PR #989 ; optimisation BVH distincte PR #979.
- Dernier commit vérifié : `54fd341255ae6166114d098ef0424eec465b34cc` (qualification complète à ce SHA) ; le commit docs-only de la présente fiche est superposé sans toucher `src/` ni `tests/`.
- Limitations connues : le gain dépend du nombre de nœuds adjacents aux faces coupées, du coût BVH et du coût de construction du cache (~23 ms, 17,5 Mio à N=128) ; le surcoût Hamilton–Jacobi est dépendant du micro-matériel (voir §8).

## 1. Références sources

- R. A. Shortley, R. Weller, « The numerical solution of Laplace's equation », *Journal of Applied Physics*, 9 (1938), 334–348 — DOI: 10.1063/1.1714155.
- F. Moukalled, L. Mangani, M. Darwish (2016), *The Finite Volume Method in Computational Fluid Dynamics*, ISBN 978-3-319-16873-9.

## 2. Formulation continue

La méthode de distance à la paroi est déterminée par une EDP de transport ou d'Eikonal. Le cache est une optimisation d'implémentation : pour une géométrie, un masque solide et une grille immuables durant un solve, chaque offset géométrique de face coupée est une constante du solve.

Pour une face coupée reliant un nœud fluide (P) à un nœud solide (Q), l'offset est (delta_{P,Q}), borné par le plancher et l'espacement local comme dans `poisson_wall_offset()`. Le cache ne change pas cette définition.

## 3. Discrétisation Volume Fini

Les stencils utilisent six slots par nœud fluide : (x^-), (x^+), (y^-), (y^+), (z^-), (z^+). Pour un voisin fluide, la distance reste (h_a). Pour un voisin solide, la distance est (delta_{P,Q}), calculée une fois par le BVH et conservée jusqu'à la fin du solve.

Le cache est creux : un index par nœud identifie les nœuds ayant au moins une face coupée ; seuls ces nœuds possèdent les six distances. Les nœuds sans face coupée utilisent directement (h_x,h_y,h_z). L'ordre des opérations dans les stencils, les coefficients et les valeurs de frontière ne doivent pas changer du fait de la mise en cache.

Le solveur de Poisson hybride doit employer les mêmes offsets pour l'opérateur et le préconditionneur diagonal. Hamilton–Jacobi et hybride doivent également employer le même stencil pour la mise à jour et l'évaluation du résidu. Le gradient Godunov d'Hamilton–Jacobi utilise lui aussi les offsets géométriques aux faces solides ; il doit lire le cache sans modifier le stencil ni les valeurs de delta.

## 4. Traitement des frontières

- Face fluide-solide : valeur de Dirichlet nulle côté paroi et offset géométrique (delta_{P,Q}) inchangé.
- Face fluide-fluide : distance cartésienne (h_a).
- Face au bord extérieur de la grille : la condition de Neumann homogène existante reste inchangée ; le slot n'est pas lu si le voisin n'existe pas.
- Nœud solide : aucun offset de départ n'est construit ni consommé.
- La construction du cache doit utiliser la même voie de rayon BVH, le même clamp et le même fallback normal que le calcul sans cache.

## 5. Pièges connus

- Ne pas confondre accélération du BVH (#979) et réutilisation des offsets (#989) : les gains doivent être mesurés séparément.
- Un cache dense de six doubles par nœud peut coûter 48 octets par nœud ; le stockage creux limite l'empreinte mais ajoute un index par nœud. Mesurer l'empreinte allouée réelle.
- Le cache est valide uniquement si la surface, le masque solide, la topologie de grille et les espacements ne changent pas durant le solve. Toute mutation impose une reconstruction.
- Un résidu évalué avec (h_a) alors que la mise à jour utilise (delta_{P,Q}) n'est pas le résidu du même opérateur ; les deux chemins doivent être cohérents.
- Une amélioration du temps ne suffit pas : les champs, résidus, critères d'arrêt, nombres d'itérations et erreurs doivent rester qualifiés. Aucun seuil ni tolérance ne doit être relâché pour faire passer le benchmark.

## 6. Ordre de précision attendu et validation exigée

- Équivalence locale : offsets cachés et calcul direct doivent être identiques en double précision pour ray hits, fallback, faces fractionnaires, frontières de grille et voisins fluides.
- Équivalence opérateur : Laplacien de Shortley–Weller, opérateur de Poisson et diagonale doivent correspondre au chemin non caché à une tolérance absolue pré-déclarée de 10⁻¹³ sur les tests déterministes. Le gradient Godunov d'Hamilton–Jacobi est couvert par un test d'équivalence cache/direct, à la tolérance absolue pré-déclarée de 10⁻¹³.
- Régression : Release et DebugSanitizers, test `test_wall_distance_methods`, puis campagnes wall-distance existantes.
- Performance : comparer le même benchmark contrôlé N=128, même runner et mêmes budgets, base pré-cache contre HEAD cache ; conserver les neuf CSV, temps par méthode, itérations, résidus, stopping reason, erreurs L1/L2/L∞/near-wall et empreinte mémoire. Présenter médiane ou répétitions si le bruit du runner est significatif.
- Critère de promotion : aucune régression numérique ; gain mesuré et reproductible sur Hamilton–Jacobi et hybride ; coût de construction et mémoire inclus dans le temps total. Si le gain est négligeable ou négatif, ne pas fusionner l'optimisation.

## 7. Statut d'implémentation CFDX

- Fichiers source concernés : `src/cfdx/physics/wall_distance.h`.
- Tests exécutables : `tests/unit/test_wall_distance_methods.cpp`, benchmark `cfdx_wall_distance_controlled_benchmark` et campagnes de validation-total.
- Écarts code ↔ fiche : aucun. Qualification Release/DebugSanitizers/`test_wall_distance_methods` et comparaison de performance exact-HEAD terminées (campagne 39374525195 au head 54fd3412 : tous les contrôles rapides verts, 9/9 fine benchmarks wall-distance N=128 verts, coarse vert, N9 Qualification/Physical verts).
- Gaps et prochaines actions : aucune bloquante. Suivi informationnel ouvert : isolement du composant responsable du surcoût HJ sur runner EPYC 7763 (cf. §8), sans impact sur la décision de fusion.

## 8. Qualification exact-HEAD et décision (2026-10-20)

**A/B pinned même runner (run 39397470405, AMD EPYC 7763, pins base `ede4f918a78802a68708b72fd03643e91b37b346` / head `54fd341255ae6166114d098ef0424eec465b34cc`)** :
- hybride Poisson–Eikonal : base 1 471 040 ms vs head 872 777 ms → **−40,7 %** (157,9 → 93,7 ms/iter ; 9 316 itérations identiques).
- Hamilton–Jacobi : base 2 364 440 ms vs head 2 708 698 ms → **+14,6 %/iter** (174,5 → 199,9 ms/iter ; 13 550 itérations identiques).
- Équivalence numérique bit-exacte : itérations, converged, stopping reason, residual_inf, l2_relative, monotonicity_violations=0 — aucun seuil, budget ou tolérance modifié.

**Campagne validation-total (run 39374525195, head 54fd3412)** : 45/45 check-runs relevés ; tous les contrôles rapides verts (Pre-commit, Validate, Audit gate, Release, DebugSanitizers, N1–N16, E2E, GUI, numerical) ; 9/9 fine benchmarks wall-distance N=128 SUCCESS ; coarse SUCCESS ; N9 Qualification + N9 Physical SUCCESS. Seuls échecs : `test_naca0012_qualification` et `test_naca0012_wall_shear_audit` — **préexistants sur master** (bisects rounds 1–3 : master rouge au merge-base aa806657 du 2026-10-14, signature de stall identique), hors du domaine wall-distance, non imputables à cette PR.

**Profiling de variantes (run 39463278884, AMD EPYC 9V45)** : 4 variantes de `wall_distance.h` construites sur l'arbre exact du head (v0=head, v2=boucle scan, v3=stencils BVH, v1=boucle HJ master) — 44/44 colonnes numériques identiques entre variantes ; sur ce runner, toutes les différences inter-variantes sont ≤ 3,1 % et l'ordre v0/v1 est inversé (v0 100,3 vs v1 103,4 ms/iter) : la régression HJ ne se reproduit pas sur ce micro-matériel.

**Caveat matériel (consigné)** : le surcoût HJ +14,6 % est dépendant de la micro-architecture — mesuré de façon cohérente sur EPYC 7763 (Zen 3, DDR4 ; 3 mesures indépendantes du head à ~2 709–2 711 s contre base à ~2 364–2 365 s), absent sur EPYC 9V45 (Zen 5, DDR5). Bilan net des deux méthodes lourdes sur le matériel le plus défavorable : **−254 s (−6,6 %)** par solve complet (−598 s hybride vs +344 s HJ). Le gain hybride est reproductible sur tous les runners mesurés.

**Décision** : le critère de promotion du §6 (« gain mesuré et reproductible sur Hamilton–Jacobi et hybride ») est satisfait par décision explicite de compromis de production (Steven, 2026-10-20) : le gain hybride −40,7 % est conservé, le surcoût HJ matériel-dépendant est accepté et documenté. Le 4e commit candidat (restauration du test `g.solid[q]` avant `offsets.delta(id,face)` dans le lambda `up` de l'hybride) reste hors périmètre de cette PR. Fusion autorisée (feu vert explicite Steven 2026-10-20) après la présente mise à jour de fiche. Aucune équation, tolérance, budget d'itérations ou critère de convergence modifié.

## Optimisation du coût par itération (2026-10-20)

- Hamilton–Jacobi et l'hybride construisent maintenant une liste des nœuds fluides, non fixes et de distance finie, dans l'ordre croissant des indices. Les boucles de mise à jour et de résidu réutilisent cette liste : l'ordre Gauss–Seidel en place est conservé, tandis que les tests solide/fixe/infini ne sont plus exécutés pour chaque nœud à chaque itération.
- Le pas advectif constant et le facteur de relaxation constant de Hamilton–Jacobi sont calculés une fois. Dans l'hybride, la somme inverse des espacements au carré et le facteur de relaxation sont aussi sortis de la boucle interne.
- Aucun budget d'itérations, seuil de convergence, coefficient PDE ni valeur de relaxation n'est modifié. Les gains et l'équivalence ont été mesurés sur le HEAD exact (A/B run 39397470405, §8).
- Hamilton–Jacobi ne répète plus les requêtes de distance exacte aux nœuds de seed : l'initialisation Eikonal applique déjà exactement le même critère de seed et calcule déjà leurs distances BVH. Le masque fixe est désormais renseigné sans recalcul géométrique.
