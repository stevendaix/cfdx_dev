# FICHE MÉTHODE 20 — Cache des offsets de faces coupées (wall-distance PDE)

**En-tête d'évidence**
- Scope : optimisation de calcul des offsets géométriques `delta` pour les stencils wall-distance Hamilton–Jacobi et hybride Poisson–Eikonal ; ne redéfinit ni la géométrie ni les équations discrètes.
- Statut : DRAFT (2026-10-19) — non qualifié tant que les benchmarks exact-HEAD et les contrôles numériques ne sont pas terminés.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json`.
- Issue / PR liées : PR #989 ; optimisation BVH distincte PR #979.
- Dernier commit vérifié : à mettre à jour après qualification.
- Limitations connues : le gain dépend du nombre de nœuds adjacents aux faces coupées, du coût BVH et du coût de construction du cache ; la qualification de performance n'est pas acquise.

## 1. Références sources

- R. A. Shortley, R. Weller, « The numerical solution of Laplace's equation », *Journal of Applied Physics*, 9 (1938), 334–348 — DOI: 10.1063/1.1714155.
- F. Moukalled, L. Mangani, M. Darwish (2016), *The Finite Volume Method in Computational Fluid Dynamics*, ISBN 978-3-319-16873-9.

## 2. Formulation continue

La méthode de distance à la paroi est déterminée par une EDP de transport ou d'Eikonal. Le cache est une optimisation d'implémentation : pour une géométrie, un masque solide et une grille immuables durant un solve, chaque offset géométrique de face coupée est une constante du solve.

Pour une face coupée reliant un nœud fluide (P) à un nœud solide (Q), l'offset est (delta_{P,Q}), borné par le plancher et l'espacement local comme dans `poisson_wall_offset()`. Le cache ne change pas cette définition.

## 3. Discrétisation Volume Fini

Les stencils utilisent six slots par nœud fluide : (x^-), (x^+), (y^-), (y^+), (z^-), (z^+). Pour un voisin fluide, la distance reste (h_a). Pour un voisin solide, la distance est (delta_{P,Q}), calculée une fois par le BVH et conservée jusqu'à la fin du solve.

Le cache est creux : un index par nœud identifie les nœuds ayant au moins une face coupée ; seuls ces nœuds possèdent les six distances. Les nœuds sans face coupée utilisent directement (h_x,h_y,h_z). L'ordre des opérations dans les stencils, les coefficients et les valeurs de frontière ne doivent pas changer du fait de la mise en cache.

Le solveur de Poisson hybride doit employer les mêmes offsets pour l'opérateur et le préconditionneur diagonal. Hamilton–Jacobi et hybride doivent également employer le même stencil géométrique dans la mise à jour et l'évaluation du résidu.

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
- Équivalence opérateur : gradient de Godunov, Laplacien de Shortley–Weller, opérateur de Poisson et diagonale doivent correspondre au chemin non caché à une tolérance absolue pré-déclarée de (10^{-13}) sur les tests déterministes.
- Régression : Release et DebugSanitizers, test `test_wall_distance_methods`, puis campagnes wall-distance existantes.
- Performance : comparer le même benchmark contrôlé N=128, même runner et mêmes budgets, base pré-cache contre HEAD cache ; conserver les neuf CSV, temps par méthode, itérations, résidus, stopping reason, erreurs L1/L2/L∞/near-wall et empreinte mémoire. Présenter médiane ou répétitions si le bruit du runner est significatif.
- Critère de promotion : aucune régression numérique ; gain mesuré et reproductible sur Hamilton–Jacobi et hybride ; coût de construction et mémoire inclus dans le temps total. Si le gain est négligeable ou négatif, ne pas fusionner l'optimisation.

## 7. Statut d'implémentation CFDX

- Fichiers source concernés : `src/cfdx/physics/wall_distance.h`.
- Tests exécutables : `tests/unit/test_wall_distance_methods.cpp`, benchmark `cfdx_wall_distance_controlled_benchmark` et campagnes de validation-total.
- Écarts code ↔ fiche : qualification Release/DebugSanitizers et comparaison de performance exact-HEAD en cours ; aucune conclusion de gain ne doit être publiée avant ces résultats.
- Gaps et prochaines actions : confirmer les contrôles CI, comparer les neuf CSV N=128 base/cache sur le même runner, examiner les changements d'itérations dus à la cohérence du résidu hybride, puis décider de la fusion.

