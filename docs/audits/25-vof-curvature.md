# FICHE MÉTHODE 25 — Courbure d’interface VOF par divergence de la normale

**En-tête d’évidence**
- Scope : contrat du helper de courbure VOF (kappa = −div(n) composante par composante) uniquement. Cette fiche ne qualifie pas la reconstruction PLIC, le calcul géométrique multidimensionnel de la normale, ni l’intégration du terme capillaire dans le solveur de production.
- Statut : PARTIEL — contrat du helper seulement.
- Issues liées : PR #1114 (contrat axis-specific de divergence de la normale).
- Limitations : aucune qualification de production ; validation analytique plans/cercles/sphères requise avant promotion.

## 1. Références sources
- S. Popinet, “An accurate adaptive solver for surface-tension-driven interfacial flows”, *Journal of Computational Physics*, 228 (2009), 5838–5866, DOI: 10.1016/j.jcp.2009.04.042.
- J. U. Brackbill, D. B. Kothe, C. Zemach, “A continuum method for modeling surface tension”, *Journal of Computational Physics*, 100 (1992), 335–354, DOI: 10.1016/0021-9991(92)90240-Y.

## 2. Formulation continue
Pour une normale unitaire à l’interface n = (n_x, n_y, n_z), la courbure signée est définie par le contrat choisi par CFDX :

kappa = −∇·n = −(∂n_x/∂x + ∂n_y/∂y + ∂n_z/∂z).

Chaque composante doit être dérivée selon sa propre coordonnée ; un stencil unidimensionnel ne peut pas réutiliser un même indice spatial pour combiner des dérivées indépendantes.

## 3. Discrétisation Volume Fini
Le helper reçoit explicitement chaque dérivée partielle axis-spécifique (∂n_i/∂i) ; aucune combinaison d’indices croisés n’est admise. En 2D, le terme ∂n_z/∂z vaut zéro et le helper reçoit explicitement zéro par défaut ; en 3D, les trois dérivées doivent être fournies. Les entrées non finies sont rejetées (fail-closed). Les tests analytiques de plans, cercles et sphères vérifient le contrat du helper, mais ne prouvent pas l’exactitude de reconstruction sur maillage arbitraire.

## 4. Traitement des frontières
Le contrat du helper ne couvre pas le traitement des interfaces proches de frontières : la normale aux cellules frontières et son impact sur la divergence restent à qualifier (interfaces coupées par une frontière du domaine, cellules dégénérées, transitions de fraction volumique). Aucune convention limite n’est revendiquée par cette fiche.

## 5. Pièges connus
- Un stencil unidimensionnel qui réutilise un même indice spatial pour combiner des dérivées indépendantes viole le contrat (défaut corrigé par la PR #1114 : divergence « axis-agnostic » rejetée au profit d’un contrat axis-spécifique).
- Le signe de la courbure dépend de la convention d’orientation de la normale ; l’appelant doit garantir une orientation cohérente.
- Les entrées non finies doivent être rejetées, jamais clampées silencieusement.
- Les tests analytiques du helper ne qualifient pas la reconstruction PLIC ni la normale multidimensionnelle.

## 6. Validation requise avant qualification de production
1. Comparer la normale et sa divergence à une interface analytique sur maillages raffinés.
2. Vérifier la convergence de la courbure sur des interfaces sphériques et circulaires, avec orientation contrôlée.
3. Tester interfaces proches de frontières, cellules dégénérées et transitions de fraction volumique.
4. Démontrer que la courbure est appelée par le solveur VOF de production et que le terme capillaire respecte le bilan de quantité de mouvement.

## 7. Statut d’implémentation CFDX
- Sources : `src/cfdx/physics/vof.h` (contrat axis-spécifique de la divergence de la normale, PR #1114).
- Tests : `tests/unit/test_vof_curvature_contract.cpp` (plans/cercles/sphères analytiques, rejet des entrées non finies).
- Écarts code ↔ fiche : aucun sur le périmètre du helper ; périmètre PLIC/production hors contrat, non couvert.
- Gaps : intégration au solveur VOF de production, bilan de quantité de mouvement du terme capillaire, interfaces proches de frontières — validation à documenter avant toute promotion.

## 8. Journal
- 2026-10-25 — création (PR #1114) ; restructurée au gabarit Rule of Gold 1 (sections 1–7) pour satisfaire la gate structurelle `check_audit_gate.py --validate-fiches` ; contenu mathématique et périmètre inchangés.
