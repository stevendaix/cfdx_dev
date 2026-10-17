# FICHE MÉTHODE 06 — Gradients : Green–Gauss et moindres carrés (linéaire/quadratique)

**En-tête d'évidence**
- Scope : opérateurs de gradient aux centres (Green–Gauss, LSQ linéaire, LSQ quadratique) pour la reconstruction, les termes sources, la diffusion non orthogonale et le post-traitement. Couvre les maillages polyédriques/tétraédriques (N2) et le calcul aux faces.
- Statut : APPROVED (2026-10-12)
- Issues liées : #461 (N2), #639 (N10) ; PRs #511 (en attente d'évaluation), #845, #852, #925 (reconstruction GAUSS_POINT aux faces)
- Dernier commit vérifié : 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : reconstruction quadratique tétraédrique en cours de qualification (PR #511) — pas de promotion sans campagne exact-HEAD.
- Limitations (2026-10-16) : reconstruction GAUSS_POINT aux faces (PR #925) en qualification — CI Release standard à confirmer, ordres de campagne à re-mesurer ; pas de promotion sans preuve.

## 1. Références sources
- F. Moukalled, L. Mangani, M. Darwish (2016), ISBN 978-3-319-16873-9, ch. 9–10 (Green–Gauss, LSQ, reconstruction).
- J. H. Ferziger, M. Perić (2002), DOI 10.1007/978-3-642-56026-2, § 8.3.
- T. Barth, D. Jespersen, « The design and application of upwind schemes on unstructured meshes », AIAA 89-0366, 1989 — DOI : 10.2514/6.1989-366. (limitation des gradients aux faces)

## 2. Formulation continue
Le gradient de cellule doit reproduire **exactement** les champs constants (gradient nul) et linéaires (gradient exact).

## 3. Discrétisation Volume Fini
**Green–Gauss** :
$$ (\\nabla\\phi)_P \\approx \\frac{1}{V_P}\\sum_f \\phi_f\\,\\mathbf{S}_f, \\qquad \\phi_f \\approx \\text{moy. pondérée de } \\phi_P,\\phi_N $$
Contrat : $\\phi_N$ doit être identifié correctement quand la cellule courante est un voisin (leçon A-B1 #419 : le bug $\\phi_f=\\phi_c$ divisait le gradient de pression par ~2 → Poiseuille faux à 44 %).
**LSQ linéaire** : minimisation sur le stencil — exact sur champs linéaires ; ordre 1 sur maillages à courbure (stencils tétraédriques).
**LSQ quadratique** : extension 9 composantes ($x,y,z$ + termes du 2e ordre), équations normales avec seuillage ; gradient = 3 premières composantes ; restaure ~$O(h^2)$ sur tétraèdres (PR #511).
**Green–Gauss point-reconstruit (GAUSS_POINT, PR #925)** :
φ_f = ½·[φ_P + (∇φ)_P·(C_f−C_P) + φ_N + (∇φ)_N·(C_f−C_N)] — moyenne symétrique des deux extrapolations affines au centre géométrique de face, gradients cellulaires (∇φ) estimés par LSQ linéaire.
Contrat : plus d'interpolation scalaire sur la ligne owner–neighbour ni de moyenne aux sommets. Exactitude linéaire conditionnelle : stencils LSQ de rang suffisant, faces planes, valeurs de frontière compatibles (voir §4).

## 4. Traitement des frontières
Aux faces frontières, $\\phi_f$ doit être la **valeur frontière** (Dirichlet) pour la somme Green–Gauss — pas la valeur de centre. Neumann : reconstruction unilatérale. Le gradient des cellules frontières doit être ré-évalué avec le traitement limite.
**GAUSS_POINT (PR #925)** : la voie point-reconstruit conserve la convention historique φ_f = φ_P (zéro-gradient) aux faces frontières, faute de valeurs limites plombées jusqu'à l'opérateur — écart assumé au contrat ci-dessus. L'exactitude linéaire n'est donc revendiquée que pour des cellules dont les faces frontières sont compatibles avec cette politique ; le raccordement de valeurs Dirichlet explicites est un suivi requis avant toute promotion.

## 5. Pièges connus
- Stencils insuffisants sur tétraèdres → LSQ linéaire d'ordre 1 (écart N2 documenté) ; stencil 3 anneaux requis pour le quadratique.
- Rang déficient : un pivot qui évite silencieusement → résultat faux ; le rang doit être **compté et diagnostiqué** (cardinalité, rang, valeurs singulières, conditionnement — diagnostics ajoutés par #511) et un stencil réellement déficient doit être **rejeté**.
- Échelle : les seuils doivent être sensibles à l'échelle du stencil, pas des constantes absolues.
- Green–Gauss aux faces pour la diffusion non orthogonale exige le même traitement limite que le gradient central.
- GAUSS_POINT (PR #925) : l'exactitude dépend du rang des stencils LSQ des deux cellules adjacentes ; un stencil déficient ne reproduit pas un champ linéaire. La reconstruction aux faces ne corrige pas à elle seule le plateau d'erreur du Laplacien polyédrique corrigé — la campagne dédiée doit re-mesurer L2/L∞ et l'ordre observé sur champ manufacturé ; un plateau reste un défaut ouvert, pas un résultat qualifié ; aucun seuil, tolérance ou test ne doit être ajusté pour le masquer.

## 6. Ordre de précision attendu et validation exigée
- Exactitude constante + linéaire : **preuve obligatoire à chaque niveau de maillage** (#845) ; LSQ quadratique : reproduction quadratique exacte sur tétraèdres intérieurs.
- Ordre mesuré sur champ lisse : $\\geq 1.8$ pour le LSQ quadratique tétraédrique, avec **les vrais ratios de raffinement** (n=4→6→8 : $h_c/h_f = 6/4$ puis $8/6$) — gate ajoutée par #511 ; l'ordre mesuré doit être documenté, jamais extrapolé.
- Campagne monômes P3 ($x^3, y^3, z^3, x^2y, …, xyz$) obligatoire avant promotion N2.
- GAUSS_POINT (PR #925) : preuve d'exactitude linéaire sur cellule intérieure (régression tests/unit/test_gradient.cpp) ; ordre mesuré sur familles tétraédriques et validation-total à documenter avant toute promotion du schéma ; la PR ne modifie aucun seuil ni tolérance.

## 7. Statut d'implémentation CFDX
- Sources : `src/cfdx/core/numerics/gradient.h’ (PR #511 : reconstruction quadratique échelonnée, seuil de pivot sensible à l'échelle, diagnostics de stencil).
- Tests : `tests/validation/test_polyhedral_gradient_campaign.cpp’ (#511), campagne Poisson polyédrique #852, contrôle d'exactitude linéaire #845.
- PR #925 : compute_gradient_gauss_point reconstruit les faces internes par extrapolation LSQ au centre géométrique (supprime la reconstruction historique point-linéaire/vertex) ; régression de champ linéaire intérieur dans tests/unit/test_gradient.cpp.
- Écarts code ↔ fiche : à auditer — PR #511 en attente : **pertinente mais à rebaser** sur le master actuel et re-auditer contre cette fiche ; le rang par valeur singulière reste un seuillage, un estimateur SVD/QR relatif est un suivi séparé assumé.
- Gaps : familles non orthogonales/skew/high-aspect-ratio du gradient (N10/#639) ; parité MPI. Raccordement des valeurs limites Dirichlet au GAUSS_POINT (§4) ; re-mesure de la campagne Laplacien polyédrique corrigé (plateau ouvert).

## 8. Journal
- 2026-10-12 — création ; synthétise #845/#852 et l'évaluation de #511.
- 2026-10-16 — PR #925 : contrat GAUSS_POINT point-reconstruit ajouté (§3–§7) ; la reconstruction historique point-linéaire/vertex est retirée du code ; qualification CI en cours (job Release standard en échec sur le head initial 98e8905b — diagnostic en cours sur la branche diag, preuves diag-evidence/logs/pr925_release_*).
