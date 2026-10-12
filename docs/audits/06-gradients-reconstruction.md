# FICHE MÉTHODE 06 — Gradients : Green–Gauss et moindres carrés (linéaire/quadratique)

**En-tête d'évidence**
- Scope : opérateurs de gradient aux centres (Green–Gauss, LSQ linéaire, LSQ quadratique) pour la reconstruction, les termes sources, la diffusion non orthogonale et le post-traitement. Couvre les maillages polyédriques/tétraédriques (N2) et le calcul aux faces.
- Statut : APPROVED (2026-10-12)
- Issues liées : #461 (N2), #639 (N10) ; PRs #511 (en attente d'évaluation), #845, #852
- Dernier commit vérifiéà: 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : reconstruction quadratique tétraédrique en cours de qualification (PR #511) — pas de promotion sans campagne exact-HEAD.

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

## 4. Traitement des frontières
Aux faces frontières, $\\phi_f$ doit être la **valeur frontière** (Dirichlet) pour la somme Green–Gauss — pas la valeur de centre. Neumann : reconstruction unilatérale. Le gradient des cellules frontières doit être ré-évalué avec le traitement limite.

## 5. Pièges connus
- Stencils insuffisants sur tétraèdres → LSQ linéaire d'ordre 1 (écart N2 documenté) ; stencil 3 anneaux requis pour le quadratique.
- Rang déficient : un pivot qui évite silencieusement → résultat faux ; le rang doit être **compté et diagnostiqué** (cardinalité, rang, valeurs singulières, conditionnement — diagnostics ajoutés par #511) et un stencil réellement déficient doit être **rejeté**.
- Échelle : les seuils doivent être sensibles à l'échelle du stencil, pas des constantes absolues.
- Green–Gauss aux faces pour la diffusion non orthogonale exige le même traitement limite que le gradient central.

## 6. Ordre de précision attendu et validation exigée
- Exactitude constante + linéaire : **preuve obligatoire à chaque niveau de maillage** (#845) ; LSQ quadratique : reproduction quadratique exacte sur tétraèdres intérieurs.
- Ordre mesuré sur champ lisse : $\\geq 1.8$ pour le LSQ quadratique tétraédrique, avec **les vrais ratios de raffinement** (n=4→6→8 : $h_c/h_f = 6/4$ puis $8/6$) — gate ajoutée par #511 ; l'ordre mesuré doit être documenté, jamais extrapolé.
- Campagne monômes P3 ($x^3, y^3, z^3, x^2y, …, xyz$) obligatoire avant promotion N2.

## 7. Statut d'implémentation CFDX
- Sources : `src/cfdx/core/numerics/gradient.h’ (PR #511 : reconstruction quadratique échelonnée, seuil de pivot sensible à l'échelle, diagnostics de stencil).
- Tests : `tests/validation/test_polyhedral_gradient_campaign.cpp’ (#511), campagne Poisson polyédrique #852, contrôle d'exactitude linéaire #845.
- Écarts code ↔ fiche : à auditer — PR #511 en attente : **pertinente mais à rebaser** sur le master actuel et re-auditer contre cette fiche ; le rang par valeur singulière reste un seuillage, un estimateur SVD/QR relatif est un suivi séparé assumé.
- Gaps : familles non orthogonales/skew/high-aspect-ratio du gradient (N10/#639) ; parité MPI.

## 8. Journal
- 2026-10-12 — création ; synthétise #845/#852 et l'évaluation de #511.
