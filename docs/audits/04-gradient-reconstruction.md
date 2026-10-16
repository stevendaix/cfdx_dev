# FICHE MÉTHODE 04 — Reconstruction du gradient cellulaire aux faces

**En-tête d'évidence**
- Scope : gradient scalaire cellulaire et reconstruction des valeurs de face sur maillages polyédriques 3D.
- Statut : DRAFT — proposition de correction en cours de qualification (2026-10-16).
- Issues liées : #922, PR #925.
- Dernier commit vérifié : aucun; la compilation et les campagnes numériques doivent encore être confirmées par CI.

## 1. Références sources
- Greenshields & Weller, *Notes on Computational Fluid Dynamics*, discussion des opérateurs de Green–Gauss et des gradients aux faces.
- La formulation ci-dessous est un contrat numérique proposé, pas une preuve expérimentale de l'ordre.

## 2. Formulation mathématique
Pour une face interne f, de centre géométrique C_f, adjacente aux cellules P et N, reconstruire deux valeurs :
\[
\phi_{f,P}=\phi_P+(\nabla\phi)_P\cdot(C_f-C_P),\qquad
\phi_{f,N}=\phi_N+(\nabla\phi)_N\cdot(C_f-C_N).
\]
Puis utiliser la moyenne symétrique
\[
\phi_f=\tfrac12(\phi_{f,P}+\phi_{f,N})
\]
dans la somme de Green–Gauss. Les gradients de cellule sont estimés par moindres carrés. Cette reconstruction est linéairement exacte sur une face plane si les deux stencils de moindres carrés sont de rang suffisant et si les valeurs de face respectent les conditions aux limites.

## 3. Discrétisation Volume Fini (forme attendue dans le code)
- Face interne : extrapolation affine depuis les deux centres de cellules jusqu'au centre géométrique de face, puis moyenne des deux valeurs.
- Face de frontière : convention historique φ_f = φ_P (gradient normal nul), faute de BoundaryField explicite. Ce n'est pas une condition de Dirichlet générale.
- L'implémentation doit rejeter les indices hors limites, volumes non positifs et valeurs reconstruites non finies.
- Aucun changement de schéma par défaut n'est couvert par cette fiche.

## 4. Limites connues et pièges
- Un stencil LS déficient en rang peut ne pas reproduire un champ linéaire; aucune exactitude universelle n'est revendiquée.
- La moyenne des gradients reconstruits aux centres ne démontre pas l'ordre deux du gradient sur toute famille de tétraèdres.
- Cette modification ne corrige pas, à elle seule, la précision du Laplacien corrigé. La campagne test_polyhedral_laplacian_campaign doit mesurer les erreurs et l'ordre observé sur champ manufacturé; un plateau reste un défaut ouvert, pas un résultat qualifié.
- Ne pas masquer un plateau en supprimant le seuil de convergence, en augmentant les tolérances ou en désactivant les tests.

## 5. Vérification exigée
1. Test unitaire sur champ affine et cellule intérieure avec faces internes skewées.
2. Vérifier les résidus constants, la finitude et l'absence de changement des tolérances.
3. Mesurer l'erreur L2, L∞ et l'ordre observé sur une famille de maillages tétraédriques.
4. Rejouer validation-total et les cas physiques pertinents avant toute promotion d'un schéma par défaut.

## 6. Statut d'implémentation CFDX
- Source : src/cfdx/core/numerics/gradient.h, compute_gradient_gauss_point.
- Test : tests/unit/test_gradient.cpp, régression de champ linéaire sur cellule intérieure.
- CI : lancée sur le HEAD de la PR; résultats en attente au moment de la rédaction.
- Laplacien polyédrique : non résolu par cette PR, validation séparée requise.

## 7. Journal
- 2026-10-16 — PR #925 remplace la reconstruction historique point/vertex par une extrapolation LS depuis les deux cellules adjacentes. La compilation et les campagnes numériques restent à confirmer par CI; cette fiche reste DRAFT jusqu'à preuve.
