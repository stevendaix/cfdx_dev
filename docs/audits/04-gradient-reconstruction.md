# Fiche méthode 04 — Reconstruction du gradient cellulaire aux faces

**Statut :** DRAFT — proposition de correction en cours de qualification (2026-10-16).  
**Portée :** gradient scalaire cellulaire et reconstruction des valeurs de face sur maillages polyédriques 3D. Cette fiche ne certifie pas à elle seule l'ordre du Laplacien diffusif.

## 1. Opérateur

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

## 2. Limites explicites

- Les faces de frontière conservent la convention historique φ_f = φ_P (gradient normal nul). Ce n'est pas une condition de Dirichlet générale.
- Un stencil LS déficient en rang peut ne pas reproduire un champ linéaire; aucune exactitude universelle n'est revendiquée.
- La moyenne de gradients reconstruits aux centres ne démontre pas l'ordre deux du gradient sur toute famille de tétraèdres.
- Cette modification ne corrige pas, à elle seule, la précision du Laplacien corrigé. La campagne test_polyhedral_laplacian_campaign doit mesurer les erreurs et l'ordre observé sur champ manufacturé; un plateau reste un défaut ouvert, pas un résultat qualifié.

## 3. Vérification exigée

1. Test unitaire sur champ affine et cellule intérieure avec faces internes skewées.
2. Vérifier les résidus constants, la finitude et l'absence de changement des tolérances.
3. Mesurer l'erreur L2, L∞ et l'ordre observé sur une famille de maillages tétraédriques.
4. Rejouer validation-total et les cas physiques pertinents avant toute promotion d'un schéma par défaut.

## 4. Références et journal

- Greenshields & Weller, Notes on Computational Fluid Dynamics, discussion des opérateurs de Green–Gauss et des gradients aux faces.
- 2026-10-16 : PR #925 remplace la reconstruction historique point/vertex par une extrapolation LS depuis les deux cellules adjacentes. La compilation et les campagnes numériques restent à confirmer par CI; cette fiche reste DRAFT jusqu'à preuve.
