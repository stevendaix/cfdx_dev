# FICHE MÉTHODE 25 — Courbure d’interface VOF par divergence de la normale

**Statut : PARTIEL — contrat du helper seulement.** Cette fiche ne qualifie pas la reconstruction PLIC, le calcul géométrique multidimensionnel de la normale, ni l’intégration de la courbure dans le solveur de production.

## 1. Contrat mathématique

Pour une normale unitaire à l’interface \(\mathbf n=(n_x,n_y,n_z)\), la courbure signée est définie par le contrat choisi par CFDX :

\[
\kappa=-\nabla\cdot\mathbf n
=-\left(\frac{\partial n_x}{\partial x}
+\frac{\partial n_y}{\partial y}
+\frac{\partial n_z}{\partial z}\right).
\]

Chaque composante doit être dérivée selon sa propre coordonnée. Un stencil unidimensionnel ne peut pas réutiliser un même indice spatial pour combiner des dérivées indépendantes.

## 2. Domaine et limites

- En 2D, le terme \(\partial n_z/\partial z\) vaut zéro et le helper reçoit explicitement zéro par défaut.
- En 3D, les trois dérivées doivent être fournies.
- Les entrées non finies sont rejetées.
- Le signe de la courbure dépend de la convention d’orientation de la normale ; l’appelant doit garantir une orientation cohérente.
- Les tests analytiques de plans, cercles et sphères vérifient le contrat du helper, mais ne prouvent pas l’exactitude de reconstruction sur maillage arbitraire.

## 3. Validation requise avant qualification de production

1. Comparer la normale et sa divergence à une interface analytique sur maillages raffinés.
2. Vérifier la convergence de la courbure sur des interfaces sphériques et circulaires, avec orientation contrôlée.
3. Tester interfaces proches de frontières, cellules dégénérées et transitions de fraction volumique.
4. Démontrer que la courbure est appelée par le solveur VOF de production et que le terme capillaire respecte le bilan de quantité de mouvement.

## 4. Références

- S. Popinet, “An accurate adaptive solver for surface-tension-driven interfacial flows”, *Journal of Computational Physics*, 228 (2009), 5838–5866, DOI: 10.1016/j.jcp.2009.04.042.
- J. U. Brackbill, D. B. Kothe, C. Zemach, “A continuum method for modeling surface tension”, *Journal of Computational Physics*, 100 (1992), 335–354, DOI: 10.1016/0021-9991(92)90240-Y.
