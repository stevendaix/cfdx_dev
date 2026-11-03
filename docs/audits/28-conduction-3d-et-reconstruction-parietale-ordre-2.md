# FICHE MÉTHODE 28 — Conduction 3D et reconstruction pariétale quadratique

**Statut :** DRAFT — implémentation partielle, non qualifiée.
**Périmètre :** diagnostic fail-closed de conduction stationnaire 3D et reconstruction de pression/traction pariétale par WLS quadratique.
**Issue/PR :** PR #1169, créée depuis le HEAD de la PR #1146.

## 1. Références sources

- F. Moukalled, L. Mangani, M. Darwish (2016), *The Finite Volume Method in Computational Fluid Dynamics*, ISBN 978-3-319-16873-9.
- J. H. Ferziger, M. Perić (2002), *Computational Methods for Fluid Dynamics*, DOI 10.1007/978-3-642-56026-2.
- J. Blazek (2015), *Computational Fluid Dynamics: Principles and Applications*, ISBN 978-0-08-099995-1.

## 2. Formulation continue

Pour la conduction stationnaire sans source, le champ manufacturé T=300+10z vérifie div(k grad(T))=0 et q=-k grad(T)=-10k e_z. Les frontières latérales sont adiabatiques et les deux plans d'extrémité imposent les températures analytiques. Le flux total attendu est constant et vaut 10kA en valeur absolue.

Pour un fluide newtonien, la traction pariétale est évaluée avec tau=mu(grad(u)+grad(u)^T), puis t_f=p_f n-tau_f n. Une reconstruction polynomiale quadratique complète doit restituer exactement un champ quadratique et donc son gradient linéaire.

## 3. Discrétisation Volume Fini

Le diagnostic 3D existant assemble l'équation sur le champ exact puis mesure le résidu A T-b. Le résidu absolu annoncé (~290) n'est pas interprétable sans sa décomposition par ligne et par face. Les contributions à auditer sont la conductance orthogonale, la correction tangente non-orthogonale, les distances normales, les gradients et les termes de frontière.

La nouvelle API wall_face_reconstruction_quadratic.h ajuste les dix monômes 3D par moindres carrés pondérés, avec coordonnées normalisées et factorisation QR de Householder à pivotement de colonnes. Le rapport des diagonales de R fournit un indicateur conservateur de rang/conditionnement ; ce n'est pas une estimation SVD du conditionnement spectral. Les gradients au centre de face sont les coefficients linéaires du polynôme centré en ce point. La contrainte visqueuse est calculée avec le tenseur symétrique complet.

## 4. Traitement des frontières

Le stencil pariétal doit contenir les valeurs des cellules voisines et la valeur de vitesse prescrite au centre de face (condition de non-glissement). Le rang insuffisant est rejeté explicitement ; aucun fallback affine silencieux n'est permis.

Un constructeur de stencil depuis l'adjacence réelle est disponible dans wall_quadratic_stencil.h : parcours BFS des cellules par faces internes, inclusion explicite de la vitesse de paroi prescrite, contrôle du patch WALL et rejet des entrées invalides. La pression au point de paroi est fournie explicitement selon le contrat de pression du patch. L'adaptateur wall_quadratic_forces.h convertit les résultats en forces::SurfaceSample, consommables par forces::integrate.

Le chemin de production est maintenant accessible par l'option explicite `--wall-forces-quadratic`. Il exige une solution convergée, une vitesse prescrite FIXED_VALUE sur chaque patch WALL et une pression FIXED_VALUE ou ZERO_GRADIENT. Il écrit les forces et moments dimensionnels dans `wall_forces.json`, avec la convention de force fluide-sur-corps et le centre de moment explicite (origine). Les coefficients non dimensionnels sont volontairement absents tant que rho_inf, U_inf, la surface et la longueur de référence ne sont pas déclarés. Tout stencil de rang déficient provoque un échec explicite ; aucun fallback affine n'est appliqué.

Les tests couvrent un cube unitaire (rejet de stencil insuffisant), un hexaèdre affine-shear 3x3x3 (géométrie non-orthogonale), et un O-grid NACA0012 extrudé en quatre couches (champ polynomial manufacturé et traction exacte sur chaque face airfoil). Ce dernier reste un test de reconstruction sur géométrie CFD réelle, pas une qualification d'un champ de solution CFD. Les stencils 2D coplanaires et les extrusions à deux couches restent volontairement rejetés.

## 5. Pièges connus

- Un résidu absolu de matrice dépend de l'échelle des coefficients ; il faut aussi rapporter résidu relatif, backward error et contributions ligne par ligne.
- Ne pas changer le signe de la correction non-orthogonale sans établir la convention du flux intégré et les signes owner/neighbour.
- Les équations normales WLS amplifient le conditionnement ; l'implémentation utilise désormais une QR de Householder pivotée et rejette les stencils dont le rapport diagonal de R est inférieur ou égal à 1e-10. Ce critère reste un proxy, pas une estimation SVD.
- Un stencil coplanaire ne peut pas identifier les dix coefficients d'un polynôme 3D.
- Le test MMS algébrique ne démontre pas à lui seul l'ordre observé sur un maillage CFD réel.

## 6. Ordre de précision attendu et validation exigée

- Conduction : exactitude du champ T=300+10z, erreur L_inf, flux analytique, fermeture globale et conservation locale.
- Paroi : pression et gradient de vitesse exacts sur champ quadratique 3D ; traction exacte ; valeur de paroi immobile explicitement incluse ; champ cubique manufacturé sur trois échelles, avec ordre observé supérieur à 1.7 sur les deux intervalles.
- Géométrie : stencil non coplanaire, faces obliques, maillages non-orthogonaux et rejet des stencils dégénérés.
- Exécution : CTest ciblé, Release, DebugSanitizers et campagnes longues pertinentes.

Les tests `test_wall_face_reconstruction_quadratic_3d`, `test_wall_quadratic_stencil` et `test_naca3d_wall_quadratic_mesh` sont enregistrés dans CTest. Le dernier utilise un O-grid NACA0012 extrudé avec quatre couches spanwise, reconstruit un champ quadratique analytique sur toutes les faces airfoil et compare pression et traction au champ exact. Le test Python de production vérifie aussi l'option CLI et l'écriture du JSON sur un volume 3D fermé. Ces cas vérifient l'algorithme et son raccordement logiciel, mais ne constituent pas une validation d'ordre sur un écoulement CFD calculé.

## 7. Statut d'implémentation

**Implémenté dans la PR #1169 :**
- API de reconstruction WLS quadratique 3D, pression et gradient de vitesse au centre de face.
- construction du stencil unilatéral depuis l'adjacence réelle des cellules, avec inclusion de la vitesse de paroi prescrite et contrôle du patch WALL.
- adaptateur de patch WALL vers forces::SurfaceSample, avec conversion explicite de la convention de normale pour forces::integrate.
- Calcul de la traction newtonienne symétrique.
- Rejet explicite des stencils 3D de rang déficient.
- MMS quadratique 3D enregistré dans CTest, avec non-glissement stationnaire et étude de convergence sur stencil synthétique ;
- résolution WLS par QR de Householder pivotée, avec rejet des stencils mal conditionnés selon un indicateur diagonal ;
- fixture thermique 3D distinct avec patches d'extrémité WALL, afin que les Dirichlet de température ne soient pas ignorés comme sur le fixture fluide à patches EMPTY ;
- test de reproduction et de traction exacte sur l'O-grid NACA0012 extrudé en quatre couches spanwise ;
- intégration dans le CLI du solveur de production via `--wall-forces-quadratic`, avec garde de convergence, contrôle des contrats de frontière, sortie des forces/moments dimensionnels et rejet fail-closed des stencils insuffisants ;
- test E2E Python de l'export des forces sur un maillage 3D wall-only.

**Non implémenté / non qualifié :**
- validation des résultats CI du HEAD courant pour le MMS de conduction et les tests de reconstruction ;
- preuve d'ordre observé 2 sur des champs de solution CFD calculés, avec étude de raffinement et comparaison GG/WLS/quadratique ;
- qualification physique des forces sur un cas expérimental 3D avec valeurs de référence ;
- support des stencils 2D coplanaires (hors contrat de l'API 3D actuelle) et estimation spectrale/SVD du conditionnement.

Aucune tolérance n'est élargie. Aucune qualification physique n'est revendiquée. La PR reste draft jusqu'à obtention de ces preuves.

## 8. Journal

- 2026-11-03 — création de la branche depuis le HEAD de #1146 ; ajout du prototype WLS quadratique 3D, du MMS CTest et de cette fiche. Le défaut de conduction 3D reste ouvert.
- 2026-11-03 — échec initial du gate documentaire : fiche complétée pour respecter la structure méthodologique exigée.
- 2026-11-03 — diagnostic de la frontière thermique : les caps front/back du fixture d'extrusion fluide sont EMPTY et sont donc exclus de l'assemblage FV ; ajout d'un fixture thermique distinct avec caps WALL et d'un CTest de conduction linéaire dans le gate Release.
- 2026-11-03 — test pariétal renforcé avec valeur de paroi nulle et mesure d'ordre sur trois échelles ; remplacement des équations normales par QR de Householder pivotée et rejet conditionné par le rapport diagonal.
- 2026-11-03 — ajout du constructeur BFS de stencil depuis les faces internes du maillage, avec contrôle WALL et échantillon de vitesse de paroi exact ; test fail-closed sur cube unitaire.
- 2026-11-03 — ajout de l'adaptateur vers forces::SurfaceSample et d'une qualification algébrique sur maillage hexaédrique 3x3x3 affine-shear, incluant reproduction quadratique, normale oblique et appel à forces::integrate.
- 2026-11-03 — ajout d'un test de reproduction quadratique sur O-grid NACA0012 extrudé en quatre couches spanwise ; test des pressions et tractions exactes sur toutes les faces airfoil.
- 2026-11-03 — raccordement opt-in au solveur de production via --wall-forces-quadratic ; sortie des forces/moments dimensionnels et test E2E sur volume 3D fermé.
