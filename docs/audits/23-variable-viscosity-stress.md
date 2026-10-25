# FICHE MÉTHODE 23 — Contrainte visqueuse transposée à viscosité variable

**En-tête d'évidence**
- Scope : contribution \(\nabla\cdot(\mu_{eff}(\nabla\mathbf U)^T)\) dans le solveur incompressible de production et dans le résidu final de quantité de mouvement. Non couvert : MMS complet Navier–Stokes, GCI et qualification polyédrique.
- Statut : PARTIEL — flux conservatif et oracle affine ajoutés ; CI exacte en attente.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json`, `NUMERICAL_METHOD_MATURITY_AUDIT.json/.md`, `NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json`.
- Issues liées : F11, #461 (N7/N12).
- Limitation : les tests affines orthogonal et hexaédrique incliné vérifient le flux discret, mais les maillages polyédriques arbitraires et la convergence globale du solveur restent à qualifier.

## 1. Références sources
- J. H. Ferziger, M. Perić, *Computational Methods for Fluid Dynamics*, 3e éd., Springer, 2002 — DOI : 10.1007/978-3-642-56026-2. Discrétisation conservative des contraintes visqueuses.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN : 978-3-319-16873-9. Forme tensorielle des contraintes et flux de face.
- P. M. Gresho, R. L. Sani, *Incompressible Flow and the Finite Element Method*, Wiley, 2000 — ISBN : 978-0-471-09850-8. Formulation incompressible et interprétation des contraintes.

## 2. Formulation continue
Pour un fluide newtonien incompressible avec viscosité dynamique effective variable :
\[
\nabla\cdot\left[\mu_{eff}\left(\nabla\mathbf U+(\nabla\mathbf U)^T\right)\right].
\]
La contribution manquante dans l'opérateur scalaire \(\nabla\cdot(\mu_{eff}\nabla U_i)\) est
\[
\left[\nabla\cdot(\mu_{eff}(\nabla\mathbf U)^T)\right]_i
=\sum_j\frac{\partial}{\partial x_j}
\left(\mu_{eff}\frac{\partial U_j}{\partial x_i}\right).
\]
Même lorsque \(\nabla\cdot\mathbf U=0\), le terme \(\sum_j(\partial_j\mu_{eff})(\partial_i U_j)\) n'est généralement pas nul si la viscosité varie spatialement.

## 3. Discrétisation Volume Fini
Le flux de la contrainte transposée pour la composante de quantité de mouvement \(i\) est
\[
F^T_{i,f}=\mu_f\sum_j\left(\frac{\partial U_j}{\partial x_i}\right)_f S_{f,j}.
\]
- Les gradients de chaque composante de vitesse sont reconstruits composante par composante. Les cellules ayant au moins trois voisins internes utilisent le gradient aux moindres carrés pour préserver l'exactitude affine sur les maillages inclinés ; les stencils sous-échantillonnés conservent le gradient Green–Gauss sensible aux conditions aux limites.
- \(\mu_f\) et les gradients sont interpolés sur les faces internes par pondération des distances normales propriétaire/voisin.
- Le flux intégré est assemblé de manière conservative : \(+F^T_{i,f}\) au second membre du propriétaire et \(-F^T_{i,f}\) au voisin.
- Le terme est explicite et recalculé à chaque itération non linéaire. Il est inclus à la fois dans les trois prédicteurs de quantité de mouvement de production et dans le résidu final réassemblé indépendamment.

## 4. Traitement des frontières
- Les valeurs imposées de vitesse sont converties en conditions scalaires fixes pour chacune des trois composantes ; les conditions de gradient nul restent des gradients nuls.
- Sur une face de frontière, le gradient cellulaire propriétaire fournit le gradient de face de premier ordre ; cette approximation n'est qualifiée ici que pour le champ affine orthogonal du test.
- Les faces internes sont strictement conservatives ; les flux de frontière sont comptabilisés dans le bilan global et ne doivent pas être annulés artificiellement.
- Les frontières polyédriques inclinées, la non-orthogonalité et les interfaces à viscosité discontinue nécessitent une campagne séparée.

## 5. Pièges connus
- Omettre le terme transposé revient à résoudre \(\nabla\cdot(\mu\nabla U_i)\) plutôt que la divergence complète de la contrainte newtonienne.
- Supposer que \(\nabla\cdot\mathbf U=0\) annule ce terme est incorrect lorsque \(\mu\) varie spatialement.
- Une contribution de contrainte ajoutée au prédicteur mais absente du résidu final peut faire accepter une solution qui ne satisfait pas l'équation effectivement visée.
- Une correction explicite peut ralentir la convergence ; cela doit être traité par la méthode non linéaire, jamais par un relâchement des critères physiques.

## 6. Ordre de précision attendu et validation exigée
- Champ manufacturé : \(\mathbf U=(a y,0,0)\), donc \(\nabla\cdot\mathbf U=0\), et \(\mu_{eff}(x)=1+x\).
- Oracle indépendant pour la composante y :
\[
\left[\nabla\cdot(\mu_{eff}(\nabla\mathbf U)^T)\right]_y
=\frac{\partial}{\partial x}(\mu_{eff}a)=a.
\]
- Sur le canal orthogonal et le cube hexaédrique incliné, le test compare le flux interne à la formule analytique et la divergence intérieure à aV_c. Les tolérances restent absolues et serrées.
- Qualification restante : MMS de quantité de mouvement sur au moins trois maillages, ordre observé, comparaison tenseur/flux sur maillage incliné ou polyédrique, et vérification d'impact sur une solution Navier–Stokes complète.

## 7. Statut d'implémentation CFDX
- Source : `src/cfdx/physics/steady_incompressible_solver.h`.
- Fonctions : `momentum_transpose_diffusion_face_fluxes`, `momentum_transpose_diffusion_rhs`, `assemble_momentum_component_with_transpose_diffusion`.
- Intégration : prédicteurs x/y/z et réassemblage du résidu final utilisent le wrapper de contrainte transposée.
- Test exécutable : `tests/unit/test_variable_viscosity_transpose_stress.cpp`, avec flux de face affine orthogonal/incliné, divergence cellulaire, conservation des flux internes, contrôle du second membre et MMS lisse sur trois raffinements (n=8/16/32) avec ordre observé.
- Registre : cas `INC-TRANSPOSE-STRESS-MMS`, statut DIAGNOSTIC ; le test unitaire ne vaut pas qualification du solveur complet.
- Gaps : CI exacte sur le HEAD final, GCI formel, maillages polyédriques arbitraires et convergence de l'équation complète.
