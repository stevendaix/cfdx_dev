# VMFL036 — cross-check du solveur principal en 3D

## Décision géométrique

Le VMFL036 de référence est l'écoulement autour d'une **sphère** traité en 2D axisymétrique. Une coupe 2D plane ne peut pas reproduire cette physique : elle représente un cylindre infini, avec une topologie d'écoulement et une traînée de référence différentes. Il ne faut donc jamais comparer son coefficient de traînée à CD=1.0895.

Le cas principal de comparaison doit être **3D autour d'une sphère**, avec :
- D=1 m, R=0.5 m ;
- rho=1 kg/m3, U=1 m/s, mu=0.01 Pa.s, donc Re=100 ;
- une paroi sphérique sans glissement ;
- un domaine lointain sphérique de rayon 50D, identique au rayon extérieur du cas axisymétrique ;
- Aref=pi D²/4 et oracle bibliographique CD=1.0895, soit FD=0.4278456495 N.

La frontière extérieure sphérique reçoit un champ de vitesse uniforme comme condition de champ lointain. C'est un premier cas de croisement numérique, pas une copie exacte des conditions de frontière du projet Fluent privé. Il faudra mesurer la sensibilité à la taille du domaine et aux conditions lointaines avant d'en faire une qualification de production.

## Comparaison du code

### Chemin général cartésien

Le fichier src/cfdx/physics/steady_incompressible_solver.h expose make_rhie_chow_mass_flux(). Il calcule le flux massique de face à partir de la vitesse interpolée et d'une correction de pression, avec des coefficients inverses de diagonale rAU directionnels. Le code commun documente explicitement l'utilisation de dAU=V/AP pour garder cohérents la correction de vitesse, Rhie–Chow et l'opérateur de continuité.

### Chemin axisymétrique VMFL036

Le solveur dédié axisymmetric_vmfl036_solver.h assemble séparément les composantes axiale et radiale, avec volumes de révolution et métriques de face 2*pi*r*ds. Son interpolation de vitesse de face et son opérateur de correction de pression sont propres à cette discrétisation. La formule de flux commune n'est pas automatiquement interchangeable avec cette formulation : la géométrie, les termes radiaux et les coefficients relaxés doivent rester cohérents.

### Hypothèse à tester, pas à présumer

L'audit a relevé dans le chemin axisymétrique une différence apparente entre le coefficient utilisé dans Rhie–Chow et celui de la correction de vitesse :
- correction de vitesse : V/(alpha_u*AP_rel) dans la forme actuellement auditée ;
- terme Rhie–Chow : V/AP_rel.

Cela correspond à un facteur alpha_u entre les deux réponses si AP_rel est bien la diagonale déjà relaxée. C'est un signal d'incohérence potentielle qui doit être dérivé à partir de l'équation de quantité de mouvement et du traitement de Majumdar avant toute modification. Le chemin général utilise un contrat rAU/dAU explicite et mérite la même vérification de bout en bout.

## Génération du maillage 3D

Le script scripts/generate_vmfl036_3d_mesh.py génère un maillage tétraédrique Gmsh 4.1 avec groupes physiques fluid, sphere_wall et farfield.

Exemple de maillage initial, volontairement grossier :

    python scripts/generate_vmfl036_3d_mesh.py --output build/validation/vmfl036-3d-coarse.msh --lc-wall 0.15 --lc-far 5.0

Affiner ensuite progressivement, par exemple lc-wall=0.10, 0.07, 0.05, en contrôlant le nombre de cellules et la capacité mémoire. Ne pas lancer immédiatement le niveau fin : le domaine 3D avec rayon extérieur 50D peut produire un maillage volumique coûteux.

## Plan de qualification

1. **Contrat de maillage** : vérifier que le maillage contient un volume fluide, une paroi sphérique et une frontière lointaine, sans faces non affectées à un groupe.
2. **Compatibilité du solveur principal** : importer le MSH via scripts/meshio_import.py, vérifier la conservation des groupes et inspecter le cas converti .cfdx.h5. Confirmer explicitement les conditions limites du solveur avant l'exécution.
3. **Cas de base 3D** : résoudre le cas laminaire stationnaire Re=100 ; enregistrer convergence momentum/continuité, flux de masse, force de pression et force visqueuse séparément.
4. **Comparaison de formulation** : sur un maillage 3D modeste, comparer le flux Rhie–Chow de production aux flux reconstruits et vérifier l'invariance de la solution convergée en faisant varier alpha_u sans modifier les autres paramètres.
5. **Raffinement** : coarse/medium/fine avec mêmes conditions limites ; comparer CD_p, CD_v, CD, bilan de masse et variation de CD.
6. **Comparaison axisymétrique–3D** : seulement après étude de convergence et sensibilité au domaine lointain. Les contributions de pression et visqueuse doivent être conservées séparément pour éviter qu'une compensation d'erreurs masque un défaut.
7. **Qualification** : n'attribuer le statut PASS qu'après exécution effective, logs retenus et critères quantitatifs pré-déclarés. Le présent ajout ne constitue pas encore un résultat CFD validé.

## Cas plan 2D complémentaire

Un cas plan 2D peut servir à vérifier le solveur principal, mais comme **benchmark distinct** (par exemple écoulement autour d'un cylindre), avec ses propres références CD, son domaine et ses critères. Il ne doit pas être présenté comme une reproduction du VMFL036 sphérique.

## Statut

Ce changement prépare la géométrie 3D et le protocole de comparaison. Il ne modifie ni le solveur, ni l'oracle bibliographique, ni les tolérances VMFL036. L'import, les conditions limites et l'exécution physique doivent encore être qualifiés sur le HEAD exact avant de déclarer le cas opérationnel.
