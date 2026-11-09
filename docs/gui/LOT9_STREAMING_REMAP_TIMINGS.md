# CFDX Workbench — Lot 9 : surface LOD, remapping conservatif, durées et conditions limites

Date : 2026-11-08 — Base : archive exacte `CFDX_Workbench_LOT8_CODE_COMPLET.zip`.

## Étendue réellement implémentée

| Domaine | Livraison effective | Limite de qualification |
|---|---|---|
| Visualisation | Aperçu de surface HDF5 **échantillonné par faces réelles** ; lecture limitée et asynchrone Qt ; rendu VTK/PyVista de la surface simplifiée | **Pas** de streaming du volume entier, pas de rendu GPU ou multi-GPU de 10–50 M cellules qualifié |
| Remapping | Intégrale du scalaire conservée par calcul exact de volumes d'intersection de boîtes hexaédriques cartésiennes, avec vérifications géométriques et refus explicites | **Pas** de maillages non orthogonaux, polyèdres généraux, ni volumes dépassant 4 096 cellules/2 M couples candidats |
| Temps | `execution.json` : `elapsed_wall_seconds`, horloge monotone `steady_clock` depuis entrée processus ; erreurs comprises | Temps total de processus, **pas** profil CPU/GPU, ni timing noyau par noyau |
| Conditions limites | Cas manufacturé indépendant de conduction Dirichlet 300–400 K, conservations vérifiées ; régressions HDF5 multi-champs et SST maintenues | Pas de nouvelle qualification physique NACA/BFS/RANS/SST/CHT complète |

## Lecture progressive de surface

`cfdx.surface_streaming.read_surface_preview(path, max_faces=5000, max_vertices=30000)` :
- N'ouvre aucun tableau volumique de cellules.
- Ne sélectionne que des **faces de bord originales** (identifiants originaux conservés).
- Contrôle les liens HDF5, les chunks, les offsets, les indices et les coordonnées finies.
- Refuse les allocations supérieures au budget déclaré ; vue 3D sans conservation de champs de volume.
- Qt : un `QRunnable` prépare le sous-ensemble ; la scène VTK n'est modifiée que par le thread graphique et uniquement si le projet n'a pas changé pendant le traitement.
- Le message de confirmation signale la **fraction de faces réellement affichées**.

Commande de diagnostic indépendant :

```sh
PYTHONPATH=src/cfdx/python python -m cfdx.surface_streaming mesh.cfdx.h5 --max-faces 5000
```

Le test de charge de métadonnées à 1 million d'entrées `patch_face_ids` est **synthétique avec des faces répétées**, et démontre seulement les bornes de lecture et d'échantillonnage ; ce n'est pas un maillage physique de 1 million de faces. La lecture d'un cube HDF5 réel CFDX a aussi été vérifiée (6 faces, 8 points).

## Remapping conservatif : champ scalaire intensif

Pour deux ensembles disjoints de boîtes fermées orthogonales à X/Y/Z décrivant un domaine identique,

\[
\bar{\phi}^{B}_{j} = \frac{1}{V^{B}_{j}} \sum_i |C^{A}_i\cap C^{B}_j|\,\bar{\phi}^{A}_{i}.
\]

L'algorithme vérifie toutes les faces et huit coins de chaque hexaèdre, la couverture complète de chaque cellule cible, celle de chaque cellule source, les volumes globaux, les unités, le temps physique et la conservation :

\[
\sum_j V^B_j\bar\phi^B_j = \sum_i V^A_i\bar\phi^A_i.
\]

Le résultat ne modifie **aucun** des fichiers. La comparaison à un champ du maillage cible retourne RMS/L2/Linf, intégrales et erreur de conservation ; `certificate=false` signale qu'il ne s'agit pas d'une qualification physique. Le Workbench affiche ce mode explicitement dans Analyze → Compare ; le comparateur maillage-identique reste inchangé.

```sh
PYTHONPATH=src/cfdx/python python -m cfdx.conservative_remap \
  source.cfdx.h5 destination.cfdx.h5 /fields/scalar/T
```

**Échecs obligatoires** : maillage cartésien déformé, trous, volumes négatifs ou incomplets, champs non finis, unités ou temps incompatibles, lien symbolique, taille trop grande, couverture non conservative. Les polyèdres non-orthogonaux demanderont un algorithme d'intersection robuste ou une bibliothèque géométrique spécialisée ; aucune interpolation non-conservative n'est substituée silencieusement.

## Temps de calcul de production

`ExecutionSummary` possède un champ facultatif `elapsed_wall_seconds`. Le code historique l'omet et écrit `null` ; le solveur de production passe une mesure monotone lors de la première écriture, de la finalisation du bilan et d'une erreur. `compare_run_reports` compare les durées seulement si elles ont été **enregistrées dans les deux exécutions** ; sinon le libellé historique `not_recorded_in_native_reports` est préservé. Une durée négative, infinie ou `NaN` est refusée. Les critères de convergence ne changent pas.

Le smoke natif avec Clang sur le cube CFDX HDF5 a terminé après 2 itérations acceptées et écrit `execution.json` avec une durée strictement positive. Il valide l'export et sa lecture, **pas** la performance de production sur N256.

## Validation et CI

- CTest (g++ Release) : `test_execution_summary`, `test_energy_independent_acceptance`, `test_case_hdf5_io`, `test_production_sst_thermal_coupling`, `test_sst_blending` : succès.
- C++ production `main.cpp` : syntaxe Clang réussie, compilation objet Clang-O0, édition de liens contre `libcfdx_core.a` construite avec GCC et exécution sur maillage natif validées.
- Python : régression ciblée 198 réussis / 3 ignorés à un point de la campagne ; tests supplémentaires pour les nouvelles fonctions, à exécuter sur le HEAD exact livré.
- VTK natif (hors Qt) : les 6 polygones et 8 points d'un cube sont reconstruits exactement.
- Les tests visuels Qt/OpenGL et captures GPU sont intégrés au workflow GUI et **non exécutés localement** faute de PySide6 et PyVistaQt.
- Pas de push GitHub, pas de CI distante, pas de qualification physique NACA/BFS/SST/rayonnement complète.

## Bibliographie et stratégie industrielle

- VTK [`vtkStreamingDemandDrivenPipeline`](https://vtk.org/doc/nightly/html/classvtkStreamingDemandDrivenPipeline.html) : support natif du traitement par morceaux (`UPDATE_PIECE`).
- [VTKHDF](https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_introduction.html) : format HDF5 partitionnable, y compris sorties parallèles.
- [ParaView Remote and Parallel Visualization](https://docs.paraview.org/en/latest/ReferenceManual/parallelDataVisualization.html) : processus de rendu distribués pour les données trop volumineuses pour un poste.
- [OpenFOAM mapFields](https://doc.openfoam.com/2312/tools/post-processing/function-objects/field/mapFields/) : transferts de champs par recouvrements volumiques, stratégie de référence à étendre aux cellules non cartésiennes.

### Verrous de qualification encore ouverts

1. Lecteur CFDX **vraiment partitionné** pour volumes polyédriques, géométrie et champs, format VTKHDF/VTM ou serveur ParaView, gestion du budget GPU et essais physiques 1/10/50 M cellules sur matériel réel.
2. Intersection de polyèdres généraux conservatrice, parallèle CPU/GPU, tests de bilans géométriques sur NACA/BFS et champs discontinus.
3. Raccordement des conditions limites multi-champs à des campagnes physiques indépendantes SST/CHT/buoyance sur géométries représentatives (avec résidus vrais et bilans).
4. Profil CPU/GPU et timing par étape, export simultané avec le commit, l'identité du cas et les preuves CI, sans remplacer les contrôles de convergence.
