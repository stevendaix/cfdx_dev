# CFDX V18 — contraintes omega des cellules de paroi et qualification BFS

**Date : 10 novembre 2026. Statut : correction logicielle ciblee verifiee ; modele de paroi experimental NON QUALIFIE physiquement.**

## Base et perimetre

Base : archive complete CFDX V17 fournie par le projet ; patch cumulatif calcule contre `cfdx_dev-master (2)(1).zip` (empreinte SHA-256 `c3a0d40206829596dca61871c10f6079724c6d0e11bbf63a47d4e894e914e62c`). Seuils numeriques et physiques inchanges. Aucun push distant.

Modifications :

1. `sst_solver.h` : contrainte atomique des lignes proprietaires de l'equation omega, avec verification des diagonales, des tailles et des donnees finies ; residualisation sur EXACTEMENT le meme operateur contraint, avec distinction entre bilan PDE et ligne contrainte.
2. `sst_wall_operator_v17.h` et `steady_incompressible_solver.h` : cible omega par cellule construite a partir du triplet paroi V17, ponderation des coins et activation exclusivement opt-in.
3. `main.cpp` et `generate_th03_openfoam_blockmesh.py` : option `v18_omega_wall_owner_constraints_experimental=enabled`, verifications des preconditions de paroi et preservation du chemin V17 inchange sans option.
4. Tests de production `test_sst_wall_owner_v18.cpp`, tests du generateur et enregistrement CTest dans CMake.

## Contrat mathematique et limites

La condition `omegaWallFunction` d'OpenFOAM modifie les champs omega et G et manipule les equations des cellules proprietaires ; une simple valeur `fixedValue` sur la face ne reproduit pas cette operation. V18 essaye une version **Dirichlet pleine** de la contrainte de cellule : pour chaque cellule proprietaire i, A[i,:]=0, A[i,i]=1, b[i]=omega_wall[i]. Les autres lignes ne sont pas modifiees et sont evaluees dans le bilan SST. Les cibles de coin sont consolidees par ponderation surfacique dans l'operateur de paroi existant. Ce traitement n'est PAS encore une reproduction exacte de tous les poids, des fonctions de lissage et des schemas temporels d'OpenFOAM ; il reste des differences de numerique et d'inlet. Le code est experimental, desactive par defaut.

L'indicateur de convergence de l'omega contraint est recalcule en reassemblant les memes lignes, pas le PDE non contraint. Il faut aussi surveiller les bilans PDE non contraints, la positivite k/omega, les valeurs y+, l'amplitude nu_t/nu, le rapport Pk/Dk et les flux conservatifs. Les residus vrais des equations et la conservation restent des portes independantes de toute comparaison physique.

## Validation logicielle realisee

- Compilation et edition de liens de `cfdx_production_solver` : PASS en Release de diagnostic (-DNDEBUG, -O0 sur l'unite lourde). Ce n'est pas un benchmark optimise.
- C++ : `test_sst_wall_owner_v18` 4/4 et `test_sst_wall_operator_v17` 4/4, soit 8/8 tests de cas.
- Python : 38/38 tests selectionnes du generateur, des contrats de paroi, des certificats globaux et des comparaisons NASA.
- La desactivation de l'option preserve la matrice et les valeurs historiques ; option sans triplet V17 correctement refusee.

## BFS reconstruit, 712 cellules — vrai solveur de production

HDF5 identique pour les deux modes a l'exception du drapeau de contrainte proprietaire. Premier essai V18 : transport SST stagne vers l'iteration interne 226 et rejette l'iteration externe. Valeurs mesurees : `omega_returned_nonlinear_raw=6.438040145440e-08`, `omega_returned_row_scaled=1.585159911771e-06`, cible imposee `1e-8`, pires cellules locales omega=77 et k=600. Aucune iteration externe acceptée, aucune contraction sur trois etats. Les cellules proprietaires de paroi contraintes sont au nombre de 69 ; le plus grand ecart omega sur ces cellules a l'iteration interne 200 est de `2.658634912223e-08`.

Essai de relaxation lineaire SST a 1.0 : stagnation FP64, omega brut `5.965398486296e-08`, pas de qualification. Les essais temporellement bornes sont des essais incomplets, pas des preuves de divergence asymptotique. La cellule 77 n'est pas une cellule de paroi. A l'iteration interne 100, son residu omega recompose est `1.687435295139e-05` et son erreur de decomposition `-2.502846833865e-14` ; les contributions convection/diffusion/production/destruction et diffusion croisee sont presentes. Il faut donc etudier les flux et la linearisation interiorieure sans inventer un gain paroi.

## BFS reconstruit, 20 540 cellules

Maillage six blocs genere a partir des parametres de geometrie OpenFOAM, **pas un fichier `blockMesh` original**. Le traceur constate un premier blocage de correction de pression SIMPLE avec residu relatif vrai `5.44746e-10` contre un objectif `1e-11`. Le calcul a ete interrompu dans le budget de diagnostic avant SST. Rien ne permet de qualifier les Schur PCD/LSC/BFBt/MGR/AMG sur cette grille a partir de cet essai SIMPLE. Les campagnes sur ces configurations restent a executer.

## Porte stricte de convergence et comparaison experimentale

Le validateur `qualify_bfs_global_contraction_v14.py` donne `strict_contraction_proved=false` (zero etat externe accepte au lieu de trois). `qualify_bfs_driver_seegmiller_v10.py` refuse de calculer le reattachement : `GLOBAL_U_P_K_OMEGA_NOT_CONVERGED_NO_QOI_COMPUTED`. La solution reste **NOT_QUALIFIED**, sans valeur calculee pour xr/h, Cf, Cp, profils et contraintes de Reynolds. La reference experimentale Driver–Seegmiller est `xr/h=6.26 +/- 0.10` (NASA TMR). Une comparaison valide exigera une grille comparable, les conditions amont, le schema SST et les fonctions de paroi de production qualifies, puis trois maillages converges et les bilans physiques. Conversion de la pression modifiee CFDX : `p_static = p_star - 2*rho*k/3` avant toute comparaison de Cp, avec une origine de pression experimentale explicite.

**Suite recommandee :** reduire le residu omega en cellule interieure 77 par une veritable coherence de la linearisation et des flux, sans faux succes par normalisation ni abaissement de tolerance ; demontrer ensuite une contraction physique U,p,k,omega sur trois etats acceptes. Sur le maillage 20k, fermer d'abord le correcteur pression puis faire tourner les 12 configurations BFS avec des preuves par Schur/preconditionneur. Refermer l'equivalence omega/G/nut au niveau des poids et du traitement des coins OpenFOAM et mesurer y+ et les contraintes de Reynolds avant la promotion de Driver–Seegmiller.

## Reproduction et preuves

```
python -m pytest -q tests/python/test_bfs_omega_owner_v18.py tests/python/test_th03_sst_numerics_v13.py tests/python/test_th03_operator_ab_v7.py tests/python/test_qualify_bfs_global_contraction_v14.py tests/python/test_qualify_bfs_nasa_fields_v16.py tests/python/test_audit_bfs_unqualified_fields_v16.py tests/python/test_bfs_kqr_v16.py
python scripts/qualify_bfs_global_contraction_v14.py --log owner_v18.log --output contraction_v18.json
```

Voir le ZIP de preuves V18 pour les HDF5 de diagnostic, journaux du solveur, tests et rapports des portes de qualification.

## Sources de methode

- Code OpenFOAM v11, `omegaWallFunctionFvPatchScalarField.C` : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html
- Documentation OpenFOAM (y compris `kqRWallFunction`) : https://doc.openfoam.com/2606/tools/processing/models/turbulence/ras/wall-functions/
- NASA TMR 2D Backward-Facing Step : https://tmbwg.github.io/turbmodels/backstep_val.html
