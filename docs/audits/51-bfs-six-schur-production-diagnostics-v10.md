# CFDX - Lot V10 : BFS turbulent, douze solveurs et verrous observes

Date : 2026-11-09. Source de verite : ZIP original `cfdx_dev-master (2)(1).zip`,
SHA256 `c3a0d40206829596dca61871c10f6079724c6d0e11bbf63a47d4e894e914e62c`.
Le lot V10 inclut integrativement les changements V8 puis V9. Il ne prouve pas
l'equivalence avec un HEAD GitHub. Aucun push n'a ete effectue.

## 1. Verdict et limites

**Verification logicielle ciblee OK ; qualification numerique/physique BFS NON ACQUISE.**
Le binaire C++ `cfdx_production_solver` a ete compile ET lie avec GCC 14.2
sur la base V9, puis recompile avec les marqueurs V10. Le gros `main.cpp`
a ete compile avec `-O0 -DNDEBUG` pour reduire les ressources de compilation ;
les calculs, seuils d'acceptation et protections restent ceux de la source V9.
Il ne s'agit pas d'un benchmark de performance en configuration Release optimisee.

Le solveur a ete execute sur **12 configurations** pour un **BFS six blocs
synthetique de 712 cellules / 2907 faces**, geometriquement inspire par le
`tutorial` OpenFOAM Driver-Seegmiller. Les premiers echec(s) sont reels mais
la grille est une **reconstruction**, pas la maille OpenFOAM issue de
`blockMesh`, et le profil d'entree est uniforme : **pas la reference
experimentale NASA**. Un second maillage du meme generateur a ete construit
(20540 cellules / 82493 faces / 113 faces de sortie), mais son premier pas
n'a pas abouti dans le budget de l'execution locale et a ete interrompu
**sans aucun etat accepte**. Le maillage complet ne fournit donc aucun
resultat de convergence ou de physique certifiable.

Le pilote des essais courts employait volontairement **12 iterations SST
internes** au lieu du budget de production **800**. Les rejets SST a 12
iterations identifient le point atteint par la chaine, mais **n'impliquent
pas** une non-convergence du SST avec 800 iterations. Aucun seuil de residu
n'a ete assoupli. Un timeout correspond a une observation incomplete,
jamais a une divergence prouvee.

## 2. Matrice reelle, premier blocage observe (712 cellules)

| Algorithme / Schur / PC demande | Execution observee | Premier obstacle | Remarques |
|---|---|---|---|
| SIMPLE | FAILURE | PRESSURE | 4.240789e-11 > 1e-11 apres 3 retries, 110 iterations |
| SIMPLEC | FAILURE | PRESSURE | plancher similaire ~4.26e-11 |
| COUPLED / PCD / auto | TIMEOUT_INCOMPLETE | KRYLOV_SCHUR | premier echec reel status=2, 52 iterations, 1.055851e-2 ; retries en cours a l'arret |
| COUPLED / Block-local / block-AMG | FAILURE | SST | COUPLED convergent, SST rejete au budget de 12 |
| COUPLED / LSC / block-AMG | TIMEOUT_INCOMPLETE | KRYLOV_SCHUR_RUNNING_TIMEOUT | pas de statut definitif avant 6 secondes |
| COUPLED / BFBt / block-AMG | TIMEOUT_INCOMPLETE | KRYLOV_SCHUR_RUNNING_TIMEOUT | pas de statut definitif avant 6 secondes |
| COUPLED / SIMPLE-Schur / block-AMG | FAILURE | SST | COUPLED convergent, SST rejete a 12 |
| COUPLED / SIMPLEC-Schur / block-AMG | FAILURE | SST | COUPLED convergent, SST rejete a 12 |
| COUPLED / PCD / pcd | TIMEOUT_INCOMPLETE | KRYLOV_SCHUR | premier echec reel a 52 iterations, retries en cours |
| COUPLED / PCD / block-AMG | TIMEOUT_INCOMPLETE | KRYLOV_SCHUR | premier echec reel a 52 iterations, retries en cours |
| COUPLED / Block-local / MGR | FAILURE | SST | COUPLED convergent, SST rejete a 12 |
| COUPLED / Block-local / Braess-Sarazin | FAILURE | SST | COUPLED convergent, SST rejete a 12 |

Les configurations sont marquees NOT_QUALIFIED par le certificat V8 pour tous
les cas. La variation du maillage **physique** est interdite entre options :
le generateur compare les empreintes de topologie/geometrie, mais les fichiers
HDF5 distincts ont naturellement des empreintes SHA256 differentes a cause
des parametres solveur. Les journaux, statuts et checksums sont persistants.

## 3. Nouvelle preuve : PC reellement resolu en COUPLED

La compilation finale du traceur C++ V10 ajoute trois champs a
`BFS_PHASE_V10 phase=COUPLED_LINEAR_END` :

- `actual_krylov` : schema exterieur effectivement resolu,
- `actual_preconditioner` : preconditionneur exterieur reel,
- `resolved_schur` : sous-modele de Schur resolu dans le solveur.

Le rejeu de PCD/auto montre :
`actual_krylov=fgmres actual_preconditioner=coupled_block_schur
resolved_schur=pcd` ; au premier echec `status=2, iterations=52,
true_relative=0.0105585`. Il ne faut pas confondre **Schur PCD** avec
**preconditionneur exterieur PCD**. Cette preuve provient du plan resolu,
non du texte de la configuration demandee. Le marqueur est observationnel
et n'autorise jamais l'acceptation d'une iteration.

## 4. Pression : les seules policies lineaires ne ferment pas le defaut

Pour SIMPLE, la correction de pression atteint 4.2407891218e-11 contre un
seuil impose de 1e-11, malgre trois retries. Les essais A/B sur cinq policies
(baseline, equilibrage, raffinement iteratif, GMRES, AMG) n'ont pas permis
de fermer la cible. Le diagnostic purement observationnel
`PRESSURE_GHOST_ROW_SUM_DIAGNOSTIC` retrouve sur la grille de 712 cellules
684 lignes quasi-Neumann FP64 et 28 lignes ancrees par la pression de sortie,
avec fuite numerique relative ~3.6e-11 a la derniere tentative. Ce resultat
est **compatible avec** un plancher d'arrondi / incoherence du residu vrai,
mais il ne demontre pas a lui seul une cause unique. Il faut examiner les
coefficients pression, les flux aux faces Rhie-Chow, les corrections
non orthogonales et la representativite du residu physique dans le solveur
sans relacher le gate 1e-11.

## 5. Convergence physique et comparaison NASA : explicitement refusees

Le validateur V8 de `convergence.json` est execute **avant toute mesure** de
reattachement, Cf ou profils. Une execution non convergee, une preuve absente,
un plan resolu incompatible, un U,p,k,omega non fini ou incomplet, ou une
acceptation d'ordre un refuse le calcul de QoI.

Exemple reel, essai SIMPLE : le fichier `convergence.json` n'existe pas car le
solveur rejette le pas de pression. Le nouveau script V10 renvoie
`GLOBAL_U_P_K_OMEGA_NOT_CONVERGED_NO_QOI_COMPUTED` et exit status=2 ;
`comparison=null`, `physical_validation=NOT_QUALIFIED` dans le JSON.

La bibliographie publique OpenFOAM et NASA TMR expose :

- `Re_h ~ 36000`, `U_entree=44.2 m/s`, `h=0.0127 m` ;
- `x_r/h = 6.26 +/- 0.10` experiment (Driver-Seegmiller) ;
- profils U et contraintes turbulentes, Cf et Cp (donnees tabulees NASA) ;
- plusieurs maillages dans NASA TMR, necessaires a la convergence spatiale.

Le postprocesseur V10 prepare une **extraction diagnostique** de Cf par
cisaillement de la paroi horizontale aval, la premiere inversion Cf negative
vers positive pour un `x_r/h` provisoire et les profils de vitesse aux
stations x/h=1,4,6,10. Les valeurs de paroi de CFDX sont seulement
`first_order_one_sided_diagnostic`, **pas le cisaillement reconstruit de
qualite qualification**. Le fichier des contraintes de Reynolds et le
Cp sur vraie face de paroi ne sont pas encore raccordes a l'export de
production. **Aucun accord chiffre avec la NASA n'est revendique**.
La condition de paroi SST CFDX n'est pas encore equivalente aux fonctions
paroi conjointes `nutUSpaldingWallFunction`, `omegaWallFunction` et la
production G d'OpenFOAM. La comparaison des methodes ne peut pas contourner
ce verrou de modele/BC.

References :

- https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/
- https://tmbwg.github.io/turbmodels/backstep_val.html
- https://tmbwg.github.io/turbmodels/backstep_grids.html

## 6. Modifications V10, base du code de production

1. `src/cfdx/physics/steady_incompressible_solver.h` : traceur par phase,
   option `CFDX_BFS_STAGE_TRACE`, incluant methode et PC resolus (pas seulement demandes).
2. `scripts/run_cfdx_coupled_campaign.py` : execution sous trace,
   empreinte SHA256 des journaux, classification du premier obstacle pour
   chaque methode et maintenance du certificat V8 fail-closed.
3. `scripts/diagnose_bfs_failure_v10.py` : parseur autonome des logs
   qui distingue phase en cours d'un echec numrique reel et lit le plan resolu.
4. `scripts/qualify_bfs_driver_seegmiller_v10.py` : refus strict de calculer
   des observables physiques lorsque le U,p,k,omega independant n'est pas ferme.
5. `scripts/generate_th03_openfoam_blockmesh.py` : metadata d'ordre de
   convection conforme au schma demande, sans fausse mention d'upwind.
6. `tests/python/test_bfs_v10_diagnostics.py` : tests des echecs et timeouts,
   de la methode reelle et des extractions de paroi, sans faux positif physique.

Ce lot est une **amelioration du diagnostic**, pas une reformulation des
operateurs de pression, du preconditionnement ou des equations SST.

## 7. Preuves testees

- `cfdx_production_solver` : compile et lie en C++ avec GCC 14.2 et `-O0`
  pour `main.cpp`, lien executable OK apres la derniere modification.
- `pytest` BFS/Schur/generator/diagnostic : **90 PASSED, 6 SKIPPED**.
- `ctest -R '^test_bfs_quick$'` : **2/2 PASS** (generation dependance + qualif).
- 12 configurations de BFS sur la meme grille : toutes lancees, aucune
  convergee physiquement ; observation des timeouts borne a 6 secondes.
- Test du gate experimental avec essai rejete : exit=2 et comparison=null.
- Maillage 20540 cellules genere : premier solve non termine ; aucune
  valeur d'essai complete, pas de validation par extrapolation.
- CTest complet, NACA N64/N128/N256, N9 long, sanitizer complet, MPI/GPU :
  **non testes dans ce lot V10**.

## 8. Reproduire

```bash
# Prerequis pour le build : dependances de CFDX, GCC, CMake, HDF5, Python h5py/numpy
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_ENABLE_MPI=OFF
cmake --build build --target cfdx_production_solver -j 1

# Les douze configurations sur un MEME cas HDF5 externe identifie
PYTHONPATH=.:src/cfdx/python python scripts/run_cfdx_coupled_campaign.py \
  --build build --only bfs --bfs-reference-case data/reference_bfs.h5 \
  --out runs/bfs_reference --iterations 1500

# Qualification refusee si la preuve globale manque
PYTHONPATH=.:src/cfdx/python python scripts/qualify_bfs_driver_seegmiller_v10.py \
  --run runs/bfs_reference/bfs_runs/bfs_pcd_auto \
  --case runs/bfs_reference/bfs_cases/bfs_pcd_auto.cfdx.h5 \
  --algorithm coupled --schur pcd --preconditioner auto \
  --exit-code 0 --density 1.2 --u-ref 44.2 \
  --output runs/bfs_reference/driver_seegmiller_pcd_auto.json
```

Le parametre `--exit-code 0` ci-dessus n'est pas une preuve a lui seul : le
validateur exige les fichiers persistants de cette meme execution et le
certificat independant V8 complet. Pour les campagnes longues, conserver un
budget SST normal et le schema d'ordre deux sur au moins trois maillages
NASA/Driver-Seegmiller adaptes ; effectuer GCI, Cf, Cp, vitesse, contraintes
Reynolds et reattachement **uniquement apres convergence globale**.

## 9. Suite priorisee pour effectivement debloquer la physique

- P0 Pression : comparer matrice CSR et divergence de flux Rhie-Chow pour la
  meme correction et le meme etat physique, avec normalisations et erreur
  d'arrondi documentees. Aucune baisse du seuil n'est acceptable.
- P0 PCD : enregistrer les normes vraies du bloc exterieur et du solve Fp
  imbrique, puis comparer sur les memes captures A,b; eviter les changements
  silencieux de Schur ou de PC.
- P0 SST : sur le budget standard 800, tracer residuals k et omega, flux,
  bilan source/puits, production G, viscosite turbulente en paroi et normes
  independantes; calibrer conjointement les fonctions de paroi.
- P1 Geometrie/conditions aux limites : importer le vrai `blockMesh` et le
  profil experimental de couche limite, conserver une identification
  exhaustive de mesh/BC et des dimensions.
- P1 Qualification : obtenir les 6 Schur converges sur reference, puis
  calculer Cf et xr/h, extraction tension de Reynolds, Cp par face et profiles,
  comparaison aux fichiers NASA, et GCI sur trois maillages.

**Critere final : aucun resultat n'est qualifie avant la preuve simultanee
linear vrai + residu global U,p,k,omega + conservations + reference physique
+ maillage.**
