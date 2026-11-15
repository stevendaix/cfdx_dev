@@ -0,0 +1,118 @@
# CFDX V64 — diagnostic causal du BFS152, sans nouveau solveur

13 novembre 2026. Base : ZIP de sources V63 livré dans cette conversation.
Cumulatif demandé : archive d'origine `cfdx_dev-physics-bfs-v6-hilo-rhie-chow-production (1).zip`.

## Décision technique

V63 a montré une résolution linéaire 4N très précise avec `selfp`, sans démontrer de convergence physique du BFS. V64 ne crée **aucun** Jacobien, préconditionneur, Schur, modèle SST ou schéma convectif nouveau. Il ajoute des observateurs facultatifs et un banc de preuve, puis utilise les solveurs existants. L'objectif est de séparer l'effet des deux mises à jour non linéaires et d'expliquer le ralentissement du COUPLED/BlockLocal.

## 1. Quatre états sur LE MÊME résidu physique

`src/cfdx/physics/bfs_causal_ablation_v64.h` évalue avec le même callback V51 :

- `F00 = Fphys(U_old,p_old,k_old,omega_old)` ;
- `F10 = Fphys(U_new,p_new,k_old,omega_old)` ;
- `F01 = Fphys(U_old,p_old,k_new,omega_new)` ;
- `F11 = Fphys(U_new,p_new,k_new,omega_new)`.

`U_old,p_old,k_old,omega_old` sont pris avant l'itération extérieure ; `U_new,p_new,k_new,omega_new` sont relevés après l'étape COUPLED et le transport SST, **avant** le certificat final. Le callback réassemble indépendamment les coefficients, la viscosité, les termes de paroi et le flux Rhie–Chow pour chacun de ces états. `F11` est recontrôlé contre l'évaluation déjà réalisée. Un résultat dépendant de l'historique est rejeté.

L'identité mixte mesurée est `I = F11 - F10 - F01 + F00`. Les six blocs `Ux,Uy,Uz,div(phi),k,omega` sont analysés dans leur unité **brute signée**. On mesure aussi les variations de `phi` et de `nu_t` pour chaque effet, ainsi que l'interaction de ces quantités, les cellules et faces les plus affectées. L'oracle exige un RC fermé à `1e-10` sur tous les états. L'identité de recomposition seule n'est évidemment pas une preuve CFD.

**Limite fondamentale :** il ne s'agit PAS de trois solveurs distincts qui auraient chacun convergé avec certains champs gelés. Ce sont des contre-factuels *réassemblés*, utilisés pour localiser les effets des mises à jour. Aucune variable ni tolérance de production n'est modifiée.

Activation bornée : `CFDX_FPHYS_CANDIDATE_V51=1 CFDX_FPHYS_CAUSAL_V64=1`. Sans `CFDX_FPHYS_CAUSAL_V64_EVERY=1`, seul le premier pas est instrumenté. Au-delà de 152 cellules, refus explicite. Le mode n'est jamais activé par défaut.

## 2. Rejeu natif à entrée identique, V64 ON/OFF

- Pilote : `TH03_debug_only`, BFS152, 634 faces, SOU, COUPLED/BlockLocal, SST expérimental. Entrée HDF5 exacte de la livraison V63.
- Le banc `scripts/run_bfs_causal_v64.py` appelle le même binaire, avec le **même SHA-256 du fichier d'entrée** ; seule la variable d'environnement qui active le diagnostic change.
- Les 9 traces physiques de référence comparables sont **textuellement identiques** (V51, V8/V13, SST et momentum post-SST). Les deux exécutions terminent non convergées (code de sortie 1).
- Le parseur `scripts/analyze_bfs_causal_v64.py` refuse 6 blocs incomplets, fermeture RC insuffisante, faux commit ou chiffre non fini ; il conserve `physical_qualification=false`.

Premier pas :

| Signal | Contrôle mesuré |
|---|---:|
| SST transport interne | 196 itérations ; k 1.71336e-11, omega 1.53552e-11 |
| COUPLED 4N | FGMRES vrai 5.20123e-13, 40 itérations |
| RC F00 / F10 / F01 / F11 | 0 / 4.60344e-12 / 0 / 2.36114e-12 |
| Interaction brute Ux | 3.36407e-2 |
| Interaction brute k | 5.36770e-1 |
| Interaction brute omega | 8.01162e+1 |
| Flux RC : effet U/p ; effet SST ; interaction | 3.53546e-3 ; 0 ; 1.05648e-3 kg/s |
| Viscosité turbulente : effet U/p ; effet SST ; interaction | 1.78372e-2 ; 1.61252e-2 ; 1.60421e-2 m2/s |
| Contrôle momentum V7 avant / après SST | 1.346709 / 1.898668 |
| Certificat global V13 | non convergé |

Le `phi_sst=0` vaut **uniquement au point U/p initial de ce premier essai**. Le flux de face varie lorsqu'on combine les deux mises à jour : aucune indépendance générale SST/pression ne peut être déduite.

## 3. Suite naturelle : 12 puis 40 itérations COUPLED existantes

Seule la limite d'itérations dans `numerics.max_iterations` des cas HDF5 de **diagnostic** est modifiée, d'abord de 1 à 3, puis à 12 et 40. Géométrie, BCs, discrétisation, paramètres SST, choix Schur et seuils restent inchangés. L'audit à quatre états a été activé sur chaque pas des 12 itérations ; à 40, les traces V13/V8/V7 déjà présentes suffisent.

| Itération | U après SST (V13) | Continuité V13 |
|---:|---:|---:|
| 1 | 1.898668 | 1.266899e-2 |
| 5 | 8.995074e-2 | 3.928894e-3 |
| 10 | 4.376415e-2 | 1.322049e-3 |
| 15 | 4.96188e-3 | 3.4559e-4 |
| 20 | 5.44345e-3 | 1.7397e-4 |
| 25 | 2.21673e-3 | 9.35e-5 |
| 30 | 2.87137e-3 | 1.0852e-4 |
| 35 | 4.66857e-3 | 5.121e-5 |
| 40 | 2.25539e-3 | 4.12874e-5 |

Observations indépendantes de la différentiation Newton :

- Les quarante résolutions COUPLED rapportent des résidus linéaires vrais `<= 9.859e-13`.
- Les quarante étapes SST rapportent `converged=1`, avec k `<=4.252e-11`, omega `<=3.716e-11`.
- Le momentum global **diminue d'un facteur ~842** entre la première et la quarantième itération, mais **remonte à 18 transitions sur 39**. Le meilleur momentum est `0.00113515` à l'itération 31, donc une partie du progrès est perdue ensuite.
- La variation relative de `nu_t` après SST, `0.964` au pas 1, vaut `0.006724` au pas 40. La forte variation de viscosité initiale ne peut pas, à elle seule, expliquer toutes les oscillations tardives.
- Le solveur ne satisfait pas les contrôles physiques finaux et ne fournit pas d'accord NASA/OpenFOAM ; il rapporte `Converged NO`.

Le parseur `scripts/audit_bfs_committed_v64.py` prend **seulement** les enregistrements `GLOBAL_NONLINEAR_COMMITTED_V13`, avec contrôles de continuité et des phases linéaires/SST, sans inventer de PASS.

## 4. Comparaison avec une globalisation DEJA présente (V13/V14/V47)

Essai indépendant, même maillage/modèle/BC/SOU et limite de 12 itérations, en activant les drapeaux expérimentaux déjà implémentés :

```
v13_joint_globalization=true
v14_monotone_global_merit=true
coupled_composite_backtracking_v47_experimental=true
coupled_composite_backtracking_retries_v47=10
```

V47 refuse l'itération extérieure 8 après 11 essais (`beta` de 1 à 9.765625e-4) : la croissance de l'équation momentum V13 reste supérieure à 1, la dernière valeur du facteur étant environ `1.91425`. La voie historique ne donne donc **pas** une convergence nouvelle en activant simplement le backtracking.

La trace dédiée `COUPLED_ZERO_STEP_FLUX_V47` documente au début du pas 8 un flux Rhie–Chow reconstitué sur l'ancien U/p avec un `rAU` nouvellement assemblé : la continuité du dernier état commis vaut `1.34992e-3`, mais celle de cette recomposition vaut `1.59056e-3`, avec saut de divergence `1.12945e-3`. La ligne de base du mérite peut donc changer **même lorsque beta approche de zéro**. Ce diagnostic correspond au problème déjà identifié en V47/V48 ; il ne suffit pas à démontrer que ce saut est l'unique cause de l'échec à l'itération 8.

La conséquence pratique : **ne pas ajouter un autre Schur/Jacobien**. Pour qualifier une globalisation, comparer le candidat `Fphys(X+beta d)` à `Fphys(X)` sous la **même définition d'opérateur état-only**, avec mise à jour transactionnelle des flux, coefficients et SST et conservation de la même discrétisation. Ne pas comparer une histoire V13 construite avec un ancien `rAU` au nouveau résidu issu de l'opérateur réassemblé.

## 5. Suite ordonnée (sans ajouter de méthode)

1. Instrumenter la fermeture des variables de face et du `rAU` **au pas beta=0** à toutes les itérations où V47 demande un rejet, en comparant aussi le véritable V51 `Fphys(old)` et `Fphys(trial)` pour éliminer le décalage d'opérateur.
2. Réutiliser les protections transactionnelles V13/V14/V47/V48 existantes ; seulement après attribution causale, corriger la référence de mérite ou le raccordement RC, sans assouplir aucun seuil.
3. Retester sur le même BFS152 et exiger plusieurs réductions du mérite global *réévalué* ; contrôle des six sous-blocs U/p/k/omega et NACA en non-régression.
4. Ensuite BFS712 et BFS20 540 avec véritables fonctions de paroi SST et références externes. Le pilote actuel utilise une diffusion moléculaire de paroi **legacy**, donc il n'est pas une qualification Driver–Seegmiller.

## 6. Reproduction

```
cmake -S . -B build-v64 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_BUILD_PYTHON=OFF -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v64 --target cfdx_production_solver \
  test_bfs_causal_ablation_v64
ctest --test-dir build-v64 -R 'test_bfs_(causal|committed).*v64' --output-on-failure
python3 scripts/run_bfs_causal_v64.py \
  --solver build-v64/cfdx_production_solver \
  --case /CHEMIN/bfs152_block_local_v63.cfdx.h5 \
  --out /tmp/cfdx-v64-replay --iterations 1
python3 scripts/audit_bfs_committed_v64.py \
  --log /CHEMIN/trace_40_iterations.log --output /tmp/committed-v64.json
```

**Qualification limitée :** tests ciblés et CPU diagnostique, sans compilation CI intégrale, sans ASAN de tout le noyau, sans GPU/MPI, sans campagne RANS industrielle ni rayon de maillage. Aucun push GitHub. Le rayonnement V57 demeure intact dans l'arbre.
