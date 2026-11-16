# CFDX V60 — fusion radiative V57 + BFS/SST V59 et audit du verrou P0

**Date : 12 novembre 2026.** Base : archive `cfdx_dev-physics-bfs-v6-hilo-rhie-chow-production (1).zip` fournie. Ce document décrit les sources effectives et les exécutions du lot V60, pas une qualification industrielle.

## Décision et état d'avancement

La branche expérimentale CPU de transport radiatif HO angulaire en `hi+lo`, correction HO–LO et SMM/AMG (V57) est maintenant **dans le même arbre source** que le solveur BFS/SST V59. La fusion conserve les fonctions de flux, la quadrature discrète, la fermeture des parois et l'oracle angulaire; elle apporte les dépendances de recherche V38–V56. Elle ne branche **pas** automatiquement l'algorithme HO–SMM dans le DOM/CHT de production. La compilation de recherche reste **OFF par défaut** (`CFDX_ENABLE_RADIATION_RESEARCH_V57`). Les trois tests ciblés de rayonnement préexistant ont réussi après la fusion; cela ne démontre pas la parité de toutes les campagnes historiques ni la qualification d'un solveur industriel.

Le travail P0 ajoute un **audit de six évaluations physiques reassemblées** sur deux véritables directions Krylov préconditionnées et leur somme. Il décompose la non-additivité entre la convection et le reste du résidu, tout en exposant les observables de Rhie–Chow, `rAU`, viscosité, taux de déformation et SST. Un nouveau garde **fail-closed** bloque toute acceptation expérimentale si un bloc du reste, le flux Rhie–Chow ou la convection dépasse le seuil *algébrique existant* de `1e-6`. L'opérateur physique et les tolérances de convergence n'ont **pas** été modifiés. Le Jacobien physique global n'est donc **pas corrigé** : il est mieux localisé et son absence de qualification est explicitement protégée.

### Matrice de qualification honnête

| Objet | Statut V60 et preuve | Ce que cela **ne** prouve pas |
|---|---|---|
| Fusion source V57 + V59 | Dépendances intégrées et tests ciblés compilés | HO/DOM/CHT en production |
| HO angulaire `hi+lo` et SMM/AMG | 5 profondeurs CPU de 640 cellules convergées individuellement | MPI/GPU, maillages complexes, validation physique indépendante |
| Correctifs préexistants du DOM | 3/3 tests production ciblés | Toute la CI, non-régression exhaustive |
| V60 différentiel | Test V60, extraction de 6 résidus physiques sur BFS152 réel | Jacobien `6N` additif ou convergent |
| Newton sur BFS152 | 0 pas accepté, 0 convergence physique | BFS712/BFS20540 |
| Six Schur | Comparaison à matrice commune `A,b,x0` vérifiée par SHA-256 | Qualification des six solveurs ou de la turbulence |

## 1. Fusion des sources : provenance et séparation des responsabilités

1. Base BFS V6 : archive utilisateur ; base V59 déjà reconstruite depuis ce ZIP et vérifiée séparément dans la livraison antérieure.
2. V57 radiatif : sept fichiers extraits depuis `CFDX_V57_RAPPORT_CODE_INTEGRAL(1).md`. **Six sont identiques octet pour octet aux empreintes SHA-256 V57 annoncées ; le septième, `scripts/run_radiation_v57_qualification.sh`, est volontairement adapté à la nouvelle option CMake `CFDX_ENABLE_RADIATION_RESEARCH_V57=ON` et n'a donc pas la même empreinte.** Les fichiers : `docs/CFDX_V57_HOLO_HILO_ANGULAR_AMG_QUALIFICATION.md`, `scripts/run_radiation_v57_qualification.sh`, `src/cfdx/physics/dom_angular_hilo_v57.h`, `src/cfdx/physics/dom_holo_hilo_v57.h`, `src/cfdx/physics/dom_moment_schur_v44.h`, `tests/validation/test_dom_angular_hilo_v57.cpp` et `tests/validation/test_dom_holo_hilo_v57.cpp`.
3. Chaîne V38–V56 requise par V57 : modules `dom_*`, `p1_marshak_v39.h`, CHT/P1, tenseur thermique, maillages de test, vérifications antérieures et scripts dédiés réintroduits depuis le patch cumulatif V57 fourni ; en cas de collision avec BFS V59, fusion ciblée des fonctions partagées de `radiation_solver.h` et `radiation_advanced.h` (pas de remplacement aveugle des fichiers).
4. La bibliothèque `cfdx_core` et le binaire `cfdx_production_solver` restent **communs**. Les études radiatives nouvelles restent un ensemble de tests CMake optionnels, sans activation automatique du DOM HO.
5. Le CTest V54 indépendant `--stress-640` est réservé à `CFDX_ENABLE_LONG_VALIDATION=ON` ; il n'est pas assimilé à une régression courte. Une tentative V60 de ce stress a été interrompue par la limite d'exécution sans verdict numérique ; la qualification des tests courts l'exclut.

## 2. Expérience radiative CPU intégrée

La qualification propre à V57 est sur 8×8×10 cellules et 64 directions ordinées, pour des profondeurs comprises entre 0,060 et 0,003 m. Rejeux V60 :

| Profondeur (m) | Cellules | Acceptations | Défaut HO final | Défaut LO final | Statut |
|---:|---:|---:|---:|---:|---|
| 0.060 | 640 | 17 | 1.81915e-10 | 2.85174e-10 | PASS individuel |
| 0.030 | 640 | 17 | 2.38761e-10 | 4.80093e-10 | PASS individuel |
| 0.012 | 640 | 18 | 1.66979e-10 | 5.22588e-10 | PASS individuel |
| 0.006 | 640 | 18 | 1.44992e-10 | 5.75007e-10 | PASS individuel |
| 0.003 | 640 | 18 | 1.02962e-10 | 5.26603e-10 | PASS individuel |

**Nuance de protocole** : la commande atomique `--family` a été interrompue par la limite de temps après les quatre premières géométries ; la cinquième (`--extreme`) a été relancée séparément avec code de retour 0. Il serait inexact de présenter un seul `--family` réussi sans interruption. La profondeur extrême a un défaut linéaire `hi+lo` de 1.13113e-14. Ces chiffres ne qualifient ni le CHT–DOM de production ni une performance industrielle.

## 3. Jacobien physique global : nouveau diagnostic V60

Pour un état initial BFS152 identique, avec BFBt et transport SST spatial, le solveur évalue le résidu physique complet six fois : `F(X+h v)`, `F(X-h v)`, `F(X+h w)`, `F(X-h w)`, `F(X+h(v+w))`, `F(X-h(v+w))`. Chaque évaluation regénère le flux de masse Rhie–Chow, les coefficients et l'opérateur SST du candidat. Les champs sont ceux de deux directions Krylov **réelles**, pas une perturbation synthétique. Pour chaque quantité Q, l'audit compare `D_hQ[v+w]` à `D_hQ[v]+D_hQ[w]` ; les six blocs de `F - convection_primal` sont comparés séparément, sans réutiliser une matrice Picard.

Résultats mesurés en option V60 (valeurs relatives, chaque observable a sa propre normalisation) :

| Observable | Non-additivité relative | Test `≤1e-6` |
|---|---:|---|
| Convection reassemblée | 9.8250454e-05 | Échec |
| Reste du résidu : globale | 2.8080728e-07 | Succès global seul |
| Reste : max des 3 blocs momentum | 2.6312070e-06 | Échec |
| Reste : continuité | 2.1968912e-06 | Échec |
| Reste : `k` | 1.3040442e-06 | Échec |
| Reste : `omega` | 2.8048853e-07 | Succès |
| Flux de faces Rhie–Chow | 1.5160250e-06 | Échec |
| `rAU` | 4.0930691e-07 | Succès |
| Viscosité effective | 5.3855663e-07 | Succès |
| Taux de déformation | 1.0530929e-06 | Échec d'observable |
| Fonction SST `F2` / contributions murales de ce cas | 0 / 0 | Non probant pour des parois SST actives |

**Important** : ce seuil `1e-6` est un garde de linéarité du diagnostic V57, non une tolérance physique ou une preuve que l'opérateur convectif est dérivable au point d'égalité du limiteur. Des valeurs proches du seuil peuvent aussi être affectées par les pas de différence finie ; cette campagne ne sépare pas encore systématiquement ces effets. Les zéros des contributions murales ne démontrent pas l'équivalence des fonctions de paroi OpenFOAM.

**Krylov global sur le même BFS152** : vrai résidu linéaire après six directions `0.957920841902`, défaut d'additivité de l'action `1.824151640401e-4`, garde `1e-6`, cible linéaire expérimentale `0.3`. Le nouveau certificat V60 indique `remainder_ok=0 additive=0`, raison `V60 freshly reassembled blocks/RC/convection fail Krylov additivity`, `0` pas Newton accepté. Le solveur finit `Converged NO`. Ce résultat démontre que figer les seules branches SOU ne suffit pas ; **le verrou P0 reste non résolu**.

### Expérience indépendante : figer `Amat` global V54 sur le BFS152 réel

Une seconde exécution sur le **même BFS152 et même état initial**, cette fois **sans tangent hybride** et en activant `CFDX_FPHYS_JFNK_V54_FROZEN=1`, a construit une matrice `Amat` explicite de 912 colonnes, obtenue à partir de **1 824 évaluations** physiques (912 différences centrées). Les six directions Arnoldi ont été recalculées avec cette matrice réellement figée :

| Indicateur V54 figé | Valeur |
|---|---:|
| Défaut d’additivité finale de `Jv` | 5.448946345473e-14 |
| Défaut maximal de fidélité des directions à la dérivée physique | **1.055273876165e-1** |
| Tolérance de fidélité V54 inchangée | **5.0e-2** |
| Défaut des dérivées latérales (diagnostic non-lisse) | 1.295767817779 |
| Inversions de signe de flux détectées | 136 |
| Résidu linéaire vrai après 6 directions | **9.722538860487e-1** |
| États Newton acceptés | **0** |

Le rapport indique `v54_frozen_amat=1`, `v54_fidelity_ok=0`, et `reason=V54 physical derivative fidelity failed for at least one Krylov direction`. **Interprétation : l’opérateur figé est bien additif sur Arnoldi, mais le lien avec le véritable résidu physique reste insuffisant.** Cela prouve qu’une simple correction de l’additivité de Krylov ne résout pas P0. Aucun test long BFS712/BFS20 540 n’a été lancé sur ce Newton rejeté. Ce calcul borné en `6N=912` est un **oracle coûteux** et non un préconditionneur généralisable industriellement sans changement d’algorithme.

### Limites et suites justifiées

- **P0-A :** construire une matrice linéaire généralisée `Amat` figée pour *toutes* les directions Arnoldi, cohérente avec le traitement des égalités triples SOU/MUSCL et les contributions non convectives ; vérifier distinctement `F(X+h v)` contre `Amat·v`. Une sélection de sous-gradient n'est pas une preuve de descente.
- **P0-B :** si nécessaire, isoler l'effet des pas finis et l'échelle absolue des traces sous `1e-6`, conserver les mêmes valeurs physiques `F(X)` dans tous les essais et l'opérateur physique primal inchangé.
- **P0-C :** démontrer une réduction acceptée du **résidu physique `U,p,k,omega`**, y compris Rhie–Chow fermé et parois SST actives, avec relecture indépendante du vrai résidu linéaire.
- **P0-D :** traiter PCD `4N` (noyau pression, conditions aux limites, nullspaces, Fp à chaque état) ; ne pas confondre sa réussite éventuelle sur matrice isolée et le Schur réel.
- **P0-E :** après BFS152 **réellement convergé**, progresser sur BFS712 puis BFS20 540 et valider le rattachement `x_r/h`, `Cf`, `Cp`, vitesse et contraintes turbulentes par rapport aux références indépendantes OpenFOAM/NASA. Aucun seuil ne doit être assoupli pour obtenir un succès nominal.

## 4. Six Schur sur matrices réelles de BFS152

Les variantes `PCD`, `LSC`, `BFBt`, `BlockLocal`, `SIMPLE-Schur` et `SIMPLEC-Schur` ont reçu le même système **COUPLED 4N** de 152 cellules (`608` lignes, `7080` non-zéros) : preuves identiques de `A`, `b` et `x0` par SHA-256 de captures Matrix Market, vérifiées pour les six méthodes. Les fichiers de capture sont également comparés entre exécutions V59 et V60 : les entrées du solveur linéaire sont identiques malgré des métadonnées de cas différentes. Ce comparateur linéaire `4N` ne qualifie pas à lui seul le Newton physique `6N`.

| Variante | Première phase linéaire | Premier obstacle de qualification | Convergence physique |
|---|---|---|---|
| PCD | FGMRES 2000 itérations, résidu vrai relatif ~6.274552e-3 | `KRYLOV_SCHUR` | NON |
| LSC | Phase linéaire terminée | `POST_MOMENTUM_NOT_CERTIFIED` | NON |
| BFBt | Phase linéaire terminée | `POST_MOMENTUM_NOT_CERTIFIED` | NON |
| BlockLocal | Phase linéaire terminée | `POST_MOMENTUM_NOT_CERTIFIED` | NON |
| SIMPLE-Schur | Phase linéaire terminée | `POST_MOMENTUM_NOT_CERTIFIED` | NON |
| SIMPLEC-Schur | Phase linéaire terminée | `POST_MOMENTUM_NOT_CERTIFIED` | NON |

La campagne est **volontairement négative** : elle rejoue une seule itération extérieure afin de contrôler l'identité des matrices et identifier les premiers obstacles. Aucune des six variantes n'est physiquement qualifiée. Le problème de PCD n'a pas été corrigé dans ce lot ; la nouvelle fusion radiative ne modifie pas les opérateurs `A,b,x0` capturés.

## 5. Tests et compilation réalisés

- Compilation ciblée des deux tests HO V57, des quatre tests radiatifs antérieurs et de l'audit V60 ; tous fonctionnels en Release CPU diagnostic. Deux configurations ont été utilisées : C++ G++ pour le noyau/tests, et un compilateur Clang++ pour l'unité principale volumineuse, avec la bibliothèque C++ commune déjà bâtie ; le binaire du solveur est une **compilation diagnostique `-O0 -g0` de l'unité principale**, pas une mesure optimisée de performance.
- **16/16 CTest ciblés** BFS/SST et radiatifs après exclusion du stress V54 long (compilation et test). **3/3 tests production radiatifs ciblés** (`test_production_radiation_contract`, `test_radiation_vv`, `test_thermal_radiation_model_matrix`). Cela fait **19/19 tests ciblés distincts**, hors stress long.
- Les fichiers journaux du BFS152 (hybride V60 **et** matrice figée V54), la comparaison six Schur et les cinq rejouages radiatifs sont conservés séparément dans la livraison ; les sorties `Converged NO` et les codes non nuls des campagnes Schur ne sont jamais changés en PASS.
- Tests **non réalisés** : totalité de la CI, sanitizers complets, MPI/GPU natifs, comparateur externe OpenFOAM, BFS712 et BFS20 540, campagnes de maillages NACA et qualification physique industrielle du DOM/SST.

## 6. Reproduction locale à partir des sources V60

```sh
cmake -S . -B build-v60 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_ENABLE_RADIATION_RESEARCH_V57=ON \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v60 -j2 --target \
  cfdx_production_solver test_physical_remainder_audit_v60 \
  test_dom_angular_hilo_v57 test_dom_holo_hilo_v57
ctest --test-dir build-v60 --output-on-failure \
  -R '^(test_physical_remainder_audit_v60|test_dom_angular_hilo_v57|test_dom_holo_hilo_v57)$'
# Les campagnes de recherche radiative sur 640 cellules sont délibérément longues :
./build-v60/test_dom_holo_hilo_v57 --extreme
./build-v60/test_dom_holo_hilo_v57 --family
# À réserver pour la recherche longue, pas à la CI courte :
#   -DCFDX_ENABLE_LONG_VALIDATION=ON
#   ctest --test-dir build-v60 -R '^test_dom_independent_angular_v54_stress$'
```

Exemple de **contrôle négatif** sur BFS152, après fourniture du vrai fichier HDF5 de cas dans `<BFS152_CASE.h5>` :

```sh
CFDX_FPHYS_CANDIDATE_V51=1 CFDX_FPHYS_JFNK_V52=1 \
CFDX_FPHYS_JFNK_V52_FREEZE_PRESTEP=1 CFDX_FPHYS_JFNK_V53_BLOCK=1 \
CFDX_FPHYS_JFNK_V52_KRYLOV=6 CFDX_FPHYS_JFNK_V57_HYBRID=1 \
CFDX_FPHYS_JFNK_V57_TRUTH=1 CFDX_FPHYS_JFNK_V59_FROZEN_SOU=1 \
CFDX_FPHYS_JFNK_V60_REMAINDER_AUDIT=1 \
./build-v60/cfdx_production_solver --mesh <BFS152_CASE.h5> \
  --output-dir build-v60/bfs152-v60 --iterations 1 --experimental-sst
```

**L'exécution ci-dessus doit normalement sortir en échec de convergence** dans l'état V60 testé : c'est le comportement fail-closed voulu, pas une réussite Newton. Pour la campagne des six variantes, utiliser `scripts/run_bfs_six_schur_v37.py` avec un fichier HDF5 commun, `--fingerprint-operator-v40` et `--capture-original-operator`. Le rapport JSON doit porter `all_six_original_operators_captured=true`, `all_six_identical_original_operator=true`, `all_six_nonlinear_converged=false` sur cette expérience.

## 7. Patchs et règle de conservation

Le ZIP utilisateur BFS V6 est **l'unique point de départ** du patch cumulatif V60. Les deux chemins vérifiés sont : `BFS V6 ZIP → V60` (patch cumulatif) et `V59 livrée → V60` (patch incrémental). La livraison doit contenir les deux `.patch`, une archive des sources finales, le rapport, les données d'expérience et un manifeste SHA-256 ; leur identité sera contrôlée après leur génération.

**Verdict V60 : intégration logicielle V57/V59 effectuée, tests radiatifs ciblés positifs et localisation du P0 améliorée. Aucun Newton physique accepté ; six Schur non qualifiés ; pas de validation BFS712/BFS20 540, parois SST OpenFOAM actives, MPI/GPU, ni CI complète. Aucune modification automatique des méthodes de production, aucun seuil assoupli et aucun push GitHub.**
