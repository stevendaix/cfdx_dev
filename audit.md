Gmail	Steven daix <steven.daix@gmail.com>
(aucun objet)
Steven daix <steven.daix@gmail.com>	25 septembre 2026 à 20:48
À : Steven daix <steven.daix@gmail.com>
# Audit détaillé — `stevendaix/cfdx_dev`

- **Date** : 2026-09-25 (v2 : codes de correction complets, revue de la PR #415)
- **Commit audité** : `master` @ `27bae126`, les branches des PRs ouvertes (`pr/<N>`), et la tête `f829e32f` de #415
- **Méthode** :
  - clone en lecture seule ;
  - build local + `ctest` + `pytest` ;
  - lecture des logs CI (`gh run`, lecture seule) ;
  - 4 revues par domaine (numérique/couplage P-V, physique, solveurs/I-O/perf, tests/CI/PRs), plus une revue dédiée de la PR #415 (§10) ;
  - **chaque correctif proposé a été appliqué dans une copie de travail** (jamais dans le clone), compilé et, si possible, testé par un test qui échoue avant et passe après. Les correctifs sont copiés dans `correctifs/` à côté de ce rapport (Annexe A).
- **Contrainte respectée** : aucun commit, aucun push, aucune écriture GitHub (ni commentaire, ni édition d'issue ou de PR). Les commandes `gh` du rapport sont **à exécuter par le propriétaire** ; aucune n'a été lancée.

**Légende de statut**
- **CONFIRMÉ** : reproduit par un test ou une exécution, ou relu ligne à ligne.
- **PLAUSIBLE** : lecture de code cohérente mais non reproduite. À vérifier avant correction.

**Statut des codes proposés** (en titre de chaque constat)
- **compilé + testé** : appliqué en copie de travail, compilé, avec un test qui prouve l'effet ;
- **compilé** : appliqué et compilé, sans test dédié ;
- **non compilé** : code écrit contre le source réel mais non construit (CUDA, MPI-HDF5, etc.).

**Légende de gravité**
- 🟥 **BLOQUANT** : résultat faux ou calcul impossible.
- 🟧 **MAJEUR** : résultat biaisé, ou robustesse compromise.
- 🟨 **MINEUR** : qualité, précision ou lisibilité.

> Les numéros de ligne renvoient à `master@27bae126`, sauf mention `pr/<N>:`. Chemins relatifs à
> `src/cfdx/` quand il n'y a pas de préfixe (`physics/`, `core/`, …).

---

## 1. Synthèse exécutive

### 1.1 Réponse courte à « pas de résultats, des décalages »

1. **Aucun écoulement piloté par la pression ne peut être juste aujourd'hui** (A-B1, §3.1).
   - Le gradient de Green-Gauss du solveur incompressible divise ∇p par ~2 : −0,428 au lieu de −0,879.
   - Un Poiseuille donne une vitesse de 0,560 au lieu de 0,996, soit un décalage de 44 %.
   - **Aucune PR ne le corrige, #415 comprise** : u_max 0,593 et ordre observé −0,06 sur #415, contre 1,87 (SIMPLE) et 1,78 (COUPLED) avec le correctif (§10.6).
   - Aucun test ne le voit, car les tests de couplage utilisent un Couette, où ∇p = 0.
   - Correctif : erreur de gradient 1,5 → 4,9e-15.
2. **Le couplage SIMPLE diverge dès que la pression initiale n'est pas nulle** (A-B2, A-B4, §3.2).
   - Cause : relaxation *explicite* de la quantité de mouvement combinée à un `rAU` calculé sur le a_P *non relaxé* ; la correction de pression est sur-amplifiée d'un facteur a_P/(a_P−Σ|a_nb|).
   - L'ancrage de pression sur une seule cellule déclenche la divergence : c'est l'échec de `test_phase9_acceptance` et du CI `master`.
   - Avec les correctifs : Couette converge avec SIMPLE, SIMPLEC, PISO et PIMPLE (382 à 387 itérations), Poiseuille en 454, ordre 2 de 16 à 32 cellules, indépendance vis-à-vis de α à 3,95e-9 (avec Majumdar).
3. **Le solveur COUPLED échoue à cause du restart GMRES, pas de la matrice** (§4 C4/C5, §10.3).
   - `gmres_solver.h:41` écrête tout restart demandé à [10, 40] ; `krylov_controls.h:32-35` *rétrécit* l'espace de Krylov quand la convergence est mauvaise. Les deux défauts sont sur `master` **et** dans #415.
   - Le « restart 512 » du commit `f829e32f` de #415 n'a donc aucun effet.
   - Les deux correctifs suffisent : COUPLED converge en 48 itérations, phase9 passe (6/6).
4. **Plusieurs modèles physiques ont des erreurs de formule directes** (§3.3 à 3.6).
   - Diffusion turbulente passée comme terme source (B1) : sur une chaîne 1D, k diverge jusqu'à 2,7e129 au lieu de rester à 1.
   - `cp` absent de la convection d'énergie (B2).
   - Signe du terme source radiatif inversé ; la correction exige une linéarisation de Newton (B3).
   - P1 sans le facteur π (B13) ; Spalart-Allmaras `fw` faux (B4) ; extrapolation d'enthalpie de signe faux (B6).
   - Nouveau : Rosseland perd la conductivité moléculaire et `cp` (B18) ; le RK2 « bas stockage » n'est que d'ordre 1 (§6).
5. **Les solveurs linéaires et les I/O mentent ou plantent** (§4).
   - NaN converti en 0 dans la norme Krylov : GMRES annonce « CONVERGED » avec x = NaN (C1).
   - Statut du CG pression jamais vérifié (C2).
   - Lecture de restart hors limites (C7) ; ordre des halos MPI faux (N2).
   - Nouveau : l'import OpenFOAM déborde la pile de `std::regex` sur un gros maillage (C20).
6. **« Pas de résultats » s'explique par la chaîne de validation** (§5, §6).
   - La matrice de validation VMFL a **0 cas PASS**.
   - `run_validation.py` rend **0 sur un `master` rouge** : il n'exécute que `ctest -L phase13`, et phase9 n'a pas ce label (D8). Le remplaçant proposé rend 1.
   - Le workflow VMFL a échoué 9 fois sur 9. Ce n'est **pas masqué** par `continue-on-error` (l'étape `Enforce validation gate` refait échouer le job), mais il ne tourne que sur `workflow_dispatch` ou sur le label `validation`.
   - 81 tests sur 93 durent moins de 0,01 s ; une partie est tautologique (D6) ; des tests ne sont jamais enregistrés (D7) ; `validation_report.py` produit un LaTeX corrompu **sans avertissement** (`\t` → tabulation).
7. **Le CI ne protège rien** (§6).
   - `master` est rouge depuis le 2026-09-23 17:18 ; 10 PRs ont été mergées rouges ; la branche n'a aucune protection.
   - Le workflow « validation » reste vert pendant que phase9 échoue, faute du label CMake `phase13` (D2).
   - Le build est cassé sans MPI (D4) ou sans HDF5 (D5).

### 1.2 Top 10

| # | Gravité | Constat | Où | Statut |
|---|---|---|---|---|
| 1 | 🟥 | Gradient Green-Gauss : φ_f = φ_c quand la cellule courante est le voisin, donc ∇p ≈ ½ (A-B1) | `physics/steady_incompressible_solver.h:137-145` | CONFIRMÉ, compilé + testé |
| 2 | 🟥 | Relaxation explicite de la QDM avec `rAU = V/a_P` non relaxé, donc divergence SIMPLE (A-B2, A-B4) | `physics/finite_volume_transport.h:567`, `steady_incompressible_solver.h:414-472` | CONFIRMÉ, compilé + testé |
| 3 | 🟥 | Restart GMRES écrêté à [10, 40] et politique adaptative inversée : COUPLED en `MAX_ITER` (C4, §10.3) | `core/linalg/gmres_solver.h:41`, `krylov_controls.h:32-35` | CONFIRMÉ, compilé + testé |
| 4 | 🟥 | Diffusion k/ε/ω passée en `extra_rhs` au lieu de `cell_diffusion` (B1) | `physics/turbulence_solver.h:96-101,181-182,231`, `sst_solver.h:96-101` | CONFIRMÉ, compilé + testé |
| 5 | 🟥 | Convection d'énergie sans `cp` (B2) | `physics/energy_solver.h:76-78,146` | CONFIRMÉ, compilé + testé |
| 6 | 🟥 | Source radiative de signe inversé, non linéarisée (B3) ; P1 sans π (B13) | `physics/radiation_solver.h:491` | CONFIRMÉ, compilé + testé |
| 7 | 🟥 | NaN converti en 0 dans la norme Krylov (C1) ; statut du CG pression ignoré (C2) | `core/linalg/krylov_reductions.h:86`, `steady_incompressible_solver.h:578-585` | CONFIRMÉ, compilé (+ testé pour C1) |
| 8 | 🟥 | `master` rouge, aucune protection, 10 PRs mergées rouges ; validation verte à tort ; pilote de campagne à 0 sur rouge (D1, D2, D8) | CI, `CMakeLists.txt:394,412-419`, `scripts/run_validation.py` | CONFIRMÉ, correctifs validés |
| 9 | 🟧 | Spalart-Allmaras `fw` faux ; enthalpie extrapolée de signe faux ; lois de paroi sans κ (B4, B6, B7) | `spalart_allmaras.h:55-57`, `thermophysical_models.h:166-167`, `turbulence_solver.h:302,309` | CONFIRMÉ, compilé + testé |
| 10 | 🟧 | ~130 Mo d'artefacts de build versionnés, 157 binaires ELF ; build cassé sans MPI ou sans HDF5 (D11, D4, D5) | `build*/`, `core/parallel/mpi_utils.h:7`, `CMakeLists.txt:87,134-136` | CONFIRMÉ, correctifs validés |

Hors Top 10 mais bloquants : lecture de restart hors limites (C7), ordre des halos (N2), débordement de pile de l'import OpenFOAM (C20).

### 1.3 Verdict sur les PRs

- **#415 (tête `f829e32f`) : NON MERGEABLE en l'état** (§10.7).
  - CI rouge sur les deux workflows (`test_gmres_solver`, `test_phase9_acceptance`).
  - Oracles faux dans `test_gmres_solver.cpp` (un système régulier attendu `DIVERGED`).
  - A-B1 et les 5 `'\\n'` toujours présents, alors que la description annonce le contraire.
  - Avec les correctifs du §10.4 : ctest 92/93 (le dernier échec est D4, préexistant) et phase9 6/6.
- **Route recommandée** (§7.2, arbitrée) :
  1. #410 + correctifs d'abord : seule branche sans conflit avec `master` et mesurée verte (ctest 96/96) ;
  2. #415 ensuite, rebasée sur ce résultat, avec ses correctifs §10.4 ;
  3. #409 et #411 fermées, car entièrement contenues dans #415.

### 1.4 Mesures réelles

| Mesure | Résultat |
|---|---|
| CI `build` sur `27bae126` (run 36036549020) | ❌ 1 échec : `test_phase9_acceptance` |
| CI `CFDX validation` sur le même commit (run 36036548902) | ✅ vert, car phase9 n'a pas le label `phase13` |
| Historique CI `build` | 1 succès sur les ~60 derniers runs ; dernier vert : `b2c69cba` (#358) |
| Build local sans HDF5 | ❌ `H5public.h` introuvable, alors que CMake annonce « HDF5-dependent features disabled » |
| Build local MPI=OFF | ❌ `test_performance_runtime` : `mpi.h` inclus sans condition |
| `ctest` local `master` (HDF5 compilé à la main, MPI OFF) | 2 échecs sur 93 (`test_performance_runtime` ne compile pas ; `test_phase9_acceptance`) |
| `ctest` #410 + correctifs | ✅ 96/96 |
| `ctest` #415 tel quel / corrigé | 90/93 → 92/93 ; phase9 6/6 après correction |
| `pytest` (après correctifs Python) | 111 passed, 2 skipped, **6 failed** : le job pytest ne doit pas encore être requis |
| `run_validation.py` sur `master` rouge | rc = 0 (original) → rc = 1 (remplaçant) |
| Durée des tests | 81 sur 93 en moins de 0,01 s ; suite `phase13` complète : 0,49 s |
| Exécution de `phase9` | Résidu BiCGStab : 1e-12 → 2e-2 au fil des itérations externes, puis échec à 2000 itérations. Message final : `Ux momentum solve did not converge (status=1, iterations=128, residual=8192, relative=8192)` |
| Logs | Chaque ligne se termine par `23662` au lieu d'un retour à la ligne : `'\\n'` est un littéral multi-caractère, 0x5C6E = 23662 (§3.6, m1) |
| Optimisations mesurées (§8) | AMG-PCG ×2 à ×2,3 en temps (645 → 44 it) ; SpMV distribué ×14 ; RCM ×160 ; BiCGStab −35 % ; espaces de travail −25 % (CG) ; OpenFOAM ×4,5. Le bounds-check (ancienne priorité 1) ne coûte que 0 à 8 % |

### 1.5 Corrections apportées à la v1 de cet audit

La v1 contenait des affirmations inexactes, corrigées dans le corps du rapport :
- VMFL n'est **pas** masqué par `continue-on-error` (§5.1) ;
- `validation/cases.json` n'existe pas (§5.1) ;
- la cible réelle de D4 est `communication_avoiding.h` ; le code D5 de la v1 était faux ;
- les noms de check de la protection de branche étaient faux (§6 D1) ;
- m6 : le « 0,75/0,25 » est un étage d'intégration temporelle ; un limiteur existe mais n'est pas câblé ;
- C16 : le « volume 0,333 » venait du maillage de test, constat rétrogradé en 🟨 ; C14 n'existe pas (trou de numérotation) ;
- §8 : l'optimisation n° 1 (×2 à ×4) est **contredite** par la mesure ; le gain OpenFOAM est ×4,5, pas ×10 ;
- B7 : y+ de croisement = 11,53 et non 11,225 ; B16 : seule la copie `m1_m4_models.h` omet A1/A2 ;
- les blocs « Code actuel » approximatifs de la v1 ont tous été remplacés par le source réel.

---

## 2. Hygiène du dépôt

> Toutes les corrections de cette section ont été appliquées **uniquement** dans une copie
> de travail (`/tmp/cfdx_scratch`, sans `.git`, base `master@27bae126`). Rien n'a été
> committé ni poussé. Les commandes `git` sont **à exécuter par le propriétaire**.

### 2.1 🟧 Artefacts de build versionnés (D11) — CONFIRMÉ

| Dossier | Fichiers suivis |
|---|---|
| `build-release/` | 456 |
| `build-debug-sanitizers/` | 456 |
| `build_test/` | 481 |
| `build-test/` | 412 |
| `poisson_mms` (ELF, mode 100755) | 1 |
| **Total** | **1806 fichiers, 99,2 Mo** (`git ls-files -z … \| xargs -0 du -cb`) |

Le `.gitignore` contient **déjà** `build*/` (commit `bcb1ea2f`, « Ignore generated CMake
build trees »). Il n'agit pas ici parce qu'un motif d'ignore **ne s'applique jamais à un
fichier déjà suivi** : les 1805 fichiers ont été ajoutés avant la règle. Ajouter des motifs
ne suffit donc pas, il faut désindexer.

Second défaut : la dernière ligne est `prompt.mdbuild_san/` (saut de ligne manquant, deux
motifs fusionnés). `build_san/` est couvert par `build*/`, mais `prompt.md` n'est pas ignoré.

**Correctif `.gitignore`** (remplace la dernière ligne ; le reste du fichier est conservé) :
```diff
 # IDE
 .vscode/
 .idea/
 *.swp

-prompt.mdbuild_san/
+# Fichiers de travail locaux
+prompt.md
+
+# Binaires compiles a la racine (poisson_mms etait versionne)
+/poisson_mms
+/_b*/
+Testing/
```

**Désindexation**, à exécuter par le propriétaire (non exécuté par l'audit) :
```bash
# Depuis une branche dediee, par ex. chore/untrack-build-artefacts
git rm -r --cached --quiet build-release build-debug-sanitizers build_test build-test poisson_mms
git status --short | grep -v '^D ' && echo "ATTENTION : autre chose que des suppressions"   # doit ne rien afficher
git commit -m "chore: stop tracking build artefacts (1806 files, 99 MB)"
```
- `--cached` retire du suivi sans toucher aux fichiers locaux.
- L'historique garde les 99 Mo. Pour les purger, il faut réécrire l'historique :
  `git filter-repo --invert-paths --path build-release --path build-debug-sanitizers --path build_test --path build-test --path poisson_mms`.
  C'est destructif : tous les clones et les PR ouvertes doivent être rebasés. **À ne faire
  qu'après** le tri des PR du §7.

**Statut** : `.gitignore` VALIDÉ en scratch (`git check-ignore` n'est pas disponible sans
`.git` ; le motif a été relu). Commandes `git` : NON EXÉCUTÉES (hors périmètre).

### 2.2 🟨 Doublons et code mort (D12) — CONFIRMÉ, suppression VALIDÉE en scratch

**Méthode.** Pour chaque candidat, on a compté les références hors du dossier lui-même
(`grep -rn` sur `src tests apps python CMakeLists.txt`, extensions `.h .cpp .txt .cmake .py`) :

| Élément | Références externes | Verdict |
|---|---|---|
| `src/cfdx/cfdx/core/**` (7 squelettes : `execution_policy.h`, `block_preconditioner.h`, `adaptive_time_step.h`, `local_time_stepping.h`, `halo_exchange.h`, `mesh_partitioner.h`, `convergence_monitor.h`) | 0 | **supprimer** |
| `src/cfdx/linear_algebra/` (`gmrès_solver.h`, `ilu0_preconditioner.h`, `rbgs_preconditioner.h`, `test_gmres.cpp`, `test_rbgs.cpp`) | 0 | **supprimer** (tests placés dans `src/`, nom de fichier accentué) |
| `src/cfdx/solvers/navier_stokes.h` | 0 | **supprimer** |
| `src/cfdx/core/visualization/mesh_visualizer.{h,cpp}` | 1 : `tests/unit/test_visualization.cpp`, jamais enregistré | **supprimer**, avec son test |
| `tests/unit/test_vtu_writer.cpp` | inclut `"vtu_writer.h"` (inexistant) et `<gtest/gtest.h>` (non fourni) | **supprimer** : ne peut pas compiler |
| `src/cfdx/thermodynamics/equation_of_state.h` | 2, **uniquement** `CMakeLists.txt:166` et `:265` (listes de sources) | **supprimer** + retirer les 2 lignes |
| `src/cfdx/thermodynamics/thermo_state.h`, `thermo_cache.h` | `physics/low_mach.h`, `physics/compressible_flux.h`, 3 tests | **CONSERVER** |
| `physics/equation_of_state.h` | `transport_models.h:12` (`"equation_of_state.h"`, même dossier) et `utils/logging_profiling.h:3` (`"physics/equation_of_state.h"`) | **CONSERVER** : c'est l'EOS utilisée |

Au total, 16 fichiers et 792 lignes supprimés.

Hors périmètre de cette suppression (encore référencés, à traiter avec le §3.2) :
- les doublons `memory_planner.h` et `temporal.h` ;
- `pressure_velocity_algorithms.h` et `pressure_velocity.h`, épinglés par `test_m1_m4_physics.cpp:19`.

**Commandes** (à exécuter par le propriétaire, dans une branche `chore/remove-dead-code`) :
```bash
git rm -r src/cfdx/cfdx src/cfdx/linear_algebra src/cfdx/solvers src/cfdx/core/visualization
git rm src/cfdx/thermodynamics/equation_of_state.h tests/unit/test_visualization.cpp tests/unit/test_vtu_writer.cpp
```

**`CMakeLists.txt`** (deux suppressions de ligne) :
```diff
 # Thermodynamics sources (P2 reorganization)
 set(CFDX_THERMODYNAMICS_SOURCES
-    src/cfdx/thermodynamics/equation_of_state.h
 )
@@
     src/cfdx/core/parallel/restart_mapping.h
-    src/cfdx/thermodynamics/equation_of_state.h
     src/cfdx/chemistry/source_terms.h
```

**Preuve** (scratch, build propre `_b2`, `-DCFDX_ENABLE_MPI=OFF`) :
- configuration et compilation OK ;
- `ctest` : 102 tests, 100 PASS ;
- les 2 échecs (`test_phase9_acceptance`, `test_poiseuille_body_force`) viennent du solveur
  de `master`, qui diverge (§5.3). Ils sont identiques avant la suppression.

**Statut** : VALIDÉ en scratch (`/tmp/cfdx_meta/ctest_d12.log`).

### 2.3 Correction d'une affirmation antérieure (D13)

`spec.md`, `pyproject.toml` et `memory_planner.md` sont bien des **fichiers** suivis, pas des
dossiers. Le plan initial disait le contraire, à cause d'un artefact d'analyse. Seul
`poisson_mms` pose problème à la racine : c'est un exécutable ELF versionné (mode `100755`,
blob `d4e3f985`). Il est traité au §2.1.

---

## 3. Physique et équations

### 3.1 Discrétisation volumes finis

> **Base des corrections.**
> - Toutes les lignes « Code actuel » sont celles de `master@27bae126`, relues dans le source.
> - Le code proposé a été appliqué dans une **surcouche hors dépôt** (`/tmp/cfdx_meta/scratch_a/fix3a/cfdx/...`, placée devant `-I /tmp/cfdx_dev/src`). Rien n'a été modifié dans le dépôt, rien n'a été commité.
> - Diffs complets :
>   - `scratch_a/fix3a_transport.diff` (`finite_volume_transport.h`, 134 lignes) ;
>   - `scratch_a/fix3a_solver.diff` (`steady_incompressible_solver.h`, 172 lignes).
> - Statut de chaque bloc :
>   - **compilé + testé** : compilé avec `g++ -std=c++20 -O2 -I fix3a -I /tmp/cfdx_dev/src`, et couvert par un test qui échoue sur `master` et passe avec la correction ;
>   - **compilé** : construit, mais sans test dédié ;
>   - **non compilé** : proposition écrite contre l'API réelle, à compiler.

#### A-B1 🟥 Gradient Green-Gauss faux côté voisin — CONFIRMÉ — *compilé + testé*
`physics/steady_incompressible_solver.h:137-145` (fonction `gauss_gradient_with_boundary`, l.118). Le code est identique dans toutes les branches `pr/*`.

La boucle parcourt les faces de la cellule `c`. Quand `c` est le **voisin** de la face, `n = own.neighbour(f)` vaut `c`. La valeur de face devient alors `0.5*(field(c)+field(c)) = field(c)` : l'autre cellule n'intervient jamais. Un Couette (dp = 0) ne peut pas le voir. **Aucune PR ne le corrige.**

**Code actuel** (l.142-145)
```cpp
if (own.neighbour(f) >= 0) {
    const std::size_t n = static_cast<std::size_t>(own.neighbour(f));
    vf = 0.5 * (field(c) + field(n));
}
```
**Code proposé** : prendre l'autre cellule, et interpoler en pondérant par la distance.
```cpp
if (own.neighbour(f) >= 0) {
    // A-B1 : l'autre cellule est le voisin si c est owner, l'owner sinon.
    const std::size_t other = owner
        ? static_cast<std::size_t>(own.neighbour(f))
        : static_cast<std::size_t>(own.owner(f));
    // Interpolation lineaire ponderee par la distance a la face.
    const double dc = (geometry.face_centres[f] - geometry.cell_centres[c]).mag();
    const double dn = (geometry.cell_centres[other] - geometry.face_centres[f]).mag();
    if (!(dc + dn > 0.0))
        throw std::runtime_error("gauss_gradient_with_boundary: degenerate face distance");
    const double w = dn / (dc + dn);
    vf = w * field(c) + (1.0 - w) * field(other);
}
```
`owner` est le booléen déjà calculé dans la boucle (`own.owner(f) == c`).

**Test** (`scratch_a/tests_fix3a.cpp`, cas `gg_linear_stretched`). Le maillage canal est étiré en x (`x_i = (i/nx)^1.4`, `channel_mesh_s.inc`). On calcule le gradient de p = 2x − 3y + 0,5, puis l'erreur max sur les cellules intérieures :

| | `master` | correction |
|---|---|---|
| erreur max ‖∇p − (2, −3, 0)‖ | **1,5** | **4,9e-15** |

À intégrer au dépôt sous la forme d'un `TEST` à tolérance 1e-12 sur les cellules intérieures d'un maillage non uniforme.

#### A-M3 🟧 Orientation des faces non contrôlée (#59) — CONFIRMÉ — *compilé + testé (C++), testé (Python)*
Correction de l'audit précédent : `build_fv_geometry` (`physics/finite_volume_transport.h:58`) **utilise déjà** `compute_cell_geometry_oriented` (l.84) et refuse les volumes non positifs (l.91). Ce qui manque, c'est la vérification que `Sf` sort de l'owner. Une face mal enroulée donne encore un volume positif, par compensation, mais un flux de signe faux.

**Code proposé** (`finite_volume_transport.h`, juste après la boucle des cellules l.91)
```cpp
// A-M3 : Sf doit sortir de l'owner (et pointer vers le voisin).
for (std::size_t f = 0; f < nf; ++f) {
    const std::size_t o = mesh.ownership().owner(f);
    const auto nb = mesh.ownership().neighbour(f);
    const Vec3 d = nb >= 0
        ? g.cell_centres[static_cast<std::size_t>(nb)] - g.cell_centres[o]
        : g.face_centres[f] - g.cell_centres[o];
    if (!(g.face_area_vectors[f].dot(d) > 0.0))
        throw std::runtime_error(
            "build_fv_geometry: face " + std::to_string(f) +
            " area vector does not point out of its owner cell " + std::to_string(o) +
            " (check the importer face winding)");
}
```
**Test C++** (`tests_fix3a.cpp`, cas `flipped_face_detected`) : on inverse l'ordre des sommets d'une face intérieure. Avec la correction, l'exception est levée (1) ; sur `master`, elle ne l'est pas (0).

**Importeur** : `scripts/meshio_import.py:9-15` (`VOLUME_FACES`). Le test du centroïde, exécuté sur chaque type, donne :

| Type | Faces rentrantes | Correction |
|---|---|---|
| tétraèdre | `(1,3,2)` | `(1,2,3)` |
| pyramide | `(0,1,2,3)` et `(0,4,1)` | `(0,3,2,1)` et `(0,1,4)` |
| hexaèdre, voxel | aucune | — |
| prisme | aucune dans l'orientation gmsh | le prisme VTK strict est enroulé à l'envers |

Les tuples proposés dans la version précédente de cet audit n'étaient pas au format du dictionnaire réel, et ils étaient en partie faux. **Plutôt que de corriger des tuples**, on propose d'orienter géométriquement dans `build()`, après la construction de `faces`/`owner`. C'est robuste pour tous les types, polyèdres compris.

**Code proposé** (`scripts/meshio_import.py`, nouvelle fonction appelée en fin de `build()`) — *testé* : `scratch_a/orient_faces.py`
```python
import numpy as np

def orient_faces(points, faces, owner, cells):
    """Retourne chaque face pour que sa normale sorte de sa cellule owner.

    Test geometrique local, independant des conventions VTK/Gmsh :
    S . (x_f - x_c) > 0, avec S la normale de Newell et x_c le centre des sommets de l'owner.

    Args:
        points: tableau (n_points, 3) des coordonnees.
        faces: liste des faces (listes d'indices de sommets), modifiee sur place.
        owner: indice de la cellule owner de chaque face.
        cells: liste (ctype, nodes, block) produite par ``build`` ; pour un polyedre,
            ``nodes`` est la liste de ses faces.

    Raises:
        ValueError: face degeneree (aire nulle) ou normale tangente au rayon owner->face.
    """
    centres = []
    for ctype, nodes, _ in cells:
        ids = sorted({v for fc in nodes for v in fc}) if ctype == "polyhedron" else nodes
        centres.append(points[ids].mean(axis=0))
    for fid, face in enumerate(faces):
        pts = points[face]
        xf = pts.mean(axis=0)
        S = 0.5 * np.cross(pts - xf, np.roll(pts, -1, axis=0) - xf).sum(axis=0)
        s = float(S @ (xf - centres[owner[fid]]))
        if not np.isfinite(s) or abs(s) <= 1e-14 * max(float(np.linalg.norm(S)), 1e-300) ** 1.5:
            raise ValueError(f"face {fid}: degenerate or tangent normal (S.d = {s})")
        if s < 0.0:
            face.reverse()
    return faces
```
**Test** : on passe par `meshio_import.build()` sur `master`, puis on appelle `orient_faces`.

| Maillage | Faces rentrantes avant | Faces rentrantes après |
|---|---|---|
| 2 tétraèdres | 2 | 0 |
| 1 pyramide | 2 | 0 |

⚠️ La structure exacte de `cells` doit être alignée sur celle de `build()` au moment de l'intégration. Le test ci-dessus reconstruit la liste `(type, nœuds)` à la main.

#### A-M4 🟧 Géométrie de cellule calculée avec des normales globales non orientées — CONFIRMÉ — *compilé + testé*
`compute_cell_geometry` est appelé avec le `Sf` global au lieu de `compute_cell_geometry_oriented(face_centres, face_Sf, face_ids, n, CellIndex cell, const FaceOwnership&)`, défini en `core/geometry/cell_geometry.h:92`. Les appels concernés :

| Fichier:ligne | Variable maillage |
|---|---|
| `runtime/dynamic_mesh.h:39` (`cell_volumes(mesh)`) | `mesh` |
| `physics/pressure_velocity.h:51` | `mesh` |
| `core/numerics/temporal.h:295` et `:386` | `mesh` |
| `core/numerics/source_term.h:213` | `mesh` |
| `core/geometry/mesh_validator.h:105` | `m` |

**Code proposé** : même remplacement aux 6 sites, un seul argument ajouté.
```cpp
// avant
const auto cg = compute_cell_geometry(face_centres, face_Sf, face_ids, n);
// après
const auto cg = compute_cell_geometry_oriented(
    face_centres, face_Sf, face_ids, n,
    static_cast<CellIndex>(c), mesh.ownership());   // m.ownership() dans mesh_validator.h
```
La fermeture de `mesh_validator.h:157-159` somme aussi `Sf` sans inversion de signe :
```cpp
// avant
sum_Sf = sum_Sf + face_Sf[f];
// après : Sf sortant de c
sum_Sf = (m.ownership().owner(f) == c) ? sum_Sf + face_Sf[f] : sum_Sf - face_Sf[f];
```
**Test** (`scratch_a/am4_check.cpp`, compilé avec `.../core/mesh/boundary.cpp`) : deux cellules cubiques de volume 0,0625 qui partagent une face.

| | `validate_mesh` V | fermeture ‖ΣSf‖ | `dynamic_mesh` V | verdict |
|---|---|---|---|---|
| `master` | [0,0208 ; 0,0625] | 0,707 | faux | **FAIL** |
| correction | [0,0625 ; 0,0625] | 0 | 0,0625 | **A-M4 OK** |

#### A-M5 🟧 Condition `FIXED_GRADIENT` fausse — CONFIRMÉ — *compilé + testé*
`physics/finite_volume_transport.h:281-285`. On a ∂φ/∂n = g, donc φ_b = φ_P + g·d. Le code actuel a trois défauts :
- il convecte `bc.value` au lieu de φ_b ;
- il ignore la convection sortante de φ_b ;
- il prend `diffusion_coefficient` même quand un Γ par cellule est fourni.

La version précédente de l'audit utilisait des noms inventés (`controls.diffusivity`, `Sf_mag`). La version ci-dessous reprend les noms réels du bloc.

**Code actuel**
```cpp
} else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
    diag[o] += std::max(F, 0.0);
    rhs[o] += (F < 0.0 ? -F * bc.value : 0.0);
    rhs[o] += diffusion_coefficient * area * bc.gradient;
    div_phi[o] += F;
```
**Code proposé**
```cpp
} else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
    // phi_b = phi_P + g*d : flux diffusif impose gamma_b*g*|Sf|,
    // flux convectif F*phi_b = F*phi_P + F*g*d (entrant ou sortant).
    const double gamma_owner = cell_diffusion
        ? (*cell_diffusion)[o] : diffusion_coefficient;
    if (gamma_owner < 0.0)
        throw std::invalid_argument("assemble_scalar_equation: negative boundary diffusion");
    diag[o] += F;
    rhs[o] -= F * bc.gradient * distance;
    rhs[o] += gamma_owner * area * bc.gradient;
    div_phi[o] += F;
```
`distance` et `area` sont les variables déjà calculées dans la branche de bord.

⚠️ Avec F < 0 et `bounded_convection = false`, `diag += F` peut rendre la diagonale négative, et l.310 lèverait alors une exception. Si ce cas est visé, appliquer la même forme différée que A-M2.

**Test** (`tests_fix3a.cpp`, cas `fixed_gradient_cell_gamma`) : diffusion pure, Γ_cell = 2, T = x attendu. On impose T = 0 à gauche et ∂T/∂n = 1 à droite.

| | `master` | correction |
|---|---|---|
| max \|T − x\| | **0,484** | **5,7e-15** |

#### A-M2 🟧 Sortie gradient nul : puits artificiel en convection non bornée — CONFIRMÉ — *compilé + testé*
`physics/finite_volume_transport.h:292`. Le commentaire du code reconnaît lui-même le problème.

**Code actuel**
```cpp
diag[o] += bounded_convection ? F : std::max(F, 0.0);
div_phi[o] += F;
```
La correction de la PR #409 (`diag[o] += F`, proposée dans la version précédente) **casse le test du dépôt** `zero_gradient_inflow_does_not_reduce_owner_diagonal` : la diagonale tombe à 0 en entrée.

**Code proposé** : correction différée. La diagonale reste positive, et le résultat est exact à convergence.
```cpp
if (bounded_convection) {
    diag[o] += F;   // annule exactement par -div(phi)*psi
} else {
    // A-M2 : phi_b = phi_P aussi en entree. On garde a_P >= 0
    // (max(F,0)) et on differe la part entrante sur l'iteré
    // courant : exact a convergence, sans puits artificiel.
    diag[o] += std::max(F, 0.0);
    if (F < 0.0 && convected_field)
        rhs[o] -= F * (*convected_field)(o);
}
div_phi[o] += F;
```
**Test** : `test_finite_volume_transport` passe 6/6, y compris le test qui avait cassé avec la version #409.

#### m6 🟨 Précision géométrique et schémas — code CONFIRMÉ, impact PLAUSIBLE — *non compilé*
**Correction de l'audit précédent** : le « mélange fixe 0,75/0,25 » n'est **pas** un schéma de convection. C'est l'étage SSP-RK3 de `low_storage_time_integration.h:24`. Le dépôt a un vrai limiteur TVD, `apply_limiter_tvd(v_owner, v_upwind, v_extrap, LimiterType)` (`core/numerics/interpolation.h:90`, MINMOD/VANLEER/SUPERBEE/VAN_ALBADA, `InterpScheme::LIMITED`). Mais il n'est **pas branché** sur le solveur de quantité de mouvement : `ConvectionScheme` (`finite_volume_transport.h:40`) n'offre que `UPWIND` et `SECOND_ORDER_UPWIND`, écrêté à [min, max] des deux cellules (l.226-244).

Les points qui restent vrais :

1. **Centre de face** : c'est la moyenne arithmétique des sommets (`core/geometry/face_geometry.h:48-70`), avec `Sf` en éventail depuis v0. Le centroïde exact d'un polygone plan est le suivant ; il est à substituer dans `compute_face_geometry`, sans changer la signature.
```cpp
// face_geometry.h, dans compute_face_geometry : remplace le calcul de centre et de Sf
Vec3 x0{0, 0, 0};
for (VertexIndex k = 0; k < n_verts; ++k) {
    const VertexIndex v = vertices[offset + k];
    x0 = x0 + Vec3{px[v], py[v], pz[v]};
}
x0 = x0 * (1.0 / static_cast<double>(n_verts));
Vec3 Sf{0, 0, 0}, weighted{0, 0, 0};
double atot = 0.0;
for (VertexIndex k = 0; k < n_verts; ++k) {
    const VertexIndex va = vertices[offset + k];
    const VertexIndex vb = vertices[offset + (k + 1) % n_verts];
    const Vec3 a{px[va], py[va], pz[va]}, b{px[vb], py[vb], pz[vb]};
    const Vec3 s = (a - x0).cross(b - x0) * 0.5;      // sous-triangle (x0, a, b)
    const double ai = s.mag();
    Sf = Sf + s;
    weighted = weighted + (x0 + a + b) * (ai / 3.0);
    atot += ai;
}
if (!(atot > 0.0))
    throw std::runtime_error("FaceGeometry: degenerate face (zero area)");
const Vec3 centre = weighted * (1.0 / atot);
```
2. **Centre de cellule** : il doit être pondéré par les volumes des pyramides (x_f, x_c0). Dans `compute_cell_geometry_oriented` (`cell_geometry.h:92`), avec `Sf_out` le vecteur orienté sortant :
```cpp
// x0 = moyenne des centres de face ; pour chaque face : Vp = Sf_out.dot(xf - x0)/3
// xc_p = x0 + 0.75*(xf - x0) ; V += Vp ; xc += Vp*xc_p ; puis xc /= V (throw si V <= 0)
```
3. **Non-orthogonalité** : la diffusion utilise `gamma*A/|d_ON|` (l.~200). La décomposition sur-relaxée ajoute une correction différée avec le gradient de cellule déjà disponible (`reconstructed_gradient`) :
```cpp
// face interieure f, e = d/|d|, Sf = |Sf| n
const Vec3 Delta = Sf * (Sf.dot(Sf) / Sf.dot(d));          // sur-relaxe
const Vec3 k     = Sf - Delta;
const double D   = gamma_f * Delta.mag() / d.mag();        // implicite (remplace gamma*A/|d|)
const Vec3 grad_f = w * gO + (1.0 - w) * gN;               // gradients de cellule interpoles
deferred_rhs[o] += gamma_f * grad_f.dot(k);
deferred_rhs[n] -= gamma_f * grad_f.dot(k);
```
4. **Option TVD** pour `assemble_scalar_equation`. On ajoute `ConvectionScheme::TVD` et un `LimiterType`, puis on reprend le bloc SOU (l.226-244) en limitant par le ratio de gradient (Darwish-Moukalled), `r = 2 d_UD·∇φ_U/(φ_D − φ_U) − 1`. La convention de `r` dans `apply_limiter_tvd` (`(v_upwind − v_owner)/(v_extrap − v_owner)`) n'est pas documentée : il faut la vérifier avant de réutiliser la fonction telle quelle.
```cpp
} else if (convection_scheme == ConvectionScheme::TVD) {
    const std::size_t U = F >= 0.0 ? o : n, D = F >= 0.0 ? n : o;
    const double phiU = (*convected_field)(U), phiD = (*convected_field)(D);
    const Vec3 dUD = geometry.cell_centres[D] - geometry.cell_centres[U];
    const Vec3 gU{gx[U], gy[U], gz[U]};
    const double dphi = phiD - phiU;
    const double r = std::abs(dphi) > 1e-30 ? 2.0 * gU.dot(dUD) / dphi - 1.0 : 0.0;
    const double psi = limiter_psi(r, limiter);          // ex. vanLeer : (r+|r|)/(1+|r|)
    const double phi_high = phiU + 0.5 * psi * dphi;
    const double correction = F * (phi_high - phiU);
    deferred_rhs[o] -= correction;
    deferred_rhs[n] += correction;
}
```
Test à ajouter : un créneau advecté en 1D, avec une variation totale non croissante et un ordre ≥ 1,8 sur un profil lisse.

### 3.2 Couplage pression-vitesse (SIMPLE / SIMPLEC / PISO)

**Mécanisme de divergence reproduit (Couette)** : pression tracée sur plusieurs itérations de `master`.

| Itération | p |
|---|---|
| N = 1 | [0 ; 0,5] |
| N = 6 | [−82,6 ; 1,94] |
| N = 12 | [−170,9 ; 8,73] |

Le signe de p oscille : la correction de pression est sur-amplifiée. Il y a trois causes cumulées :
- la relaxation explicite avec un `rAU` non relaxé (A-B2) ;
- le gradient Green-Gauss faux (A-B1) ;
- l'ancrage sur une seule cellule (A-B4).

**Bilan avant / après sur le jeu complet de corrections** (`scratch_a/couette.cpp`, MAXIT = 1500, maillage canal) :

| Cas | `master` | correction (A-B1, A-B2, A-B4, A-M1…A-M7, Majumdar, C2) |
|---|---|---|
| Couette SIMPLE | ❌ `Ux momentum solve did not converge (iterations=128, residual=8192, relative=8192)`, p ~ 1e10 | ✅ 382 it, L∞(u) = 2,51e-7, continuité 4,5e-10 |
| Couette SIMPLEC | ❌ `invalid momentum diagonal` à l'itération 1 | ✅ 387 it |
| Couette PISO / PIMPLE | ❌ même divergence que SIMPLE | ✅ 382 it |
| Couette SOU / non borné | ❌ | ✅ 382 / 383 it |
| Poiseuille piloté en pression | ❌ diverge | ✅ 454 it, L∞ = 0,0039, continuité 6,9e-18 |

| Test unitaire (`scratch_a/tests_fix3a.cpp`) | `master` | correction |
|---|---|---|
| Poiseuille 16 cellules | ❌ `Ux momentum solve did not converge` | ✅ 593 it, erreur 0,00390625 |
| Poiseuille 32 cellules | — | ✅ erreur 0,000976562 (ratio 4 : **ordre 2**) |
| Majumdar : \|U(α=0,7) − U(α=0,5)\| | — | 3,95e-9 |
| `relax_implicit` : point fixe inchangé | — | écart 0 |

Tests du dépôt compilés contre la surcouche :

| Test | Résultat |
|---|---|
| `test_steady_incompressible_solver` | 11/11 (identique à `master`) |
| `test_finite_volume_transport` | 6/6 |
| `test_analytical_benchmarks` | exit 0 |
| `test_discretization_verification` | ❌ `RK2 temporal refinement insufficient` (ordre observé ~1,0), **identique sur `master`** |

Le dernier échec est préexistant et sans lien avec §3.1-3.2 (voir §3.6, B10).

L'erreur de Poiseuille (0,00390625 = 1/256 = Δy² à 16 cellules, divisée par 4 à 32 cellules) est l'erreur de discrétisation d'ordre 2 attendue, pas un défaut de couplage.

#### A-B2 🟥 Relaxation explicite + `rAU` non relaxé (cause racine) — CONFIRMÉ — *compilé + testé*
Les trois sites en cause :
- `physics/finite_volume_transport.h:567` (et le cas 1×1 l.373) : `solution(i) += controls.relaxation * (candidate(i) - solution(i));` **après** le solve ;
- `physics/steady_incompressible_solver.h:414-422` : l'appel `solve_scalar_equation(ex, ux, {..., controls.coupling.alpha_u})` ;
- `steady_incompressible_solver.h:445-472` : `rAU = V/a_P`, avec le a_P **non relaxé**.

SIMPLE suppose u* = (H + b − ∇p)/(a_P/α). Avec une relaxation a posteriori, p′ est calculé avec un `rAU` qui ne correspond pas à la matrice réellement résolue.

**Code actuel** (`steady_incompressible_solver.h:414-422`)
```cpp
const auto rx = solve_scalar_equation(ex, ux, {controls.linear_max_iterations,
                                                 controls.linear_tolerance,
                                                 controls.coupling.alpha_u});
// idem ry (ey, uy), rz (ez, uz)
```
**Code proposé 1/2** : nouvelle fonction dans `finite_volume_transport.h`, avant `scalar_equation_residual_inf`. Elle utilise l'API CSR réelle : `values_data()` non const, `row_offsets_data()`, `columns_data()`.
```cpp
// A-B2 : sous-relaxation implicite de Patankar.
//   (a_P/alpha) phi_P - sum a_nb phi_nb = b + (1-alpha)/alpha a_P phi_P^old
// A appeler AVANT solve_scalar_equation (avec controls.relaxation = 1).
// eq.diagonal recoit a_P/alpha : rAU = V/(a_P/alpha), coherent avec SIMPLE.
inline void relax_implicit(ScalarEquation& eq, const cfdx::core::Vector& phi_old, double alpha)
{
    if (!(alpha > 0.0 && alpha <= 1.0) || !std::isfinite(alpha))
        throw std::invalid_argument("relax_implicit: alpha must be in (0,1]");
    if (phi_old.size() != eq.diagonal.size() || eq.rhs.size() != eq.diagonal.size())
        throw std::invalid_argument("relax_implicit: dimension mismatch");
    if (alpha == 1.0) return;
    const auto* rows = eq.matrix.row_offsets_data();
    const auto* cols = eq.matrix.columns_data();
    auto* vals = eq.matrix.values_data();
    for (std::size_t r = 0; r < eq.diagonal.size(); ++r) {
        const double aP = eq.diagonal[r];
        const double aP_relaxed = aP / alpha;
        bool found = false;
        for (auto k = rows[r]; k < rows[r + 1]; ++k) {
            if (cols[k] == r) { vals[k] = aP_relaxed; found = true; }
        }
        if (!found)
            throw std::runtime_error("relax_implicit: missing diagonal entry at row " + std::to_string(r));
        eq.rhs(r) += (1.0 - alpha) / alpha * aP * phi_old(r);
        eq.diagonal[r] = aP_relaxed;
    }
}
```
**Code proposé 2/2** : l'appel dans `steady_incompressible_solver.h`, qui remplace l.414-422.
```cpp
// A-B2 : relaxation implicite AVANT le solve ; ex/ey/ez.diagonal = a_P/alpha
// alimente rAU plus bas. Plus aucune relaxation a posteriori.
relax_implicit(ex, ux, controls.coupling.alpha_u);
relax_implicit(ey, uy, controls.coupling.alpha_u);
relax_implicit(ez, uz, controls.coupling.alpha_u);
const ScalarSolveControls momentum_solve{controls.linear_max_iterations,
                                         controls.linear_tolerance, 1.0};
const auto rx = solve_scalar_equation(ex, ux, momentum_solve);
const auto ry = solve_scalar_equation(ey, uy, momentum_solve);
const auto rz = solve_scalar_equation(ez, uz, momentum_solve);
```
Le bloc `rAU` (l.445-472) n'est pas modifié : il lit `ex.diagonal`, qui vaut désormais a_P/α.

**Terme de Majumdar** (α-indépendance de la solution convergée) — *compilé + testé*. `make_rhie_chow_mass_flux` (l.200) reçoit des paramètres **par défaut** : les autres appelants ne changent pas.
```cpp
inline cfdx::core::Field<double,cfdx::core::Location::FACE> make_rhie_chow_mass_flux(
    /* ... mesh, geometry, U, p, rAU, rho, */
    const VelocityBoundaryConditions& bcs,
    const ScalarBoundaryConditions& pressure_bcs = {},
    const cfdx::core::Field<double,cfdx::core::Location::FACE>* previous_flux = nullptr,
    const cfdx::core::Field<double,cfdx::core::Location::CELL>* previous_U = nullptr,
    double alpha_u = 1.0)
{
    // ... validations existantes ...
    if ((previous_flux == nullptr) != (previous_U == nullptr))
        throw std::invalid_argument("make_rhie_chow_mass_flux: previous_flux and previous_U go together");
    if (!(alpha_u > 0.0 && alpha_u <= 1.0))
        throw std::invalid_argument("make_rhie_chow_mass_flux: alpha_u must be in (0,1]");
    auto flux=make_mass_flux(mesh,geometry,U,rho,bcs);
    auto gradp=gauss_gradient_with_boundary(p,mesh,geometry,pressure_bcs);   // A-M6
    // Majumdar : phi_f += (1-alpha)*(phi_f^old - rho*U_f^old.Sf)
    std::optional<Field<double, Location::FACE>> interp_old;
    if (previous_flux && alpha_u < 1.0)
        interp_old = make_mass_flux(mesh,geometry,*previous_U,rho,bcs);
    auto majumdar = [&](std::size_t f) {
        return interp_old ? (1.0-alpha_u)*((*previous_flux)(f)-(*interp_old)(f)) : 0.0;
    };
    // ... boucle des faces : flux(f) += majumdar(f); sur faces internes ET de bord (A-M7)
}
```
Ajouter `#include <optional>`. L'appel l.474 devient :
```cpp
mass_flux=make_rhie_chow_mass_flux(
    mesh,geometry,U,p,rAU,controls.density,velocity_bcs,
    pressure_bcs,&mass_flux_old,&U_old,controls.coupling.alpha_u);
```
**Tests** :
- `relax_implicit` laisse inchangé un point fixe déjà convergé (écart 0) ;
- Couette et Poiseuille convergent (tableau ci-dessus) ;
- \|U(α=0,7) − U(α=0,5)\|∞ = 3,95e-9.

#### A-B3 🟥 SIMPLEC lève une exception dès l'itération 1 — CONFIRMÉ — *corrigé par A-B2 (compilé + testé), garde-fou compilé*
`steady_incompressible_solver.h:452-468`. Le lambda `corrected_diagonal` renvoie a_P + Σ_{k≠r} a_k, et les a_nb sont stockés négatifs. Sans relaxation, cette somme vaut ~0, car la matrice est à diagonale faiblement dominante. La l.470 lève alors `"invalid momentum diagonal"`.

Avec A-B2, `eq.diagonal` et l'entrée diagonale CSR valent a_P/α, et a_P/α − Σ|a_nb| > 0 pour α < 1 : SIMPLEC converge en 387 it. Il reste le cas α_u = 1, où le défaut réapparaît. Garde-fou à placer en tête de `solve_steady_incompressible` :
```cpp
if (controls.algorithm == PressureVelocityAlgorithm::SIMPLEC &&
    !(controls.coupling.alpha_u < 1.0))
    throw std::invalid_argument(
        "solve_steady_incompressible: SIMPLEC requires alpha_u < 1 (a_P/alpha - sum|a_nb| must be > 0)");
```

#### A-B4 🟧 Ancrage de pression sur une seule cellule (déclencheur de phase9) — CONFIRMÉ — *compilé + testé*
`steady_incompressible_solver.h:590-591`. Écraser une seule cellule crée une marche de pression locale, que le gradient réinjecte dans u au correcteur suivant.

**Code actuel**
```cpp
if (!has_fixed_pressure_boundary)
    p(controls.pressure_reference_cell) = controls.pressure_reference_value;
```
**Code proposé** (même principe que `pr/415:1561-1566`)
```cpp
if (!has_fixed_pressure_boundary) {
    // A-B4 : jauge par decalage uniforme (grad p inchange), pas par
    // ecrasement d'une seule cellule.
    const double shift =
        controls.pressure_reference_value - p(controls.pressure_reference_cell);
    for (std::size_t c = 0; c < nc; ++c) p(c) += shift;
}
```
Test : Couette tout-Neumann en pression (tableau ci-dessus). Il converge avec la correction et diverge sur `master`.

**C2** (renvoi §4.1) : le statut du CG de pression (l.579) est ignoré. La vérification appliquée dans la surcouche est la suivante :
```cpp
if (rp.status != cfdx::core::SolverStatus::CONVERGED)
    throw std::runtime_error(
        "solve_steady_incompressible: pressure-correction CG did not converge (status=" +
        std::to_string(static_cast<int>(rp.status)) +
        ", iterations=" + std::to_string(rp.iterations) +
        ", relative=" + std::to_string(rp.residual_relative) + ")");
```

#### A-M1 🟧 Flux de masse recalculé à chaque itération — CONFIRMÉ — *compilé + testé*
`steady_incompressible_solver.h:364` : `auto mass_flux = make_mass_flux(...)` est reconstruit depuis u **dans** la boucle externe. Le flux conservatif corrigé par p′ à l'itération précédente est jeté, et la convection de l'itération suivante n'est donc pas à divergence nulle. La PR `pr/415:1201` corrige ce point.

**Code actuel**
```cpp
for (std::size_t iter = 1; iter <= controls.convergence.max_iterations; ++iter) {
    const auto U_old = U;
    const auto p_old = p;
    auto mass_flux = make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);
```
**Code proposé**
```cpp
// A-M1 : le flux conservatif corrige par p' est porte d'une iteration a l'autre.
auto mass_flux = make_mass_flux(mesh, geometry, U, controls.density, velocity_bcs);
for (std::size_t iter = 1; iter <= controls.convergence.max_iterations; ++iter) {
    const auto U_old = U;
    const auto p_old = p;
    const auto mass_flux_old = mass_flux;   // sert aussi au terme de Majumdar
```
Test : la continuité finale du Poiseuille vaut 6,9e-18 (tableau ci-dessus).

#### A-M6 / A-M7 🟧 Rhie-Chow — CONFIRMÉ — *compilé + testé*
- **A-M6** (l.217) : `gradp = gauss_gradient_with_boundary(p, mesh, geometry, {})`. Avec des CL vides, la valeur de face au bord est p_owner au lieu de p_b, donc ∇p est faux dans toutes les cellules de bord. **Proposé** : passer `pressure_bcs`, nouveau paramètre défauté (voir le bloc Majumdar).
- **A-M7** (l.219-220) : `if(nr<0) continue;`. Les faces de bord à pression imposée ne reçoivent pas de terme Rhie-Chow, alors que la correction p′ leur est appliquée (l.~625-635). Le flux prédit et le flux corrigé sont alors incohérents : c'est ce qui faisait stagner le Poiseuille.

**Code actuel** (l.219-220)
```cpp
const auto nr=mesh.ownership().neighbour(f);
if(nr<0) continue;
```
**Code proposé**
```cpp
const auto nr=mesh.ownership().neighbour(f);
if(nr<0) {
    // A-M7 : Rhie-Chow de bord sur les faces a pression imposee.
    const std::size_t patch=geometry.face_patch[f];
    if(patch>=mesh.boundary().n_patches()) continue;
    const auto it=pressure_bcs.find(mesh.boundary().patch(patch).name);
    if(it==pressure_bcs.end() || it->second.type!=ScalarBoundaryType::FIXED_VALUE) continue;
    const std::size_t o=mesh.ownership().owner(f);
    const Vec3 Sf=geometry.face_area_vectors[f];
    const double d=(geometry.face_centres[f]-geometry.cell_centres[o]).mag();
    if(!(d>0.0)) throw std::runtime_error("make_rhie_chow_mass_flux: degenerate boundary face");
    const Vec3 gpo{gradp.component_data(0)[o],gradp.component_data(1)[o],gradp.component_data(2)[o]};
    const double compact=(it->second.value-p(o))/d*Sf.mag();
    flux(f)-=rho*rAU[o]*(compact-gpo.dot(Sf));
    flux(f)+=majumdar(f);
    continue;
}
```
Sur les faces internes, on ajoute seulement `flux(f)+=majumdar(f);` après `flux(f)-=rho*rface*(orth-gradface);`.

Test : Poiseuille piloté en pression, L∞ = 0,0039 avec l'ordre 2 confirmé ; il divergeait sur `master`.

#### m3 🟨 PISO = SIMPLE — CONFIRMÉ — *non compilé*
La boucle des correcteurs (l.488, `pcorr` défini l.346) réutilise le même `rAU` et le même U sans recalculer H(u). Le 2ᵉ correcteur n'est donc qu'une deuxième passe de p′ sur le même opérateur. Ce n'est pas faux, mais PISO et PIMPLE ne diffèrent de SIMPLE que par le coût, ce que confirment les nombres d'itérations identiques (382 it).

**Proposé** : à partir du 2ᵉ correcteur, recalculer HbyA avec la vitesse corrigée, puis le flux prédit. Ces lignes sont à insérer en tête de la boucle `for (corr ...)` quand `corr > 0` :
```cpp
if (corr > 0) {
    // H(u) = b - sum_{nb} a_nb u_nb  (matrice ex/ey/ez deja relaxee, diagonale exclue)
    auto H = [&](const ScalarEquation& eq, std::size_t comp) {
        std::vector<double> h(nc);
        const auto* rows = eq.matrix.row_offsets_data();
        const auto* cols = eq.matrix.columns_data();
        const auto* vals = eq.matrix.values_data();
        for (std::size_t r = 0; r < nc; ++r) {
            double s = eq.rhs(r) + geometry.cell_volumes[r] * grad_p.component_data(comp)[r];
            for (auto k = rows[r]; k < rows[r + 1]; ++k)
                if (cols[k] != r) s -= vals[k] * U.component_data(comp)[cols[k]];
            h[r] = s;
        }
        return h;
    };
    const auto hx = H(ex, 0), hy = H(ey, 1), hz = H(ez, 2);
    for (std::size_t c = 0; c < nc; ++c) {
        U.component_data(0)[c] = rAU[c] / geometry.cell_volumes[c] * hx[c] - rAU[c] * grad_p.component_data(0)[c];
        U.component_data(1)[c] = rAU[c] / geometry.cell_volumes[c] * hy[c] - rAU[c] * grad_p.component_data(1)[c];
        U.component_data(2)[c] = rAU[c] / geometry.cell_volumes[c] * hz[c] - rAU[c] * grad_p.component_data(2)[c];
    }
    mass_flux = make_rhie_chow_mass_flux(mesh, geometry, U, p, rAU, controls.density,
                                         velocity_bcs, pressure_bcs);
}
```
⚠️ Le second membre de quantité de mouvement contient déjà −V·∇p. La reconstruction exacte de H dépend de la façon dont `body_x`/`grad_p` sont injectés vers l.380-410 : c'est à vérifier avant de compiler. `grad_p` doit aussi être recalculé après chaque mise à jour de p. Si PISO n'est pas prioritaire, **l'alternative honnête** est de retirer `PISO` de `PressureVelocityAlgorithm` pour le stationnaire, ou de lever une exception.

Test à ajouter : un transitoire (Taylor-Green 2D) sans relaxation. PISO à 2 correcteurs doit y donner une erreur de divergence au moins 10× plus faible que 1 correcteur, au même pas de temps.

#### m4 🟨 Critères de convergence — CONFIRMÉ — *non compilé*
`steady_incompressible_solver.h:668-757`. Quatre défauts :
- `velocity_scale = 1.0` et `pressure_scale = 1.0` servent de **planchers** (l.~690). Pour U ~ 1e-3 m/s, la variation relative est sous-estimée d'un facteur 1000.
- `momentum_rhs_scale = 1.0` sert aussi de plancher.
- `h.continuity_normalized` est calculé mais **pas utilisé** : le test l.757 porte sur `h.continuity_linf`, qui est absolu.
- `h.momentum_residual` et `h.pressure_residual` sont des résidus **du solveur linéaire**, pas du système non linéaire.

**Code actuel** (extrait l.~688-690, l.757)
```cpp
double velocity_scale = 1.0;
double pressure_scale = 1.0;
// ...
h.continuity_linf <= controls.convergence.continuity_tolerance &&
```
**Code proposé**
```cpp
// Echelles physiques, sans plancher arbitraire : tiny seulement contre 0/0.
constexpr double tiny = 1e-300;
double velocity_scale = 0.0, pressure_range_lo = +HUGE_VAL, pressure_range_hi = -HUGE_VAL;
for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
    for (std::size_t d = 0; d < 3; ++d)
        velocity_scale = std::max(velocity_scale, std::abs(U.component_data(d)[c]));
    pressure_range_lo = std::min(pressure_range_lo, p(c));
    pressure_range_hi = std::max(pressure_range_hi, p(c));
}
const double pressure_scale = std::max(pressure_range_hi - pressure_range_lo, tiny); // p est a une jauge pres
velocity_change_inf /= std::max(velocity_scale, tiny);
pressure_change_inf /= pressure_scale;

// Continuite : normalisee par le flux de reference sum|F| (definition Fluent-like).
double flux_ref = 0.0;
for (std::size_t f = 0; f < mesh.n_faces(); ++f) flux_ref += std::abs(mass_flux(f));
h.continuity_normalized = l1 / std::max(flux_ref, tiny);

// Residu non lineaire : r_it / r_it0 (premiere iteration), comme un solveur industriel.
if (iter == 1) result.reference_momentum_residual = std::max(final_momentum_residual, tiny);
h.momentum_equation_residual_relative = final_momentum_residual / result.reference_momentum_residual;
// ...
h.continuity_normalized <= controls.convergence.continuity_tolerance &&
```
`result.reference_momentum_residual` est un nouveau champ de `IncompressibleResult` (défaut 0). Si l'on veut garder la compatibilité avec les tolérances existantes de phase9, il faut ajouter un nouveau critère plutôt que modifier `continuity_tolerance` (voir §5.3).

Test à ajouter : le même Couette mis à l'échelle ×1e-3 en vitesse doit converger en autant d'itérations que l'original.

#### m2 🟨 Repli LU dense mal renseigné — CONFIRMÉ — *compilé + testé*
`physics/finite_volume_transport.h:~541-552`, dans la cascade bicgstab → gmres(64) → gauss_seidel → cg → LU dense (n ≤ 256). La version précédente de l'audit proposait des champs inexistants (`result.method`, `relative_residual`) : `SolverResult` n'a que `{status, iterations, residual, residual_relative}` (`core/linalg/cg_solver.h:40`).

Le repli rapporte :
- `iterations = n` alors qu'il n'a fait aucune itération ;
- un résidu **absolu** sous le nom `residual_relative` ;
- un statut CONVERGED testé sur le résidu absolu.

D'où le message `iterations=128, relative=8192`.

**Code actuel**
```cpp
result = {
    residual <= controls.tolerance
        ? cfdx::core::SolverStatus::CONVERGED
        : cfdx::core::SolverStatus::MAX_ITER_REACHED,
    n, residual, residual
};
```
**Code proposé**
```cpp
double rhs_norm = 0.0;
for (std::size_t i = 0; i < n; ++i)
    rhs_norm = std::max(rhs_norm, std::abs(equation.rhs(i)));
const double relative = residual / std::max(rhs_norm, 1e-300);
if (std::isfinite(residual)) {
    result = {
        relative <= controls.tolerance
            ? cfdx::core::SolverStatus::CONVERGED
            : cfdx::core::SolverStatus::MAX_ITER_REACHED,
        0,          // solveur direct : aucune iteration
        residual,   // norme infinie absolue
        relative    // normalisee par ||b||_inf
    };
}
```
Le littéral `'\\n'` des traces de la cascade (l.388 et l.405) est corrigé en `'\n'` dans la même surcouche (voir m1, §3.6).

#### m5 🟨 Code P-V mort avec formules fausses — CONFIRMÉ — *non compilé*
Correction de l'audit précédent : les deux fonctions n'ont **pas** le même défaut.

1. `physics/pressure_velocity_algorithms.h:42-57`, `rhie_chow_face_flux`. Ce fichier est inclus par `steady_incompressible_solver.h:11`, mais la fonction n'y est pas appelée. Trois défauts :
   - le coefficient `d = 0.5(1/aP_o + 1/aP_n)` n'a pas de volume, alors qu'il faudrait V/a_P ;
   - le signe est faux : le code renvoie `interp − d(p_o − p_n)A/d_f`, alors que le terme compact correct est `interp − d(p_n − p_o)A/d_f` ;
   - il n'y a pas de terme de gradient interpolé.

   Le test `test_m1_m4_physics.cpp:19` **épingle la mauvaise réponse** : `rhie_chow_face_flux(2,5,3,1,2,2,1)` attend 1,0, alors que la forme compacte correcte donne 3,0.
2. `physics/m1_m4_models.h:66-81`, `rhie_chow_mass_flux`. Le signe est correct, mais le coefficient est `d_f` (une distance, en m) à la place de `V/a_P` (en m³·s/kg) : c'est dimensionnellement faux. Elle est utilisée par 3 tests de validation.

**Proposé (préféré)** : supprimer les deux fonctions et leurs tests. Le solveur réel utilise `make_rhie_chow_mass_flux`. **Sinon**, les corriger ainsi :
```cpp
// pressure_velocity_algorithms.h
inline double rhie_chow_face_flux(double interpolated_flux,
                                  double pressure_owner, double pressure_neighbour,
                                  double grad_p_face_dot_Sf,   // (grad p)_f interpole . Sf
                                  double d_f, double rAU_face,  // rAU_f = moyenne de V/a_P
                                  double area)
{
    if (!(d_f > 0.0) || !(rAU_face > 0.0) || !(area >= 0.0) ||
        !std::isfinite(interpolated_flux + pressure_owner + pressure_neighbour + grad_p_face_dot_Sf))
        throw std::invalid_argument("rhie_chow_face_flux: invalid input");
    const double compact = (pressure_neighbour - pressure_owner) / d_f * area;
    return interpolated_flux - rAU_face * (compact - grad_p_face_dot_Sf);
}

// m1_m4_models.h : remplacer d_f par rAU_face (V/a_P) comme coefficient ; d_f reste au denominateur.
const double correction = rho * rAU_face * (dp_over_d - grad_p_face) * area;
```
Test à remplacer : avec un champ de pression linéaire, le terme Rhie-Chow doit s'annuler (compact = interpolé), donc `rhie_chow_face_flux(F, p_o, p_n, (p_n-p_o)/d*A, d, r, A) == F` à 1e-14 près. C'est une propriété physique, pas une valeur épinglée.

#### m7 🟨 Tolérances de phase9 trop lâches — CONFIRMÉ
Traité en §5.3 (tolérances des tests d'acceptation) : non répété ici. Avec les corrections ci-dessus, Couette atteint L∞ = 2,5e-7 et Poiseuille atteint l'erreur de discrétisation d'ordre 2. Les tolérances de phase9 peuvent donc être resserrées à ces niveaux.

**Artefacts de vérification** (hors dépôt, rien de commité) :

| Artefact | Rôle |
|---|---|
| `/tmp/cfdx_meta/scratch_a/fix3a/cfdx/` | surcouche corrigée (2 en-têtes physics + 6 sites A-M4) |
| `fix3a_transport.diff`, `fix3a_solver.diff` | diffs appliquables |
| `tests_fix3a.cpp` + `channel_mesh_s.inc` | A-B1, A-M3, A-M5, A-B2, Poiseuille 16/32, Majumdar |
| `couette.cpp` | Couette et Poiseuille, tous algorithmes, `master` contre correction |
| `am4_check.cpp` | A-M4 (volumes, fermeture, `dynamic_mesh`) |
| `orient_faces.py` | orientation géométrique pour l'importeur |

Commande type : `g++ -std=c++20 -O2 -I fix3a -I /tmp/cfdx_dev/src [-I /tmp/cfdx_dev/tests] <fichier>.cpp`.

---

### 3.3 Turbulence

> **Base de cette section.** Tout ce qui suit a été vérifié sur le source réel, `stevendaix/cfdx_dev@27bae126`. Chaque correctif existe sous forme de patch appliqué à une **copie** des en-têtes (`/tmp/cfdx_meta/scratch_b/s3b/overlay/src`). Le dépôt n'a pas été modifié : rien n'a été commité ni poussé.
> Les numéros de ligne « Code actuel » renvoient au fichier d'origine (côté `-` des diffs).
> Les tests de non-régression comparent chaque résultat à une **référence analytique**. Ils ont tourné sur l'original, où ils doivent échouer, puis sur la copie corrigée, où ils doivent passer.
> Chaque item porte l'un de trois statuts :
> - **compilé + testé** : le patch compile en `-Wall -Wextra -Wpedantic`, le test dédié passe sur la copie corrigée et échoue sur l'original, et les suites existantes restent vertes ;
> - **compilé** : le patch compile et les suites existantes passent, mais aucun test analytique dédié n'existe ;
> - **non compilé** : proposition de code, non construite.
>
> Récapitulatif chiffré :
> - `test_s3b_behaviour.cpp` : 0/10 sur l'original, 10/10 sur la copie corrigée ;
> - `test_s3b_new_api.cpp` : 11/11 sur la copie corrigée (il teste les fonctions ajoutées, qui n'existent pas dans l'original) ;
> - `test_s3b_b17.cpp` (sur la PR #398) : 0/2 sur l'original, 2/2 sur la copie corrigée ;
> - suites du dépôt : 19/19 sur l'original ; sur la copie corrigée, toutes passent sauf 3, qui échouent **par changement d'oracle voulu** (§3.6.3).

#### B1 🟥 Diffusion turbulente passée comme terme source — CONFIRMÉ, **compilé + testé**

**Code actuel.** Voici l'ordre réel des paramètres de `assemble_scalar_equation` (`finite_volume_transport.h:137`) :
`(mesh, geometry, mass_flux, diffusivity, su, sp, bcs, bounded, face_values*, extra_diagonal*, extra_rhs*, cell_diffusion* [const std::vector<double>*], scheme, convected*)`.
Les appels passent `nullptr, nullptr, &gamma_k` : Γ atterrit donc en **8ᵉ position pointeur**, dans `extra_rhs`. La diffusion vaut 0 et Γ·V s'ajoute au second membre comme une source :
- k-ε : `turbulence_solver.h:96-101` ;
- RNG : `turbulence_solver.h:181-182` ;
- k de k-ω : `turbulence_solver.h:231`.

Les appels ε et ω ont le même défaut.

**Code proposé.** Ajouter un `nullptr` pour décaler Γ dans `cell_diffusion`. Ce hunk est complété par ceux de B8 pour k-ω.
`src/cfdx/physics/turbulence_solver.h`
```diff
@@ -93,12 +94,13 @@
                 (controls.molecular_viscosity+nut/controls.sigma_epsilon);
         }
 
+        // face_values, extra_diagonal, extra_rhs, cell_diffusion
         auto eqk=assemble_scalar_equation(
             mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
-            nullptr,nullptr,&gamma_k);
+            nullptr,nullptr,nullptr,&gamma_k);
         auto eqe=assemble_scalar_equation(
             mesh,geometry,mass_flux,0.0,se,spe,epsilon_bcs,true,
-            nullptr,nullptr,&gamma_e);
+            nullptr,nullptr,nullptr,&gamma_e);
 
         ScalarSolveControls sc;
         sc.max_iterations=2000;
@@ -178,8 +180,8 @@
             gamma_k[i]=controls.density*(controls.molecular_viscosity+nut/controls.rng_sigma_k);
             gamma_e[i]=controls.density*(controls.molecular_viscosity+nut/controls.rng_sigma_epsilon);
         }
-        auto eqk=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,nullptr,nullptr,&gamma_k);
-        auto eqe=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,se,spe,epsilon_bcs,true,nullptr,nullptr,&gamma_e);
+        auto eqk=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,nullptr,nullptr,nullptr,&gamma_k);
+        auto eqe=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,se,spe,epsilon_bcs,true,nullptr,nullptr,nullptr,&gamma_e);
         ScalarSolveControls sc{2000,tolerance,0.7};
         cfdx::core::Vector ks(n,0.0), es(n,0.0);
         for(std::size_t i=0;i<n;++i){ks(i)=k(i);es(i)=epsilon(i);}
```

> 💡 **Recommandation.** Supprimer la classe entière de bugs en remplaçant les 4 pointeurs positionnels par une structure nommée, passée par référence :
> ```cpp
> struct ScalarExtraTerms { const std::vector<double>* diagonal=nullptr; const std::vector<double>* rhs=nullptr;
>                           const std::vector<double>* cell_diffusion=nullptr; };
> ```
> Statut : **non compilé**. Toutes les signatures appelantes changent ; le changement est à faire dans une PR dédiée.

**Test** (`test_s3b_behaviour.cpp`, cas `B1`).
- **Montage :** chaîne 1D de 10 cellules, k fixé à 1 aux deux bords, ε = 1e-3, pas de production.
- **Référence :** k ≡ 1, à 1e-3 près.
- **Résultats :** l'original donne **k = 2,716·10¹²⁹** (divergence, puisque Γ est injecté en source) ; la copie corrigée donne k = 1,000000 et 0,999999.

#### B4 🟥 Spalart-Allmaras : `fw` et S̃ faux — CONFIRMÉ, **compilé + testé**

**Code actuel.** `spalart_allmaras.h:55-57` présente trois défauts :
- `g = r + cw2 (r³ − r)` utilise **r³ au lieu de r⁶** ;
- l'exposant **^(1/6)** manque ;
- r n'est **pas plafonné à 10**.

S̃ (`turbulence_solver.h:270-275`) peut par ailleurs devenir négatif dès que fv2 < 0. L'écart relevé dans l'audit initial est confirmé et complété sur deux points :
- le plafond r ≤ 10 manque aussi ;
- ν_t de SA (`turbulence_transport.h:84`, `return second;`) **oublie fv1**.

**Code proposé.** Réécrire `fw` avec r⁶, l'exposant 1/6 et le plafond r ≤ 10 :
`src/cfdx/physics/spalart_allmaras.h`
```diff
@@ -52,9 +52,12 @@
     double destruction_coefficient(double r) const {
         if(!std::isfinite(r)||r<0.0)
             throw std::invalid_argument("invalid SA destruction ratio");
-        const double g=r+cw2*(std::pow(r,3.0)-r);
-        const double g6=std::pow(g,6.0);
-        return g*(1.0+std::pow(cw3,6.0))/(g6+std::pow(cw3,6.0));
+        // Spalart & Allmaras (1992) : r plafonne a 10, g = r + cw2 (r^6 - r),
+        // fw = g [(1 + cw3^6)/(g^6 + cw3^6)]^(1/6). fw(1) = 1, fw -> ~2 quand r -> 10.
+        const double rr=std::min(r,10.0);
+        const double g=rr+cw2*(std::pow(rr,6.0)-rr);
+        const double cw36=std::pow(cw3,6.0);
+        return g*std::pow((1.0+cw36)/(std::pow(g,6.0)+cw36),1.0/6.0);
     }
 
     double production_coefficient(double chi,double vorticity) const {
```
Ajouter le S̃ d'Allmaras 2012 (note c), qui n'est jamais négatif, et le brancher dans le solveur :
`src/cfdx/physics/turbulence_solver.h`
```diff
@@ -242,6 +290,21 @@
     return result;
 }
 
+// S-tilde d'Allmaras et al. (2012, note c) : jamais negatif, meme quand
+// fv2 < 0 (chi ~ 1..10) ou pres de la paroi. cv2 = 0.7, cv3 = 0.9.
+inline double sa_modified_vorticity(
+    double vorticity, double nu_tilde, double chi, double wall_distance,
+    const SpalartAllmarasModel& sa)
+{
+    if(!(wall_distance>0.0) || !std::isfinite(vorticity) || vorticity<0.0)
+        throw std::invalid_argument("sa_modified_vorticity: invalid inputs");
+    constexpr double cv2=0.7, cv3=0.9;
+    const double sbar=nu_tilde*sa.fv2(chi)/(sa.kappa*sa.kappa*wall_distance*wall_distance);
+    if(sbar>=-cv2*vorticity) return vorticity+sbar;
+    return vorticity+vorticity*(cv2*cv2*vorticity+cv3*sbar)/
+                     ((cv3-2.0*cv2)*vorticity-sbar);
+}
+
 inline TurbulenceTransportResult solve_spalart_allmaras_transport(
     const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
     const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
@@ -270,8 +333,10 @@
         for(std::size_t i=0;i<n;++i){
             const double wt=std::max(nu_tilde(i),0.0), nu=controls.molecular_viscosity, d=wall_distance(i), vort=std::max(vorticity(i),1e-20);
             if(!(d>0)||!std::isfinite(d)||!std::isfinite(vort)) throw std::invalid_argument("invalid SA wall/vorticity field");
-            const double chi=wt/nu, st=std::max(vort+wt*sa.fv2(chi)/(sa.kappa*sa.kappa*d*d),1e-20);
-            const double r=wt/(st*sa.kappa*sa.kappa*d*d), fw=sa.destruction_coefficient(r), ft2=sa.ft2(chi);
+            const double chi=wt/nu;
+            const double st=std::max(sa_modified_vorticity(vort,wt,chi,d,sa),1e-20);
+            const double r=std::min(wt/(st*sa.kappa*sa.kappa*d*d),10.0);
+            const double fw=sa.destruction_coefficient(r), ft2=sa.ft2(chi);
             const double prod=sa.cb1*(1-ft2)*st;
             const double destr=std::max(sa.cw1*fw-sa.cb1*ft2/(sa.kappa*sa.kappa),0.0)*wt/(d*d);
             const double grad2=gx[i]*gx[i]+gy[i]*gy[i]+gz[i]*gz[i];
```
Pour ν_t = ν̃·fv1, voir le hunk `turbulence_nu_t` de B16 (case `SPALART_ALLMARAS`).

**Test.**
- `B4_sa_fw` (behaviour) :
  - **Références :** fw(1) = 1 ; fw(10) = fw(1000) = g·(65/(g⁶+64))^(1/6), avec g = 10 + 0,3·(10⁶ − 10).
  - **Résultats :** l'original s'écarte de **2,005** ; la copie corrigée passe.
- `B4_sa_modified_vorticity` (new_api) :
  - **Références :** S̃ ≥ 0 pour χ ∈ {0,5 ; 2 ; 5 ; 10} et d ∈ {1e-4 ; 1e-2} ; S̃(ν̃ = 0) = Ω.
  - **Résultat :** passe.

#### B7 🟧 Lois de paroi — CONFIRMÉ, **compilé + testé**

**Code actuel.** Quatre défauts :
- `turbulence_solver.h:291-303` (`wall_k_epsilon`) : ε_w = Cμ^¾k^1,5/y et ω_w = √k/(Cμ^¼y), **sans κ**.
- `wall_functions.h:30` : `omega_log` utilise **β₁ = 0,075** au lieu de β* et omet κ.
- Le raccord y⁺ est codé en dur à 11 dans deux endroits :
  - `wall_functions.h:28` ;
  - `m1_m4_models.h:155`.

  L'audit initial citait 11,225, ce qui est inexact : la racine de y = ln(Ey)/κ vaut **11,53** (E = 9,793) et 11,54 (E = 9,8).

**Code proposé.**
`src/cfdx/physics/wall_functions.h`
```diff
@@ -23,17 +23,32 @@
     return std::log(E*yplus)/kappa;
 }
 
+// y+ de raccord sous-couche visqueuse / couche log : racine de y = ln(E y)/kappa.
+// kappa = 0.41 : 11.53 (E = 9.793) et 11.54 (E = 9.8), pas 11.
+inline double log_layer_crossover(double kappa=0.41,double E=9.8) {
+    if(kappa<=0.0 || E<=1.0) throw std::invalid_argument("log_layer_crossover: invalid inputs");
+    double y=11.0;
+    for(int it=0; it<200; ++it) {
+        const double next=std::log(E*y)/kappa;
+        if(std::abs(next-y)<1e-13*y) return next;
+        y=next;
+    }
+    throw std::runtime_error("log_layer_crossover: no convergence");
+}
+
 inline double turbulent_viscosity_log(double wall_distance,double friction_velocity,double nu,double kappa=0.41,double E=9.8) {
     const double yp=y_plus(wall_distance,friction_velocity,nu);
-    if(yp<=11.0) return 0.0;
+    if(yp<=log_layer_crossover(kappa,E)) return 0.0;
     const double up=u_plus_log(yp,kappa,E);
     if(up<=0.0) return 0.0;
     return std::max(0.0, wall_distance*friction_velocity/up-nu);
 }
 
-inline double omega_log(double wall_distance,double friction_velocity,double beta1=0.075) {
-    if(wall_distance<=0.0 || friction_velocity<=0.0 || beta1<=0.0) throw std::invalid_argument("omega_log: invalid inputs");
-    return friction_velocity/(std::sqrt(beta1)*wall_distance);
+// Couche log : omega = u_tau / (sqrt(beta*) kappa y). La constante est beta*
+// (0.09), pas beta1 (0.075), et kappa est au denominateur.
+inline double omega_log(double wall_distance,double friction_velocity,double beta_star=0.09,double kappa=0.41) {
+    if(wall_distance<=0.0 || friction_velocity<=0.0 || beta_star<=0.0 || kappa<=0.0) throw std::invalid_argument("omega_log: invalid inputs");
+    return friction_velocity/(std::sqrt(beta_star)*kappa*wall_distance);
 }
 
 } // namespace cfdx::physics::wall
```
`src/cfdx/physics/turbulence_solver.h`
```diff
@@ -288,25 +353,37 @@
     return result;
 }
 
+// Condition de paroi sur k seule. epsilon/omega se posent en valeur de
+// cellule via wall_epsilon_from_k / wall_omega_from_k.
+inline ScalarBoundaryCondition wall_k_fixed_value(double k_value)
+{
+    if(!std::isfinite(k_value) || k_value<0.0)
+        throw std::invalid_argument("invalid wall k value");
+    return {ScalarBoundaryType::FIXED_VALUE,k_value,0.0};
+}
+
+[[deprecated("epsilon_value est ignore : utiliser wall_k_fixed_value")]]
 inline ScalarBoundaryCondition wall_k_epsilon(double k_value, double epsilon_value)
 {
-    if(k_value<0.0 || epsilon_value<=0.0)
+    if(epsilon_value<=0.0)
         throw std::invalid_argument("invalid wall turbulence values");
-    return {ScalarBoundaryType::FIXED_VALUE,k_value,0.0};
+    return wall_k_fixed_value(k_value);
 }
 
-inline double wall_epsilon_from_k(double k, double y, double Cmu=0.09)
+// Couche log : epsilon_P = Cmu^(3/4) k^(3/2) / (kappa y).
+inline double wall_epsilon_from_k(double k, double y, double Cmu=0.09, double kappa=0.41)
 {
-    if(k<0.0 || y<=0.0 || Cmu<=0.0)
+    if(k<0.0 || y<=0.0 || Cmu<=0.0 || kappa<=0.0)
         throw std::invalid_argument("invalid wall epsilon input");
-    return std::pow(Cmu,0.75)*std::pow(std::max(k,0.0),1.5)/y;
+    return std::pow(Cmu,0.75)*std::pow(k,1.5)/(kappa*y);
 }
 
-inline double wall_omega_from_k(double k, double y, double betaStar=0.09)
+// Couche log : omega_P = sqrt(k) / (beta*^(1/4) kappa y), soit epsilon/(beta* k).
+inline double wall_omega_from_k(double k, double y, double betaStar=0.09, double kappa=0.41)
 {
-    if(k<0.0 || y<=0.0 || betaStar<=0.0)
+    if(k<0.0 || y<=0.0 || betaStar<=0.0 || kappa<=0.0)
         throw std::invalid_argument("invalid wall omega input");
-    return std::sqrt(std::max(k,0.0))/(std::sqrt(betaStar)*y);
+    return std::sqrt(k)/(std::pow(betaStar,0.25)*kappa*y);
 }
 
 } // namespace cfdx::physics
```
`src/cfdx/physics/m1_m4_models.h`
```diff
@@ -1,5 +1,6 @@
 #pragma once
 
+#include "cfdx/physics/wall_functions.h"
 #include <algorithm>
 #include <cmath>
 #include <cstddef>
@@ -152,7 +153,7 @@
 inline double wall_function_u_plus(double y_plus, double E = 9.793)
 {
     if (y_plus <= 0.0 || E <= 0.0) throw std::invalid_argument("invalid wall-function input");
-    if (y_plus < 11.0) return y_plus;
+    if (y_plus < cfdx::physics::wall::log_layer_crossover(kappa, E)) return y_plus;
     return std::log(E * y_plus) / kappa;
 }
 
```

**Test.**
- `B7_wall_values_include_kappa` (behaviour) :
  - **Références :** ε = Cμ^¾k^1,5/(κy) et ω = √k/(Cμ^¼κy), à 1e-12.
  - **Résultats :** l'original s'écarte d'un facteur 1/κ ; la copie corrigée passe.
- `B7_log_layer_crossover_and_omega` (new_api) :
  - **Références :** y_c vérifie y = ln(9,8y)/0,41 à 1e-10, et y_c = 11,53 ± 5e-3 ; `omega_log(0,01 ; 0,5) = 0,5/(0,3·0,41·0,01)`.
  - **Résultat :** passe.

#### B8 🟧 k-ω — CONFIRMÉ, **compilé + testé**

**Code actuel.**
- `turbulence_solver.h:204-240` réutilise les constantes k-ε (σ_k, σ_ε, C1, C2).
- ν_t = `a1·k/ω` (`turbulence_transport.h:80`) est la forme SST, pas celle de Wilcox.
- Le limiteur de contrainte est absent, ainsi que le terme croisé σ_d.

**Code proposé.**
- Les constantes Wilcox (2006) sont ajoutées à `TurbulenceTransportControls` et validées.
- Un noyau cellule pur est isolé pour pouvoir être testé, puis branché dans le solveur.
- `cell_gradient_dot` calcule ∇k·∇ω par Green-Gauss.

`src/cfdx/physics/turbulence_transport.h`
```diff
@@ -44,6 +44,20 @@
     double sa_cb1 = 0.1355, sa_cb2 = 0.622, sa_sigma = 2.0/3.0;
     double sa_kappa = 0.41, sa_cw2 = 0.3, sa_cw3 = 2.0, sa_cv1 = 7.1;
     double sa_ct3 = 1.2, sa_ct4 = 0.5;
+    // Wilcox (2006) k-omega. Les sigma du k-epsilon ne s'appliquent pas ici.
+    double komega_alpha = 13.0/25.0;
+    double komega_beta0 = 0.0708;
+    double komega_sigma_k = 0.6;       // sigma*
+    double komega_sigma_w = 0.5;       // sigma
+    double komega_sigma_d0 = 1.0/8.0;  // diffusion croisee si grad k . grad omega > 0
+    double komega_clim = 7.0/8.0;      // limiteur de contrainte
+    // Menter (2003) SST : sigma melanges par F1, limiteur de production.
+    double sst_sigma_k1 = 0.85, sst_sigma_k2 = 1.0;
+    double sst_sigma_w1 = 0.5, sst_sigma_w2 = 0.856;
+    double sst_production_limiter = 10.0;
+    // LES / DES algebriques.
+    double smagorinsky_Cs = 0.17;
+    double des_Cdes = 0.65;
 };
 
 inline void enforce_turbulence_bounds(
@@ -106,7 +137,11 @@
     const double values[] = {c.C_mu,c.C1,c.C2,c.beta_star,c.beta1,c.beta2,
         c.gamma1,c.gamma2,c.a1,c.rng_C_mu,c.rng_C1,c.rng_C2,c.rng_sigma_k,
         c.rng_sigma_epsilon,c.rng_eta0,c.rng_beta,c.sa_cb1,c.sa_cb2,c.sa_sigma,
-        c.sa_kappa,c.sa_cw2,c.sa_cw3,c.sa_cv1,c.sa_ct3,c.sa_ct4};
+        c.sa_kappa,c.sa_cw2,c.sa_cw3,c.sa_cv1,c.sa_ct3,c.sa_ct4,
+        c.komega_alpha,c.komega_beta0,c.komega_sigma_k,c.komega_sigma_w,
+        c.komega_sigma_d0,c.komega_clim,c.sst_sigma_k1,c.sst_sigma_k2,
+        c.sst_sigma_w1,c.sst_sigma_w2,c.sst_production_limiter,
+        c.smagorinsky_Cs,c.des_Cdes};
     for (double v : values)
         if (!std::isfinite(v)) throw std::invalid_argument("non-finite turbulence coefficient");
     if(c.C_mu<=0.0 || c.C1<0.0 || c.C2<0.0 || c.beta_star<=0.0 ||
@@ -115,7 +150,12 @@
        c.rng_C2<=0.0 || c.rng_sigma_k<=0.0 || c.rng_sigma_epsilon<=0.0 ||
        c.rng_eta0<=0.0 || c.rng_beta<=0.0 || c.sa_cb1<=0.0 || c.sa_cb2<0.0 ||
        c.sa_sigma<=0.0 || c.sa_kappa<=0.0 || c.sa_cw2<0.0 || c.sa_cw3<=0.0 ||
-       c.sa_cv1<=0.0 || c.sa_ct3<0.0 || c.sa_ct4<0.0)
+       c.sa_cv1<=0.0 || c.sa_ct3<0.0 || c.sa_ct4<0.0 ||
+       c.komega_alpha<=0.0 || c.komega_beta0<=0.0 || c.komega_sigma_k<=0.0 ||
+       c.komega_sigma_w<=0.0 || c.komega_sigma_d0<0.0 || c.komega_clim<0.0 ||
+       c.sst_sigma_k1<=0.0 || c.sst_sigma_k2<=0.0 || c.sst_sigma_w1<=0.0 ||
+       c.sst_sigma_w2<=0.0 || c.sst_production_limiter<=0.0 ||
+       c.smagorinsky_Cs<0.0 || c.des_Cdes<=0.0)
         throw std::invalid_argument("invalid turbulence coefficients");
 }
 
```
`src/cfdx/physics/turbulence_solver.h`
```diff
@@ -9,6 +9,7 @@
 #include <cmath>
 #include <cstddef>
 #include <stdexcept>
+#include <vector>
 
 namespace cfdx::physics {
 
@@ -200,6 +202,54 @@
 }
 
 
+// Termes sources et diffusivites Wilcox (2006) pour une cellule :
+//   nu_t = k/omega_t,  omega_t = max(omega, Clim S/sqrt(beta*))
+//   k : Pk - beta* rho k omega,          Gamma_k = rho(nu + sigma* k/omega)
+//   w : alpha (w/k) Pk - beta0 rho w^2 + rho sigma_d/w grad k . grad w,
+//       Gamma_w = rho(nu + sigma k/omega), sigma_d = 1/8 si grad k . grad w > 0
+// f_beta (etirement tourbillonnaire) = 1 : exact en 2D, conservatif en 3D.
+struct KOmegaCellSources {
+    double nut=0.0, sk=0.0, spk=0.0, sw=0.0, spw=0.0, gamma_k=0.0, gamma_w=0.0;
+};
+
+inline KOmegaCellSources komega2006_cell_sources(
+    double k, double omega, double strain, double grad_k_dot_grad_omega,
+    const TurbulenceTransportControls& c)
+{
+    if(!std::isfinite(k) || !std::isfinite(omega) || !std::isfinite(strain) ||
+       !std::isfinite(grad_k_dot_grad_omega))
+        throw std::invalid_argument("komega2006_cell_sources: non-finite input");
+    const double ki=std::max(k,c.k_min), wi=std::max(omega,c.omega_min);
+    const double S=std::max(strain,0.0);
+    const double omega_t=std::max(wi,c.komega_clim*S/std::sqrt(c.beta_star));
+    KOmegaCellSources s;
+    s.nut=ki/omega_t;
+    const double Pk=c.density*s.nut*S*S;
+    s.sk=Pk;
+    s.spk=-c.density*c.beta_star*wi;
+    const double sigma_d=grad_k_dot_grad_omega>0.0 ? c.komega_sigma_d0 : 0.0;
+    s.sw=c.komega_alpha*(wi/ki)*Pk + c.density*sigma_d/wi*grad_k_dot_grad_omega;
+    s.spw=-c.density*c.komega_beta0*wi;
+    s.gamma_k=c.density*(c.molecular_viscosity+c.komega_sigma_k*ki/wi);
+    s.gamma_w=c.density*(c.molecular_viscosity+c.komega_sigma_w*ki/wi);
+    return s;
+}
+
+inline std::vector<double> cell_gradient_dot(
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& a,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& b,
+    const cfdx::core::Mesh& mesh)
+{
+    const auto ga=cfdx::core::compute_gradient_gauss(a,mesh);
+    const auto gb=cfdx::core::compute_gradient_gauss(b,mesh);
+    std::vector<double> dot(mesh.n_cells());
+    for(std::size_t i=0;i<mesh.n_cells();++i)
+        dot[i]=ga.component_data(0)[i]*gb.component_data(0)[i]+
+               ga.component_data(1)[i]*gb.component_data(1)[i]+
+               ga.component_data(2)[i]*gb.component_data(2)[i];
+    return dot;
+}
+
 inline TurbulenceTransportResult solve_komega_transport(
     const cfdx::core::Mesh& mesh,const FvGeometry& geometry,
     const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
@@ -220,15 +270,13 @@
         cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),sw(n,"Sw","W/m3",1);
         cfdx::core::Field<double,cfdx::core::Location::CELL> spk(n,"Spk","kg/m3/s",1),spw(n,"Spw","kg/m3/s",1);
         std::vector<double> gk(n),gw(n);
+        const auto cross=cell_gradient_dot(k,omega,mesh);
         for(std::size_t i=0;i<n;++i){
-            const double ki=std::max(k(i),controls.k_min), wi=std::max(omega(i),controls.omega_min);
-            const double nut=controls.a1*ki/wi, P=2.0*controls.density*nut*strain_rate(i)*strain_rate(i);
-            sk(i)=P; spk(i)=-controls.density*controls.beta_star*wi;
-            sw(i)=controls.gamma1*P/std::max(nut,1e-20); spw(i)=-controls.density*controls.beta1*wi;
-            gk[i]=controls.density*(controls.molecular_viscosity+controls.sigma_k*nut);
-            gw[i]=controls.density*(controls.molecular_viscosity+controls.sigma_epsilon*nut);
+            const auto s=komega2006_cell_sources(k(i),omega(i),strain_rate(i),cross[i],controls);
+            sk(i)=s.sk; spk(i)=s.spk; sw(i)=s.sw; spw(i)=s.spw;
+            gk[i]=s.gamma_k; gw[i]=s.gamma_w;
         }
-        auto eqk=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,nullptr,nullptr,&gk);
+        auto eqk=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,nullptr,nullptr,nullptr,&gk);
         auto eqw=assemble_scalar_equation(mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,nullptr,nullptr,nullptr,&gw);
         ScalarSolveControls sc{2000,tolerance,0.7}; cfdx::core::Vector ks(n,0.0),ws(n,0.0);
         for(std::size_t i=0;i<n;++i){ks(i)=k(i);ws(i)=omega(i);}
```

**Test** (`B8_komega2006_kernel`, new_api). Montage : k = 0,1, ω = 10, ν = 1e-5.
- **S = 1** (limiteur inactif) :
  - ν_t = 0,01 ;
  - P_k = ν_tS² = 0,01 ;
  - S_ω = α·(ω/k)·P_k = (13/25)·100·0,01 ;
  - Sp_ω = −β₀ω = −0,708 ;
  - Γ_k = ν + 0,6ν_t ;
  - Γ_ω = ν + 0,5ν_t.
- **S = 100** (limiteur actif) : ν_t = k/(C_lim·S/√β*) = 0,1/291,67.
- **Terme croisé :** +σ_d/ω·(∇k·∇ω) = 0,125/10·2 lorsque le produit est positif, et 0 lorsqu'il est négatif.
- **Résultat :** tous les points passent, à 1e-12.

#### B9 🟧 SST incomplet — CONFIRMÉ, **compilé + testé**

**Code actuel.**
- `sst_solver.h:35-58` : la surcharge à 4 arguments de `compute_sst_blending` **n'a pas CD_kω**. F1 reste donc ≈ 1 hors de la couche limite.
- `sst_solver.h:61-104` : ρ manque dans la production et dans le terme croisé.
- Le limiteur P_k ≤ 10β*kω est absent.
- ν_t est calculé avec F2 = 1 en dur (`turbulence_transport.h:83-86`).

**Code proposé.** Ajouter une surcharge `compute_sst_blending(..., CDkw, sigma_w2)` et un noyau `sst_cell_sources`, puis les brancher dans les deux boucles :
`src/cfdx/physics/sst_solver.h`
```diff
@@ -35,6 +35,59 @@
 }
 
 
+// Menter (1994/2003) avec CD_kw : c'est lui qui fait basculer F1 vers 0
+// hors couche limite. Sans ce terme (surcharge ci-dessus), F1 reste proche de 1
+// dans le sillage et le modele degenere en k-omega pur, sensible a omega libre.
+inline std::pair<double,double> compute_sst_blending(
+    double k, double omega, double wall_distance,
+    double molecular_viscosity, double beta_star,
+    double grad_k_dot_grad_omega, double density, double sigma_w2)
+{
+    if (!std::isfinite(grad_k_dot_grad_omega) || !(density>0.0) || !(sigma_w2>0.0))
+        throw std::invalid_argument("compute_sst_blending: invalid CDkw inputs");
+    const auto base = compute_sst_blending(k, omega, wall_distance, molecular_viscosity, beta_star);
+    const double ki = std::max(k, 0.0);
+    const double wi = std::max(omega, 1e-20);
+    const double y = wall_distance;
+    const double cd_kw = std::max(2.0*density*sigma_w2/wi*grad_k_dot_grad_omega, 1e-10);
+    const double arg1 = std::min(
+        std::max(std::sqrt(ki)/(beta_star*wi*y), 500.0*molecular_viscosity/(y*y*wi)),
+        4.0*density*sigma_w2*ki/(cd_kw*y*y));
+    return {std::clamp(std::tanh(std::pow(arg1, 4.0)), 0.0, 1.0), base.second};
+}
+
+// Sources SST-2003 d'une cellule (NASA TMR "SST-2003") :
+//   nu_t = a1 k / max(a1 w, S F2),  Pk = min(rho nu_t S^2, 10 beta* rho k w)
+//   k : Pk - beta* rho k w                     Gamma_k = rho(nu + sigma_k nu_t)
+//   w : gamma rho S^2 - beta rho w^2 + 2(1-F1) rho sigma_w2 / w grad k . grad w
+//                                             Gamma_w = rho(nu + sigma_w nu_t)
+// sigma_k, sigma_w, beta, gamma melanges par F1. La part negative du terme
+// croise est implicitee (spw), la positive reste explicite (sw).
+inline KOmegaCellSources sst_cell_sources(
+    double k, double omega, double strain, double F1, double F2,
+    double grad_k_dot_grad_omega, const TurbulenceTransportControls& c)
+{
+    if(!std::isfinite(k) || !std::isfinite(omega) || !std::isfinite(strain) ||
+       !std::isfinite(F1) || !std::isfinite(F2) || !std::isfinite(grad_k_dot_grad_omega))
+        throw std::invalid_argument("sst_cell_sources: non-finite input");
+    const double ki=std::max(k,c.k_min), wi=std::max(omega,c.omega_min);
+    const double S=std::max(strain,0.0);
+    const double f1=std::clamp(F1,0.0,1.0), f2=std::clamp(F2,0.0,1.0);
+    auto blend=[f1](double a, double b){ return f1*a+(1.0-f1)*b; };
+    KOmegaCellSources s;
+    s.nut=c.a1*ki/std::max(c.a1*wi,S*f2);
+    const double Pk=std::min(c.density*s.nut*S*S,
+                             c.sst_production_limiter*c.beta_star*c.density*ki*wi);
+    s.sk=Pk;
+    s.spk=-c.density*c.beta_star*wi;
+    const double cross=2.0*(1.0-f1)*c.density*c.sst_sigma_w2/wi*grad_k_dot_grad_omega;
+    s.sw=blend(c.gamma1,c.gamma2)*c.density*S*S + std::max(cross,0.0);
+    s.spw=-c.density*blend(c.beta1,c.beta2)*wi + std::min(cross,0.0)/wi;
+    s.gamma_k=c.density*(c.molecular_viscosity+blend(c.sst_sigma_k1,c.sst_sigma_k2)*s.nut);
+    s.gamma_w=c.density*(c.molecular_viscosity+blend(c.sst_sigma_w1,c.sst_sigma_w2)*s.nut);
+    return s;
+}
+
 inline TurbulenceTransportResult solve_sst_transport(
     const cfdx::core::Mesh& mesh,
     const FvGeometry& geometry,
@@ -61,8 +114,6 @@
     for(std::size_t iter=1;iter<=max_iterations;++iter) {
         auto oldk=k;
         auto oldw=omega;
-        cfdx::core::Field<double,cfdx::core::Location::CELL> P(
-            n,"Pk","W/m3",1);
         cfdx::core::Field<double,cfdx::core::Location::CELL> sk(
             n,"Sk","W/m3",1);
         cfdx::core::Field<double,cfdx::core::Location::CELL> sw(
@@ -73,32 +124,19 @@
             n,"Spw","kg/m3/s",1);
         std::vector<double> gamma_k(n),gamma_w(n);
 
+        const auto cross=cell_gradient_dot(k,omega,mesh);
         for(std::size_t i=0;i<n;++i) {
-            const double ki=std::max(k(i),controls.k_min);
-            const double wi=std::max(omega(i),controls.omega_min);
-            const double f1=std::clamp(F1(i),0.0,1.0);
-            const double f2=std::clamp(F2(i),0.0,1.0);
-            const double beta=f1*controls.beta1+(1.0-f1)*controls.beta2;
-            const double gamma=f1*(5.0/9.0)+(1.0-f1)*0.44;
-            const double nut=controls.a1*ki/
-                std::max(controls.a1*wi,strain_rate(i)*f2);
-            P(i)=2.0*nut*strain_rate(i)*strain_rate(i);
-            gamma_k[i]=controls.density*
-                (controls.molecular_viscosity+controls.sigma_k*nut);
-            gamma_w[i]=controls.density*
-                (controls.molecular_viscosity+controls.sigma_epsilon*nut);
-            sk(i)=P(i);
-            spk(i)=-controls.density*controls.beta_star*wi;
-            sw(i)=gamma*P(i)/std::max(nut,1e-20);
-            spw(i)=-controls.density*beta*wi;
+            const auto s=sst_cell_sources(k(i),omega(i),strain_rate(i),F1(i),F2(i),cross[i],controls);
+            sk(i)=s.sk; spk(i)=s.spk; sw(i)=s.sw; spw(i)=s.spw;
+            gamma_k[i]=s.gamma_k; gamma_w[i]=s.gamma_w;
         }
 
         auto eqk=assemble_scalar_equation(
             mesh,geometry,mass_flux,0.0,sk,spk,k_bcs,true,
-            nullptr,nullptr,&gamma_k);
+            nullptr,nullptr,nullptr,&gamma_k);
         auto eqw=assemble_scalar_equation(
             mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
-            nullptr,nullptr,&gamma_w);
+            nullptr,nullptr,nullptr,&gamma_w);
 
         ScalarSolveControls sc{2000,tolerance,0.7};
         cfdx::core::Vector k_solution(n, 0.0);
@@ -151,8 +189,10 @@
     for(std::size_t iter=1;iter<=max_iterations;++iter){
         auto old_k=k,old_w=omega;
         cfdx::core::Field<double,cfdx::core::Location::CELL> F1(n,"F1","1",1),F2(n,"F2","1",1);
+        const auto cross=cell_gradient_dot(k,omega,mesh);
         for(std::size_t i=0;i<n;++i){
-            const auto b=compute_sst_blending(k(i),omega(i),wall_distance(i),controls.molecular_viscosity,controls.beta_star);
+            const auto b=compute_sst_blending(k(i),omega(i),wall_distance(i),controls.molecular_viscosity,
+                controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
             F1(i)=b.first; F2(i)=b.second;
         }
         const auto inner=solve_sst_transport(mesh,geometry,mass_flux,k,omega,strain_rate,F1,F2,controls,k_bcs,omega_bcs,1,tolerance);
```

**Test** (`B9_sst_kernel`, new_api). Montage : ρ = 2.
- F1 = 1, F2 = 0 : ν_t = k/ω = 0,01 et P_k = ρν_tS² = 0,02.
- S = 1000 : P_k est écrêté à 10β*ρkω = 1,8.
- F1 = 0 : le terme croisé vaut 2ρσ_ω2/ω·(∇k·∇ω) = 2·2·0,856/10·3.
- F2 = 1, S grand : ν_t = a1k/S.
- **Contrôle de sens :** un grand CD_kω fait baisser F1.
- **Résultat :** passe.

#### B15 🟨 Contrat sur S — CONFIRMÉ (facteur 2 selon l'appelant), **compilé + testé**

**Code actuel.**
- Aucune documentation ne dit si `strain_rate` vaut √(2SijSij) ou SijSij.
- `turbulence_solver.h:33` et `turbulence_solver.h:170` n'ont pas ρ dans P_k = ν_tS².

**Code proposé.**
- Le contrat S = √(2SijSij) est documenté.
- Un calculateur de référence `strain_rate_magnitude(grad_u)` est ajouté.
- La production s'écrit ρν_tS².

`src/cfdx/physics/turbulence.h`
```diff
@@ -2,6 +2,7 @@
 
 #include "cfdx/core/field/field.h"
 #include <algorithm>
+#include <array>
 #include <cmath>
 #include <cstddef>
 #include <stdexcept>
@@ -52,6 +53,23 @@
     return std::sqrt(2.0 * (s.x*s.x + s.y*s.y + s.z*s.z));
 }
 
+// Contrat CFDX sur S (tous les modeles a viscosite turbulente) :
+//   Sij = (dui/dxj + duj/dxi)/2,  S = sqrt(2 Sij Sij),  P_k = rho nu_t S^2.
+// grad_u[i][j] = dui/dxj. Cisaillement pur du/dy = g : S = |g|.
+// La surcharge Vec3 ci-dessus ne voit que la diagonale (pas de cisaillement).
+inline double strain_rate_magnitude(const std::array<std::array<double,3>,3>& grad_u)
+{
+    double sum = 0.0;
+    for (std::size_t i = 0; i < 3; ++i)
+        for (std::size_t j = 0; j < 3; ++j) {
+            const double sij = 0.5 * (grad_u[i][j] + grad_u[j][i]);
+            if (!std::isfinite(sij))
+                throw std::invalid_argument("strain_rate_magnitude: non-finite gradient");
+            sum += sij * sij;
+        }
+    return std::sqrt(2.0 * sum);
+}
+
 inline double smagorinsky_eddy_viscosity(double delta, double strain,
                                          double Cs = 0.17)
 {
```
`src/cfdx/physics/turbulence_solver.h`
```diff
@@ -30,7 +31,7 @@
     for(std::size_t i=0;i<mesh.n_cells();++i) {
         const double nut=c.C_mu*std::max(k(i),c.k_min)*std::max(k(i),c.k_min)/
                          std::max(epsilon(i),c.epsilon_min);
-        production(i)=c.density*2.0*nut*strain_rate(i)*strain_rate(i);
+        production(i)=c.density*nut*strain_rate(i)*strain_rate(i);
     }
 }
 
@@ -167,7 +169,7 @@
             const double ki=std::max(k(i),controls.k_min), ei=std::max(epsilon(i),controls.epsilon_min);
             const double nut=controls.rng_C_mu*ki*ki/ei;
             const double S=std::max(strain_rate(i),0.0);
-            const double P=2.0*controls.density*nut*S*S;
+            const double P=controls.density*nut*S*S;
             const double eta=S*ki/ei;
             const double C1star=controls.rng_C1-
                 eta*(1.0-eta/controls.rng_eta0)/(1.0+controls.rng_beta*eta*eta*eta);
```

**Test** (`B15_strain_contract_pure_shear`, new_api).
- **Référence :** pour un cisaillement pur ∂u/∂y = 7, S = 7 à 1e-14, et P_k = ρC_μ(k²/ε)S² avec ρ = 1,2.
- **Résultat :** passe.

#### B16 🟨 LES et DES — CONFIRMÉ, **compilé + testé**

**Code actuel.** `turbulence_transport.h:85` (Smagorinsky) et `turbulence_transport.h:90` (DES) calculent `std::cbrt(1.0)`, soit **Δ = 1 m quel que soit le maillage**. La proposition de l'audit initial (`geometry.cell_volumes[c]`) ne compile pas : la fonction ne reçoit pas la géométrie.

**Code proposé.** `turbulence_nu_t` reçoit `cell_volume`, qui devient obligatoire en LES et DES. Le même hunk corrige ν_t de k-ω (B8), de SA (fv1, B4) et de SST (F2, B9).
`src/cfdx/physics/turbulence_transport.h`
```diff
@@ -59,37 +73,54 @@
     }
 }
 
+// Contrat : strain = S = sqrt(2 Sij Sij). second = epsilon (k-epsilon, RNG),
+// omega (k-omega, SST) ou nu_tilde (SA). cell_volume est obligatoire en
+// LES/DES (Delta = V^(1/3)) ; F2 est le second melange SST.
 inline double turbulence_nu_t(
     double k, double second, double strain, double wall_distance,
-    const TurbulenceTransportControls& c)
+    const TurbulenceTransportControls& c,
+    double cell_volume = -1.0, double F2 = 1.0)
 {
+    const bool omega_based =
+        c.model==TurbulenceModel::SST || c.model==TurbulenceModel::KOMEGA;
     k=std::max(k,c.k_min);
-    second=std::max(second,
-        c.model==TurbulenceModel::SST ? c.omega_min : c.epsilon_min);
+    if(c.model!=TurbulenceModel::SPALART_ALLMARAS)
+        second=std::max(second, omega_based ? c.omega_min : c.epsilon_min);
+    auto les_delta=[&]() {
+        if(wall_distance<=0.0) throw std::invalid_argument("wall distance must be positive");
+        if(!(cell_volume>0.0) || !std::isfinite(cell_volume))
+            throw std::invalid_argument("LES/DES eddy viscosity requires the cell volume");
+        return std::cbrt(cell_volume);
+    };
     switch(c.model) {
     case TurbulenceModel::LAMINAR: return 0.0;
     case TurbulenceModel::KEPSILON:
         return c.C_mu*k*k/second;
     case TurbulenceModel::RNG_KEPSILON:
         return c.rng_C_mu*k*k/second;
-    case TurbulenceModel::KOMEGA:
-        return c.a1*k/second;
-    case TurbulenceModel::SPALART_ALLMARAS:
-        return second;
-    case TurbulenceModel::SST: {
-        const double F2=1.0;
-        return c.a1*k/std::max(c.a1*second,strain*F2);
+    case TurbulenceModel::KOMEGA: {
+        // Wilcox 2006 : nu_t = k / omega_tilde, omega_tilde >= Clim S / sqrt(beta*).
+        const double omega_tilde=std::max(
+            second, c.komega_clim*std::max(strain,0.0)/std::sqrt(c.beta_star));
+        return k/omega_tilde;
     }
-    case TurbulenceModel::SMAGORINSKY: {
-        if(wall_distance<=0.0) throw std::invalid_argument("wall distance must be positive");
-        const double delta=std::cbrt(1.0);
-        return smagorinsky_eddy_viscosity(delta,strain);
+    case TurbulenceModel::SPALART_ALLMARAS: {
+        // nu_t = nu_tilde fv1(chi), chi = nu_tilde/nu.
+        const double nt=std::max(second,0.0);
+        if(!(c.molecular_viscosity>0.0))
+            throw std::invalid_argument("SA eddy viscosity requires nu > 0");
+        const double chi3=std::pow(nt/c.molecular_viscosity,3.0);
+        return nt*chi3/(chi3+std::pow(c.sa_cv1,3.0));
     }
-    case TurbulenceModel::DES: {
-        if(wall_distance<=0.0) throw std::invalid_argument("wall distance must be positive");
-        const double length=std::min(wall_distance,0.65*std::cbrt(1.0));
-        return length*length*strain;
+    case TurbulenceModel::SST: {
+        if(F2<0.0 || F2>1.0) throw std::invalid_argument("SST F2 must lie in [0,1]");
+        return c.a1*k/std::max(c.a1*second,std::max(strain,0.0)*F2);
     }
+    case TurbulenceModel::SMAGORINSKY:
+        return smagorinsky_eddy_viscosity(les_delta(),strain,c.smagorinsky_Cs);
+    case TurbulenceModel::DES:
+        return des_eddy_viscosity(les_delta(),wall_distance,strain,
+                                  c.smagorinsky_Cs,c.des_Cdes);
     }
     throw std::invalid_argument("unknown turbulence model");
 }
```

**Test** (`B16_les_filter_width`, new_api).
- **Référence :** V = 1e-6 m³ donne Δ = 1e-2 et ν_t = (C_s·Δ)²·S.
- **Contrôle d'erreur :** l'appel sans volume lève `std::invalid_argument`.
- **Résultat :** passe.
- **Oracle existant modifié :** `test_phase10_model_registry:43` (§3.6.3).

### 3.4 Énergie, thermophysique, CHT

#### B2 🟥 `cp` absent de la convection d'énergie — CONFIRMÉ, **compilé + testé**

**Code actuel.** L'équation est écrite en T, mais `assemble_energy_equation` (`energy_solver.h:77-78`) passe `mass_flux` [kg/s] tel quel. Le terme convectif vaut donc ṁ·T au lieu de ṁ·cp·T. Deux autres endroits sont incohérents :
- le bilan `energy_balance_relative` (`energy_solver.h:146`) ;
- Rosseland (`radiation_advanced.h:590`), qui est un **nouveau constat** (voir B18).

**Code proposé.** Ajouter `enthalpy_face_flux`, l'appliquer dans l'assemblage et dans le bilan, et ouvrir `source_implicit` (Sp ≤ 0), dont B3 a besoin :
`src/cfdx/physics/energy_solver.h`
```diff
@@ -43,6 +43,31 @@
     if (c.dt < 0.0) throw std::invalid_argument("energy time step must be >= 0");
 }
 
+// L'equation est ecrite en temperature : div(m_dot cp T) - div(k grad T) = S.
+// Le flux convectif transporte de l'enthalpie, donc m_dot [kg/s] x cp [J/kg/K].
+inline cfdx::core::Field<double,cfdx::core::Location::FACE> enthalpy_face_flux(
+    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux, double cp)
+{
+    if(!(cp>0.0) || !std::isfinite(cp))
+        throw std::invalid_argument("enthalpy_face_flux: cp must be positive");
+    auto flux=mass_flux;
+    for(std::size_t f=0;f<flux.size();++f) flux(f)=mass_flux(f)*cp;
+    return flux;
+}
+
+inline void validate_implicit_source(
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit,
+    std::size_t n_cells)
+{
+    if(!source_implicit) return;
+    if(source_implicit->size()!=n_cells)
+        throw std::invalid_argument("energy implicit source size mismatch");
+    for(std::size_t i=0;i<n_cells;++i)
+        if(!std::isfinite((*source_implicit)(i)) || (*source_implicit)(i)>0.0)
+            throw std::invalid_argument("energy implicit source must be finite and <= 0");
+}
+
+// source_implicit = Sp [W/m3/K] : S = Su + Sp T, Sp <= 0 (diagonale renforcee).
 inline ScalarEquation assemble_energy_equation(
     const cfdx::core::Mesh& mesh,
     const FvGeometry& geometry,
@@ -51,11 +76,13 @@
     const cfdx::core::Field<double,cfdx::core::Location::CELL>& old_temperature,
     const EnergySolverControls& c,
     const ScalarBoundaryConditions& bcs = {},
-    const ScalarBoundaryFaceValues* face_values = nullptr)
+    const ScalarBoundaryFaceValues* face_values = nullptr,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr)
 {
     validate_energy_controls(c);
     if (source.size()!=mesh.n_cells() || old_temperature.size()!=mesh.n_cells())
         throw std::invalid_argument("energy field size mismatch");
+    validate_implicit_source(source_implicit,mesh.n_cells());
 
     cfdx::core::Field<double,cfdx::core::Location::CELL> su(
         mesh.n_cells(),"energy_source","W/m3",1);
@@ -66,15 +93,16 @@
 
     for(std::size_t i=0;i<mesh.n_cells();++i) {
         su(i)=source(i);
-        sp(i)=0.0;
+        sp(i)=source_implicit ? (*source_implicit)(i) : 0.0;
         if(c.dt>0.0) {
             transient_diag[i]=c.density*c.cp*geometry.cell_volumes[i]/c.dt;
             transient_rhs[i]=transient_diag[i]*old_temperature(i);
         }
     }
 
+    const auto enthalpy_flux=enthalpy_face_flux(mass_flux,c.cp);
     return assemble_scalar_equation(
-        mesh,geometry,mass_flux,c.conductivity,su,sp,bcs,true,
+        mesh,geometry,enthalpy_flux,c.conductivity,su,sp,bcs,true,
         face_values,&transient_diag,&transient_rhs);
 }
 
@@ -87,15 +115,21 @@
     const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
     const EnergySolverControls& controls,
     const ScalarBoundaryConditions& bcs,
-    const ScalarBoundaryFaceValues* face_values = nullptr)
+    const ScalarBoundaryFaceValues* face_values = nullptr,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr,
+    const std::vector<double>* cell_conductivity = nullptr)
 {
+    validate_implicit_source(source_implicit,mesh.n_cells());
+    if(cell_conductivity && cell_conductivity->size()!=mesh.n_cells())
+        throw std::invalid_argument("energy balance conductivity size mismatch");
     double net_flux = 0.0;
     double source_total = 0.0;
     double accumulation = 0.0;
     const auto& own = mesh.ownership();
 
     for(std::size_t c=0;c<mesh.n_cells();++c) {
-        source_total += source(c) * geometry.cell_volumes[c];
+        const double sp = source_implicit ? (*source_implicit)(c) : 0.0;
+        source_total += (source(c) + sp*temperature(c)) * geometry.cell_volumes[c];
         if(controls.dt > 0.0)
             accumulation += controls.density * controls.cp *
                 geometry.cell_volumes[c] *
@@ -139,11 +173,11 @@
                 Tf=temperature(o)+bc.gradient*d;
         }
 
-        const double k=controls.conductivity;
+        const double k=cell_conductivity ? (*cell_conductivity)[o] : controls.conductivity;
         double conductive=-k*bc.gradient*area;
         if(d>0.0 && (bc.type==ScalarBoundaryType::FIXED_VALUE || face_override_valid))
             conductive=-k*(Tf-temperature(o))/d*area;
-        const double convective=F*(F>=0.0 ? temperature(o) : Tf);
+        const double convective=controls.cp*F*(F>=0.0 ? temperature(o) : Tf);
         net_flux += convective + conductive;
     }
 
@@ -161,7 +195,8 @@
     const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
     const EnergySolverControls& controls = {},
     const ScalarBoundaryConditions& bcs = {},
-    const ScalarBoundaryFaceValues* face_values = nullptr)
+    const ScalarBoundaryFaceValues* face_values = nullptr,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr)
 {
     validate_energy_controls(controls);
     if(temperature.size()!=mesh.n_cells() || source.size()!=mesh.n_cells())
@@ -171,7 +206,7 @@
     auto old=temperature;
     for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
         auto eq=assemble_energy_equation(
-            mesh,geometry,mass_flux,source,old,controls,bcs,face_values);
+            mesh,geometry,mass_flux,source,old,controls,bcs,face_values,source_implicit);
         cfdx::core::Vector candidate(temperature.size(),0.0);
         for(std::size_t i=0;i<temperature.size();++i)
             candidate(i)=temperature(i);
@@ -185,7 +220,7 @@
             temperature(i)=candidate(i);
 
         const double imbalance=energy_balance_relative(
-            mesh,geometry,mass_flux,temperature,old,source,controls,bcs,face_values);
+            mesh,geometry,mass_flux,temperature,old,source,controls,bcs,face_values,source_implicit);
         result.history.push_back({iter,res,imbalance});
         result.iterations=iter;
 
```

**Test** (`B2_energy_convection_uses_cp`, behaviour).
- **Montage :** chaîne de 40 cellules, ρ = 1, cp = 1000, k = 10, u = 0,01. Le Péclet réel est Pe = ρcpuL/k = 1.
- **Référence :** solution exacte de convection-diffusion 1D, (e^{Pe·x} − 1)/(e^{Pe} − 1), à 5e-3.
- **Résultats :** l'original donne T(19) = **0,4874**, ce qui correspond à Pe = 0,001, c'est-à-dire de la diffusion pure ; la copie corrigée donne 0,3670.
- **Complément :** `B2_enthalpy_face_flux` (new_api) vérifie −2·1005 = −2010 et le rejet de cp = 0.

#### B6 🟧 Enthalpie tabulée et bornée — CONFIRMÉ, étendu, **compilé + testé**

**Code actuel.** `thermophysical_models.h:162-178` a deux défauts d'extrapolation. Hors de la table, l'intégrale est prise sur [lo − (lo − a), lo], ce qui est faux :
- en extrapolation basse, quand b < lo ;
- en extrapolation haute, quand a > hi.

Nouveau constat : avec `enforce_bounds`, les formes fermées (`thermophysical_models.h:119-139`) intègrent la loi cp **sans son écrêtage**. Or `minimum` et `maximum` bornent la **valeur** de cp, pas T. Enfin, `thermophysical_models.h:135` a une indentation trompeuse : `integral=0.0` est placé sous un `throw`, dans la branche POLYNOMIAL.

**Code proposé.**
`src/cfdx/physics/thermophysical_models.h`
```diff
@@ -116,7 +116,15 @@
             throw std::invalid_argument("enthalpy: invalid inputs");
         const double a=std::min(T,Tref), b=std::max(T,Tref);
         double integral=0.0;
-        if (heat_capacity.model == ScalarPropertyModel::CONSTANT) {
+        // Les formes fermees ignorent l'ecretage [minimum, maximum] : avec
+        // enforce_bounds, on integre cp(T) tel qu'il est reellement evalue.
+        if (heat_capacity.enforce_bounds) {
+            const std::size_t n=static_cast<std::size_t>(std::ceil(intervals));
+            const double h=(b-a)/static_cast<double>(n);
+            double s=cp(a)+cp(b);
+            for(std::size_t i=1;i<n;++i) s+=2*cp(a+h*static_cast<double>(i));
+            integral=0.5*h*s;
+        } else if (heat_capacity.model == ScalarPropertyModel::CONSTANT) {
             integral=heat_capacity.reference_value*(b-a);
         } else if (heat_capacity.model == ScalarPropertyModel::LINEAR) {
             const double c0=heat_capacity.reference_value;
@@ -140,7 +148,6 @@
             auto integrate_in_range = [&](double lo,double hi) {
                 if (hi<=lo) return 0.0;
                 double result=0.0;
-                double x0=lo, y0=tab.evaluate(lo,"heat capacity table");
                 for (std::size_t i=0;i+1<tab.temperature.size();++i) {
                     const double l=std::max(lo,tab.temperature[i]);
                     const double r=std::min(hi,tab.temperature[i+1]);
@@ -150,7 +157,6 @@
                     const double yr=tab.value[i] + (r-tab.temperature[i]) *
                         (tab.value[i+1]-tab.value[i])/(tab.temperature[i+1]-tab.temperature[i]);
                     result += 0.5*(yl+yr)*(r-l);
-                    x0=r; y0=yr;
                 }
                 return result;
             };
@@ -159,23 +165,22 @@
                     throw std::out_of_range("enthalpy: heat-capacity table outside range");
                 const double lo=tab.temperature.front(), hi=tab.temperature.back();
                 integral=integrate_in_range(std::max(a,lo),std::min(b,hi));
+                // Sous lo : cp(x) = v0 + s0 (x - lo), integre sur [a, min(b,lo)].
                 if (a<lo) {
-                    if (tab.extrapolation == ExtrapolationPolicy::CLAMP)
-                        integral += tab.value.front()*(lo-a);
-                    else {
-                        const double slope=(tab.value[1]-tab.value[0])/(tab.temperature[1]-tab.temperature[0]);
-                        integral += tab.value.front()*(lo-a) + 0.5*slope*((lo-a)*(lo-a));
-                    }
+                    const double c=std::min(b,lo);
+                    const double s0=tab.extrapolation == ExtrapolationPolicy::CLAMP ? 0.0 :
+                        (tab.value[1]-tab.value[0])/(tab.temperature[1]-tab.temperature[0]);
+                    integral += tab.value.front()*(c-a) +
+                        0.5*s0*((c-lo)*(c-lo)-(a-lo)*(a-lo));
                 }
+                // Au-dessus de hi : cp(x) = vN + sN (x - hi), integre sur [max(a,hi), b].
                 if (b>hi) {
-                    if (tab.extrapolation == ExtrapolationPolicy::CLAMP)
-                        integral += tab.value.back()*(b-hi);
-                    else {
-                        const double slope=(tab.value.back()-tab.value[tab.value.size()-2])/
-                                           (tab.temperature.back()-tab.temperature[tab.temperature.size()-2]);
-                        const double d=b-hi;
-                        integral += tab.value.back()*d + 0.5*slope*d*d;
-                    }
+                    const double d=std::max(a,hi);
+                    const std::size_t m=tab.value.size();
+                    const double sN=tab.extrapolation == ExtrapolationPolicy::CLAMP ? 0.0 :
+                        (tab.value[m-1]-tab.value[m-2])/(tab.temperature[m-1]-tab.temperature[m-2]);
+                    integral += tab.value.back()*(b-d) +
+                        0.5*sN*((b-hi)*(b-hi)-(d-hi)*(d-hi));
                 }
             } else {
                 integral=integrate_in_range(a,b);
```

**Test** (`B6_*`, behaviour).
- **Cas table linéaire :**
  - h(250 → 300) = −48 750 ; l'original donne **−51 250**, soit un écart de 2 500 ;
  - h(200 → 250) = −46 250 ;
  - h(450 → 400) = 1100·50 + ½·50².
- **Cas borné** (min 100, max 800, CONSTANT 1000) :
  - cp(350) = 800 ;
  - h(400 ← 320) = **64 000**, alors que la forme fermée donnait 80 000.

#### B14 🟧 Équation d'état incompressible — CONFIRMÉ, **compilé + testé**

**Code actuel.** Il y a deux copies de l'EOS incompressible :
- `physics/equation_of_state.h:122-143` renvoie c = √(ρ/β), dp/dρ|s = 1/(βρ) et dp/dT|ρ = ρβ ;
- `thermodynamics/equation_of_state.h:103-122` fait la même chose et ne valide pas ses paramètres.

Ces formules traitent β, coefficient de dilatation **thermique**, comme une compressibilité **isotherme**. À ρ constant, (∂ρ/∂p)_T = 0 : c est infini et (∂p/∂T)_ρ n'est pas défini.

**Code proposé.**
`src/cfdx/physics/equation_of_state.h`
```diff
@@ -119,11 +119,11 @@
         return params_.Cp * std::log(T / params_.T_ref);
     }
 
+    // Densite constante (beta ne sert qu'a la flottabilite de Boussinesq) :
+    // (drho/dp)_T = 0, donc c = sqrt((dp/drho)_s) est infinie, et
+    // (dp/dT)_rho = -(drho/dT)_p / (drho/dp)_T n'est pas defini.
     double speed_of_sound(double /*p*/, double /*T*/) const override {
-        if (params_.beta <= 0) {
-            return std::numeric_limits<double>::infinity();
-        }
-        return std::sqrt(params_.rho / params_.beta);
+        return std::numeric_limits<double>::infinity();
     }
 
     double temperature_from_enthalpy(double /*p*/, double h) const override {
@@ -135,11 +135,11 @@
     }
 
     double dp_drho_s(double /*p*/, double /*T*/) const override {
-        return 1.0 / (params_.beta * params_.rho);
+        return std::numeric_limits<double>::infinity();
     }
 
     double dp_dT_rho(double /*p*/, double /*T*/) const override {
-        return params_.rho * params_.beta;
+        throw std::domain_error("IncompressibleEOS: (dp/dT)_rho undefined at constant density");
     }
 
     double drho_dp_T(double /*p*/, double /*T*/) const override {
```
`src/cfdx/thermodynamics/equation_of_state.h`
```diff
@@ -25,6 +25,7 @@
 #include <memory>
 #include <stdexcept>
 #include <cmath>
+#include <limits>
 
 namespace cfdx {
 namespace thermodynamics {
@@ -88,7 +89,15 @@
     IncompressibleParams params_;
 
 public:
-    IncompressibleEOS(const IncompressibleParams& params = {}) : params_(params) {}
+    IncompressibleEOS(const IncompressibleParams& params = {}) : params_(params) { validate_params(); }
+
+    void validate_params() const {
+        if (!std::isfinite(params_.rho) || !(params_.rho > 0.0) ||
+            !std::isfinite(params_.Cp) || !(params_.Cp > 0.0) ||
+            !std::isfinite(params_.T_ref) || !(params_.T_ref > 0.0) ||
+            !std::isfinite(params_.beta) || params_.beta < 0.0)
+            throw std::invalid_argument("invalid incompressible EOS parameters");
+    }
     EquationOfStateType type() const override { return EquationOfStateType::INCOMPRESSIBLE; }
 
     double density(double /*p*/, double /*T*/) const override {
@@ -103,8 +112,11 @@
         return params_.Cp * std::log(T / params_.T_ref);
     }
 
+    // Densite constante (beta ne sert qu'a la flottabilite de Boussinesq) :
+    // (drho/dp)_T = 0, donc c = sqrt((dp/drho)_s) est infinie, et
+    // (dp/dT)_rho = -(drho/dT)_p / (drho/dp)_T n'est pas defini.
     double speed_of_sound(double /*p*/, double /*T*/) const override {
-        return 1.0 / params_.beta;
+        return std::numeric_limits<double>::infinity();
     }
 
     double temperature_from_enthalpy(double /*p*/, double h) const override {
@@ -116,11 +128,11 @@
     }
 
     double dp_drho_s(double /*p*/, double /*T*/) const override {
-        return 1.0 / (params_.beta * params_.rho);
+        return std::numeric_limits<double>::infinity();
     }
 
     double dp_dT_rho(double /*p*/, double /*T*/) const override {
-        return params_.rho * params_.beta;
+        throw std::domain_error("IncompressibleEOS: (dp/dT)_rho undefined at constant density");
     }
 
     double drho_dp_T(double /*p*/, double /*T*/) const override {
@@ -139,7 +151,7 @@
         return internal_energy(0, T) + 0.5 * u_mag2;
     }
 
-    void set_params(const IncompressibleParams& params) { params_ = params; }
+    void set_params(const IncompressibleParams& params) { params_ = params; validate_params(); }
     const IncompressibleParams& params() const { return params_; }
 };
 
```

**Test** (`B14_*`, behaviour).
- **Montage :** ρ = 1000, β = 2e-4.
- **Références :**
  - c = ∞ pour les deux copies, alors que l'original donnait √(ρ/β) ;
  - dp/dρ|s = ∞, alors que l'original donnait 1/β ;
  - dp/dT|ρ lève `std::domain_error`.

#### B12 🟧 Boussinesq non couplé — CONFIRMÉ, fonction **compilé + testé**, câblage solveur **compilé**

**Code actuel.**
- `boussinesq.h` fournit une accélération scalaire qui n'est lue nulle part.
- `steady_incompressible_solver.h:366-372` n'ajoute que `density·body_force` : aucune flottabilité.

**Code proposé.** Ajouter une force vectorielle f = −ρ_ref·β(T − T_ref)·**g** et la brancher sur `IncompressibleSolverControls::temperature`. Le pointeur n'est pas possédé et vaut `nullptr` par défaut, ce qui préserve le comportement actuel.
`src/cfdx/physics/boussinesq.h`
```diff
@@ -1,8 +1,16 @@
 #pragma once
+#include "cfdx/core/field/field.h"
 #include <cmath>
 #include <stdexcept>
 namespace cfdx::physics {
 struct BoussinesqModel{double rho_ref=1.0,beta=0.0,T_ref=300.0,gravity=9.81;
 double density(double T)const{if(rho_ref<=0.0||beta<0.0)throw std::invalid_argument("invalid Boussinesq parameters");return rho_ref*(1.0-beta*(T-T_ref));}
 double buoyancy_acceleration(double T)const{return gravity*beta*(T-T_ref);}};
+// Force volumique de flottabilite (N/m3) a ajouter a rho_ref*g :
+//   f = (rho(T) - rho_ref) g = -rho_ref beta (T - T_ref) g
+// g est le VECTEUR gravite (ex. {0,0,-9.81}) : f_z > 0 pour T > T_ref.
+inline cfdx::core::Vec3 boussinesq_body_force(const BoussinesqModel& m,double T,const cfdx::core::Vec3& g){
+    if(m.rho_ref<=0.0||m.beta<0.0||!std::isfinite(T))throw std::invalid_argument("invalid Boussinesq parameters");
+    const double s=-m.rho_ref*m.beta*(T-m.T_ref);
+    return cfdx::core::Vec3{s*g.x,s*g.y,s*g.z};}
 }
```
`src/cfdx/physics/steady_incompressible_solver.h`
```diff
@@ -7,6 +7,7 @@
 #include "cfdx/core/linalg/vector.h"
 #include "cfdx/core/numerics/gradient.h"
 #include "cfdx/io/restart/dat_restart.h"
+#include "cfdx/physics/boussinesq.h"
 #include "cfdx/physics/finite_volume_transport.h"
 #include "cfdx/physics/pressure_velocity_algorithms.h"
 #include "cfdx/physics/solver_control.h"
@@ -62,6 +63,11 @@
     double kinematic_viscosity = 1.0e-3;
     double turbulent_viscosity = 0.0;
     cfdx::core::Vec3 body_force{0.0, 0.0, 0.0};
+    // Flottabilite de Boussinesq, active si temperature != nullptr. Le champ doit
+    // survivre au solve (pointeur non possede). gravity est le vecteur g.
+    const cfdx::core::Field<double, cfdx::core::Location::CELL>* temperature = nullptr;
+    BoussinesqModel boussinesq;
+    cfdx::core::Vec3 gravity{0.0, 0.0, -9.81};
     std::size_t pressure_reference_cell = 0;
     double pressure_reference_value = 0.0;
     bool use_bounded_convection = true;
@@ -367,10 +373,19 @@
         Field<double, Location::CELL> body_x(mesh.n_cells(), "body_x", "N/m3", 1);
         Field<double, Location::CELL> body_y(mesh.n_cells(), "body_y", "N/m3", 1);
         Field<double, Location::CELL> body_z(mesh.n_cells(), "body_z", "N/m3", 1);
+        if (controls.temperature != nullptr && controls.temperature->size() != mesh.n_cells())
+            throw std::invalid_argument("Boussinesq temperature field size must match mesh cells");
         for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
             body_x(c) = controls.density * controls.body_force.x;
             body_y(c) = controls.density * controls.body_force.y;
             body_z(c) = controls.density * controls.body_force.z;
+            if (controls.temperature != nullptr) {
+                const auto fb = boussinesq_body_force(
+                    controls.boussinesq, (*controls.temperature)(c), controls.gravity);
+                body_x(c) += fb.x;
+                body_y(c) += fb.y;
+                body_z(c) += fb.z;
+            }
         }
 
         ScalarBoundaryConditions ubc_x, ubc_y, ubc_z;
```

**Test** (`B12_boussinesq_body_force`, new_api).
- **Référence :** ρ = 1,2, β = 3e-3, ΔT = 10 et g = (0, 0, −9,81) donnent f_z = +0,35316 N/m³, à 1e-12.
- **Résultat :** passe.
- **Limite :** le câblage dans le solveur compile, et les suites « steady incompressible » restent vertes avec `temperature = nullptr`. Il n'est **pas** exercé par un cas de convection naturelle.
- **À ajouter en CI :** une cavité différentiellement chauffée de de Vahl Davis, Ra = 10³, dont la référence est Nu = 1,118.

#### B11 🟧 CHT — CONFIRMÉ, **compilé + testé**

**Code actuel.** `cht_solver.h:125-131` (échange) et `cht_solver.h:152-158` (contrôle de flux) :
- `d = |x_f − x_P|` est une distance **euclidienne** ; elle doit être **normale**.
- Il y a deux sources de vérité pour k, `controls.conductivity1/2` et `energy1/2.conductivity`, qui peuvent diverger en silence.

Le même `.mag()` apparaît dans deux autres endroits :
- `finite_volume_transport.h:264`, pour toutes les conditions de bord de valeur fixée ;
- `energy_solver.h:123`, dans le bilan.

**Code proposé.**
`src/cfdx/physics/finite_volume_transport.h`
```diff
@@ -133,6 +134,20 @@
 // with first-order upwind convection and two-point orthogonal diffusion.
 // For steady problems, bounded convection subtracts div(phi)*psi, matching
 // the standard bounded finite-volume treatment.
+// Distance normale centre -> face de bord, |(x_f - x_P).n|. C'est la distance
+// que suppose le flux diffusif Gamma*A*(phi_f - phi_P)/d. La norme euclidienne
+// |x_f - x_P| la surestime des qu'un centre est decale tangentiellement.
+inline double boundary_normal_distance(
+    const FvGeometry& geometry, std::size_t face, std::size_t cell)
+{
+    const auto& Sf = geometry.face_area_vectors[face];
+    const double area = Sf.mag();
+    if (!(area > 0.0) || !std::isfinite(area))
+        throw std::runtime_error("boundary_normal_distance: degenerate face");
+    const auto d = geometry.face_centres[face] - geometry.cell_centres[cell];
+    return std::abs(d.dot(Sf)) / area;
+}
+
 inline ScalarEquation assemble_scalar_equation(
     const cfdx::core::Mesh& mesh,
     const FvGeometry& geometry,
@@ -261,7 +276,7 @@
             }
 
             const double area = geometry.face_area_vectors[f].mag();
-            const double distance = (geometry.face_centres[f] - geometry.cell_centres[o]).mag();
+            const double distance = boundary_normal_distance(geometry, f, o);
             min_face_area = std::min(min_face_area, area);
             max_face_area = std::max(max_face_area, area);
             if (!(distance > 0.0) || !(area > 0.0))
```
`src/cfdx/physics/cht_solver.h`
```diff
@@ -95,6 +95,20 @@
     double interface_temperature_change=0.0;
 };
 
+// Temperature d'interface a flux normal continu entre deux cellules :
+//   k1 (T1 - Ti)/d1 = k2 (Ti - T2)/d2  =>  Ti = (h1 T1 + h2 T2)/(h1 + h2), h = k/d.
+// d1, d2 sont les distances NORMALES centre-face (boundary_normal_distance).
+inline double cht_interface_temperature(double k1, double d1, double T1,
+                                        double k2, double d2, double T2)
+{
+    if(!(k1>0.0) || !(k2>0.0) || !(d1>0.0) || !(d2>0.0) ||
+       !std::isfinite(k1) || !std::isfinite(k2) || !std::isfinite(d1) ||
+       !std::isfinite(d2) || !std::isfinite(T1) || !std::isfinite(T2))
+        throw std::invalid_argument("cht_interface_temperature: invalid inputs");
+    const double h1=k1/d1, h2=k2/d2;
+    return (h1*T1+h2*T2)/(h1+h2);
+}
+
 inline ChtSolveResult solve_two_region_cht(
     const cfdx::core::Mesh& mesh1, const FvGeometry& g1,
     const cfdx::core::Mesh& mesh2, const FvGeometry& g2,
@@ -110,6 +124,14 @@
     const ScalarBoundaryConditions& bcs1 = {},
     const ScalarBoundaryConditions& bcs2 = {})
 {
+    // Une seule source de verite pour k : celle des equations d'energie. Les
+    // champs conductivity1/2 sont conserves pour compatibilite mais doivent
+    // concorder, sinon l'interface et le volume utiliseraient deux k differents.
+    if(controls.conductivity1!=energy1.conductivity ||
+       controls.conductivity2!=energy2.conductivity)
+        throw std::invalid_argument(
+            "CHT: ChtInterfaceControls::conductivity1/2 must equal EnergySolverControls::conductivity");
+    const double k1=energy1.conductivity, k2=energy2.conductivity;
     const auto pairs=match_cht_interface(mesh1,g1,mesh2,g2,controls);
     ChtSolveResult result;
     std::vector<double> previous_interface_temperature(pairs.size(), std::numeric_limits<double>::quiet_NaN());
@@ -122,13 +144,11 @@
         // sides in the fixed-point limit.
         for(std::size_t i=0;i<pairs.size();++i) {
             const auto& p=pairs[i];
-            const double d1=(g1.face_centres[p.face1]-g1.cell_centres[p.cell1]).mag();
-            const double d2=(g2.face_centres[p.face2]-g2.cell_centres[p.cell2]).mag();
+            const double d1=boundary_normal_distance(g1,p.face1,p.cell1);
+            const double d2=boundary_normal_distance(g2,p.face2,p.cell2);
             if(!(d1>0.0) || !(d2>0.0) || !std::isfinite(d1) || !std::isfinite(d2))
                 throw std::runtime_error("CHT interface face-to-cell distance is invalid");
-            const double h1=controls.conductivity1/d1;
-            const double h2=controls.conductivity2/d2;
-            const double Tint_new=(h1*T1(p.cell1)+h2*T2(p.cell2))/(h1+h2);
+            const double Tint_new=cht_interface_temperature(k1,d1,T1(p.cell1),k2,d2,T2(p.cell2));
             double Tint=Tint_new;
             if(std::isfinite(previous_interface_temperature[i]))
                 Tint=controls.relaxation*Tint_new +
@@ -149,13 +169,11 @@
         double qimb=0.0, qscale=1e-30, dtint=0.0;
         for(std::size_t i=0;i<pairs.size();++i) {
             const auto& p=pairs[i];
-            const double d1=(g1.face_centres[p.face1]-g1.cell_centres[p.cell1]).mag();
-            const double d2=(g2.face_centres[p.face2]-g2.cell_centres[p.cell2]).mag();
+            const double d1=boundary_normal_distance(g1,p.face1,p.cell1);
+            const double d2=boundary_normal_distance(g2,p.face2,p.cell2);
             const double Tint=fv1.values.at(controls.region1_patch)[p.face1];
-            const double qflux1=controls.conductivity1*
-                (T1(p.cell1)-Tint)/d1;
-            const double qflux2=controls.conductivity2*
-                (Tint-T2(p.cell2))/d2;
+            const double qflux1=k1*(T1(p.cell1)-Tint)/d1;
+            const double qflux2=k2*(Tint-T2(p.cell2))/d2;
             qimb=std::max(qimb,std::abs(qflux1-qflux2));
             qscale=std::max(qscale,std::abs(qflux1));
             qscale=std::max(qscale,std::abs(qflux2));
```
`src/cfdx/physics/energy_solver.h`
```diff
@@ -120,7 +154,7 @@
         const double sign=(Sf0.dot(dvec)>=0.0)?1.0:-1.0;
         const auto Sf=Sf0*sign;
         const double area=Sf.mag();
-        const double d=dvec.mag();
+        const double d=boundary_normal_distance(geometry,f,o);
         const double F=mass_flux(f);
         double Tf=temperature(o);
         bool face_override_valid=false;
```

**Test** (`B11_cht_normal_distance`, new_api).
- **Montage :** face (0, 0, 0), normale x ; cellule 1 en (−0,5 ; 0,5 ; 0), donc décalée tangentiellement ; cellule 2 en (0,5 ; 0 ; 0) ; k = 1 ; T1 = 1, T2 = 0.
- **Référence :** d₁ = 0,5, donc T_i = **0,5**. L'ancien `.mag()` donnait 0,414.
- **Cas conductivités inégales :** T_i = (100·400 + 10·300)/110.
- **Résultat :** passe.

### 3.5 Rayonnement

#### B3 🟥 Couplage rayonnement–énergie : signe et linéarisation — CONFIRMÉ, **compilé + testé**

**Code actuel.** `radiation_solver.h:490-491` : `source = S_nr + qrad`. Or qrad = κ(4πI_b − G) est l'énergie **émise nette** par le milieu, donc un puits pour l'énergie. La correction de signe seule (`S_nr − qrad`) proposée par l'audit initial ne suffit pas : en explicite, une cellule adiabatique donne une diagonale nulle, et le solveur lève `non-positive diagonal`.

**Code proposé.** Linéarisation de Newton :
- Sp = −16κσT*³ ;
- Su = S_nr − qrad* − Sp·T*.

Sp est transmis par `source_implicit` (B2).
`src/cfdx/physics/radiation_solver.h`
```diff
@@ -487,12 +487,23 @@
         if(!rr.converged)
             throw std::runtime_error("radiation inner solve did not converge");
 
-        for(std::size_t c=0;c<nc;++c)
-            source(c)=non_radiative_source(c)+qrad(c);
+        // qrad = kappa (4 pi Ib - G) est l'energie EMISE nette par le milieu :
+        // c'est un puits pour l'equation d'energie. Linearisation de Newton autour
+        // de T* : S(T) = -qrad* + Sp (T - T*), Sp = d(-kappa 4 sigma T^4)/dT
+        // = -16 kappa sigma T*^3 <= 0 (diagonale renforcee, pas de relaxation
+        // explicite necessaire).
+        cfdx::core::Field<double,cfdx::core::Location::CELL> sp(
+            nc,"radiation_sp","W/m3/K",1);
+        for(std::size_t c=0;c<nc;++c) {
+            const double Tstar=temperature(c);
+            sp(c)=-16.0*controls.radiation.absorption*
+                  STEFAN_BOLTZMANN*Tstar*Tstar*Tstar;
+            source(c)=non_radiative_source(c)-qrad(c)-sp(c)*Tstar;
+        }
 
         auto er=solve_energy(
             mesh,geometry,mass_flux,temperature,source,
-            controls.energy,thermal_bcs);
+            controls.energy,thermal_bcs,nullptr,&sp);
         if(!er.converged)
             throw std::runtime_error("energy inner solve did not converge");
 
```

**Test** (`B3_radiative_equilibrium`, behaviour).
- **Montage :** cube dont toutes les parois sont à I_b(1000 K), κ = 1, T₀ = 800, sans autre source.
- **Référence :** T = 1000 K.
- **Résultats :** l'original lève « assemble_scalar_equation: non-positive diagonal at cell 0 » ; la copie corrigée donne T = 1000.

#### B18 🟧 (nouveau, ex-N1) Rosseland : k moléculaire et cp perdus, bilan incohérent — **compilé**

**Code actuel.** `radiation_advanced.h:586-592` et `radiation_advanced.h:615` :
- la conductivité passée vaut `k_R` seule, sans k_mol ;
- le flux convectif n'est pas multiplié par cp ;
- le bilan relit `energy_controls.conductivity`, et non la conductivité réellement utilisée.

**Code proposé.**
`src/cfdx/physics/radiation_advanced.h`
```diff
@@ -583,11 +611,13 @@
             }
         }
 
+        // k_eff = k_moleculaire + k_Rosseland ; flux convectif = m_dot cp.
         std::vector<double> conductivity_values(mesh.n_cells(),0.0);
         for(std::size_t c=0;c<mesh.n_cells();++c)
-            conductivity_values[c]=conductivity(c);
+            conductivity_values[c]=energy_controls.conductivity+conductivity(c);
+        const auto enthalpy_flux=enthalpy_face_flux(mass_flux,energy_controls.cp);
         auto eq=assemble_scalar_equation(
-            mesh,geometry,mass_flux,0.0,su,sp,bcs,true,nullptr,
+            mesh,geometry,enthalpy_flux,0.0,su,sp,bcs,true,nullptr,
             &transient_diag,&transient_rhs,&conductivity_values);
 
         cfdx::core::Vector candidate(mesh.n_cells(),0.0);
@@ -612,7 +642,8 @@
         const double rel=max_delta/scale;
         result.temperature_residuals.push_back(rel);
         result.energy_balance_residuals.push_back(energy_balance_relative(
-            mesh,geometry,mass_flux,temperature,old,source,energy_controls,bcs));
+            mesh,geometry,mass_flux,temperature,old,source,energy_controls,bcs,
+            nullptr,nullptr,&conductivity_values));
         result.iterations=iter;
         if(rel<=controls.tolerance &&
            result.energy_balance_residuals.back()<=controls.tolerance) {
```

**Test.** Aucun test analytique dédié. `test_phase12_radiation` et les suites Rosseland passent sur la copie corrigée.
**Test à ajouter :** plaque 1D d'épaisseur L, sans convection. On doit retrouver T linéaire et un flux q = (k_mol + 16σn²T̄³/(3β_R))·ΔT/L lorsque k_R est pris constant.

#### B5 🟧 Facteurs de forme par lancer de rayons : échantillonnage corrélé — CONFIRMÉ, **compilé + testé**

**Code actuel.** `radiation_advanced.h:475` : `r1 = (s+0.5)/N` reprend le même u que le choix du triangle source (`radiation_advanced.h:484`). Le point d'émission et l'angle polaire sont donc liés, ce qui biaise F dès que la source compte plusieurs triangles.

**Code proposé.** Une suite de Kronecker indépendante :
`src/cfdx/physics/radiation_advanced.h`
```diff
@@ -472,7 +497,10 @@
         tangent=radiation_scale(tangent,1.0/radiation_norm(tangent));
         auto bitangent=radiation_cross(n,tangent);
 
-        const double r1=(static_cast<double>(s)+0.5)/static_cast<double>(samples);
+        // r1 ne doit pas reprendre u = (s+0.5)/N : u choisit le triangle source,
+        // r1 l'angle polaire. Identiques, ils lient le point d'emission a la
+        // direction et biaisent F des que la source a plusieurs triangles.
+        const double r1=std::fmod(static_cast<double>(s)*0.8566748838545029+0.5,1.0);
         const double r2=std::fmod((static_cast<double>(s)*0.569840296)+0.25,1.0);
         const double phi=2.0*M_PI*r2;
         const double z=std::sqrt(1.0-r1);
```

**Test** (`B5_ray_traced_view_factor_offset_plates`, behaviour).
- **Montage :** carrés unitaires parallèles, distants de 0,5 et décalés de 1 selon x ; 65 536 rayons.
- **Référence :** F = 0,09367 (Hottel, aires croisées), à 2e-3.
- **Résultats :** l'original donne 0,1082 ou 0,0792 selon l'ordre des triangles (écart 0,0143) ; la copie corrigée donne 0,0938 dans les deux ordres.

#### B16 🟨 Autres constats rayonnement — CONFIRMÉ ; A1/A2 et réciprocité **compilé + testé**, DOM **non compilé**

**Rapport des surfaces A1/A2.** `m1_m4_models.h:231-240` n'a pas de rapport A1/A2. Précision par rapport à l'audit initial : `radiation.h::two_surface_net_exchange` a **déjà** une surcharge avec aires correcte ; seule la copie m1_m4 est fausse.
`src/cfdx/physics/m1_m4_models.h`
```diff
@@ -228,14 +229,18 @@
     return emissivity * blackbody(T);
 }
 
+// Flux net par unite de surface de 1 [W/m2] :
+//   q1 = sigma (T1^4 - T2^4) / [ (1-e1)/e1 + 1/F12 + (A1/A2)(1-e2)/e2 ]
+// A1_over_A2 = 1 : plaques paralleles infinies ; r1/r2 : cylindres coaxiaux.
 inline double two_surface_exchange(double e1, double e2, double T1, double T2,
-                                   double F12)
+                                   double F12, double A1_over_A2 = 1.0)
 {
     if (e1 <= 0.0 || e1 > 1.0 || e2 <= 0.0 || e2 > 1.0 ||
-        T1 < 0.0 || T2 < 0.0 || F12 < 0.0 || F12 > 1.0)
+        T1 < 0.0 || T2 < 0.0 || F12 < 0.0 || F12 > 1.0 ||
+        !(A1_over_A2 > 0.0) || !std::isfinite(A1_over_A2))
         throw std::invalid_argument("invalid radiation exchange input");
     if (F12 == 0.0) return 0.0;
-    const double resistance = (1.0-e1)/e1 + 1.0/F12 + (1.0-e2)/e2;
+    const double resistance = (1.0-e1)/e1 + 1.0/F12 + A1_over_A2*(1.0-e2)/e2;
     return sigma_sb * (std::pow(T1,4)-std::pow(T2,4)) / resistance;
 }
 
```

**Réciprocité des facteurs de forme.** `radiation_advanced.h:246-254` normalise les lignes, et seulement si leur somme dépasse 1, ce qui **casse** la réciprocité. La symétrisation proposée par l'audit initial, suivie d'une normalisation ligne par ligne, la casse aussi.

Corrigé : on symétrise G_ij = ½(A_iF_ij + A_jF_ji), puis on applique une mise à l'échelle **globale** si une ligne dépasse 1. Il ne faut pas forcer Σ_j F_ij = 1 : l'enceinte peut être ouverte.
`src/cfdx/physics/radiation_advanced.h`
```diff
@@ -231,6 +231,38 @@
                                         (M_PI * r2)));
 }
 
+// Reciprocite A_i F_ij = A_j F_ji imposee sur la moyenne G_ij, puis mise a
+// l'echelle globale si une ligne depasse 1 (la mise a l'echelle ligne par ligne
+// casse la reciprocite).
+inline std::vector<double> symmetrize_view_factors(
+    const std::vector<double>& F, const std::vector<double>& areas)
+{
+    const std::size_t n=areas.size();
+    if (n==0 || F.size()!=n*n)
+        throw std::invalid_argument("symmetrize_view_factors: dimension mismatch");
+    for (double a : areas)
+        if (!(a>0.0) || !std::isfinite(a))
+            throw std::invalid_argument("symmetrize_view_factors: areas must be positive");
+    std::vector<double> G(n*n,0.0);
+    for (std::size_t i=0;i<n;++i)
+        for (std::size_t j=0;j<n;++j) {
+            if (!std::isfinite(F[i*n+j]) || F[i*n+j]<0.0)
+                throw std::invalid_argument("symmetrize_view_factors: invalid entry");
+            G[i*n+j]=0.5*(areas[i]*F[i*n+j]+areas[j]*F[j*n+i]);
+        }
+    double worst=1.0;
+    for (std::size_t i=0;i<n;++i) {
+        double row=0.0;
+        for (std::size_t j=0;j<n;++j) row+=G[i*n+j]/areas[i];
+        worst=std::max(worst,row);
+    }
+    std::vector<double> out(n*n,0.0);
+    for (std::size_t i=0;i<n;++i)
+        for (std::size_t j=0;j<n;++j)
+            out[i*n+j]=G[i*n+j]/(worst*areas[i]);
+    return out;
+}
+
 inline std::vector<double> estimate_view_factor_matrix(
     const std::vector<ViewFactorPatch>& patches)
 {
@@ -243,16 +275,9 @@
             if (i != j)
                 F[i*n+j] = patch_pair_view_factor(patches[i], patches[j]);
 
-    // Enforce enclosure closure while preserving the physically required
-    // reciprocity as far as the centroid approximation permits.
-    for (std::size_t i=0; i<n; ++i) {
-        double sum = 0.0;
-        for (std::size_t j=0; j<n; ++j) sum += F[i*n+j];
-        if (sum > 1.0) {
-            for (std::size_t j=0; j<n; ++j) F[i*n+j] /= sum;
-        }
-    }
-    return F;
+    std::vector<double> areas(n);
+    for (std::size_t i=0; i<n; ++i) areas[i]=patches[i].area;
+    return symmetrize_view_factors(F, areas);
 }
 
 // -----------------------------------------------------------------------------
```

**Test** (`B16_view_factor_reciprocity`, new_api).
- `symmetrize_view_factors({0 ; 0,3 ; 0,2 ; 0}, {1 ; 2})` donne F12 = 0,35 et F21 = 0,175, soit A1F12 = A2F21.
- Cylindres gris coaxiaux, avec A1/A2 = 0,5, ε = 0,8 et 0,6 : on retrouve la formule de référence à 1e-3.

**Condition de paroi DOM** (`radiation_solver.h:377-379`). Une seule `wall_intensity_bcs` est appliquée à **toutes** les directions. Le bon opérateur, `apply_diffuse_gray_wall` (`radiation_advanced.h:161`), existe mais n'est appelé que par un test unitaire.

⚠️ Il faut aussi vérifier sa convention de normale. `hemispherical_irradiation` compte comme incidentes les directions telles que `d·n < 0`, alors que le paramètre s'appelle `outward_normal`. Avec une vraie normale sortante, ce sont les directions `d·n > 0` qui arrivent sur la paroi. Le test `test_m1_m4_physics.cpp:89` passe `{1,0,0}` sans que le sens soit fixé.

**Proposition** (**non compilée**). Construire, direction par direction, des `ScalarBoundaryFaceValues` à partir de l'itéré précédent :
```cpp
// Dans solve_participating_radiation (radiation_solver.h), avant l'assemblage de la direction m.
// Hypothèse : chaque patch de paroi est déclaré FIXED_VALUE dans wall_intensity_bcs, pour que la
// valeur par face soit prise (finite_volume_transport.h:247-256).
// wall_emissivity et wall_temperature sont deux nouveaux std::map<std::string,double>
// de RadiationTransportControls.
ScalarBoundaryFaceValues fv;
const auto& own = mesh.ownership();
for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
    if (!own.is_boundary(f)) continue;
    const auto& name = mesh.boundary().patch(geometry.face_patch[f]).name;
    const auto eps = controls.wall_emissivity.find(name);
    if (eps == controls.wall_emissivity.end()) continue;
    const std::size_t o = own.owner(f);
    const auto& Sf = geometry.face_area_vectors[f];
    const double A = Sf.mag();
    const auto n_out = (Sf.dot(geometry.face_centres[f] - geometry.cell_centres[o]) >= 0.0) ? Sf / A : Sf / (-A);
    double G = 0.0;                                        // éclairement incident sur la paroi
    for (std::size_t k = 0; k < directions.size(); ++k) {
        const double mu = directions[k].dx*n_out.x + directions[k].dy*n_out.y + directions[k].dz*n_out.z;
        if (mu > 0.0) G += directions[k].weight * mu * intensities[k](o);
    }
    const double mu_m = directions[m].dx*n_out.x + directions[m].dy*n_out.y + directions[m].dz*n_out.z;
    auto& v = fv.values[name];
    if (v.size() != mesh.n_faces()) v.resize(mesh.n_faces(), 0.0);
    v[f] = (mu_m < 0.0)                                    // direction quittant la paroi
        ? gray_diffuse_wall_intensity(eps->second, controls.wall_temperature.at(name), G)
        : intensities[m](o);                               // direction sortante : ignorée par l'upwind
}
auto eq = assemble_scalar_equation(mesh, geometry, directional_flux, 0.0, source, sp,
                                   wall_intensity_bcs, true, &fv);
```
Les noms `own.is_boundary` et `own.owner` sont à confirmer sur `core/mesh` : c'est pour cela que ce code n'est pas compilé.

**Test à ajouter :** deux plaques grises parallèles infinies, simulées par une tranche 1D en milieu transparent. Le flux doit être q = σ(T1⁴ − T2⁴)/(1/ε1 + 1/ε2 − 1) et ne pas dépendre de la quadrature.

#### B13 🟥 P1 : facteur π manquant — CONFIRMÉ (était PLAUSIBLE), **compilé + testé**

**Code actuel.** Deux fonctions renvoient `4κ(I_b − J)` alors que J = G/(4π) :
- `m1_m4_models.h:279` (`p1_source`) ;
- `radiation.h:143` (`p1_radiative_source`).

Le puits radiatif est donc **sous-estimé d'un facteur π** : avec J = 0, S = 4κσT⁴/π au lieu de 4κσT⁴.

**Code proposé.**
`src/cfdx/physics/m1_m4_models.h`
```diff
@@ -275,8 +280,8 @@
 {
     if (absorption < 0.0 || mean_intensity < 0.0 || T < 0.0)
         throw std::invalid_argument("invalid P1 source input");
-    // P1 source convention: S_r = 4*kappa*(I_b - J).
-    return 4.0 * absorption * (p1_blackbody_intensity(T) - mean_intensity);
+    // S_r = kappa (4 pi I_b - G) = 4 pi kappa (I_b - J), J = G/(4 pi).
+    return 4.0 * M_PI * absorption * (p1_blackbody_intensity(T) - mean_intensity);
 }
 
 struct Direction {
```
`src/cfdx/physics/radiation.h`
```diff
@@ -140,7 +140,8 @@
     validate_finite(temperature, "temperature");
     if (absorption < 0.0 || mean_intensity < 0.0 || temperature < 0.0)
         throw std::invalid_argument("p1_radiative_source: invalid input");
-    return 4.0 * absorption *
+    // S_r = kappa (4 pi I_b - G) = 4 pi kappa (I_b - J), J = G/(4 pi) [W/m2/sr].
+    return 4.0 * M_PI * absorption *
            (blackbody_intensity(temperature) - mean_intensity);
 }
 
```

**Test** (`B13_p1_source_4pi`, behaviour).
- **Montage :** κ = 1, J = 0, T = 1000.
- **Référence :** 4σT⁴ = 226 814,98 W/m³.
- **Résultats :** l'original s'écarte de 0,68 en relatif (= 1 − 1/π) ; la copie corrigée passe.

#### B17 🟧 S2S dans la PR #398 — CONFIRMÉ (était PLAUSIBLE), **compilé + testé** (sur une copie de `radiation_s2s.h` de `pr/398`)

**Code actuel.** `radiation_s2s.h`, à la tête de `refs/remotes/pr/398`, a deux défauts :
- l. 95-96 et 140-143 : `Gext = ambient_irradiation` est ajouté **en entier** à chaque surface, même en enceinte fermée (ΣF = 1). Il faut le pondérer par F_open = 1 − Σ_j F_ij.
- l. 158-175 : `estimate_s2s_view_factors_monte_carlo` estime F_ij et F_ji indépendamment. La matrice obtenue n'est pas réciproque, et **`validate_s2s_view_factors` du même fichier la rejette**.

**Code proposé.** La sémantique de `external_irradiation[i]` est conservée : c'est un éclairement déjà rapporté à la surface, ajouté tel quel.
`src/cfdx/physics/radiation_s2s.h` (pr/398)
```diff
@@ -13,12 +13,27 @@
 
 namespace cfdx::physics {
 
+// ambient_irradiation : eclairement d'un environnement noir hors enceinte
+// (sigma T_amb^4). Il n'atteint la surface i que par sa fraction ouverte
+// F_open,i = 1 - sum_j F_ij. external_irradiation[i], s'il est fourni, est un
+// eclairement deja rapporte a la surface i et s'ajoute tel quel.
 struct S2SControls {
     double tolerance = 1e-10;
     std::size_t max_iterations = 1;
     double ambient_irradiation = 0.0;
 };
 
+inline double s2s_external_irradiation(
+    const std::vector<double>& F, std::size_t n, std::size_t i,
+    const S2SControls& controls, const std::vector<double>& external_irradiation)
+{
+    double row_sum=0.0;
+    for (std::size_t j=0;j<n;++j) row_sum+=F[i*n+j];
+    const double F_open=std::max(0.0,1.0-row_sum);
+    return F_open*controls.ambient_irradiation +
+        (external_irradiation.empty() ? 0.0 : external_irradiation[i]);
+}
+
 struct S2SResult {
     bool converged = false;
     std::vector<double> radiosity;
@@ -92,8 +107,8 @@
             !std::isfinite(temperatures[i]) || temperatures[i]<=0.0)
             throw std::invalid_argument("invalid S2S surface state");
         const double E=STEFAN_BOLTZMANN*std::pow(temperatures[i],4);
-        const double Gext=external_irradiation.empty()
-            ? controls.ambient_irradiation : external_irradiation[i];
+        const double Gext=s2s_external_irradiation(
+            view_factors,n,i,controls,external_irradiation);
         if (!std::isfinite(Gext) || Gext<0.0)
             throw std::invalid_argument("invalid S2S external irradiation");
         const double rho=1.0-emissivities[i];
@@ -137,8 +152,8 @@
     r.net_flux.resize(n,0.0);
     double absorbed=0.0, emitted=0.0;
     for (std::size_t i=0;i<n;++i) {
-        const double Gext=external_irradiation.empty()
-            ? controls.ambient_irradiation : external_irradiation[i];
+        const double Gext=s2s_external_irradiation(
+            view_factors,n,i,controls,external_irradiation);
         for (std::size_t j=0;j<n;++j) r.irradiation[i]+=view_factors[i*n+j]*x[j];
         r.irradiation[i]+=Gext;
         r.net_flux[i]=x[i]-r.irradiation[i];
@@ -171,7 +186,12 @@
                 surfaces[i], surfaces[j], blockers, controls.samples);
         }
     }
-    return F;
+    // L'estimation Monte-Carlo de F_ij et F_ji est independante : sans cette
+    // etape validate_s2s_view_factors rejette la matrice (reciprocite).
+    std::vector<double> areas(n,0.0);
+    for (std::size_t i=0;i<n;++i)
+        for (const auto& t : surfaces[i]) areas[i]+=radiation_triangle_area(t);
+    return symmetrize_view_factors(F,areas);
 }
 
 } // namespace cfdx::physics
```

**Test** (`test_s3b_b17.cpp`).
- **Cas 1 :** plaques noires 1×1 en regard, F = 0,4152, T = 1000 et 600 K, T_amb = 300 K.
  - **Référence :** G1 = F·σT2⁴ + (1 − F)·σT_amb⁴.
  - **Cas fermé (F = 1) :** G1 = σT2⁴, l'ambiance ne doit rien apporter.
  - **Résultats :** l'original s'écarte de 5,7 % ; la copie corrigée passe à 1e-12.
- **Cas 2 :** Monte-Carlo sur des plaques d'aires 1 et 2.
  - **Référence :** A1F12 = A2F21 à 1e-12, et `validate_s2s_view_factors` ne doit pas lever.
  - **Résultats :** l'original s'écarte de 3,7e-4 et la validation rejette la matrice ; la copie corrigée passe.

### 3.6 Schémas temporels et divers

#### B10 🟧 Schémas temporels — CONFIRMÉ ; scalaire et BDF2 variable **compilé + testé**, Newton champ **non compilé**

**Code actuel.**
- `physics/temporal.h:24-45` : `detail::fixed_point` fait une substitution de Picard, x ← c + βΔt·f(x). Elle ne converge que si |βΔt·f′| < 1, ce qui retire aux schémas implicites toute stabilité en raideur.
- `core/numerics/temporal.h:137-173` : la variante champ a le même défaut.
- BDF2 utilise les coefficients 3, −4, 1 à pas constant (`core/numerics/temporal.h:160-161`), alors que Δt est adaptatif.

**Code proposé (scalaire).** Newton avec dérivée par différences centrées, et `bdf2_step_variable` :
`src/cfdx/physics/temporal.h`
```diff
@@ -30,18 +30,28 @@
     int max_iterations = 100,
     double tolerance = 1e-12)
 {
+    // Newton sur R(x) = x - constant - rhs_scale f(x). La substitution
+    // x <- constant + rhs_scale f(x) ne converge que si |rhs_scale f'(x)| < 1,
+    // soit dt < 1/|lambda| : elle retirait aux schemas implicites leur
+    // stabilite en raideur. R'(x) = 1 - rhs_scale f'(x) (differences centrees).
+    (void)dt;
     double x = initial;
     for (int iteration = 0; iteration < max_iterations; ++iteration) {
-        const double x_new = constant + rhs_scale * f(x);
-        if (!std::isfinite(x_new)) {
-            throw std::runtime_error("temporal fixed-point iteration produced a non-finite value");
+        const double fx = f(x);
+        const double residual = x - constant - rhs_scale * fx;
+        const double h = 1e-7 * std::max(1.0, std::abs(x));
+        const double dfdx = (f(x + h) - f(x - h)) / (2.0 * h);
+        const double jacobian = 1.0 - rhs_scale * dfdx;
+        if (!std::isfinite(residual) || !std::isfinite(jacobian) || std::abs(jacobian) < 1e-300) {
+            throw std::runtime_error("temporal Newton iteration produced a non-finite or singular step");
         }
+        const double x_new = x - residual / jacobian;
         if (std::abs(x_new - x) <= tolerance * std::max(1.0, std::abs(x_new))) {
             return x_new;
         }
         x = x_new;
     }
-    throw std::runtime_error("temporal fixed-point iteration did not converge");
+    throw std::runtime_error("temporal Newton iteration did not converge");
 }
 
 } // namespace detail
@@ -91,6 +101,26 @@
         constant);
 }
 
+// BDF2 a pas variable (omega = dt/dt_prev) :
+//   a0 phi^{n+1} = (1+omega) phi^n - omega^2/(1+omega) phi^{n-1} + dt f(phi^{n+1}),
+//   a0 = (1+2 omega)/(1+omega). omega = 1 redonne bdf2_step.
+template<class Function>
+inline double bdf2_step_variable(
+    double phi_n,
+    double phi_prev,
+    double dt,
+    double dt_prev,
+    Function&& f)
+{
+    detail::validate_dt(dt);
+    detail::validate_dt(dt_prev);
+    const double omega = dt / dt_prev;
+    const double a0 = (1.0 + 2.0 * omega) / (1.0 + omega);
+    const double constant =
+        ((1.0 + omega) * phi_n - omega * omega / (1.0 + omega) * phi_prev) / a0;
+    return detail::fixed_point(phi_n, std::forward<Function>(f), dt, dt / a0, constant);
+}
+
 template<class Function>
 inline double derivative(Function&& f, double phi, double h) {
     if (!(h > 0.0) || !std::isfinite(h)) {
```

**Code proposé (champ).** BDF2 à pas variable, avec `dt_prev` mémorisé dans `TimeIntegrationContext` :
`src/cfdx/core/numerics/temporal.h`
```diff
@@ -70,6 +70,7 @@
     Field<double, Location::CELL> phi_curr;  // φ^n
 
     bool has_prev = false;  // true after first step
+    double dt_prev = 0.0;   // pas de temps ayant produit phi_curr (0 = inconnu -> pas constant)
 
     TimeIntegrationContext() = default;
 
@@ -90,6 +91,11 @@
         phi_curr = phi_new;
         has_prev = true;
     }
+
+    void shift(const Field<double, Location::CELL>& phi_new, double dt) {
+        shift(phi_new);
+        dt_prev = dt;
+    }
 };
 
 // Advance field in time by one step
@@ -157,9 +163,13 @@
                         candidate = phi_data[cell] +
                             0.5 * dt * (rhs_old[cell] + rhs_data[cell]);
                     } else {
+                        // BDF2 a pas variable ; omega = 1 si dt_prev inconnu.
                         const double* prev = ctx->phi_prev.component_data(d);
-                        candidate = (4.0 * phi_data[cell] - prev[cell] +
-                            2.0 * dt * rhs_data[cell]) / 3.0;
+                        const double omega = ctx->dt_prev > 0.0 ? dt / ctx->dt_prev : 1.0;
+                        const double a0 = (1.0 + 2.0 * omega) / (1.0 + omega);
+                        candidate = ((1.0 + omega) * phi_data[cell] -
+                            omega * omega / (1.0 + omega) * prev[cell] +
+                            dt * rhs_data[cell]) / a0;
                     }
                     max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                     new_data[cell] = candidate;
@@ -244,7 +254,7 @@
     }
 
     if (ctx) {
-        ctx->shift(phi_new);
+        ctx->shift(phi_new, dt);
     }
     return phi_new;
 }
```

**Newton champ** (**non compilé**). Pour un second membre local à la cellule, un Jacobien diagonal suffit. Dans `core/numerics/temporal.h`, cela remplace la boucle de point fixe des schémas implicites :
```cpp
// R_c = phi_c - phi_c^n - beta*dt*rhs_c(phi)  ;  J_c ≈ 1 - beta*dt * d rhs_c / d phi_c
Field<double, Location::CELL> pert = phi_new, rhs_pert = rhs_new;
for (std::size_t c = 0; c < n; ++c) pert(c) += 1e-7 * std::max(1.0, std::abs(phi_new(c)));
rhs(pert, rhs_pert);
for (std::size_t c = 0; c < n; ++c) {
    const double h = pert(c) - phi_new(c);
    const double J = 1.0 - beta * dt * (rhs_pert(c) - rhs_new(c)) / h;
    const double R = phi_new(c) - phi_n(c) - beta * dt * rhs_new(c);
    if (!std::isfinite(J) || std::abs(J) < 1e-300) throw std::runtime_error("advance_time: singular Newton step");
    phi_new(c) -= R / J;
}
```
Ce Jacobien diagonal n'est **valide que si rhs_c ne dépend que de φ_c**. Pour un second membre diffusif, il faut passer par la matrice assemblée (`assemble_scalar_equation` avec `extra_diagonal` = ρV/Δt) : c'est le chemin qu'emprunte déjà `solve_energy`.

**Test.**
- `B10_implicit_euler_stiff` (behaviour) :
  - **Référence :** `implicit_euler_step(1 ; Δt = 1 ; f = −1e4·x)` = 1/(1 + 1e4).
  - **Résultats :** l'original lève « temporal fixed-point iteration produced a non-finite value » ; la copie corrigée passe.
- `B10_bdf2_variable_step` (new_api) :
  - **Montage :** φ = t², f = 2√φ, t = 1 → 1,5 → 1,8. BDF2 est exact pour un polynôme de degré 2, d'où la **référence φ = 3,24**.
  - **Vérifié dans trois configurations :**
    - `bdf2_step_variable` ;
    - `advance_time(BDF2)` avec `dt_prev = 0,5` ;
    - ω = 1, qui doit redonner `bdf2_step`.
  - **Résultat :** passe dans les trois cas.
- **Constat annexe :** `tests/numerical/test_temporal.cpp` **ne compile ni sur l'original ni sur la copie corrigée**, car l'ordre des arguments d'`advance_time` y est périmé. Il n'est pas enregistré dans CMake, ce qui explique qu'il ait pourri sans être vu. Seul `tests/unit/test_temporal.cpp` (`add_cfdx_test(test_temporal_physics …)`, `CMakeLists.txt:379`) tourne : 5/5 sur l'original comme sur la copie corrigée.

#### m1 🟨 Littéral `'\\n'` : pollution des logs par « 23662 » — CONFIRMÉ, **compilé + testé**

**Code actuel.** `finite_volume_transport.h:388,405` : `'\\n'` est une constante multi-caractère, affichée comme l'entier 23662 (0x5C6E). La trace est en outre inconditionnelle dès que N ≤ 256.

**Code proposé.** La trace est activée par `CFDX_SOLVER_TRACE`, et le `'\n'` est corrigé.
`src/cfdx/physics/finite_volume_transport.h`
```diff
@@ -20,6 +20,7 @@
 #include <vector>
 #include <limits>
 #include <iostream>
+#include <cstdlib>
 
 namespace cfdx::physics {
 
@@ -383,9 +398,10 @@
     auto result = cfdx::core::solve_bicgstab(
         equation.matrix, equation.rhs, candidate,
         controls.max_iterations, controls.tolerance);
-    if (solution.size() <= 256)
+    static const bool trace_cascade = std::getenv("CFDX_SOLVER_TRACE") != nullptr;
+    if (trace_cascade && solution.size() <= 256)
         std::cerr << "CFDX solver cascade: bicgstab status=" << static_cast<int>(result.status)
-                  << " iter=" << result.iterations << " residual=" << result.residual << '\\n';
+                  << " iter=" << result.iterations << " residual=" << result.residual << '\n';
 
     // Keep all retries anchored to the same nonlinear iterate; the accepted
     // predictor is updated only after a solver reports convergence.
@@ -400,9 +416,9 @@
         result = cfdx::core::solve_gmres(
             equation.matrix, equation.rhs, candidate,
             64, controls.max_iterations, controls.tolerance);
-        if (solution.size() <= 256)
+        if (trace_cascade && solution.size() <= 256)
             std::cerr << "CFDX solver cascade: gmres status=" << static_cast<int>(result.status)
-                      << " iter=" << result.iterations << " residual=" << result.residual << '\\n';
+                      << " iter=" << result.iterations << " residual=" << result.residual << '\n';
     }
 
     if (result.status != cfdx::core::SolverStatus::CONVERGED) {
```

**Test.**
- Sur l'original, la sortie de `test_s3b_behaviour` contient « …residual=2.75e-14**23662**CFDX solver cascade… », et `-Wall` émet 2 avertissements `-Wmultichar`, aux lignes 388:90 et 405:94.
- Sur la copie corrigée : 0 avertissement, aucune trace par défaut.
- À ajouter dans `CMakeLists.txt` : `add_compile_options(-Werror=multichar)`.

#### 3.6.3 Oracles existants à mettre à jour (changement voulu, pas régression)

Trois suites du dépôt échouent sur la copie corrigée, et c'est attendu : leurs oracles encodent l'ancien comportement.

| Test | Ligne | Ancien oracle | Nouvel oracle |
|---|---|---|---|
| `tests/unit/test_phase10_turbulence_hardening.cpp` | 40 | ε sans κ | `std::pow(0.09,0.75)*std::pow(k,1.5)/(0.41*y), 1e-12` |
| idem | 42 | ω sans κ | `std::sqrt(k)/(std::pow(0.09,0.25)*0.41*y), 1e-12` |
| `tests/unit/test_m2_m4_solver_validation.cpp` | 44 | ε_w sans κ | `std::pow(0.09,0.75)/(0.41*0.1), 1e-12` |
| `tests/unit/test_phase10_model_registry.cpp` | 43 | `turbulence_nu_t(0.1, 0.02, 10.0, 0.01, c)` en LES | `turbulence_nu_t(0.1, 0.02, 10.0, 0.01, c, 1.0e-6)` : le volume est désormais obligatoire |

Les versions mises à jour se trouvent dans `scratch_b/s3b/tests_updated/`. Elles passent toutes sur la copie corrigée.

Suites vertes sur la copie corrigée :
- level_c_coupled, sst_limiter, sst_blending, cht_validation, core_temporal, temporal_physics ;
- benchmark_matrix (38/38) ;
- m1_m4_validation, m1_m4_physics ;
- level_b, scalar_diffusion, field_wall ;
- phase12_radiation, thermophysical, EOS consistency, transport_models.

#### 3.6.4 Reproduire

Le dépôt n'est jamais modifié : tout se passe dans `/tmp/cfdx_meta/scratch_b/s3b/`.

```bash
cd /tmp/cfdx_meta/scratch_b/s3b
rm -rf overlay/src && cp -r orig_src overlay/src          # orig_src = copie intacte de src/ @27bae126
cd patches && for p in p_fvt p_turb_transport p_turb_solver p_sst p_wall_rad p_energy p_coupled \
  p_rosseland p_vf p_thermo p_eos p_limits p_cht p_boussinesq p_temporal; do python3 $p.py; done; cd ..
build/build_one.sh orig tests/test_s3b_behaviour.cpp      # attendu : 0 passed, 10 failed
build/build_one.sh ovl  tests/test_s3b_behaviour.cpp      # attendu : 10 passed
build/build_one.sh ovl  tests/test_s3b_new_api.cpp        # attendu : 11 passed
build/run_suite.sh overlay/src ovl /tmp/cfdx_dev/tests/unit/*.cpp   # attendu : 3 échecs d'oracle, voir 3.6.3
# B17 : cp radiation_s2s.h de pr/398 dans b17/{orig,ovl}/src, puis python3 patches/p_b17.py b17/ovl/src
```
Le diff complet de la copie corrigée par rapport à l'original (20 fichiers, +548/−158) est dans `build/s3b.diff`, et celui de B17 dans `build/b17.diff`.

#### 3.6.5 Tableau de synthèse

| Item | Gravité | Statut | Test / référence |
|---|---|---|---|
| B1 Γ en source | 🟥 | compilé + testé | chaîne 1D, k ≡ 1 (orig : 2,7e129) |
| B4 SA fw, S̃, fv1 | 🟥 | compilé + testé | fw(1) = 1, plafond r = 10, S̃ ≥ 0 |
| B7 lois de paroi | 🟧 | compilé + testé | κ dans ε/ω ; y⁺_c = 11,53 |
| B8 k-ω 2006 | 🟧 | compilé + testé | noyau cellule, formules Wilcox |
| B9 SST | 🟧 | compilé + testé | limiteur 10β*kω, CD_kω, F2 |
| B15 contrat S | 🟨 | compilé + testé | cisaillement pur S = ∂u/∂y |
| B16 LES/DES Δ | 🟨 | compilé + testé | Δ = V^(1/3) |
| B2 cp convection | 🟥 | compilé + testé | Pe = 1 exact (orig ≈ diffusion pure) |
| B6 enthalpie | 🟧 | compilé + testé | intégrales exactes et cp borné |
| B14 EOS | 🟧 | compilé + testé | c = ∞, domain_error |
| B12 Boussinesq | 🟧 | fonction compilé + testé ; câblage compilé | f_z = 0,35316 |
| B11 CHT | 🟧 | compilé + testé | T_i = 0,5 (orig 0,414) |
| B3 couplage radiatif | 🟥 | compilé + testé | équilibre T = 1000 K (orig : diagonale nulle) |
| B18 Rosseland | 🟧 | compilé | suites existantes vertes |
| B5 lancer de rayons | 🟧 | compilé + testé | Hottel 0,09367 |
| B16 A1/A2, réciprocité | 🟨 | compilé + testé | A1F12 = A2F21 |
| B16 paroi DOM | 🟨 | non compilé | plaques grises, à ajouter |
| B13 P1 ×π | 🟥 | compilé + testé | 4σT⁴ |
| B17 S2S #398 | 🟧 | compilé + testé | F_open ; réciprocité MC |
| B10 Newton scalaire, BDF2 variable | 🟧 | compilé + testé | 1/(1 + 1e4) ; φ = 3,24 |
| B10 Newton champ | 🟧 | non compilé | — |
| B1-bis `ScalarExtraTerms` | — | non compilé | — |
| m1 `'\\n'` | 🟨 | compilé + testé | 0 -Wmultichar |

---

## 4. Solveurs linéaires, I/O, parallélisme

Référence : `/tmp/cfdx_dev@27bae126`. Les numéros de ligne « Code actuel » renvoient à ce commit.

Le « Code proposé » est extrait des fichiers patchés `scratch_c/psrc/cfdx/...`. Ils ont été compilés en scratch et n'ont **jamais** été écrits dans le dépôt.

Le patch complet (diff -ru, 2237 lignes) est dans `sections/s4_s8.patch`. Pour l'appliquer :

```
cd cfdx_dev/src && patch -p0 --dry-run < s4_s8.patch
```

Les chemins du patch sont absolus : adapter `-p`.

Statuts employés :
- **compilé+testé** : patch compilé, et le cas de reproduction bascule de l'ancien au nouveau comportement ;
- **compilé** : patch compilé, sans test dédié ;
- **non compilé** : code à relire avant intégration.

Numérotation : C1 à C20. **C14 n'existe pas** : c'est un trou de numérotation, pas un constat retiré. C19, C20 et les items N1 à N4 sont nouveaux.

Tests existants à mettre à jour si le patch est appliqué :
- `tests/test_bicgstab_solver.cpp:132-142` attend un refus sans diagonale : incompatible avec C12.
- `tests/test_vector.cpp:36` attend `out_of_range` en Release : incompatible avec la proposition §8-1.
- `test_ghia_cavity` échoue **avant et après** le patch. L'échec est pré-existant et hors périmètre. Sa durée passe de 3,60 s à 2,11 s.

### 4.1 Krylov

#### C1 🟥 NaN transformé en 0 dans la norme — CONFIRMÉ, compilé+testé
**Code actuel** : `core/linalg/krylov_reductions.h:86`
```cpp
return std::sqrt(std::max(0.0, krylov_sum(local, policy)));   // max(0, NaN) == 0
```
**Code proposé** (`krylov_norm2`) :
```cpp
const double sum = krylov_sum(local, policy);
// NaN/Inf must propagate: std::max(0.0, NaN) returns 0.0 and would turn a
// corrupted vector into an apparently converged one.
if (!std::isfinite(sum)) return std::numeric_limits<double>::quiet_NaN();
return std::sqrt(std::max(0.0, sum));
```
Le patch ajoute aussi une surcharge `krylov_dot(const double*, const double*, n, ...)`. Elle sert au MGS de GMRES (§8-3). Les trois solveurs testent ensuite `isfinite` sur β, h, ρ, pAp et rz.

**Résultat** : second membre contenant un NaN.
- GMRES : `CONVERGED` avant le patch, `DIVERGED` après.
- CG et BiCGStab : `DIVERGED` après le patch.

#### C2 🟥 Statut du CG pression ignoré — CONFIRMÉ, compilé
**Code actuel** : `physics/steady_incompressible_solver.h:578-585`. `solve_cg(...)` est appelé, puis `p_corr` est utilisé sans lire `status`.

**Code proposé** :
```cpp
const auto rp = solve_cg(A, b, p_corr, controls_.linear_max_iterations, controls_.linear_tolerance);
if (rp.status == SolverStatus::DIVERGED || rp.status == SolverStatus::NOT_APPLICABLE ||
    !std::isfinite(rp.residual_relative))
    throw std::runtime_error(std::string("pressure correction failed: status=") +
                             to_string(rp.status) + " it=" + std::to_string(rp.iterations) +
                             " rel=" + std::to_string(rp.residual_relative));
// MAX_ITER_REACHED toléré : une correction partielle reste utile en SIMPLE.
```
Le site d'appel fixe déjà `pressure_reference_cell` : la ligne est vidée, la diagonale vaut 1 et la colonne est effacée. Le cas Neumann pur n'arrive donc pas ici.

Défauts à rappeler : `linear_max_iterations = 1000` (l.59), `linear_tolerance = 1e-10` (l.60).

#### C3 🟧 CG : critère sur la norme préconditionnée — CONFIRMÉ, compilé+testé
**Code actuel** : `cg_solver.h:139,196-201`.
- `tol_abs = tolerance*max(b_norm,1e-15)` ;
- critère `sqrt(|r·z|)`, c'est-à-dire la norme préconditionnée ;
- `Ap`, `p_vector` et `Ap_vector` sont alloués à chaque itération.

**Code proposé** : nouvelle fonction `solve_pcg`. `solve_cg` garde sa signature et lui délègue.
```cpp
inline SolverResult solve_pcg(const SparseMatrix& A, const Vector& b, Vector& x,
                              const Preconditioner* M,            // nullptr => Jacobi (diag sommée)
                              std::size_t max_iter = 1000, double tolerance = 1e-12,
                              PrecisionPolicy precision = {}, KrylovReductionPolicy reduction = {});
...
Vector r(n), z(n), p(n), Ap(n);                 // alloués une fois
const double tol_abs = tolerance * b_norm;      // relatif ; b_norm == 0 => x = 0, CONVERGED
...
rnorm = krylov_norm2(r, redp, reduction);       // ||r||_2, pas sqrt(r·z)
if (!std::isfinite(rnorm)) return finish(SolverStatus::DIVERGED, iter);
if (rnorm <= tol_abs) {
    const SolverResult f = finish(SolverStatus::CONVERGED, iter);   // recalcule b - A x
    if (f.residual <= 10.0 * tol_abs) return f;                      // récurrence fiable
    result.status = SolverStatus::NOT_APPLICABLE; rnorm = f.residual; // sinon on continue
}
```
Les cas `pAp <= 0` et `rz <= 0` renvoient `NOT_APPLICABLE` (matrice non SPD).

**Résultat** (matrice mise à l'échelle 1e6) :

| | Résidu rapporté | Résidu vrai |
|---|---|---|
| Avant | 9,6e-9 | 1,36e-5 |
| Après | — | 2,1e-13 |

#### C4 / C5 🟧 GMRES : tolérance, restart, contrôles perdus — CONFIRMÉ, compilé+testé
**Code actuel** : `gmres_solver.h`
- l.53 : `tol = tolerance*max(b_norm,1.0)`. Le critère devient absolu dès que ‖b‖ < 1.
- l.41 : le restart est clampé à `[restart_min, restart_max]` = [10, 40] (`krylov_controls.h:10-11`).
- l.249 : la surcharge `SparseMatrix` passe `{}` au lieu des `controls` de l'appelant.
- `krylov_controls.h:32-35` (identique sur `master` et `pr/415`) : la politique adaptative est **inversée**. Une mauvaise réduction par cycle *rétrécit* l'espace de Krylov, jusqu'à 10. Correctif et preuve : §10.4.2 (COUPLED passe de `MAX_ITER` à 48 itérations).

**Code proposé** :
```cpp
controls.restart_min = std::max(1, std::min(controls.restart_min, restart));  // élargir, ne pas clamper
controls.restart_max = std::max(controls.restart_max, restart);
int current_restart = std::min(restart, n_int);
const double tol = tolerance * b_norm;          // relatif ; b_norm == 0 => x = 0
...
inline SolverResult solve_gmres(const SparseMatrix& A, const Vector& b, Vector& x, int restart = 30,
    std::size_t max_iter = 1000, double tolerance = 1e-12, Preconditioner* preconditioner = nullptr,
    PrecisionPolicy precision = {}, KrylovControls controls = {}) {   // C5
    ...
    return solve_gmres(op, b, x, restart, max_iter, tolerance, preconditioner, controls);
}
```
Arguments ajoutés à `SparseMatrix` : `precision` était déjà là. `controls` est nouveau, en dernière position, avec une valeur par défaut. L'ABI source est donc compatible.

**Résultats** :
- Pour ‖b‖ de 1e-4 à 1e-12 : avant, `CONVERGED` à l'itération 0 avec un résidu relatif vrai de 1 ; après, convergence réelle.
- Balayage du restart (326/451/719/1189 itérations) : identique à scipy.

#### C12 🟨 BiCGStab exige une diagonale sans préconditionneur — CONFIRMÉ, compilé+testé
**Code actuel** : `bicgstab_solver.h:46-62`. Toute ligne à diagonale nulle donne `NOT_APPLICABLE`, même quand `preconditioner == nullptr`.

**Code proposé** : bloc supprimé. Le refus revient au préconditionneur, qui l'exprime dans `setup()` :
```cpp
if(preconditioner && !preconditioner->setup(A)) return result;
```
Mettre à jour `test_bicgstab_solver.cpp:132-142`.

#### C13 🟨 GMRES : résidu vrai sur un x périmé — CONFIRMÉ, compilé
**Code actuel** : `gmres_solver.h`, bloc « residual replacement » dans la boucle interne. Il évalue `b - A x` avant la mise à jour de x du cycle. Le test est mort et coûte une matvec.

**Code proposé** : bloc supprimé. Le résidu vrai est calculé une seule fois par cycle, **après** `x += Z y` :
```cpp
beta = true_residual();
if (!std::isfinite(beta)) return finish(SolverStatus::DIVERGED, iterations, beta);
if (beta <= tol)          return finish(SolverStatus::CONVERGED, iterations, beta);
```
Ce β et ce `w.r` servent directement au départ du cycle suivant, sans matvec supplémentaire.

#### C17 🟨 Seuils absolus — CONFIRMÉ, compilé+testé
**Code actuel** :
- seuils 1e-30 dans BiCGStab (ρ, r̂·v) et dans la remontée GMRES ;
- pivot LU 1e-14 dans l'AMG.

**Code proposé** :
```cpp
// BiCGStab : breakdown relatif + redémarrage (au plus 20)
if (std::abs(rho_new) <= eps * rhat_norm * res) { if (++restarts > max_restarts) return finish(DIVERGED, ...); restart(); ... }
if (std::abs(rv) <= eps * rhat_norm * nrm(w.v)) { ... }
// GMRES : happy breakdown et remontée relatifs
breakdown = hnext <= 1e3 * eps * wnorm0;
if (std::abs(diag) <= eps * rmax) return finish(SolverStatus::DIVERGED, ...);
// AMG : pivot relatif (voir N-AMG)
const double pivot_tol = 1e3 * eps * amax * n;
```
**Résultat** : BiCGStab à ‖b‖ ≈ 1e-8 renvoie `DIVERGED` avant le patch et converge après.

#### C18 🟨 try/catch comme flux de contrôle — CONFIRMÉ, compilé+testé (np=2)
**Code actuel** : `distributed_poisson.h:74-87`. Un `try { local_index(gid) } catch (...)` est exécuté pour chaque voisin, à chaque matvec.

**Code proposé** (`distributed_execution.h`) :
```cpp
inline std::optional<std::size_t> find_local_index(const DistributedMeshView& v, std::size_t gid) noexcept;
// distributed_poisson.h
if (const auto li = find_local_index(view, gid)) { /* cellule locale */ } else { /* halo */ }
```
**Résultat** : matvec Poisson de 10,75 ms à 0,76 ms (×14).

#### N-DUP 🟧 Doublons CSR non fusionnés — nouveau, CONFIRMÉ, compilé+testé
**Code actuel** :
- `sparse_matrix.h:82` (`finalize`) garde les doublons (i, j) ;
- `preconditioner.h:32-45,61-70` (Jacobi) prend la **première** entrée diagonale ;
- `cg_solver.h` fait de même (boucle « premier trouvé + break »).

**Code proposé** :
- `finalize` trie de façon stable et fusionne les doublons (i, j) par somme, en conservant les zéros explicites ;
- Jacobi et `summed_diagonal(A)` somment toutes les entrées diagonales.

```cpp
for (std::size_t k = Ar[i]; k < Ar[i + 1]; ++k)
    if (Ac[k] == static_cast<std::uint32_t>(i)) d[i] += Av[k];
```
**Résultat** : sur une ligne {(i,i,−1), (i,i,+3)}, le CG renvoie `NOT_APPLICABLE` avant le patch (diagonale lue = −1) et `CONVERGED` après.

#### N-AMG 🟧 AMG : LU grossier refait à chaque application, échec en Neumann — nouveau, CONFIRMÉ, compilé+testé
**Code actuel** : `amg_preconditioner.h`, `smooth_coarsest`.
- Le LU dense est reconstruit **à chaque** `apply`.
- Le pivot absolu vaut 1e-14.
- Sur un Neumann pur, le noyau constant fait échouer le LU.

**Code proposé** : factorisation au `setup`. L'inconnue du noyau est fixée à 0 quand seul le dernier pivot est nul.
```cpp
if (!factor_coarsest()) return false;             // dans setup()
...
const double pivot_tol = 1e3 * eps * amax * static_cast<double>(n);
if (!(std::abs(lu_[p * n + k]) > pivot_tol)) {
    if (k + 1 == n) { singular_last_ = true; break; }   // noyau 1-D (Neumann)
    return false;                                       // rang < n-1 : non géré
}
// smooth_coarsest : descentes/remontées triangulaires seulement
if (singular_last_) x(n - 1) = 0.0;
```
Les tampons `ax`, `xc`, `rc` et `coarse_rhs_` deviennent des membres `mutable` par niveau.

⚠️ Conséquence : `apply` n'est plus réentrant ni thread-safe. Il faut une instance par thread.

**Résultat** (grille 256², `solve_pcg(A,b,y,&amg,20000,1e-8)`, `AgglomeratedAMGPreconditioner amg(op,0.7,2,2)`) :

| | Jacobi-CG | AMG-PCG |
|---|---|---|
| Dirichlet | 645 it, 0,375 s | 44 it, 0,166 s |
| Neumann | 952 it, 0,567 s | 66 it, 0,278 s |

Le branchement dans la pression est traité en §8-2 (non appliqué).

#### C6 🟧 Préconditionneur de la PR #415
Hors de cette section. Voir `sections/s10_pr415.md`.

#### C15 🟨 / C19 🟧 CUDA : critères d'arrêt et breakdown maquillé — code seul, **non compilé**
Pas de nvcc dans le scratch.

**Code actuel** :
- `cuda_poisson.cu:265` : critère d'arrêt sur un résidu absolu ;
- `cuda_poisson.cu:335` : un breakdown (pAp ≤ 0 ou NaN) sort de la boucle et est rapporté `MAX_ITER` (C19, nouveau).

**Code proposé** (aligner sur `solve_pcg`) :
```cpp
const double tol_abs = tolerance * b_norm;                    // relatif, b_norm calculé une fois
...
if (!isfinite(pAp)) { status = SolverStatus::DIVERGED; break; }
if (!(pAp > 0.0))   { status = SolverStatus::NOT_APPLICABLE; break; }   // C19 : plus de MAX_ITER
...
if (rnorm <= tol_abs) { /* recalcul b - A x sur GPU, garde 10×tol comme CPU */ }
```
Le pool de tampons GPU est traité en §8-7.

### 4.2 I/O

#### C7 🟥 Lecture de restart hors limites — CONFIRMÉ, compilé+testé
**Code actuel** : `io/restart/dat_restart.h:81-172`. Les tailles de U et p ne sont jamais vérifiées. U est écrit via `component_data` (≈l.131), p vers la ligne 139.

**Code proposé** (au début de `read_dat_restart_fields`) :
```cpp
if (U.dimension() != 3 || U.size() != mesh.n_cells())
    throw std::invalid_argument("read_dat_restart_fields: U must be 3-D with n_cells entries");
if (p.dimension() != 1 || p.size() != mesh.n_cells())
    throw std::invalid_argument("read_dat_restart_fields: p must be 1-D with n_cells entries");
...
if (!(in >> value) || !std::isfinite(value)) throw std::runtime_error("read_dat_restart_fields: invalid field value");
```
**Résultat** : sur des `Fields` vides, ASan détecte un SEGV avant le patch ; après, une exception est levée.

Reste à faire : lire dans des temporaires puis faire `swap`. Aujourd'hui, une exception levée en cours de lecture laisse un U partiellement écrasé.

#### C10 🟧 Écriture de restart non atomique, NaN écrits — CONFIRMÉ, compilé+testé
**Code actuel** : `dat_restart.h:36-80`. `std::ofstream out(path)` tronque le seul restart. Aucune vérification du flux, aucune vérification `isfinite`.

**Code proposé** :
```cpp
for (std::size_t c = 0; c < mesh.n_cells(); ++c)
    if (!std::isfinite(U.component_data(0)[c]) || !std::isfinite(U.component_data(1)[c]) ||
        !std::isfinite(U.component_data(2)[c]) || !std::isfinite(p(c)))
        throw std::invalid_argument("write_dat_restart_fields: non-finite value at cell " + std::to_string(c));
// Ecriture atomique : fichier temporaire puis rename(), l'ancien restart reste intact si on echoue.
const std::string tmp_path = path + ".tmp";
std::ofstream out(tmp_path, std::ios::trunc);
if (!out) throw std::runtime_error("write_dat_restart_fields: cannot open " + tmp_path);
/* ... écriture ... */
out.flush();
if (!out) { out.close(); std::remove(tmp_path.c_str());
    throw std::runtime_error("write_dat_restart_fields: write failed for " + tmp_path); }
out.close();
if (std::rename(tmp_path.c_str(), path.c_str()) != 0) { std::remove(tmp_path.c_str());
    throw std::runtime_error("write_dat_restart_fields: cannot rename " + tmp_path + " -> " + path); }
```
**Résultat** : avec un NaN dans U, l'ancien code écrit un fichier corrompu. Le nouveau lève une exception et **conserve** le restart précédent.

#### C8 🟧 HDF5 : chemin fixe `fields/values` — CONFIRMÉ, compilé+testé
**Code actuel** : `io/hdf5/hdf5_writer.cpp:366-403`. Tous les champs vont dans `fields/values` : le 2ᵉ écrase le 1ᵉʳ. Le lecteur ne vérifie pas le nom.

**Code proposé** (writer) :
```cpp
if (fname.empty() || fname.find('/') != std::string::npos || fname == "values")
    return false;  // le nom sert de chemin HDF5 : vide, '/' et "values" (layout legacy) interdits
const std::string path = "fields/" + fname;                    // layout v2 : un groupe par champ
if (H5Lexists(file, path.c_str(), H5P_DEFAULT) > 0 &&
    H5Ldelete(file, path.c_str(), H5P_DEFAULT) < 0) { H5Fclose(file); return false; }
hid_t grp = H5Gcreate2(file, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
bool ok = grp >= 0 && write_dataset(grp, "values", flat_values.data(), dims, 1) >= 0 && /* attrs */;
```
**Code proposé** (reader) :
- nouvelle surcharge `read_field_hdf5(filename, name, field)` ;
- elle lit d'abord `fields/<name>/values`, sinon l'ancien `fields/values`, en contrôlant l'attribut `name` ;
- l'ancienne signature lui délègue avec `field.name()`.

**Résultat** : écriture de p puis de T, relecture de T.
- Avant : on relit p (101325).
- Après : T = 300.

⚠️ Relire un fichier v2 dans un `Field` sans nom échoue : il faut passer `name`.

#### C9 🟧 `herr_t` ignorés — CONFIRMÉ, compilé
**Code actuel** : `hdf5_writer.cpp:288-362` (`write_mesh_hdf5`). Les retours de `write_dataset*` et `write_attr_str` sont ignorés, puis la fonction renvoie `true`.

**Code proposé** :
```cpp
if (write_dataset(file, "mesh/points", ...) < 0) { H5Fclose(file); return false; }
/* idem pour chaque write_dataset* / write_attr_str */
return H5Fclose(file) >= 0;
```
Une macro `CFDX_H5_CHECK` qui lèverait une exception est équivalente. Le retour `bool` a été retenu pour respecter l'API existante.

#### C16 🟨 Écrivain VTU — CONFIRMÉ (lecture du code), **non compilé**
**Code actuel** : `io/vtu/vtu_writer.cpp`
- **Flux non vérifié** : l.53 `return true;` sans test de `os`. Un disque plein donne un VTU tronqué signalé comme réussi.
- **Champs de face ignorés** : le paramètre `fields_face` (l.17) n'est **jamais lu** dans `write()`.
- **Tailles en 32 bits** :
  - `std::vector<std::uint32_t> cell_face_offsets` (l.33) ;
  - boucles `for (std::uint32_t i = begin; ...)` et `for (std::uint32_t v = fv_begin; ...)` (≈l.191-214).

  Au-delà de 2³² entrées face→sommet, la troncature est silencieuse, alors que les DataArray sont déclarés `UInt64`.
- **Orientation** : le flux `faces` écrit chaque face avec l'ordre de sommets du maillage, c'est-à-dire avec la normale sortante du propriétaire, **sans inversion pour la cellule voisine** (≈l.186-200). Pour le voisin, les faces partagées sont donc rentrantes.

  Mesure : `vtkCellSizeFilter` et `IsInside` / probe tolèrent cette orientation (volume 1,0 pour les deux cubes, T sondé correct). Le « volume 0,333 » de l'ancienne version de cet audit venait de mon maillage de test, dont la face 5 rentrait dans son propriétaire : il est **retiré**.

  Gravité ramenée à 🟨. Le risque subsiste pour les filtres qui exigent des faces sortantes.
- Tous les DataArray sont en ASCII (l.123-269) : voir §8-10.

**Code proposé** :
```cpp
// l.33 et boucles : 64 bits
std::vector<std::uint64_t> cell_face_offsets;
for (std::uint64_t i = begin; i < end; ++i) { ...
    for (std::uint64_t v = fv_begin; v < fv_end; ++v) ... }

// flux faces : inverser l'ordre pour la cellule voisine
const bool owner = mesh.faces().owner(face) == orig_c;          // accesseur à confirmer dans Mesh
if (owner) for (auto v = fv_begin; v < fv_end; ++v) os << fv[v] << ' ';
else       for (auto v = fv_end; v-- > fv_begin; )  os << fv[v] << ' ';

// champs de face : pas de support VTK direct -> avertir au lieu d'ignorer
if (!fields_face.empty())
    std::cerr << "VtuWriter: " << fields_face.size() << " face field(s) not exported (unsupported)\n";

// fin de write()
os.flush();
if (!os) { std::cerr << "VtuWriter: write failed for " << filename << '\n'; return false; }
return true;
```
Le nom de l'accesseur owner (`mesh.faces().owner`) n'a pas été vérifié : le tracer avant d'intégrer. L'issue #74 elle-même est résolue (POLYHEDRON natif depuis `fef59b99`).

#### C20 🟥 Import OpenFOAM : débordement de pile de `std::regex` — nouveau, CONFIRMÉ, compilé+testé
**Code actuel** : `io/openfoam/openfoam_importer.cpp`
- regex dans `strip_comments` (l.17-21), `read_points` (l.41-58) et `read_faces` (l.73-97) ;
- `read_label_list` (l.60-70) ignore le nombre déclaré par `read_declared_count` (l.33-39).

**Code proposé** : `strip_comments` linéaire, curseur à la main, contrôle du nombre déclaré et de la fin de liste.
```cpp
struct Cursor {
    const char* p; const char* end;
    void skip_ws() { while (p < end && std::isspace(static_cast<unsigned char>(*p))) ++p; }
    bool expect(char c) { skip_ws(); if (p < end && *p == c) { ++p; return true; } return false; }
    bool read_double(double& v) { skip_ws(); char* e; v = std::strtod(p, &e); if (e == p) return false; p = e; return true; }
    bool read_size(std::size_t& v) { skip_ws(); char* e; v = std::strtoull(p, &e, 10); if (e == p) return false; p = e; return true; }
};
// open_list : saute l'en-tête FoamFile (après le premier '}'), lit "N (", borne à la dernière ')'
// read_points / read_label_list / read_faces : lisent exactement N entrées puis exigent ')' final
if (n_read != declared) throw std::runtime_error("OpenFOAM: " + file + " declares " +
                                                 std::to_string(declared) + " entries, got " + std::to_string(n_read));
```
`read_boundary` garde sa regex : le fichier est petit. Les commentaires périmés des l.28 et 117 sont corrigés. La compilation avec `-Wall -Wextra -Wconversion -Wsign-conversion` ne donne aucun avertissement.

**Résultats** (banc `scratch_c/foam/`, `gen.py N` + `t_foam.cpp`, hash du maillage identique avant/après) :

| Cas | Avant | Après |
|---|---|---|
| Commentaire `/* */` de 200 kB dans `points` | **SEGV** rc=139 (pile de 8 Mo) | OK |
| Commentaire de 20 kB | OK | OK |
| bad1 : `owner` déclaré 3301 pour 3300 | accepté | **rejeté** |
| bad2 : jeton parasite | rejeté | rejeté |
| bad3 : faces sur une ligne + commentaire inline | OK | OK |

Performances : voir §8-9.

### 4.3 Parallélisme, OOC, GPU

#### C11 🟧 Budget OOC — CONFIRMÉ, compilé+testé
**Code actuel** : `runtime/ooc/ooc_executor.h:76-88,118-127`.
- Le budget compte un working set alors que le pipeline en garde deux en mémoire.
- Il lève une exception quand la charge **tient** en mémoire.
- Avec `staging=1`, le pool est épuisé.

**Code proposé** :
```cpp
enum class BudgetVerdict { UNCONSTRAINED, TILE_TOO_LARGE, FITS_IN_DEVICE, OOC_REQUIRED };
// 1 working set en séquentiel, 2 en pipeline (tuile i+1 chargée pendant le calcul de i)
BudgetVerdict assess_memory_budget(std::size_t field_components = 1,
                                   std::size_t resident_working_sets = 1) const;
```
- `for_each_tile` ne lève plus que sur `TILE_TOO_LARGE`.
- `for_each_tile_pipelined` repasse en séquentiel dans trois cas : moins de 2 tampons, une seule tuile, ou 2 working sets qui ne tiennent pas.
- `Pending::release` attend `load_future`, ce qui supprime la course du double buffer.
- `validate_memory_budget` garde son contrat (appelé par `runtime_pipeline.h:74`). `validate_memory_budget_if_configured` est supprimé.

**Résultats** :
- `staging=1` : l'ancien code lève « pool exhausted », le nouveau calcule 3/3 tuiles.
- Budget de 48 B : l'ancien pic résident vaut 2×32 B (hors budget), le nouveau garde une seule tuile résidente.
- Tests OOC du dépôt : OK.

Le pic n'est observable qu'avec 5 ms de `sleep` dans le calcul.

#### N1 🟧 Partition Morton : rangs vides et 10 bits — nouveau, CONFIRMÉ, compilé+testé
**Code actuel** : `mesh_partitioner.h:113-144`.
- Quantification sur 10 bits par axe.
- Découpe par `ceil(n/parts)`, qui laisse les derniers rangs vides.

**Code proposé** :
```cpp
const auto quant = [qmax](double t) { /* 21 bits : qmax = (1u << 21) - 1 */ };
const std::uint64_t morton = morton3d(quant(nx), quant(ny), quant(nz));   // sfc_ordering.h
std::sort(keys.begin(), keys.end(), [](auto& a, auto& b){ return a.key != b.key ? a.key < b.key : a.index < b.index; });
part.cell_rank[keys[i].index] = static_cast<int>((i * parts) / n_cells);  // découpe équilibrée
```
**Résultats** :

| Cas | Avant | Après |
|---|---|---|
| N=27, P=8 | 4,4,4,4,4,4,3,**0** | 4,3,4,3,3,4,3,3 |
| N=27, P=10 | dernier rang vide | 2 ou 3 cellules par rang |

Sur un maillage gradué 64³ (262 144 cellules), on obtient 176 128 clés distinctes à 10 bits contre 241 664 à 21 bits. `test_mpi_partition` à np=2 passe 11/11.

#### N2 🟥 Ordre des halos — nouveau, CONFIRMÉ (auparavant PLAUSIBLE), compilé+testé
**Code actuel** : `distributed_execution.h`.
- ≈l.159 : la liste d'envoi suit l'ordre local, alors que la réception suit l'ordre des gid.
- ≈l.191-250 : 3 `MPI_Sendrecv` pour les comptes et 1 pour les données, avec **tous** les rangs.

**Code proposé** :
- liste d'envoi triée par gid ;
- échange limité aux voisins, non bloquant :

```cpp
MPI_Irecv(recv_buffers[k].data(), (int)recv_buffers[k].size(), MPI_DOUBLE, peer, 612, comm, &requests.back());
MPI_Isend(send_buffers[k].data(), (int)send_buffers[k].size(), MPI_DOUBLE, peer, 612, comm, &requests.back());
...
MPI_Waitall((int)requests.size(), requests.data(), statuses.data());
int received_count = 0;
MPI_Get_count(&statuses[2 * k], MPI_DOUBLE, &received_count);
if ((std::size_t)received_count != recv_buffers[k].size())
    throw std::runtime_error("exchange_distributed_cell_halo: peer sent an unexpected halo size");
```
Un plan asymétrique (envoi vide mais réception non vide) lève une exception. ⚠️ Un message plus long que prévu produit `MPI_ERR_TRUNCATE`, pas l'exception ci-dessus.

**Résultats** (gids inversés) :
- avant : 16/16 valeurs de halo fausses par rang à np=2 ;
- après : 0 fausse à np=2 et à np=4 ;
- tests du dépôt à np=2 : OK.

#### N3 🟨 Réductions : Gather + Bcast — compilé+testé
**Code actuel** : `mpi_utils.h:101-113`. Somme déterministe par `MPI_Gather` vers le rang 0, puis `MPI_Bcast`.

**Code proposé** : `MPI_Allgather` des contributions, puis une somme dans l'ordre des rangs sur chaque rang. C'est toujours déterministe, en un seul collectif. `MPI_Allreduce` a été écarté : son ordre de sommation dépend de l'implémentation, ce qui casse la reproductibilité bit à bit.

#### N4 🟨 RCM en O(Nc·Nf) — nouveau, compilé+testé
Détail en §8-8.

#### Architecture (constats sans patch)
- La production utilise uniquement CG-Jacobi. L'AMG existe mais n'est pas branché : voir §8-2.
- Le chemin GPU n'est exercé que par les tests et n'est jamais compilé en CI.

---

## 5. Validation et tests

> Correctifs appliqués **uniquement** dans des copies de travail :
> - `/tmp/cfdx_scratch` : `master@27bae126` ;
> - `/tmp/cfdx_pr410` : tête de #410 ;
> - `/tmp/cfdx_pytest_scratch` : copie de #410 pour la partie Python.
>
> Aucun fichier suivi de `/tmp/cfdx_dev` n'a été modifié. Rien n'a été committé ni poussé, rien n'a été écrit sur GitHub.
>
> Statut donné pour chaque bloc :
> - **VALIDÉ** : compilé et exécuté en scratch ;
> - **NON COMPILÉ** : l'environnement nécessaire n'existe pas ici ;
> - **NON EXÉCUTÉ** : commande à lancer par le propriétaire.
>
> Builds utilisés :
> - `_b2` : `-DCFDX_ENABLE_MPI=OFF -DHDF5_ROOT=/tmp/hdf5_inst` ;
> - `_b3` : MPI ON, HDF5 série ;
> - `/tmp/cfdx_pr410/_b` : même configuration que `_b2`.

### 5.1 🟥 Le pilote de campagne rend 0 sur un `master` rouge (D8) — CONFIRMÉ

**Constat.** `scripts/run_validation.py` agrège tout en un seul `ctest -L phase13`.
- Un test de campagne **sans** LABEL `phase13` n'est donc jamais exécuté. C'est le cas de `test_phase9_acceptance`.
- Un test **absent** du registre n'est jamais signalé.
- Mesure : sur `master` (`_b2`), le script d'origine rend **0**, alors que `test_phase9_acceptance` échoue.

**Correction d'une affirmation antérieure.** La proposition précédente de ce §, basée sur `validation/cases.json`, `exe`, `metric` et `reference`, **n'existe pas dans le dépôt** : elle est retirée.

Même chose pour « le workflow VMFL est masqué par `continue-on-error` ». La vérification sur les fichiers a donné ceci :
- `cfdx-vmfl-validation.yml:57` met bien `continue-on-error` sur l'étape du rapport ;
- mais l'étape finale `Enforce validation gate` (`if: steps.validation_report.outcome == 'failure'` puis `exit 1`) refait échouer le job ;
- le 9/9 rouge est donc **visible**, pas masqué ;
- ce workflow ne se déclenche que sur `workflow_dispatch` ou sur le label `validation`.

**Remplacement complet de `scripts/run_validation.py`**, à mettre à la place du fichier existant :

Changements :
- résolution de chaque test dans le registre CTest (`ctest --show-only=json-v1`) ;
- exécution unitaire par nom exact (`^nom$`, `--no-tests=error`) ;
- statuts par test : `PASS`, `FAIL`, `MISSING`, `UNLABELLED_PASS`, `UNLABELLED_FAIL` ;
- statuts par campagne : `PASS`, `FAIL`, `MISSING_TEST`, `PASS_UNLABELLED`, `SKIPPED_BACKEND_OFF` (13.11 sans MPI), `BLOCKED_WITHOUT_HARDWARE` (13.12, 13.13) ;
- JUnit avec un `testcase` par test et un par campagne, et des `skipped` explicites.

```python
#!/usr/bin/env python3
"""Run CFDX numerical-model verification campaigns.

Phase 13 deliberately separates executable evidence from the closure decision:
software tests must pass; hardware-only evidence (CUDA/OOC) is reported as
BLOCKED when the required device is unavailable rather than being counted as
PASS.

Each campaign test is resolved against the CTest registry of the build
directory (``ctest --show-only=json-v1``) and executed individually, so that
a test that is missing, not labelled ``phase13`` or failing is attributed to
its campaign instead of being hidden behind a single aggregated exit code.
Use --strict to also fail while the documented closure matrix still contains
pending campaigns.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path

P13_LABEL = "phase13"

# Campaigns with executable evidence already present in the repository.
SOFTWARE_CAMPAIGNS: dict[str, list[str]] = {
    "13.1": ["test_mms_scalar_diffusion"],
    "13.2": ["test_mesh_refinement_order", "test_mms_scalar_diffusion"],
    "13.3": ["test_core_temporal", "test_temporal_physics"],
    "13.4": ["test_m1_m4_validation", "test_m2_m4_solver_validation", "test_cht_validation"],
    "13.5": ["test_matrix_free_fv_operator"],
    "13.6": ["test_analytical_benchmarks", "test_benchmark_matrix"],
    "13.7": ["test_steady_incompressible_solver", "test_phase9_acceptance",
             "test_poiseuille_body_force"],
    "13.8": ["test_cht_validation", "test_cht_two_slab_conduction",
             "test_level_c_coupled_verification"],
    "13.9": ["test_phase10_turbulence_hardening", "test_thermophysical_models_vv"],
    "13.10": ["test_level_c_coupled_verification", "test_analytical_benchmarks"],
    "13.11": ["test_phase8_mpi_poisson"],
}

# Campaigns that only exist when the build enables an optional backend.
# A missing test there is reported SKIPPED_BACKEND_OFF, never PASS.
OPTIONAL_BACKEND_CAMPAIGNS: dict[str, str] = {
    "13.11": "MPI (configure with -DCFDX_ENABLE_MPI=ON)",
}

# These require capabilities not guaranteed on a normal CPU runner.
HARDWARE_CAMPAIGNS: dict[str, str] = {
    "13.12": "real CUDA CPU/GPU equivalence for production M1-M4 solvers",
    "13.13": "real CUDA out-of-core CFD equivalence beyond device VRAM",
}

PENDING_CLOSURE: list[str] = [
    "13.1 full NS/transient/turbulence/thermal/radiation MMS",
    "13.2 production-scheme systematic order matrix",
    "13.3 transient PDE order for every production integrator",
    "13.4 independent integral conservation across all M1-M4 physics",
    "13.5 production assembled/matrix-free solution equivalence",
    "13.6 production conditioning/scale robustness",
    "13.7 pressure-velocity manufactured solution",
    "13.8 thermal/CHT MMS and interface convergence",
    "13.9 turbulence transport MMS for k-epsilon/SST/SA",
    "13.10 radiation transport and angular convergence campaign",
    "13.11 serial/MPI equivalence for production physics",
    "13.12 real CUDA production CPU/GPU equivalence",
    "13.13 real CUDA beyond-VRAM OOC equivalence",
]


@dataclass(frozen=True)
class TestResult:
    """Outcome of one CTest test."""

    name: str
    status: str  # PASS | FAIL | MISSING | UNLABELLED_PASS | UNLABELLED_FAIL
    returncode: int | None
    seconds: float
    output: str = ""


@dataclass
class CampaignResult:
    """Aggregated outcome of one Phase-13 campaign."""

    item: str
    status: str
    tests: list[str] = field(default_factory=list)
    detail: str = ""


def discover_tests(build_dir: Path) -> dict[str, list[str]]:
    """Return {test name: labels} from the CTest registry of ``build_dir``.

    Raises:
        RuntimeError: if ctest cannot list the tests or returns invalid JSON.
    """
    proc = subprocess.run(
        ["ctest", "--test-dir", str(build_dir), "--show-only=json-v1"],
        capture_output=True, text=True, check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"ctest --show-only failed ({proc.returncode}): {proc.stderr.strip()}")
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"ctest returned invalid JSON: {exc}") from exc

    registry: dict[str, list[str]] = {}
    for test in payload.get("tests", []):
        labels: list[str] = []
        for prop in test.get("properties", []):
            if prop.get("name") == "LABELS":
                labels = list(prop.get("value", []))
        registry[test["name"]] = labels
    return registry


def run_one(build_dir: Path, name: str, labels: list[str], timeout: int) -> TestResult:
    """Run a single registered test by exact name."""
    start = time.monotonic()
    try:
        proc = subprocess.run(
            ["ctest", "--test-dir", str(build_dir), "-R", f"^{name}$",
             "--output-on-failure", "--no-tests=error", "--timeout", str(timeout)],
            capture_output=True, text=True, check=False,
        )
        rc, out = proc.returncode, proc.stdout + proc.stderr
    except OSError as exc:
        rc, out = 127, f"cannot run ctest: {exc}"
    seconds = time.monotonic() - start
    ok = rc == 0
    if P13_LABEL in labels:
        status = "PASS" if ok else "FAIL"
    else:
        # The test runs, but `ctest -L phase13` in CI would never select it.
        status = "UNLABELLED_PASS" if ok else "UNLABELLED_FAIL"
    return TestResult(name, status, rc, seconds, out if not ok else "")


def evaluate_campaigns(
    registry: dict[str, list[str]], results: dict[str, TestResult]
) -> list[CampaignResult]:
    """Derive one status per campaign from the per-test results."""
    campaigns: list[CampaignResult] = []
    for item, tests in SOFTWARE_CAMPAIGNS.items():
        missing = [t for t in tests if t not in registry]
        if missing and item in OPTIONAL_BACKEND_CAMPAIGNS and len(missing) == len(tests):
            campaigns.append(CampaignResult(item, "SKIPPED_BACKEND_OFF", tests,
                                            OPTIONAL_BACKEND_CAMPAIGNS[item]))
            continue
        statuses = {t: results[t].status for t in tests if t in results}
        if missing:
            status, detail = "MISSING_TEST", "not registered: " + ", ".join(missing)
        elif any(s in ("FAIL", "UNLABELLED_FAIL") for s in statuses.values()):
            status = "FAIL"
            detail = "failed: " + ", ".join(t for t, s in statuses.items() if s.endswith("FAIL"))
        elif any(s == "UNLABELLED_PASS" for s in statuses.values()):
            status = "PASS_UNLABELLED"
            detail = "missing LABEL phase13: " + ", ".join(
                t for t, s in statuses.items() if s == "UNLABELLED_PASS")
        else:
            status, detail = "PASS", ""
        campaigns.append(CampaignResult(item, status, tests, detail))
    for item, description in HARDWARE_CAMPAIGNS.items():
        campaigns.append(CampaignResult(item, "BLOCKED_WITHOUT_HARDWARE", [], description))
    return campaigns


def write_junit(path: Path, results: list[TestResult], campaigns: list[CampaignResult]) -> None:
    """Write one <testcase> per test and per campaign."""
    suite = ET.Element("testsuite", {"name": "cfdx-phase13"})
    failures = skipped = 0
    for res in results:
        case = ET.SubElement(suite, "testcase", {
            "classname": "ctest", "name": res.name, "time": f"{res.seconds:.3f}"})
        if res.status == "MISSING":
            failures += 1
            ET.SubElement(case, "failure", {"message": "test not registered in CTest"})
        elif res.status.endswith("FAIL"):
            failures += 1
            node = ET.SubElement(case, "failure", {"message": f"exit code {res.returncode}"})
            node.text = res.output[-20000:]
    for camp in campaigns:
        case = ET.SubElement(suite, "testcase", {"classname": "campaign", "name": camp.item})
        if camp.status in ("FAIL", "MISSING_TEST", "PASS_UNLABELLED"):
            failures += 1
            ET.SubElement(case, "failure", {"message": f"{camp.status}: {camp.detail}"})
        elif camp.status in ("BLOCKED_WITHOUT_HARDWARE", "SKIPPED_BACKEND_OFF"):
            skipped += 1
            ET.SubElement(case, "skipped", {"message": camp.detail})
    suite.set("tests", str(len(results) + len(campaigns)))
    suite.set("failures", str(failures))
    suite.set("skipped", str(skipped))
    suite.set("errors", "0")
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.ElementTree(suite).write(path, encoding="utf-8", xml_declaration=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=Path("build"))
    parser.add_argument("--junit", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--timeout", type=int, default=1800, help="per-test timeout [s]")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="also fail while any documented Phase-13 campaign lacks closure evidence",
    )
    args = parser.parse_args()

    if not (args.build_dir / "CTestTestfile.cmake").exists():
        print(f"Build directory is not configured: {args.build_dir}", file=sys.stderr)
        return 2
    try:
        registry = discover_tests(args.build_dir)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    unique_tests = sorted({t for tests in SOFTWARE_CAMPAIGNS.values() for t in tests})
    results: dict[str, TestResult] = {}
    for name in unique_tests:
        if name not in registry:
            results[name] = TestResult(name, "MISSING", None, 0.0)
            continue
        print(f"$ ctest -R ^{name}$", flush=True)
        results[name] = run_one(args.build_dir, name, registry[name], args.timeout)
        print(f"  -> {results[name].status} ({results[name].seconds:.1f} s)", flush=True)

    campaigns = evaluate_campaigns(registry, results)

    print("\nPHASE 13 CAMPAIGN MATRIX")
    print("========================")
    for camp in campaigns:
        tail = f" -> {camp.detail}" if camp.detail else ""
        print(f"{camp.item:>6}: {camp.status}{tail}")

    # Optional-backend campaigns may legitimately have MISSING tests.
    optional_missing = {t for item in OPTIONAL_BACKEND_CAMPAIGNS
                        for t in SOFTWARE_CAMPAIGNS[item] if t not in registry}
    hard_fail = any(c.status in ("FAIL", "MISSING_TEST") for c in campaigns)
    unlabelled = any(c.status == "PASS_UNLABELLED" for c in campaigns)

    report = {
        "phase": 13,
        "tests": {r.name: {"status": r.status, "returncode": r.returncode,
                           "seconds": round(r.seconds, 3)}
                  for r in results.values() if r.name not in optional_missing},
        "campaigns": {c.item: {"status": c.status, "tests": c.tests, "detail": c.detail}
                      for c in campaigns},
        "pending_closure": PENDING_CLOSURE,
        "closure": "VALIDATION_IN_PROGRESS",
    }
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.junit:
        write_junit(args.junit,
                    [r for r in results.values() if r.name not in optional_missing], campaigns)

    print("\nCLOSURE: VALIDATION_IN_PROGRESS")
    print("The driver never promotes skipped/unavailable hardware into PASS.")
    if hard_fail:
        print("FAIL: at least one software campaign failed or references a missing test.")
        return 1
    if unlabelled:
        print("FAIL: campaign tests without LABEL phase13 are invisible to `ctest -L phase13`.")
        return 1
    if args.strict:
        print("STRICT: FAIL - Phase 13 still has pending closure campaigns.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

**Mesures** (scratch, après les correctifs des §5.2 à §5.4) :

| Build | Script d'origine | Script corrigé | Détail |
|---|---|---|---|
| `master` `_b2` | rc = 0 | **rc = 1** | 13.7 FAIL : `test_phase9_acceptance`, `test_poiseuille_body_force` ; 13.11 SKIPPED_BACKEND_OFF ; 13.12 et 13.13 BLOCKED. JUnit : tests = 30, failures = 4, skipped = 3 |
| #410 + correctifs | — | **rc = 0** (`--strict` : 1, clôture encore en attente) | 13.1 à 13.10 PASS |

**Statut** : VALIDÉ en scratch. Sorties : `/tmp/cfdx_meta/rv/{master_r.json,master_j.xml,pr410_r.json}`.

### 5.2 🟧 Tests tautologiques (D6) — CONFIRMÉ, 3 remplaçants VALIDÉS

| Test existant | Problème |
|---|---|
| `test_level_b_reference_benchmarks.cpp:13-23` | compare des constantes à elles-mêmes |
| `test_benchmark_matrix.cpp:65-67,77-78,108,111,131-132` | idem, aucun appel au solveur |
| `test_fluent_vmfl_reference.cpp:64,73` | aucun appel au solveur, affiche PASS |
| `test_phase10_sst_limiter.cpp:12` | teste une constante |
| `test_cht_validation` | cubes d'une seule cellule |
| `test_nonorthogonal_skew_campaign.cpp:43-60` | modifie `face_Sf` à la main, ce qui rend le cache géométrique incohérent ; seul critère : `correction > 0` |
| `test_analytical_benchmarks.cpp:274-297` | oracle auto-référencé |

Les vrais tests sont à préserver :
- MMS diffusion : ordre ≥ 1,90 (`test_mms_scalar_diffusion.cpp:357`) et ≥ 1,80 (`:381`) ;
- `test_mesh_refinement_order` ;
- `test_analytical_benchmarks.cpp:222-273`.

**Aide partagée `tests/common/channel_mesh.h`** (nouveau fichier). C'est le constructeur de canal extrait de `test_phase9_acceptance.cpp`, avec les en-têtes standard explicites :

```cpp
#pragma once
// Canal hexaedrique [0,1]^3 decoupe en nx x ny x 1, patches inlet/outlet/bottom/top/front/back.
// Extrait tel quel de tests/validation/test_phase9_acceptance.cpp (audit D6) pour que les
// tests solveur partagent un seul constructeur de maillage.

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::test {

using namespace cfdx::core;

inline Mesh make_channel_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 2 || ny < 4) throw std::invalid_argument("channel mesh too small");

    Mesh mesh;
    const std::size_t plane = (nx + 1) * (ny + 1);
    mesh.points().resize(2 * plane);
    const auto id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };

    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            const double x = static_cast<double>(i) / static_cast<double>(nx);
            const double y = static_cast<double>(j) / static_cast<double>(ny);
            mesh.points().set(id(i, j, 0), x, y, 0.0);
            mesh.points().set(id(i, j, 1), x, y, 1.0);
        }
    }

    std::map<std::vector<std::size_t>, std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx * ny);

    auto add_face = [&](std::initializer_list<std::size_t> vertices,
                        std::size_t cell) {
        std::vector<std::size_t> key(vertices);
        std::sort(key.begin(), key.end());
        const auto found = face_map.find(key);
        if (found != face_map.end()) {
            mesh.ownership().set_neighbour(found->second, static_cast<int>(cell));
            return found->second;
        }

        const std::size_t f = mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        mesh.ownership().set_owner(f, cell);
        mesh.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        face_map.emplace(std::move(key), f);
        return f;
    };

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t c = j * nx + i;
            const auto a = id(i, j, 0);
            const auto b = id(i + 1, j, 0);
            const auto c0 = id(i + 1, j + 1, 0);
            const auto d = id(i, j + 1, 0);
            const auto e = id(i, j, 1);
            const auto f = id(i + 1, j, 1);
            const auto g = id(i + 1, j + 1, 1);
            const auto h = id(i, j + 1, 1);
            cell_faces[c] = {
                add_face({a, d, c0, b}, c),
                add_face({e, f, g, h}, c),
                add_face({a, b, f, e}, c),
                add_face({d, h, g, c0}, c),
                add_face({a, e, h, d}, c),
                add_face({b, c0, g, f}, c)};
        }
    }

    for (const auto& faces : cell_faces) mesh.cells().push_cell(faces);

    Patch inlet{"inlet", PatchType::INLET, {}};
    Patch outlet{"outlet", PatchType::OUTLET, {}};
    Patch bottom{"bottom", PatchType::WALL, {}};
    Patch top{"top", PatchType::WALL, {}};
    Patch front{"front", PatchType::EMPTY, {}};
    Patch back{"back", PatchType::EMPTY, {}};

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f) >= 0) continue;

        const auto& vertices = mesh.faces().vertices();
        const auto begin = vertices.begin() +
            static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end = begin +
            static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));

        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        for (auto it = begin; it != end; ++it) {
            x += mesh.points().x(*it);
            y += mesh.points().y(*it);
            z += mesh.points().z(*it);
        }
        const double n = static_cast<double>(mesh.faces().face_size(f));
        x /= n; y /= n; z /= n;

        constexpr double tol = 1e-12;
        if (std::abs(x) < tol) inlet.face_ids.push_back(f);
        else if (std::abs(x - 1.0) < tol) outlet.face_ids.push_back(f);
        else if (std::abs(y) < tol) bottom.face_ids.push_back(f);
        else if (std::abs(y - 1.0) < tol) top.face_ids.push_back(f);
        else if (std::abs(z) < tol) front.face_ids.push_back(f);
        else if (std::abs(z - 1.0) < tol) back.face_ids.push_back(f);
        else throw std::runtime_error("channel boundary face is unclassified");
    }

    mesh.boundary().add_patch(inlet);
    mesh.boundary().add_patch(outlet);
    mesh.boundary().add_patch(bottom);
    mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(front);
    mesh.boundary().add_patch(back);
    return mesh;
}

} // namespace cfdx::test
```

**`tests/validation/test_poiseuille_body_force.cpp`** (nouveau). Le solveur SIMPLE produit `u`, et le test vérifie l'ordre observé sur 8, 16 et 32 mailles :

```cpp
// Poiseuille plan entraine par une force de volume : u(y) = f/(2 nu) y (1 - y).
// Remplace les controles formule-contre-formule de test_benchmark_matrix et de
// test_fluent_vmfl_reference (audit D6) : ici c'est le solveur qui produit u, et on
// verifie l'ordre spatial observe sur trois maillages.
#include "common/channel_mesh.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

constexpr double kNu = 0.1;
constexpr double kFx = 1.0;

double exact_u(double y) { return kFx / (2.0 * kNu) * y * (1.0 - y); }

struct PoiseuilleError {
    double l2 = 0.0;
    double linf = 0.0;
    double max_abs_uy = 0.0;
};

PoiseuilleError run_poiseuille(std::size_t ny)
{
    const std::size_t nx = 4;
    Mesh mesh = cfdx::test::make_channel_mesh(nx, ny);
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(0.0);

    VelocityBoundaryConditions ubc;
    ubc["inlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["bottom"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};
    ubc["top"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};
    ubc["front"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["back"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};

    ScalarBoundaryConditions pbc;
    for (const char* patch : {"inlet", "outlet", "bottom", "top", "front", "back"})
        pbc[patch] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.coupling.alpha_u = 0.7;
    c.coupling.alpha_p = 0.3;
    c.convergence.max_iterations = 3000;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-8;
    c.linear_max_iterations = 2000;
    c.linear_tolerance = 1e-9;   // GMRES(64) plafonne a 1,7e-10 sur le Poisson pur Neumann ny=32
    c.density = 1.0;
    c.kinematic_viscosity = kNu;
    c.body_force = {kFx, 0.0, 0.0};

    const auto solve = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!solve.converged)
        throw std::runtime_error("Poiseuille ny=" + std::to_string(ny) + " did not converge");

    PoiseuilleError e;
    for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell) {
        const std::size_t j = cell / nx;
        const double y = (static_cast<double>(j) + 0.5) / static_cast<double>(ny);
        const double err = U.component_data(0)[cell] - exact_u(y);
        e.l2 += err * err;
        e.linf = std::max(e.linf, std::abs(err));
        e.max_abs_uy = std::max(e.max_abs_uy, std::abs(U.component_data(1)[cell]));
    }
    e.l2 = std::sqrt(e.l2 / static_cast<double>(mesh.n_cells()));
    std::cout << "POISEUILLE ny=" << ny << " iterations=" << solve.iterations
              << " L2=" << e.l2 << " Linf=" << e.linf
              << " |Uy|max=" << e.max_abs_uy << "\n";
    return e;
}

} // namespace

int main()
{
    run_case("poiseuille_body_force_second_order", [] {
        const auto e8 = run_poiseuille(8);
        const auto e16 = run_poiseuille(16);
        const auto e32 = run_poiseuille(32);

        const double order_1 = std::log2(e8.l2 / e16.l2);
        const double order_2 = std::log2(e16.l2 / e32.l2);
        std::cout << "POISEUILLE observed order: " << order_1 << " " << order_2 << "\n";

        // Diffusion a deux points + demi-maille paroi : ordre 2 attendu.
        EXPECT_TRUE(order_1 > 1.8 && order_1 < 2.2);
        EXPECT_TRUE(order_2 > 1.8 && order_2 < 2.2);
        // Borne absolue sur le maillage fin : u_max = 1.25, erreur relative < 0.5 %.
        EXPECT_TRUE(e32.linf < 5e-3 * exact_u(0.5));
        // Ecoulement parallele : aucune vitesse transverse.
        EXPECT_TRUE(e32.max_abs_uy < 1e-7);
    });
    return run_all();
}
```

**`tests/validation/test_cht_two_slab_conduction.cpp`** (nouveau). Deux dalles `k1 = 5` et `k2 = 1`, couplées par `solve_two_region_cht`, comparées au profil exact affine par morceaux :

```cpp
// Conduction 1D a travers deux dalles de conductivites differentes, couplees par
// solve_two_region_cht. Remplace le cas mono-cellule de test_cht_validation et l'oracle
// auto-reference de test_analytical_benchmarks:274-297 (audit D6).
//
// Solution exacte (L1 = L2 = 1, pas de source) :
//   q      = (T_hot - T_cold) / (1/k1 + 1/k2)
//   T_int  = T_hot - q / k1
//   T1(x)  = T_hot - q x / k1          x dans [0, 1]
//   T2(x)  = T_int - q (x - 1) / k2    x dans [1, 2]
// Le profil est affine par morceaux : le schema a deux points est exact, l'ecart
// residuel ne vient que de l arret du point fixe (mesure : 1,8e-6 K sur 100 K).
#include "cfdx/physics/cht_solver.h"
#include "common/channel_mesh.h"
#include "common/test_harness.h"

#include <cmath>
#include <cstddef>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {

Mesh shifted_slab(std::size_t nx, double dx)
{
    Mesh m = cfdx::test::make_channel_mesh(nx, 4);
    for (std::size_t i = 0; i < m.points().size(); ++i)
        m.points().set(i, m.points().x(i) + dx, m.points().y(i), m.points().z(i));
    return m;
}

ScalarBoundaryConditions slab_bcs(const char* fixed_patch, double value)
{
    ScalarBoundaryConditions bcs;
    for (const char* patch : {"bottom", "top", "front", "back"})
        bcs[patch] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    bcs[fixed_patch] = {ScalarBoundaryType::FIXED_VALUE, value, 0.0};
    return bcs;
}

} // namespace

int main()
{
    run_case("two_slab_series_conduction_matches_exact_profile", [] {
        constexpr std::size_t nx = 10;
        constexpr double T_hot = 400.0;
        constexpr double T_cold = 300.0;
        constexpr double k1 = 5.0;
        constexpr double k2 = 1.0;
        const double q = (T_hot - T_cold) / (1.0 / k1 + 1.0 / k2);
        const double T_int = T_hot - q / k1;

        Mesh m1 = shifted_slab(nx, 0.0);
        Mesh m2 = shifted_slab(nx, 1.0);
        const auto g1 = build_fv_geometry(m1);
        const auto g2 = build_fv_geometry(m2);

        Field<double, Location::FACE> phi1(m1.n_faces(), "phi1", "kg/s", 1);
        Field<double, Location::FACE> phi2(m2.n_faces(), "phi2", "kg/s", 1);
        phi1.fill(0.0);
        phi2.fill(0.0);
        Field<double, Location::CELL> T1(m1.n_cells(), "T1", "K", 1);
        Field<double, Location::CELL> T2(m2.n_cells(), "T2", "K", 1);
        Field<double, Location::CELL> s1(m1.n_cells(), "s1", "W/m3", 1);
        Field<double, Location::CELL> s2(m2.n_cells(), "s2", "W/m3", 1);
        T1.fill(350.0);
        T2.fill(350.0);
        s1.fill(0.0);
        s2.fill(0.0);

        EnergySolverControls e1;
        e1.conductivity = k1;
        e1.relaxation = 1.0;
        e1.max_iterations = 2000;
    e1.tolerance = 1e-10;
        EnergySolverControls e2 = e1;
        e2.conductivity = k2;

        ChtInterfaceControls c;
        c.region1_patch = "outlet";   // x = 1 cote dalle 1
        c.region2_patch = "inlet";    // x = 1 cote dalle 2
        c.conductivity1 = k1;
        c.conductivity2 = k2;
        c.max_iterations = 500;
        c.relaxation = 1.0;
        c.tolerance = 1e-10;
        c.temperature_tolerance = 1e-10;
        c.matching_tolerance = 1e-12;

        const auto r = solve_two_region_cht(
            m1, g1, m2, g2, phi1, phi2, T1, T2, s1, s2, e1, e2, c,
            slab_bcs("inlet", T_hot), slab_bcs("outlet", T_cold));
        std::cout << "CHT converged=" << r.converged << " iterations=" << r.iterations
                  << " imbalance=" << r.interface_imbalance << "\n";
        EXPECT_TRUE(r.converged);

        double err = 0.0;
        for (std::size_t cell = 0; cell < m1.n_cells(); ++cell) {
            const double x = g1.cell_centres[cell].x;
            err = std::max(err, std::abs(T1(cell) - (T_hot - q * x / k1)));
        }
        for (std::size_t cell = 0; cell < m2.n_cells(); ++cell) {
            const double x = g2.cell_centres[cell].x;
            err = std::max(err, std::abs(T2(cell) - (T_int - q * (x - 1.0) / k2)));
        }
        std::cout << "CHT max |T - T_exact| = " << err << " K (T_int=" << T_int << ")\n";
        // Mesure : 1,8e-6 K apres 356 iterations de point fixe (arret a c.tolerance).
        EXPECT_TRUE(err < 1e-5);
        EXPECT_NEAR(r.interface_imbalance, 0.0, 1e-8);
    });
    return run_all();
}
```

**`tests/validation/test_nonorthogonal_laplacian_mesh.cpp`** (nouveau). Déforme les sommets puis **recalcule** la géométrie. Il remplace la campagne « skew », qui modifiait `face_Sf` à la main :

```cpp
// Laplacien sur un maillage reellement non orthogonal : les sommets interieurs sont
// deplaces, puis la geometrie est RECALCULEE (make_geometry_cache). Remplace
// test_nonorthogonal_skew_campaign, qui modifiait face_Sf a la main et rendait le cache
// geometrique incoherent avec le maillage, avec pour seul critere "correction > 0" (audit D6).
//
// Champ quadratique phi = x^2 + y^2 : laplacien exact = 4 (le maillage est extrude en z).
// On mesure l'erreur sur les cellules dont aucune face n'est une frontiere (l'operateur
// traite les frontieres en gradient nul, sans valeur de bord).
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/numerics/laplacian.h"
#include "common/channel_mesh.h"
#include "common/test_harness.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

constexpr double kPi = 3.14159265358979323846;

// Deplacement lisse, nul sur les bords du carre unite : les faces laterales restent
// planes (meme deplacement sur les deux plans z) et le maillage reste valide.
Mesh make_distorted_mesh(std::size_t n, double amplitude)
{
    Mesh m = cfdx::test::make_channel_mesh(n, n);
    const double h = 1.0 / static_cast<double>(n);
    for (std::size_t i = 0; i < m.points().size(); ++i) {
        const double x = m.points().x(i);
        const double y = m.points().y(i);
        const double bump = std::sin(kPi * x) * std::sin(kPi * y);
        m.points().set(i,
            x + amplitude * h * bump * std::sin(2.0 * kPi * y),
            y + amplitude * h * bump * std::sin(2.0 * kPi * x),
            m.points().z(i));
    }
    return m;
}

struct LaplacianError {
    double orthogonal = 0.0;
    double corrected = 0.0;
};

LaplacianError interior_error(std::size_t n, double amplitude)
{
    const Mesh mesh = make_distorted_mesh(n, amplitude);
    const GeometryCache geometry = make_geometry_cache(mesh);
    ScalarCellField phi(mesh.n_cells(), "phi", "1", 1);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const Vec3 x = geometry.cell_centres[c];
        phi(c) = x.x * x.x + x.y * x.y;
    }
    const auto orth = compute_laplacian(phi, mesh, geometry, LaplacianScheme::ORTHOGONAL);
    const auto corr = compute_laplacian(phi, mesh, geometry, LaplacianScheme::CORRECTED);

    // Cellules "coeur" : ni elles ni leurs voisines ne touchent une frontiere, pour que
    // le gradient de Gauss utilise par la correction soit lui aussi interieur.
    const auto* faces = mesh.cells().faces_data();
    const auto* offsets = mesh.cells().offsets_data();
    std::vector<bool> touches_boundary(mesh.n_cells(), false);
    for (std::size_t c = 0; c < mesh.n_cells(); ++c)
        for (Offset k = offsets[c]; k < offsets[c + 1]; ++k) {
            const auto f = faces[k];
            const bool lateral = std::abs(geometry.face_Sf[f].z) < 1e-12;
            if (lateral && mesh.ownership().neighbour(f) < 0) touches_boundary[c] = true;
        }

    LaplacianError e;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        bool core = !touches_boundary[c];
        for (Offset k = offsets[c]; core && k < offsets[c + 1]; ++k) {
            const auto nb = mesh.ownership().neighbour(faces[k]);
            const auto ow = mesh.ownership().owner(faces[k]);
            const std::size_t other = ow == c ? static_cast<std::size_t>(nb) : ow;
            if (nb >= 0 && touches_boundary[other]) core = false;
        }
        if (!core) continue;
        e.orthogonal = std::max(e.orthogonal, std::abs(orth(c) - 4.0));
        e.corrected = std::max(e.corrected, std::abs(corr(c) - 4.0));
    }
    std::cout << "NONORTH n=" << n << " amplitude=" << amplitude
              << " err_orthogonal=" << e.orthogonal
              << " err_corrected=" << e.corrected << "\n";
    return e;
}

} // namespace

int main()
{
    run_case("orthogonal_mesh_is_exact_for_both_schemes", [] {
        const auto e = interior_error(16, 0.0);
        EXPECT_TRUE(e.orthogonal < 1e-9);
        EXPECT_TRUE(e.corrected < 1e-9);
    });
    run_case("correction_reduces_error_on_distorted_mesh", [] {
        const auto e = interior_error(16, 0.3);
        EXPECT_TRUE(e.orthogonal > 1e-3);          // la distorsion est bien visible
        EXPECT_TRUE(e.corrected < 0.5 * e.orthogonal);
    });
    run_case("corrected_scheme_converges_under_refinement", [] {
        const auto e16 = interior_error(16, 0.3);
        const auto e32 = interior_error(32, 0.3);
        const auto e64 = interior_error(64, 0.3);
        std::cout << "NONORTH corrected ratios: " << e16.corrected / e32.corrected
                  << " " << e32.corrected / e64.corrected << "\n";
        // Mesure : rapports 16,1 et 16,2 (superconvergence sur deformation lisse).
        // Seuil a 4 (ordre 2) : on exige au moins l'ordre nominal du schema corrige.
        EXPECT_TRUE(e16.corrected / e32.corrected > 4.0);
        EXPECT_TRUE(e32.corrected / e64.corrected > 4.0);
        // Le schema orthogonal, lui, ne converge qu'a l'ordre 1 sur ce maillage.
        EXPECT_TRUE(e32.orthogonal / e64.orthogonal < 2.5);
    });
    return run_all();
}
```

**Enregistrement** (`CMakeLists.txt`, après `add_cfdx_test(test_phase9_acceptance …)`) :
```cmake
    add_cfdx_test(test_poiseuille_body_force tests/validation/test_poiseuille_body_force.cpp)
    add_cfdx_test(test_cht_two_slab_conduction tests/validation/test_cht_two_slab_conduction.cpp)
    add_cfdx_test(test_nonorthogonal_laplacian_mesh tests/validation/test_nonorthogonal_laplacian_mesh.cpp)
    set_tests_properties(test_poiseuille_body_force test_cht_two_slab_conduction
        test_nonorthogonal_laplacian_mesh PROPERTIES LABELS "phase13;validation")
```

**Mesures** (`/tmp/cfdx_pr410/_b`, sortie `ctest -V`) :

| Test | Résultat #410 + correctifs | `master` |
|---|---|---|
| Poiseuille | L2 = 1,95e-2 / 4,88e-3 / 1,22e-3 ; **ordre 2,000 et 2,001** ; \|Uy\|max = 3,1e-9 | **FAIL** : le solveur diverge (§5.3) |
| CHT deux dalles | convergé en 356 it ; **max \|T − T_exact\| = 1,77e-6 K** ; déséquilibre d'interface 7,0e-14 | PASS |
| Laplacien non orthogonal | orthogonal exact (0) sans déformation. Avec déformation, err orth / corr = 2,05 / 1,96e-2 (n = 16), 1,08 / 1,22e-3 (32), 0,54 / 7,50e-5 (64) ; **rapports corrigés 16,1 et 16,2** ; l'orthogonal est d'ordre 1 | PASS |

**Statut** : VALIDÉ (compilés, 4/4 PASS sur #410 avec `test_phase9_acceptance`). Le test Poiseuille échoue volontairement sur `master` : c'est le signal attendu.

**Correctif associé : norme du résidu externe de l'énergie** (`src/cfdx/physics/energy_solver.h`, ~l.183, #410) :
- BiCGStab s'arrête sur une norme L2 **relative** ;
- la boucle externe testait un résidu **absolu** en norme infinie, inatteignable dès que ‖b‖ > 1 (T ≈ 400 K) ;
- la boucle tournait donc `max_iterations` fois sans progrès.

```diff
--- a/src/cfdx/physics/energy_solver.h
+++ b/src/cfdx/physics/energy_solver.h
@@ -180,7 +180,19 @@
         sc.tolerance=controls.tolerance;
         sc.relaxation=controls.relaxation;
         const auto linear=solve_scalar_equation(eq,candidate,sc);
-        double res=scalar_equation_residual_inf(eq,candidate);
+        // BiCGStab arrete sur ||b - A x||_2 <= tol*||b||_2 (bicgstab_solver.h:83). Le test
+        // externe doit employer la MEME norme relative : l'ancien residu absolu en norme
+        // infinie etait inatteignable des que ||b|| > 1 (T ~ 400 K) et la boucle tournait
+        // max_iterations fois sans progres.
+        double r2=0.0,b2=0.0;
+        for(std::size_t i=0;i<eq.rhs.size();++i) {
+            double r=-eq.rhs(i);
+            for(auto k=eq.matrix.row_offsets_data()[i];k<eq.matrix.row_offsets_data()[i+1];++k)
+                r+=eq.matrix.values_data()[k]*candidate(eq.matrix.columns_data()[k]);
+            r2+=r*r;
+            b2+=eq.rhs(i)*eq.rhs(i);
+        }
+        const double res=std::sqrt(r2)/std::max(std::sqrt(b2),1e-300);
         for(std::size_t i=0;i<temperature.size();++i)
             temperature(i)=candidate(i);
 
```
**Statut** : VALIDÉ (#410 + correctifs : ctest 96/96).

**Tests encore manquants** (backlog, **aucun code proposé**), chacun relié à un bug du §3 :
1. Poiseuille piloté en **pression** : u_max = 1,5·ū (A-B1, A-M7).
2. Gradient d'un champ linéaire sur maillage non uniforme (A-B1).
3. Advection-diffusion 1D, Pe = 10 (B2).
4. Couette avec p0 ≠ 0 et α_u = 0,7 (A-B2, A-B4).
5. Facteur de forme entre plaques décalées (B5).
6. Canal turbulent Re_τ = 395 contre la DNS de Moser (B1, B7, B9).
7. Cavité chauffée Ra = 1e4, Nu = 2,243 (B12).

### 5.3 🟧 Tolérances de phase9 trop lâches (D9) — CONFIRMÉ, correctif VALIDÉ

**Constat** (#410, `tests/validation/test_phase9_acceptance.cpp`) :
- tolérances de profil L2 = 5e-2 et L∞ = 1e-1 sur un Couette **linéaire** ;
- porte `Umax` en `0,90 < max_u < 1,10`.

Un solveur faux de 5 % passe ce test. Mesure : l'erreur réelle est L2 ≤ 2e-7 et L∞ ≤ 3e-7.

La proposition antérieure `kCouetteTol = 1e-8` était **trop serrée**. Sa valeur `kExpectedMaxU = 0,96875 = 15,5/16` était juste ; elle est reprise ci-dessous.

Le bloc `ENVIRONMENT "CFDX_DEBUG_CELL=33;CFDX_FREEZE_STATE=1;CFDX_DEBUG_COUPLED=1"` pose deux problèmes :
- il injecte des variables de **débogage lues par le solveur de production** (`steady_incompressible_solver.h:1135`) dans le test d'acceptation ;
- `CFDX_FREEZE_STATE` modifie le chemin de calcul.

Il est remplacé par un LABEL qui rend le test visible pour la campagne phase13.

```diff
--- a/tests/validation/test_phase9_acceptance.cpp
+++ b/tests/validation/test_phase9_acceptance.cpp
@@ -352,8 +352,13 @@
              ConvectionScheme::UPWIND, true},
         };
 
-        constexpr double profile_l2_tolerance = 5.0e-2;
-        constexpr double profile_linf_tolerance = 1.0e-1;
+        // Couette lineaire : upwind et SOU sont exacts sur un profil affine, la
+        // seule erreur est l'arret des iterations (mesure : L2 <= 2e-7, Linf <= 3e-7).
+        // Une tolerance de 5e-2 laissait passer un solveur faux a 5 %.
+        constexpr double profile_l2_tolerance = 1.0e-6;
+        constexpr double profile_linf_tolerance = 1.0e-6;
+        // Maximum discret aux centres de cellules : u(y_{ny-1}) = (ny - 1/2)/ny.
+        constexpr double couette_umax_exact = (16.0 - 0.5) / 16.0;
         constexpr double transverse_velocity_tolerance = 1.0e-7;
         constexpr double pressure_uniformity_tolerance = 1.0e-7;
         constexpr double boundary_velocity_tolerance = 1.0e-8;
@@ -474,7 +479,7 @@
                     gates.push_back("Uy");
                 if (!(max_abs_uz < transverse_velocity_tolerance))
                     gates.push_back("Uz");
-                if (!(max_u > 0.90 && max_u < 1.10))
+                if (!(std::abs(max_u - couette_umax_exact) < profile_linf_tolerance))
                     gates.push_back("Umax");
                 if (!(min_u > -boundary_velocity_tolerance))
                     gates.push_back("Umin");
@@ -588,7 +593,7 @@
                       << " momentum_eq_rel=" << h.momentum_equation_residual_relative
                       << " corrected_flux_continuity=" << h.corrected_flux_continuity_linf
                       << " reconstructed_velocity_continuity=" << h.reconstructed_velocity_continuity_linf
-                      << "\\n";
+                      << "\n";
             if (!(error.l2 < profile_l2_tolerance &&
                   error.linf < profile_linf_tolerance))
                 throw std::runtime_error(
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -392,12 +392,9 @@
     add_executable(test_ghia_cavity tests/validation/test_ghia_cavity.cpp)
     target_link_libraries(test_ghia_cavity PRIVATE cfdx_core)
     add_cfdx_test(test_phase9_acceptance tests/validation/test_phase9_acceptance.cpp)
-    # Phase 9 diagnostics must be deterministic under CTest. The validation
-    # workflow may set these variables at the shell level, but CTest test
-    # properties are the authoritative per-test environment and avoid losing
-    # the microscope when the test is launched through CTest.
+    # Porte d'acceptation Phase 9 : selectionnee par la campagne phase13 (audit D9).
     set_tests_properties(test_phase9_acceptance PROPERTIES
-        ENVIRONMENT "CFDX_DEBUG_CELL=33;CFDX_FREEZE_STATE=1;CFDX_DEBUG_COUPLED=1")
+        LABELS "phase9;phase13;validation")
     add_cfdx_test(test_m2_m4_solver_validation tests/validation/test_m2_m4_solver_validation.cpp)
     add_cfdx_test(test_phase10_turbulence_hardening tests/validation/test_phase10_turbulence_hardening.cpp)
     add_cfdx_test(test_phase10_model_registry tests/validation/test_phase10_model_registry.cpp)
```

**Messages de non-convergence illisibles.** `std::to_string` imprime 6 décimales fixes : un résidu de 3e-9 s'affichait `0.000000`.

```diff
--- a/src/cfdx/physics/steady_incompressible_solver.h
+++ b/src/cfdx/physics/steady_incompressible_solver.h
@@ -20,6 +20,7 @@
 #include <limits>
 #include <map>
 #include <numeric>
+#include <sstream>
 #include <stdexcept>
 #include <string>
 #include <vector>
@@ -27,6 +28,16 @@
 
 namespace cfdx::physics {
 
+// std::to_string imprime 6 decimales fixes : un residu de 3e-9 s'affichait "0.000000"
+// dans les messages de non-convergence, ce qui les rendait inexploitables.
+inline std::string format_residual(double value)
+{
+    std::ostringstream os;
+    os.precision(3);
+    os << std::scientific << value;
+    return os.str();
+}
+
 struct VelocityBoundaryCondition {
     enum class Type { FIXED_VALUE, ZERO_GRADIENT };
     Type type = Type::ZERO_GRADIENT;
@@ -1225,8 +1236,8 @@
                     "solve_steady_incompressible: coupled momentum-continuity solve did not converge "
                     "(status=" + std::to_string(static_cast<int>(coupled_result.status)) +
                     ", iterations=" + std::to_string(coupled_result.iterations) +
-                    ", residual=" + std::to_string(coupled_result.residual) +
-                    ", relative=" + std::to_string(coupled_result.residual_relative) + ")");
+                    ", residual=" + format_residual(coupled_result.residual) +
+                    ", relative=" + format_residual(coupled_result.residual_relative) + ")");
             }
 
             // The coupled linear system is a Picard/Newton linearization of the
@@ -1280,8 +1291,8 @@
                     " momentum solve did not converge (status=" +
                     std::to_string(static_cast<int>(solve.status)) +
                     ", iterations=" + std::to_string(solve.iterations) +
-                    ", residual=" + std::to_string(solve.residual) +
-                    ", relative=" + std::to_string(solve.residual_relative) + ")");
+                    ", residual=" + format_residual(solve.residual) +
+                    ", relative=" + format_residual(solve.residual_relative) + ")");
         };
         require_linear_convergence("Ux", rx);
         require_linear_convergence("Uy", ry);
@@ -1481,8 +1492,8 @@
                     "solve_steady_incompressible: pressure-correction solve did not converge "
                     "(status=" + std::to_string(static_cast<int>(rp.status)) +
                     ", iterations=" + std::to_string(rp.iterations) +
-                    ", residual=" + std::to_string(rp.residual) +
-                    ", relative=" + std::to_string(rp.residual_relative) + ")");
+                    ", residual=" + format_residual(rp.residual) +
+                    ", relative=" + format_residual(rp.residual_relative) + ")");
 
             // Relax the physical pressure exactly once. The reference is a
             // gauge: enforce it by a uniform shift, never by overwriting one
```

**Ghia Re = 100** : `add_executable` sans `add_test`. Le test n'était jamais lancé et **échoue sur les deux branches** :
- `master` : diverge, résidu 3,2e13 ;
- #410 : la correction de pression ne converge pas après 5000 it.

Il est enregistré mais `DISABLED`, pour qu'il apparaisse dans ctest :
```cmake
    # Ghia Re=100 : diverge sur master (residu 3,2e13) et ne converge pas sur #410
    # (correction de pression apres 5000 it). Enregistre mais DISABLED : ctest
    # l'affiche "Not Run (Disabled)" au lieu de l'oublier. Reactiver avec l'issue.
    add_cfdx_test(test_ghia_cavity tests/validation/test_ghia_cavity.cpp)
    set_tests_properties(test_ghia_cavity PROPERTIES DISABLED TRUE LABELS "validation;known-failure")
```

**Statut** :
- d9 + format_residual : VALIDÉ, `test_phase9_acceptance` PASS en 11,8 s sur #410 ; `ctest -L phase13` sélectionne 19 tests ;
- Ghia : VALIDÉ (« Not Run (Disabled) ») ;
- `CFDX_FREEZE_STATE` et `CFDX_DEBUG_*` restent lus par le code de production : à retirer avec #415 (§7).

### 5.4 🟧 Tests jamais exécutés (D7) — CONFIRMÉ, correctifs VALIDÉS

**C++ : 8 sources `tests/**/test_*.cpp` compilées par aucune cible.** `test_ghia_cavity` en plus était compilé mais non exécuté.
- Deux sont supprimées au §2.2 : `test_visualization` et `test_vtu_writer`.
- Les autres sont raccrochées dans `CMakeLists.txt`, après `add_cfdx_test(test_thermophysical_models_vv …)` :

```cmake
    # Sources de test qui n'etaient compilees par aucune cible (audit D7).
    add_cfdx_test(test_numerical_temporal tests/numerical/test_temporal.cpp)
    add_cfdx_test(test_application_layer tests/unit/test_application_layer.cpp)
    add_cfdx_test(test_memory_budget tests/unit/test_memory_budget.cpp)
    add_cfdx_test(test_memory_planner tests/unit/test_memory_planner.cpp)
    add_cfdx_test(test_discretization_verification tests/validation/test_discretization_verification.cpp)
    add_cfdx_test(test_numerical_model_verification tests/validation/test_numerical_model_verification.cpp)
```

Dans le bloc `if(CFDX_ENABLE_MPI AND MPI_FOUND)` :
```cmake
        # Orphelins raccroches (audit D7).
        add_executable(test_mpi_deterministic_fused tests/unit/test_mpi_deterministic_fused.cpp)
        target_link_libraries(test_mpi_deterministic_fused PRIVATE cfdx_core MPI::MPI_CXX)
        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND)
            add_executable(test_parallel_hdf5 tests/unit/test_parallel_hdf5.cpp)
            target_link_libraries(test_parallel_hdf5 PRIVATE cfdx_core MPI::MPI_CXX)
            # ... (cibles existantes)
# puis, a cote des add_test MPI existants :
        add_test(NAME test_mpi_deterministic_fused
                 COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                         ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mpi_deterministic_fused> ${MPIEXEC_POSTFLAGS})
        if(CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND)
            add_test(NAME test_parallel_hdf5
                     COMMAND ${MPIEXEC_EXECUTABLE} ${CFDX_MPIEXEC_NUMPROCS_FLAG} 2
                             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_parallel_hdf5> ${MPIEXEC_POSTFLAGS})
```

**`tests/numerical/test_temporal.cpp` ne compilait plus.**
- L'ordre des arguments était faux : la signature réelle est `advance_time(phi, dt, rhs, scheme[, &ctx])`.
- Les oracles CN et BDF2 supposaient des schémas explicites, alors que `advance_time` les résout implicitement.

```diff
--- a/tests/numerical/test_temporal.cpp
+++ b/tests/numerical/test_temporal.cpp
@@ -76,7 +76,7 @@
         phi.fill(1.0);
 
         double dt = 0.1;
-        auto phi_new = advance_time(TimeScheme::EULER_EXPLICIT, phi, dt, decay_rhs);
+        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::EULER_EXPLICIT);
 
         // Exact: φ(t) = φ0 * exp(-t)
         // Euler: φ^{n+1} = φ^n * (1 - dt)
@@ -89,15 +89,12 @@
         phi.fill(1.0);
 
         double dt = 0.1;
-        auto phi_new = advance_time(TimeScheme::CRANK_NICOLSON, phi, dt, decay_rhs);
+        auto phi_new = advance_time(phi, dt, decay_rhs, TimeScheme::CRANK_NICOLSON);
 
-        // CN predictor-corrector:
-        // Predictor: φ* = φ^n + Δt * RHS(φ^n) = 1.0 + 0.1 * (-1.0) = 0.9
-        // Corrector: φ^{n+1} = φ^n + 0.5*Δt * (RHS(φ^n) + RHS(φ*))
-        // = 1.0 + 0.5*0.1 * (-1.0 + -0.9) = 1.0 - 0.05 * 1.9 = 1.0 - 0.095 = 0.905
-        double phi_star = 1.0 + dt * (-1.0);
-        double expected = 1.0 + 0.5 * dt * (-1.0 + -phi_star);
-        EXPECT_NEAR(phi_new(0), expected, 1e-12);
+        // advance_time resout CN implicitement (point fixe converge) :
+        // phi^{n+1} = phi^n + dt/2 (RHS(phi^n) + RHS(phi^{n+1}))  =>  (1 - dt/2)/(1 + dt/2)
+        const double expected = (1.0 - 0.5 * dt) / (1.0 + 0.5 * dt);
+        EXPECT_NEAR(phi_new(0), expected, 1e-10);
     });
 
     run_case("bdf2_requires_history", [&]() {
@@ -109,16 +106,14 @@
         ctx.initialize(phi);
 
         double dt = 0.1;
-        // First step (falls back to Euler)
-        auto phi1 = advance_time(TimeScheme::BDF2, phi, dt, decay_rhs, &ctx);
-        EXPECT_NEAR(phi1(0), 1.0 * (1.0 - dt), 1e-12);
-
-        // Second step (uses BDF2 explicit):
-        // φ2 = (4φ1 - φ0 + 2Δt * RHS(φ1)) / 3
-        // RHS(φ1) = -φ1
-        auto phi2 = advance_time(TimeScheme::BDF2, phi1, dt, decay_rhs, &ctx);
-        double expected = (4.0 * phi1(0) - 1.0 + 2.0 * dt * (-phi1(0))) / 3.0;
-        EXPECT_NEAR(phi2(0), expected, 1e-12);
+        // Premier pas : demarrage par Euler implicite, phi1 = phi0 / (1 + dt)
+        auto phi1 = advance_time(phi, dt, decay_rhs, TimeScheme::BDF2, &ctx);
+        EXPECT_NEAR(phi1(0), 1.0 / (1.0 + dt), 1e-10);
+
+        // Second pas, BDF2 implicite : 3 phi2 - 4 phi1 + phi0 = 2 dt RHS(phi2) = -2 dt phi2
+        auto phi2 = advance_time(phi1, dt, decay_rhs, TimeScheme::BDF2, &ctx);
+        const double expected = (4.0 * phi1(0) - 1.0) / (3.0 + 2.0 * dt);
+        EXPECT_NEAR(phi2(0), expected, 1e-10);
     });
 
     run_case("compute_cfl_time_step", [&]() {
```

**`test_nonorthogonal_skew_campaign.cpp:62`** : le littéral `"\\n"` imprimait un antislash suivi de `n`.
```diff
--- a/tests/validation/test_nonorthogonal_skew_campaign.cpp
+++ b/tests/validation/test_nonorthogonal_skew_campaign.cpp
@@ -59,7 +59,7 @@
             if(skew==0.0) EXPECT_NEAR(correction,0.0,1e-12);
             else EXPECT_TRUE(correction>0.0);
             std::cout<<"PHASE3_6 skew="<<skew<<" orth="<<m.orth<<" corrected="<<m.corrected
-                     <<" correction="<<correction<<" conservation="<<m.conservation<<"\\n";
+                     <<" correction="<<correction<<" conservation="<<m.conservation<<"\n";
         });
     }
     return run_all();
```

**Garde-fou de configuration, nouveau module `cmake/CfdxCheckTestRegistration.cmake`.** Toute source de test orpheline fait échouer `cmake`.

La proposition antérieure (`if(NOT TARGET ${name})`) ne fonctionnait pas, pour deux raisons :
- les noms de cible ne suivent pas les noms de fichier (`test_numerical_temporal`) ;
- les tests MPI n'existent pas sans MPI.

```cmake
# Garde-fou de configuration (audit D7) : toute source tests/**/test_*.cpp doit
# etre compilee par une cible, sinon elle pourrit sans que personne ne le voie
# (constat : 8 sources orphelines, dont test_temporal.cpp qui ne compilait plus).
#
# Usage, a la fin du bloc if(CFDX_BUILD_TESTS) :
#   include(CfdxCheckTestRegistration)
#   cfdx_check_test_registration(EXCLUDE <chemins relatifs a la racine>...)

function(_cfdx_collect_targets dir out_var)
    get_property(targets DIRECTORY "${dir}" PROPERTY BUILDSYSTEM_TARGETS)
    get_property(subdirs DIRECTORY "${dir}" PROPERTY SUBDIRECTORIES)
    foreach(sub IN LISTS subdirs)
        _cfdx_collect_targets("${sub}" sub_targets)
        list(APPEND targets ${sub_targets})
    endforeach()
    set(${out_var} ${targets} PARENT_SCOPE)
endfunction()

function(cfdx_check_test_registration)
    cmake_parse_arguments(ARG "" "" "EXCLUDE" ${ARGN})

    file(GLOB_RECURSE test_sources RELATIVE "${CMAKE_SOURCE_DIR}"
         CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/tests/*/test_*.cpp")

    _cfdx_collect_targets("${CMAKE_SOURCE_DIR}" targets)
    set(compiled "")
    foreach(target IN LISTS targets)
        get_target_property(type ${target} TYPE)
        if(type STREQUAL "INTERFACE_LIBRARY")
            continue()
        endif()
        get_target_property(srcs ${target} SOURCES)
        get_target_property(src_dir ${target} SOURCE_DIR)
        foreach(src IN LISTS srcs)
            if(NOT IS_ABSOLUTE "${src}")
                set(src "${src_dir}/${src}")
            endif()
            file(RELATIVE_PATH rel "${CMAKE_SOURCE_DIR}" "${src}")
            list(APPEND compiled "${rel}")
        endforeach()
    endforeach()

    set(orphans "")
    foreach(src IN LISTS test_sources)
        list(FIND compiled "${src}" idx_compiled)
        list(FIND ARG_EXCLUDE "${src}" idx_excluded)
        if(idx_compiled EQUAL -1 AND idx_excluded EQUAL -1)
            list(APPEND orphans "${src}")
        endif()
    endforeach()

    if(orphans)
        list(JOIN orphans "\n  " orphan_text)
        message(FATAL_ERROR
            "Sources de test compilees par aucune cible :\n  ${orphan_text}\n"
            "Les enregistrer (add_cfdx_test) ou les ajouter a EXCLUDE avec une justification.")
    endif()
endfunction()
```
Appel, à la fin du bloc `if(CFDX_BUILD_TESTS)`, juste avant son `endif()` :
```cmake
    # Tests MPI : compiles uniquement si CFDX_ENABLE_MPI, donc exclus du garde-fou sinon.
    set(_cfdx_backend_only_tests "")
    if(NOT (CFDX_ENABLE_MPI AND MPI_FOUND))
        file(GLOB _cfdx_backend_only_tests RELATIVE "${CMAKE_SOURCE_DIR}"
             "${CMAKE_SOURCE_DIR}/tests/*/test_mpi_*.cpp"
             "${CMAKE_SOURCE_DIR}/tests/*/test_phase8_mpi_*.cpp"
             "${CMAKE_SOURCE_DIR}/tests/*/test_phase5_distributed*.cpp"
             "${CMAKE_SOURCE_DIR}/tests/*/test_parallel_hdf5*.cpp"
             "${CMAKE_SOURCE_DIR}/tests/*/test_phase5_restart_*.cpp")
    elseif(NOT (CFDX_ENABLE_PARALLEL_HDF5 AND HDF5_FOUND))
        file(GLOB _cfdx_backend_only_tests RELATIVE "${CMAKE_SOURCE_DIR}"
             "${CMAKE_SOURCE_DIR}/tests/*/test_phase5_distributed*.cpp"
             "${CMAKE_SOURCE_DIR}/tests/*/test_parallel_hdf5*.cpp"
             "${CMAKE_SOURCE_DIR}/tests/*/test_phase5_restart_*.cpp")
    endif()
    if(NOT CFDX_ENABLE_GPU)
        file(GLOB _cfdx_cuda_only_tests RELATIVE "${CMAKE_SOURCE_DIR}"
             "${CMAKE_SOURCE_DIR}/tests/*/test_*cuda*.cpp")
        list(APPEND _cfdx_backend_only_tests ${_cfdx_cuda_only_tests})
    endif()
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
    include(CfdxCheckTestRegistration)
    cfdx_check_test_registration(EXCLUDE ${_cfdx_backend_only_tests})
```

**Validation du garde-fou**, trois configurations :
1. MPI OFF : la configuration passe.
2. Un fichier sonde orphelin : `FATAL_ERROR`, avec le chemin du fichier.
3. MPI ON (`_b3`) : il a d'abord signalé `test_mpi_deterministic_fused`, `test_parallel_hdf5` et `test_phase5_distributed`. Une fois ceux-ci enregistrés, la configuration passe. `test_mpi_deterministic_fused` est **PASS** à 2 rangs.

**`tests/unit/test_parallel_hdf5.cpp` : `assert` avec effets de bord.** Tous les appels HDF5 (`H5Pset_fapl_mpio`, `H5Dwrite`, `H5Dread`…) sont **dans** des `assert`. Sous `NDEBUG` (Release, le mode de la CI), ils disparaissent : le test « passe » sans rien écrire ni relire. Remplacement :
```diff
# diff normal (sans en-tete) : original -> corrige, tests/unit/test_parallel_hdf5.cpp
10d9
< #include <cassert>
14a14,24
> // CFDX_CHECK() disparait sous NDEBUG (Release) : les appels HDF5 places dans CFDX_CHECK()
> // n'etaient donc jamais executes en CI Release et le test passait sans rien ecrire.
> // CFDX_CHECK evalue toujours son argument et interrompt tous les rangs en cas d'echec.
> #define CFDX_CHECK(expr)                                                              \
>     do {                                                                              \
>         if (!(expr)) {                                                                \
>             std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #expr); \
>             MPI_Abort(MPI_COMM_WORLD, 1);                                             \
>         }                                                                             \
>     } while (0)
>
23c33
<     assert(size == 2);
---
>     CFDX_CHECK(size == 2);
30,31c40,41
<     assert(fapl >= 0);
<     assert(H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);
---
>     CFDX_CHECK(fapl >= 0);
>     CFDX_CHECK(H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);
34,35c44,45
<     assert(file >= 0);
<     H5Pclose(fapl);
---
>     CFDX_CHECK(file >= 0);
>     CFDX_CHECK(H5Pclose(fapl) >= 0);
38c48
<     assert(filespace >= 0);
---
>     CFDX_CHECK(filespace >= 0);
41c51
<     assert(dataset >= 0);
---
>     CFDX_CHECK(dataset >= 0);
45c55
<     assert(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr,
---
>     CFDX_CHECK(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr,
49c59
<     assert(memspace >= 0);
---
>     CFDX_CHECK(memspace >= 0);
55,56c65,66
<     assert(dxpl >= 0);
<     assert(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE) >= 0);
---
>     CFDX_CHECK(dxpl >= 0);
>     CFDX_CHECK(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE) >= 0);
58c68
<     assert(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
---
>     CFDX_CHECK(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
60,64c70,74
<     assert(H5Pclose(dxpl) >= 0);
<     assert(H5Sclose(memspace) >= 0);
<     assert(H5Sclose(filespace) >= 0);
<     assert(H5Dclose(dataset) >= 0);
<     assert(H5Fclose(file) >= 0);
---
>     CFDX_CHECK(H5Pclose(dxpl) >= 0);
>     CFDX_CHECK(H5Sclose(memspace) >= 0);
>     CFDX_CHECK(H5Sclose(filespace) >= 0);
>     CFDX_CHECK(H5Dclose(dataset) >= 0);
>     CFDX_CHECK(H5Fclose(file) >= 0);
69,70c79,80
<     assert(fapl >= 0);
<     assert(H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);
---
>     CFDX_CHECK(fapl >= 0);
>     CFDX_CHECK(H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL) >= 0);
72,73c82,83
<     assert(file >= 0);
<     H5Pclose(fapl);
---
>     CFDX_CHECK(file >= 0);
>     CFDX_CHECK(H5Pclose(fapl) >= 0);
76c86
<     assert(dataset >= 0);
---
>     CFDX_CHECK(dataset >= 0);
78,79c88,89
<     assert(filespace >= 0);
<     assert(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr,
---
>     CFDX_CHECK(filespace >= 0);
>     CFDX_CHECK(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, nullptr,
82c92
<     assert(memspace >= 0);
---
>     CFDX_CHECK(memspace >= 0);
86,88c96,98
<     assert(dxpl >= 0);
<     assert(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE) >= 0);
<     assert(H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
---
>     CFDX_CHECK(dxpl >= 0);
>     CFDX_CHECK(H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE) >= 0);
>     CFDX_CHECK(H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
92c102
<         assert(read_values[i] == static_cast<double>(offset + i));
---
>         CFDX_CHECK(read_values[i] == static_cast<double>(offset + i));
94,98c104,108
<     H5Pclose(dxpl);
<     H5Sclose(memspace);
<     H5Sclose(filespace);
<     H5Dclose(dataset);
<     H5Fclose(file);
---
>     CFDX_CHECK(H5Pclose(dxpl) >= 0);
>     CFDX_CHECK(H5Sclose(memspace) >= 0);
>     CFDX_CHECK(H5Sclose(filespace) >= 0);
>     CFDX_CHECK(H5Dclose(dataset) >= 0);
>     CFDX_CHECK(H5Fclose(file) >= 0);
```
**Statut** : NON COMPILÉ, pas de HDF5 parallèle ici. `cfdx-build.yml` configure `-DCFDX_ENABLE_PARALLEL_HDF5=ON` : il y sera compilé.

**Statut global C++** :
- `_b2` (MPI OFF), `master` + correctifs : ctest 104 tests. Seuls échouent `test_phase9_acceptance` et `test_poiseuille_body_force`, à cause du solveur de `master`. Les deux « Not Run » (`test_visualization`, `test_vtu_writer`) disparaissent avec le §2.2.
- #410 + correctifs : **96/96**.

**Python : la CI n'exécute qu'un seul des 40 fichiers `tests/python/*.py`** (`test_application_e2e.py`). Cinq bugs ont été trouvés et corrigés. Ils sont aussi présents sur `master` :
- `checkpoint.py` depuis `5c63f552` (#134) ;
- `test_output.py` depuis `18362077` (#119) ;
- `validation_report.py` depuis `70b7e139` (#142).

1. `python/cfdx/checkpoint.py:11` : saut de ligne littéral dans la chaîne, **SyntaxError**. Le module ne s'importe pas.
```diff
--- a/python/cfdx/checkpoint.py
+++ b/python/cfdx/checkpoint.py
@@ -8,8 +8,7 @@
 
 
 def save_checkpoint(checkpoint: Checkpoint, path: Path) -> None:
-    path.write_text(json.dumps(checkpoint.__dict__, indent=2, sort_keys=True) + "
-", encoding="utf-8")
+    path.write_text(json.dumps(checkpoint.__dict__, indent=2, sort_keys=True) + "\n", encoding="utf-8")
 
 
 def load_checkpoint(path: Path) -> Checkpoint:
```
2. `tests/python/test_output.py` : le fichier entier tient en **0 ligne**, avec des `\n` littéraux. De plus, le puits attendu par `OutputThrottle` est `Callable[[str, bool], None]` (`output.py:13`), pas `list.append`. Fichier complet de remplacement :
```python
from cfdx.output import OutputThrottle


def test_throttle_batches_and_preserves_stream() -> None:
    emitted = []
    throttle = OutputThrottle(lambda line, is_stderr: emitted.append((line, is_stderr)), interval=0.05, max_lines=2)
    throttle.push("a", False)
    throttle.push("b", True)
    assert throttle.flush(now=1.0) == 2
    assert emitted == [("a", False), ("b", True)]
    throttle.push("c", False)
    assert throttle.flush(now=1.01) == 0
    assert throttle.pending() == 1
    assert throttle.flush(now=1.06) == 1
    assert emitted[-1] == ("c", False)


def test_throttle_validates_limits() -> None:
    emitted = []
    for interval, max_lines in ((0, 10), (0.1, 0)):
        try:
            OutputThrottle(emitted.append, interval=interval, max_lines=max_lines)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid limits must raise ValueError")
```
3. `tests/python/test_validation_report.py` : `@dataclass` exige que le module soit dans `sys.modules`.
```diff
--- a/tests/python/test_validation_report.py
+++ b/tests/python/test_validation_report.py
@@ -1,3 +1,4 @@
+import sys
 from importlib.util import module_from_spec, spec_from_file_location
 from pathlib import Path
 
@@ -9,6 +10,9 @@
 )
 assert SPEC is not None and SPEC.loader is not None
 MODULE = module_from_spec(SPEC)
+# @dataclass resout les annotations differees via sys.modules[cls.__module__] :
+# le module doit y etre enregistre avant exec_module, sinon AttributeError a l'import.
+sys.modules[SPEC.name] = MODULE
 SPEC.loader.exec_module(MODULE)
 
 
```
4. `scripts/validation_report.py` : LaTeX corrompu, quatre défauts.
   - l.206 : `\i` est une séquence invalide (SyntaxError sous `-W error`) et `\t` devient une **tabulation**.
   - l.251, 262 et 277 : `f"\texttt…"` produit une tabulation suivie de `exttt`. Ce défaut est **silencieux**, aucun avertissement.
   - `…\\"` en f-string ne produit qu'**un** antislash : les lignes de `longtable` ne sont pas terminées.
   - `c.ident` n'était pas échappé.
```diff
--- a/scripts/validation_report.py
+++ b/scripts/validation_report.py
@@ -203,7 +203,7 @@
         lines += [
             r"\begin{figure}[h]",
             r"\centering",
-            f"\includegraphics[width=0.92\textwidth]{{{latex_escape(plot_name)}}}",
+            rf"\includegraphics[width=0.92\textwidth]{{{latex_escape(plot_name)}}}",
             r"\caption{Ghia Re=100 centreline RMS error reported by the CFDX solver at three meshes.}",
             r"\end{figure}",
         ]
@@ -248,7 +248,7 @@
         r"CFDX geometry, boundary conditions and coupled solver are executed.",
     ]
     for ident, (rc, _) in logs.items():
-        lines.append(f"\texttt{{{latex_escape(ident)}}}: return code {rc}.")
+        lines.append(rf"\texttt{{{latex_escape(ident)}}}: return code {rc}.")
 
     lines += [
         r"\section{Numerical model verification}",
@@ -259,8 +259,8 @@
     ]
     for m in model_results:
         lines.append(
-            f"\texttt{{{latex_escape(m['name'])}}} & {latex_escape(m['metric'])} & "
-            f"{latex_escape(m['value'])} & {latex_escape(m['reference'])}\\"
+            rf"\texttt{{{latex_escape(m['name'])}}} & {latex_escape(m['metric'])} & "
+            rf"{latex_escape(m['value'])} & {latex_escape(m['reference'])}\\"
         )
     lines += [
         r"\bottomrule",
@@ -274,8 +274,8 @@
     ]
     for c in cases:
         lines.append(
-            f"\texttt{{{c.ident}}} & {latex_escape(c.name)} & {latex_escape(c.status)} & "
-            f"{latex_escape(c.comparison)}\\"
+            rf"\texttt{{{latex_escape(c.ident)}}} & {latex_escape(c.name)} & {latex_escape(c.status)} & "
+            rf"{latex_escape(c.comparison)}\\"
         )
     lines += [
         r"\bottomrule",
```
5. `python/cfdx/case_io.py:192` : `validate_case_bundle` utilise `data` **avant** `data = json.loads(raw)`. D'où un `UnboundLocalError` systématique : un bundle valide est toujours refusé.
```diff
--- a/python/cfdx/case_io.py
+++ b/python/cfdx/case_io.py
@@ -186,6 +186,12 @@
             raise ValueError("mesh data is missing from the self-contained case")
         if "case_hash" not in h5.attrs:
             raise ValueError("case integrity hash is missing")
+        # La configuration doit etre decodee AVANT le controle d'integrite :
+        # `data` etait lu apres son premier usage (UnboundLocalError systematique).
+        raw = h5[_CASE_DATASET][()]
+        if isinstance(raw, bytes):
+            raw = raw.decode("utf-8")
+        data = json.loads(raw)
         mesh_hash = h5.attrs.get("mesh_hash")
         if isinstance(mesh_hash, bytes):
             mesh_hash = mesh_hash.decode("utf-8")
@@ -197,10 +203,6 @@
             raise ValueError("case integrity hash mismatch")
         if "fields" not in h5 or "values" not in h5["fields"]:
             raise ValueError("fields data is missing from the self-contained case")
-        raw = h5[_CASE_DATASET][()]
-        if isinstance(raw, bytes):
-            raw = raw.decode("utf-8")
-        data = json.loads(raw)
         if not isinstance(data.get("physics"), dict):
             raise ValueError("physics configuration is missing")
         if not isinstance(data.get("numerics"), dict):
```
6. `tests/python/test_production_solver_e2e.py` exige `CFDX_PRODUCTION_SOLVER` et `CFDX_PRODUCTION_MESH`, fournis par ctest. Sous un `pytest` nu, il échoue sur `assert None`. On le fait sauter explicitement :
```diff
--- a/tests/python/test_production_solver_e2e.py
+++ b/tests/python/test_production_solver_e2e.py
@@ -4,6 +4,8 @@
 import tempfile
 from pathlib import Path
 
+import pytest
+
 from cfdx import CFDXSession, ExecutionController, SolverRunner
 from cfdx.dat_io import read_dat_restart
 
@@ -27,8 +29,11 @@
 def test_production_solver_full_application_e2e(tmp_path: Path) -> None:
     solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
     mesh = os.environ.get("CFDX_PRODUCTION_MESH")
-    assert solver and Path(solver).is_file()
-    assert mesh and Path(mesh).is_file()
+    if not solver or not mesh:
+        # Fournis par ctest (test_production_solver_e2e) ; un pytest nu ne connait pas le binaire.
+        pytest.skip("CFDX_PRODUCTION_SOLVER / CFDX_PRODUCTION_MESH non definis (lance par ctest)")
+    assert Path(solver).is_file(), solver
+    assert Path(mesh).is_file(), mesh
 
     first_dir = tmp_path / "first"
     session = CFDXSession()
```

**Mesure pytest** (`/tmp/cfdx_pytest_scratch`, Python 3.14, `QT_QPA_PLATFORM=offscreen`, `--ignore` sur `test_gui.py` et `test_gui_3d.py`) :

| | Avant | Après correctifs |
|---|---|---|
| Collecte | 2 fichiers en erreur d'import ou de syntaxe | 0 |
| Résultat | — | **111 passed, 2 skipped, 6 failed** |

Échecs restants, sans correctif proposé (cause non tranchée) :

| Fichier | Échec | Piste |
|---|---|---|
| `test_execution.py` (3) | `'FAILED' != 'CONVERGED'` ; 1 échantillon de moniteur au lieu de 2 | parseur de métriques / redémarrage du contrôleur |
| `test_results_series.py` (2) | une trame VTU illisible est marquée `complete=True` ; `ValueError` non levée | VTK journalise l'erreur de parsing sans lever d'exception : il faut vérifier la sortie du lecteur |
| `test_renderer_pyvista.py` (1) | `.vtu` sauvegardé sur un `PolyData` | le type de maillage n'est pas converti |
| `test_gui.py` | bloque (attente sur une socket) | job xvfb dédié |
| `test_gui_3d.py` | core dump sans OpenGL | job xvfb dédié |
| `test_cfdx_cli.py`, `test_python_api.py` | « no tests ran » | ce sont des scripts sans fonction `test_*` |

**Statut** : correctifs 1 à 6 VALIDÉS (`compileall -W error::SyntaxWarning` OK, fichiers ciblés verts).

---

## 6. CI

### D1 🟥 `master` rouge sans garde-fou — CONFIRMÉ
- Rouge depuis le 2026-09-23 17:18 ; dernier vert : `b2c69cba` (#358).
- Régression introduite par **#365** (`29798aa0`), mergée rouge (D3).
- PRs mergées rouges : #363, #365, #360, #376, #367, #385, #384, #387, #395, #396.
- `master` n'a aucune protection de branche. L'issue #32 a été fermée sans commentaire.

**Protection de branche**, **à exécuter par le propriétaire** (NON EXÉCUTÉ). Les contextes sont les **noms de check réels** lus dans `.github/workflows/` :
- `ci.yml` : job `Build and test (${{ matrix.build_type }})` avec `Release` et `DebugSanitizers` ;
- `cfdx-build.yml` : `build-and-test` ;
- `cfdx-validation.yml` : `M1-M4 validation`.

La proposition antérieure (`build`, `CFDX validation`) ne correspondait à aucun check.

```bash
gh api -X PUT repos/stevendaix/cfdx_dev/branches/master/protection \
  -H "Accept: application/vnd.github+json" --input - <<'JSON'
{
  "required_status_checks": {
    "strict": true,
    "contexts": [
      "Build and test (Release)",
      "Build and test (DebugSanitizers)",
      "build-and-test",
      "M1-M4 validation"
    ]
  },
  "enforce_admins": false,
  "required_pull_request_reviews": null,
  "restrictions": null
}
JSON
```
À n'appliquer **qu'après** le merge de la PR d'hygiène. Sinon, plus aucune PR n'est mergeable tant que `master` est rouge. Le job `Python tests (tests/python)` (D10) sera ajouté aux contextes une fois les 6 échecs du §5.4 corrigés.

### D2 🟥 Validation verte pendant que phase9 échoue — CONFIRMÉ, correctif VALIDÉ
Cause : `test_phase9_acceptance` n'a pas le LABEL `phase13`, et le workflow filtre par label. La correction CMake est au §5.3 (`LABELS "phase9;phase13;validation"`). Le label « acceptance;validation » proposé précédemment ne rejoignait pas la campagne.

Workflow `cfdx-validation.yml` (base #410) : fichier complet dans `/tmp/cfdx_meta/workflows/cfdx-validation.yml`. Changements :
- ajout des 3 nouveaux tests à la liste de build ;
- suppression de l'étape `ctest -L phase13`, redondante avec le pilote ;
- suppression de l'étape « Couette residual » (`CFDX_DEBUG_*`) ;
- le pilote n'est plus en `if: always()`.

```diff
--- a/.github/workflows/cfdx-validation.yml
+++ b/.github/workflows/cfdx-validation.yml
@@ -48,26 +48,20 @@
             test_analytical_benchmarks test_benchmark_matrix \
             test_level_b_reference_benchmarks test_level_c_coupled_verification \
             test_phase9_acceptance test_steady_incompressible_solver \
+            test_poiseuille_body_force test_cht_two_slab_conduction \
+            test_nonorthogonal_laplacian_mesh \
             test_mms_scalar_diffusion test_mesh_refinement_order \
             test_m1_m4_validation test_m2_m4_solver_validation \
             test_cht_validation test_phase10_turbulence_hardening \
             test_thermophysical_models_vv test_matrix_free_fv_operator \
             test_temporal_physics test_core_temporal --parallel 2
 
-      - name: Run Phase 13 validation tests
-        run: |
-          ctest --test-dir build-validation -L phase13 --output-on-failure
-
-      - name: Run Couette residual validation
-        if: always()
-        run: |
-          set -o pipefail
-          CFDX_DEBUG_CELL=auto CFDX_DEBUG_COUPLED=1 ctest --test-dir build-validation \
-            -R 'test_(phase9_acceptance|steady_incompressible_solver)' \
-            -V 2>&1 | tee build-validation/couette-residuals.log
-
+      # Le pilote resout chaque test de campagne dans le registre CTest, l execute
+      # seul et echoue si un test manque, echoue ou n a pas le LABEL phase13.
+      # Plus d etape "Couette residual" : CFDX_DEBUG_* ne sont pas des options de
+      # validation (lues par le solveur) et test_phase9_acceptance est desormais
+      # labellise phase13, donc execute par le pilote.
       - name: Run campaign driver
-        if: always()
         run: |
           .venv-validation/bin/python scripts/run_validation.py \
             --build-dir build-validation \
```
**Statut** : YAML valide (PyYAML). `actionlint` non exécuté (absent ici). Le pilote est VALIDÉ (§5.1).

### D4 🟧 Build cassé avec MPI = OFF — CONFIRMÉ, correctif VALIDÉ
**Correction d'une affirmation antérieure** : le défaut n'est pas dans `mpi_utils.h`. `src/cfdx/core/linalg/communication_avoiding.h` inclut `mpi_utils.h` et définit `mpi_fused_reduction` sans garde, et il est inclus par du code non MPI. Le `#error` proposé avant aurait cassé le build **davantage**.
```diff
--- a/src/cfdx/core/linalg/communication_avoiding.h
+++ b/src/cfdx/core/linalg/communication_avoiding.h
@@ -1,6 +1,8 @@
 #pragma once
 #include "cfdx/core/linalg/vector.h"
+#if CFDX_HAS_MPI
 #include "cfdx/core/parallel/mpi_utils.h"
+#endif
 #include <algorithm>
 #include <cmath>
 #include <cstddef>
@@ -25,6 +27,7 @@
     return r;
 }
 
+#if CFDX_HAS_MPI
 // Batch the two additive Krylov reductions into one MPI collective.
 // max_abs is reduced separately because it requires MPI_MAX rather than SUM.
 // This API is deliberately independent of a solver so distributed Krylov
@@ -72,4 +75,6 @@
     return mpi_fused_reduction(fused_reduction(a, b), comm);
 }
 
+#endif // CFDX_HAS_MPI
+
 } // namespace cfdx::core
```
**Statut** : VALIDÉ (`_b2` MPI OFF : configuration et build complets ; `_b3` MPI ON : `test_mpi_*` compilés).

### D5 🟧 « HDF5-dependent features disabled » est faux — CONFIRMÉ, correctif VALIDÉ
**Correction d'une affirmation antérieure** : HDF5 n'est pas optionnel. `hdf5_writer.cpp`, `hdf5_reader.cpp`, `mesh_importer.cpp` et `gmsh_importer.cpp` l'incluent sans garde. Rendre HDF5 optionnel demande un chantier (§3). Le correctif honnête est d'échouer tôt, avec un message utile :
```diff
@@ CMakeLists.txt ~l.86 @@
 else()
-    message(WARNING "HDF5 not found. HDF5-dependent features disabled.")
+    # hdf5_writer.cpp, hdf5_reader.cpp, mesh_importer.cpp et gmsh_importer.cpp
+    # incluent HDF5 sans garde : cfdx_core ne peut pas etre construit sans lui.
+    message(FATAL_ERROR
+        "HDF5 introuvable : cfdx_core en depend inconditionnellement "
+        "(io/hdf5, io/mesh, io/gmsh). Installer libhdf5-dev, passer "
+        "-DHDF5_ROOT=<prefixe>, ou construire third_party/hdf5.")
 endif()
```
**Statut** : VALIDÉ. Sans `HDF5_ROOT`, la configuration s'arrête avec le message ; avec, elle passe.

### D10 🟧 Couverture CI
Constats :
- GPU jamais compilé ;
- MPI à 2 rangs seulement ;
- LSAN désactivé (`CMakeLists.txt:501-502`) ;
- `ci.yml:6-7` filtre `pull_request` par branche : les PRs empilées (#415→#411, #409→#406, #403→#397) n'ont **aucun** job `CFDX CI` ;
- aucun job MPI = OFF ;
- pytest limité à 1 fichier sur 40.

**Correctif `ci.yml`** : filtre de branche retiré, nouveau job `python-tests`. Fichier complet dans `/tmp/cfdx_meta/workflows/ci.yml`.
```diff
--- a/.github/workflows/ci.yml
+++ b/.github/workflows/ci.yml
@@ -3,8 +3,9 @@
 on:
   push:
     branches: [master, main, develop]
+  # Sans filtre de branche : les PR empilees (#415->#411, #409->#406, #403->#397)
+  # visent une branche de fonctionnalite et ne declenchaient aucun job CI.
   pull_request:
-    branches: [master, main, develop]
   workflow_dispatch:
 
 concurrency:
@@ -86,3 +87,43 @@
           name: compile-commands
           path: build/compile_commands.json
           if-no-files-found: error
+
+  python-tests:
+    name: Python tests (tests/python)
+    runs-on: ubuntu-24.04
+    timeout-minutes: 20
+    steps:
+      - name: Checkout
+        uses: actions/checkout@v4
+
+      - name: Install Python test stack
+        run: |
+          python3 -m venv .venv
+          .venv/bin/python -m pip install --upgrade pip
+          .venv/bin/python -m pip install -e . pytest h5py meshio
+
+      - name: Byte-compile every Python source
+        # Echoue sur une erreur de syntaxe (python/cfdx/checkpoint.py:11 ne compilait pas)
+        # et, avec -W error, sur les sequences d echappement invalides (validation_report.py:206).
+        run: .venv/bin/python -W error::SyntaxWarning -m compileall -q python scripts tests/python
+
+      - name: Run pytest
+        # test_gui.py (bloque sur une socket) et test_gui_3d.py (core dump sans GL) exigent
+        # un affichage : a sortir dans un job xvfb dedie. Mesure de l audit apres correctifs :
+        # 111 passed, 2 skipped, 6 failed (test_execution x3, test_results_series x2,
+        # test_renderer_pyvista x1) : ce job restera ROUGE tant qu ils ne sont pas corriges,
+        # ne le rendre obligatoire qu ensuite.
+        run: >-
+          .venv/bin/python -m pytest -q tests/python
+          --ignore=tests/python/test_gui.py --ignore=tests/python/test_gui_3d.py
+          --junitxml=pytest-results.xml
+        env:
+          QT_QPA_PLATFORM: offscreen
+
+      - name: Upload pytest results
+        if: always()
+        uses: actions/upload-artifact@v4
+        with:
+          name: pytest-results
+          path: pytest-results.xml
+          if-no-files-found: warn
```
Vérification locale :
- `compileall -W error::SyntaxWarning` sur la version d'origine détecte bien `checkpoint.py:11` et `validation_report.py:206` (rc = 1) ;
- les `\t` silencieux des l.251, 262 et 277 ne sont **pas** détectables ainsi, d'où le correctif explicite du §5.4.

**Statut** : YAML valide ; `actionlint` non exécuté ; le job restera rouge tant que les 6 échecs pytest restent ouverts.

**Job MPI = OFF supplémentaire** (NON EXÉCUTÉ en CI ; la même configuration est VALIDÉE localement dans `_b2`). À ajouter sous `jobs:` dans `ci.yml` :
```yaml
  build-no-mpi:
    name: Build and test (MPI OFF)
    runs-on: ubuntu-24.04
    timeout-minutes: 60
    steps:
      - uses: actions/checkout@v4
      - name: Install dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y cmake ninja-build g++ zlib1g-dev libhdf5-dev python3-venv
      - name: Configure (MPI OFF)
        run: cmake -S . -B build-nompi -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF
      - name: Build
        run: cmake --build build-nompi --parallel 4
      - name: Test
        run: ctest --test-dir build-nompi --output-on-failure
```
Il n'y a pas de job « sans HDF5 » : HDF5 est obligatoire (D5). La compilation CUDA (`-DCFDX_ENABLE_GPU=ON`, et non `CFDX_ENABLE_CUDA`) dans `nvidia/cuda:12.4.1-devel-ubuntu22.04` reste à écrire après la vérification de la chaîne CUDA du §3 : NON PROPOSÉ en code.

### Traces de débogage en production — CONFIRMÉ, correctif VALIDÉ
`src/cfdx/physics/finite_volume_transport.h:386` et `:403` (`master`), soit l.418 et 435 sur #410, contiennent deux défauts :
- un `std::cerr` à **chaque** résolution linéaire dès que `n ≤ 256` ;
- un terminateur `'\\n'`, caractère multi-octets qui vaut `23662` et que la sortie ctest affiche collé au nombre suivant (`residual=6.55533e-1723662`).

```diff
--- a/src/cfdx/physics/finite_volume_transport.h
+++ b/src/cfdx/physics/finite_volume_transport.h
@@ -383,9 +383,6 @@
     auto result = cfdx::core::solve_bicgstab(
         equation.matrix, equation.rhs, candidate,
         controls.max_iterations, controls.tolerance);
-    if (solution.size() <= 256)
-        std::cerr << "CFDX solver cascade: bicgstab status=" << static_cast<int>(result.status)
-                  << " iter=" << result.iterations << " residual=" << result.residual << '\\n';
 
     // Keep all retries anchored to the same nonlinear iterate; the accepted
     // predictor is updated only after a solver reports convergence.
@@ -400,9 +397,6 @@
         result = cfdx::core::solve_gmres(
             equation.matrix, equation.rhs, candidate,
             64, controls.max_iterations, controls.tolerance);
-        if (solution.size() <= 256)
-            std::cerr << "CFDX solver cascade: gmres status=" << static_cast<int>(result.status)
-                      << " iter=" << result.iterations << " residual=" << result.residual << '\\n';
     }
 
     if (result.status != cfdx::core::SolverStatus::CONVERGED) {
```
Garde-fou suggéré contre la récidive (NON COMPILÉ), à placer après `project()` :
```cmake
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Werror=multichar>)
endif()
```
**Statut** : dbg.diff VALIDÉ sur `master` (`_b2`). À appliquer aussi sur #410.

### RK2 à faible stockage d'ordre 1 — CONFIRMÉ, correctif VALIDÉ
`src/cfdx/physics/low_storage_time_integration.h` enchaîne deux demi-pas d'Euler explicite, ce qui ne donne que l'ordre 1. Remplacement par le RK2 2N de Williamson (point milieu) :
```diff
--- a/src/cfdx/physics/low_storage_time_integration.h
+++ b/src/cfdx/physics/low_storage_time_integration.h
@@ -60,11 +60,23 @@
     if (!(dt > 0.0) || !std::isfinite(dt))
         throw std::invalid_argument("invalid dt");
     const std::size_t n = u.size();
+    // Williamson 2N-storage RK2 (A = {0, -1/2}, B = {1/2, 1}) : identique au
+    // point milieu u^{n+1} = u^n + dt f(u^n + dt/2 f(u^n)), ordre 2, deux
+    // registres (u, dq). L'ancienne version enchainait deux demi-pas d'Euler
+    // explicite et n'etait que d'ordre 1.
     cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs");
+    cfdx::core::Field<double,cfdx::core::Location::CELL> dq(n, "rk2_register");
+    dq.fill(0.0);
     rhs(u, k);
-    for (std::size_t i = 0; i < n; ++i) u(i) += 0.5 * dt * k(i);
+    for (std::size_t i = 0; i < n; ++i) {
+        dq(i) = dt * k(i);
+        u(i) += 0.5 * dq(i);
+    }
     rhs(u, k);
-    for (std::size_t i = 0; i < n; ++i) u(i) += 0.5 * dt * k(i);
+    for (std::size_t i = 0; i < n; ++i) {
+        dq(i) = -0.5 * dq(i) + dt * k(i);
+        u(i) += dq(i);
+    }
 }
 
 } // namespace cfdx::physics
```
**Statut** : VALIDÉ. Ordre observé sur `du/dt = −u` : **2,05 / 2,03 / 2,01**, contre ≈ 1 avant.

---

## 7. Workflow PR / issues

### 7.1 Carte des PRs ouvertes (révisée)

| PR | Relation | Couette borné | Couette non borné | Poiseuille | phase9 locale |
|---|---|---|---|---|---|
| `master` | — | diverge | diverge | diverge | ❌ |
| #406 | ⊂ #410 | converge | échec (0,33) | L∞ 0,42 | ❌ |
| #409 | ⊂ #411 | converge | converge | L∞ 0,42 | ❌ |
| #410 | CLEAN, seule mergeable | converge | converge | — | ❌ tel quel |
| **#410 + correctifs D4, D9, énergie, format_residual** | scratch | converge | converge | **ordre 2,00** (nouveau test) | ✅ **ctest 96/96** |
| #411 | ⊂ #415 | converge | converge | L∞ 0,42 | ❌ (COUPLED) |
| #415 | tête de chaîne | converge | converge | L∞ 0,42 | ❌ (COUPLED) ; lit `CFDX_FREEZE_STATE` (l.1204-1210) |

Le L∞ de 0,42 sur Poiseuille des PRs de la chaîne est la signature d'A-B1 (Green-Gauss, §3.1).

Le nouveau test Poiseuille (force de volume, parois FIXED_VALUE) atteint l'ordre 2 sur #410 + correctifs. Il ne couvre **pas** le Poiseuille piloté en pression (backlog §5.2, n° 1).

### 7.2 Recommandation (révisée)
La recommandation précédente (« garder #409 → #411 → #415 ») est **remplacée**. Motif : #410 est la seule PR sans conflit, et avec les correctifs de cet audit elle est la seule branche **mesurée verte** (96/96).

Base retenue : **#410 + correctifs** :
- D4 (`d4.diff`) ;
- D9 (`d9.diff`) ;
- `energy.diff` ;
- `format_residual.diff` ;
- `dbg.diff` ;
- les 3 tests D6 et `channel_mesh.h` ;
- le garde-fou D7 ;
- le pilote D8.

Puis **#415**, rebasée sur `master` une fois #410 + correctifs mergée :
- avec ses correctifs du §10.4 (restart GMRES, oracles, A-B1, `'\\n'`, C6) ; avec eux, COUPLED converge et **n'a plus à être sorti des gates** ;
- sans `CFDX_FREEZE_STATE` ni `CFDX_DEBUG_*` dans le code de production (§10.4.9).

#409 et #411 sont **entièrement contenues** dans #415 (§10.1) : à fermer, pas à merger séparément.

**Arbitrage avec la variante « #415 seule » (examinée au §10.1).** Merger #415 seule et fermer #410 est
possible, mais fait perdre les 24 commits de #410 absents de #415 (dont le workflow
`cfdx-validation.yml` et les tests associés), et part d'une branche mesurée rouge. L'ordre
« #410 d'abord » part de la seule base mesurée verte ; le coût est un rebase de #415 avec conflit
sur 5 fichiers, dont un seul conflit de fond :
- `finite_volume_transport.h` : les deux correctifs `diag += F` sont équivalents (`d9d865a1` côté #410), garder celui de #410 ;
- `tests/unit/test_finite_volume_transport.cpp:72-76` : garder **les deux** cas, le `EXPECT_THROW` de #415 et le cas `extra_diagonal{10}` / `EXPECT_NEAR(diag, 10)` de #410 (`b01f0498`) ;
- `cfdx-validation.yml`, `preconditioner.h`, `steady_incompressible_solver.h` : prendre #410 puis réappliquer les ajouts de #415.

### 7.3 Issues (D15, D16)
- #32 et #112 ont été fermées sans preuve : les rouvrir, ou les commenter avec le lien vers le commit et le test.
- #382 est déclarée corrigée par #383, qui est encore ouverte.
- #388 (enthalpie) et #380 (E2E) ont été corrigées par #396 mais restent ouvertes. B6 montre que l'enthalpie reste fausse en extrapolation basse.
- #59 (orientation) : le contrôle `Sf·(xf−xc)` (§3.1 A-M3) n'existe toujours pas.
- #74 (VTU) : résolue par `fef59b99`, à fermer.
- La feuille de route #34 est trop optimiste au regard des §3 et §5.

### 7.4 Triage et ordre de merge (à faire par le propriétaire)
1. PR « hygiène » sur la base #410 : §2 (D11, D12), D1 (après merge), D2, D4, D5, D7, D8, D9, D10, dbg, RK2, Python.
2. Fermetures : #373 (après récupération de `test_phase11_thermal_acceptance.cpp`), #374, #375, #381, #391, #397, #400, #402, #403, #406, #408, #409.
3. #415 rebasée sur le résultat de l'étape 1, avec les correctifs §10.8 (restart GMRES, oracles, A-B1, m1, C6), sans `CFDX_FREEZE_STATE`. #409 et #411 fermées (contenues dans #415).
4. #383.
5. #398 puis #414, avec le correctif B5 appliqué aux deux.
6. Protection de branche (§6 D1).

**Commandes**, **à exécuter par le propriétaire** (NON EXÉCUTÉES par l'audit) :
```bash
R=stevendaix/cfdx_dev
# Doublons / PR contenues ailleurs
gh pr close 374 -R $R -c "Doublon de #375, integre par #385. Ferme suite a l'audit (AUDIT_cfdx_dev.md §7)."
gh pr close 375 -R $R -c "Integre par #385. Ferme suite a l'audit (§7)."
gh pr close 391 -R $R -c "Doublon de #397. Ferme suite a l'audit (§7)."
gh pr close 397 -R $R -c "Perimee (base obsolete), doublon de #391. Ferme suite a l'audit (§7)."
gh pr close 403 -R $R -c "Perimee, empilee sur #397 fermee. Ferme suite a l'audit (§7)."
gh pr close 406 -R $R -c "Contenue dans #410 (retenue comme base, §7.2)."
gh pr close 408 -R $R -c "Contenue dans #406/#410."
gh pr close 409 -R $R -c "Entierement contenue dans #415 (audit §10.1)."
gh pr close 411 -R $R -c "Ancetre de #415, qui la contient (audit §10.1) ; #415 sera rebasee sur #410 corrigee (§7.2)."
gh pr close 381 -R $R -c "Ferme suite a l'audit (§7.4)."
gh pr close 400 -R $R -c "Ferme suite a l'audit (§7.4)."
gh pr close 402 -R $R -c "Ferme suite a l'audit (§7.4)."
# #373 : recuperer d'abord le test, puis fermer
git fetch origin pull/373/head:pr373
git switch chore/hygiene   # branche de la PR d hygiene
git checkout pr373 -- tests/validation/test_phase11_thermal_acceptance.cpp   # puis l enregistrer (add_cfdx_test)
gh pr close 373 -R $R -c "test_phase11_thermal_acceptance.cpp recupere dans la PR d'hygiene ; le reste est perime."
# Issues
gh issue close 74 -R $R -c "Resolue par fef59b99."
gh issue reopen 32 -R $R && gh issue comment 32 -R $R -b "Rouverte : fermee sans preuve ; master rouge depuis #365 (audit §6 D1)."
gh issue reopen 112 -R $R && gh issue comment 112 -R $R -b "Rouverte : fermee sans commit ni test de preuve (audit §7.3)."
gh issue comment 388 -R $R -b "Corrigee par #396 pour le domaine nominal ; l'extrapolation basse reste fausse (audit B6). A garder ouverte ou a scinder."
```
Pour #380 et #382, il faut attendre le merge de #396 et de #383 avant de commenter ou de fermer.

---

## 8. Optimisations

Chaque gain est qualifié de l'une de trois façons :
- **mesuré** : banc en scratch, ancien et nouveau code compilés à l'identique ;
- **estimé** : aucune mesure ;
- **contredit** : la mesure infirme l'estimation de la version précédente de cet audit.

Le code des items patchés se trouve dans `sections/s4_s8.patch`. Pour les correctifs de justesse associés, voir §4.

| # | Optimisation | Où | Ancienne estimation | Résultat | Statut |
|---|---|---|---|---|---|
| 1 | Bounds-check dans `Vector::operator()` ; `-O2` forcé ; pas d'OpenMP | `vector.h:37,41`, `CMakeLists.txt:53-55` | ×2 à ×4 | bounds-check 0-8 % (bruit) ; `-O3` : 0 ; OpenMP SpMV ×2 au mieux | **contredit**, non compilé |
| 2 | AMG en préconditionneur du CG pression | `steady_incompressible_solver.h:578` | ×5-×20 it | 645→44 it, ×2,3 en temps (Dirichlet) ; 952→66 it, ×2,0 (Neumann) | **mesuré** (AMG-PCG), branchement non compilé |
| 3 | Espaces de travail Krylov préalloués | CG, GMRES, BiCGStab | 10-30 % | CG −25 % ; GMRES(30) −10 % | **mesuré** |
| 4 | BiCGStab : 2 matvecs par itération au lieu de 3 | `bicgstab_solver.h` | ~30 % | −35 % de temps (et 469→366 it) | **mesuré** |
| 5 | Halo non bloquant, voisins seulement | `distributed_execution.h:191-250` | scalabilité | justesse vérifiée np=2/4 ; gain de temps non mesuré | compilé+testé, gain **estimé** |
| 6 | LU grossier AMG factorisé une fois | `amg_preconditioner.h` | important | inclus dans le 2 ; rend aussi le Neumann soluble | **mesuré** (via 2) |
| 7 | Pool de tampons GPU | `cuda_poisson.cu` | latence | — | **estimé**, non compilé |
| 8 | RCM par CSR | `memory_planner.cpp:336-460` | pré-traitement | 27k cellules : 1,958 s → 12,1 ms (×160) | **mesuré** |
| 9 | Import OpenFOAM sans regex | `openfoam_importer.cpp` | ×10 | ×4,5 (512k cellules) | **mesuré**, ×10 corrigé |
| 10 | VTU binaire (base64 ou appended) | `vtu_writer.cpp:123-269` | taille ÷3, temps ÷5 | — | **estimé**, non compilé |
| 11 | Try/catch de `distributed_poisson` | `distributed_poisson.h:74-87` | — | matvec 10,75 → 0,76 ms (×14) | **mesuré** (nouveau) |

Priorisation après mesure :
1. items 2 et 6 (AMG) ;
2. item 8 (RCM, le seul à devenir un blocage à 1M cellules) ;
3. item 11 ;
4. item 9.

L'item 1, qui était classé premier, tombe en dernier.

### 8-1. Noyaux, bounds-check, options de compilation — **contredit**, code non compilé
**Code actuel** :
- `vector.h:37,41` : `operator()` appelle `check_index` (qui lève `out_of_range`, l.131). `data()` est défini aux l.53-54.
- `CMakeLists.txt:53-55` :

```cmake
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug" AND NOT CMAKE_BUILD_TYPE STREQUAL "DebugSanitizers")
    add_compile_options(-O2 -DNDEBUG)      # écrase le -O3 de Release et les flags utilisateur
endif()
```

**Mesures** (SpMV, n = 1e6) :

| Configuration | Temps |
|---|---|
| Bounds-check | 0-8 % d'écart, dans le bruit |
| `-O3` au lieu de `-O2` | aucun gain |
| OpenMP, 1 thread | 5,2 ms |
| OpenMP, 4 threads | 2,7 ms |
| OpenMP, 8 threads | 2,4 ms |
| OpenMP, 16 threads | 7,6 ms (dégradé) |

Le noyau est limité par la bande passante mémoire. Le « ×2 à ×4 » est donc faux. Les solveurs patchés utilisent déjà `data()` dans les boucles chaudes.

**Code proposé** (intérêt : hygiène, et ×2 sur la SpMV avec OpenMP) :
```cmake
# supprimer les l.53-55 : laisser CMAKE_BUILD_TYPE décider (-O3 -DNDEBUG en Release)
option(CFDX_NATIVE "Optimise for host CPU" OFF)
option(CFDX_OPENMP "Enable OpenMP kernels" OFF)
if(CFDX_NATIVE) add_compile_options(-march=native) endif()
if(CFDX_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(cfdx_core PUBLIC OpenMP::OpenMP_CXX)   # nom de cible à confirmer
endif()
```
```cpp
// vector.h : vérification en Debug seulement
double& operator()(std::size_t i) {
#ifndef NDEBUG
    check_index(i);
#endif
    return data_[i];
}
```
```cpp
// SpMV CSR (si CFDX_OPENMP), à plafonner au nombre de cœurs physiques d'un socket
#pragma omp parallel for schedule(static)
for (std::size_t r = 0; r < n; ++r) { double s = 0.0;
    for (auto k = rp[r]; k < rp[r + 1]; ++k) s += val[k] * x[col[k]]; y[r] = s; }
```
- `tests/test_vector.cpp:36`, qui attend `out_of_range` en Release, doit être restreint à Debug.
- Le nom de la cible CMake n'a pas été vérifié.

### 8-2. AMG en préconditionneur de la pression — AMG-PCG **mesuré**, branchement **non compilé**
**Code actuel** : `steady_incompressible_solver.h:578` utilise `solve_cg(A, b, p_corr, ...)`, c'est-à-dire Jacobi-CG. `AgglomeratedAMGPreconditioner` n'est utilisé nulle part en production.

**Code proposé** (repose sur `solve_pcg` et N-AMG du §4) :
```cpp
// IncompressibleSolverControls (l.55) : nouvelle option
enum class PressurePreconditioner { JACOBI, AMG };
PressurePreconditioner pressure_preconditioner = PressurePreconditioner::JACOBI;  // défaut inchangé

// l.578
SolverResult rp;
if (controls_.pressure_preconditioner == PressurePreconditioner::AMG) {
    FunctionalLinearOperator op(A.n_rows(), [&A](const Vector& in, Vector& out) { A.multiply(in, out); });
    AgglomeratedAMGPreconditioner amg(op, 0.7, 2, 2);
    if (!amg.setup(A)) throw std::runtime_error("pressure AMG setup failed");
    rp = solve_pcg(A, b, p_corr, &amg, controls_.linear_max_iterations, controls_.linear_tolerance);
} else {
    rp = solve_cg(A, b, p_corr, controls_.linear_max_iterations, controls_.linear_tolerance);
}
// + contrôle de statut C2
```
Points à vérifier avant intégration :
- le nom exact de la matvec de `SparseMatrix`, écrit ici `A.multiply` ;
- la réutilisation du `setup` d'une itération externe à l'autre quand la matrice de pression ne change que par ses coefficients.

**Mesure** (256², tolérance 1e-8) :

| | Jacobi-CG | AMG-PCG |
|---|---|---|
| Dirichlet | 645 it, 0,375 s | 44 it, 0,166 s |
| Neumann | 952 it, 0,567 s | 66 it, 0,278 s |

Soit ×15 en itérations et ×2 à ×2,3 en temps. Le gain en temps croît avec la taille du maillage (estimé).

### 8-3. Espaces de travail préalloués — **mesuré**, compilé+testé
**Code actuel** :
- CG : allocations dans la boucle (`cg_solver.h:47-220`).
- GMRES :
  - `Vector vi_vector(n)` pour chaque produit scalaire du MGS ;
  - réallocation et remise à zéro de (m+1)·n + m·n doubles quand le restart change.

**Code proposé** :
```cpp
Vector r(n), z(n), p(n), Ap(n);                                  // CG, une fois
GmresWorkspace w; w.resize(n, static_cast<std::size_t>(restart_cap));   // stride fixe = restart_cap
const double h = krylov_dot(vi, w.w.data(), n, SolverPrecision::FP64, controls.reduction);  // MGS sans copie
```
**Mesure** (b_krylov, n = 90k) :

| Solveur | Avant | Après |
|---|---|---|
| CG | 0,277 s | 0,207 s |
| GMRES(30) | 1,116 s | 1,005 s |

### 8-4. BiCGStab à 2 matvecs — **mesuré**, compilé+testé
**Code actuel** : `bicgstab_solver.h` fait une 3ᵉ matvec par itération pour recalculer le résidu vrai.

**Code proposé** : le résidu vrai est recalculé seulement sur convergence apparente et toutes les `true_residual_period = 50` itérations. Ce nouveau paramètre est placé en dernière position, avec une valeur par défaut.
```cpp
res = nrm(w.r);                                       // résidu récursif
const bool periodic = true_residual_period > 0 && iter % true_residual_period == 0;
if (res <= tol || periodic) {
    res = true_residual();                            // matvec de contrôle, rare
    if (res <= tol) return finish(SolverStatus::CONVERGED, iter, res);
}
```
**Mesure** (n = 90k) : 469 it en 1,79 s avant, 366 it en 1,17 s après. La baisse du nombre d'itérations vient des redémarrages sur breakdown relatif (C17).

### 8-5. Halo non bloquant — compilé+testé, gain **estimé**
Le code est donné en §4, N2. Deux changements :
- les pairs passent de P−1 à la liste des vrais voisins ;
- 4 `Sendrecv` bloquants par pair sont remplacés par 1 `Irecv` et 1 `Isend`, suivis d'un seul `Waitall`.

Le gain attendu croît avec P. Il n'a pas été mesuré à np > 4.

### 8-6. LU grossier factorisé une fois — **mesuré** via 8-2
Le code est donné en §4, N-AMG. Le coût O(n³) du grossier sort de chaque `apply` et passe au `setup`.

### 8-7. Pool de tampons GPU — **estimé**, non compilé (pas de nvcc)
**Code actuel** : `cuda_poisson.cu` fait des `cudaMalloc` et `cudaFree` à chaque appel du solveur.

**Code proposé** :
```cpp
struct DeviceBuffers {
    double *x = nullptr, *r = nullptr, *p = nullptr, *Ap = nullptr; std::size_t n = 0;
    void ensure(std::size_t m) {
        if (m <= n) return;
        release();
        for (double** b : {&x, &r, &p, &Ap})
            if (cudaMalloc(b, m * sizeof(double)) != cudaSuccess) { release(); throw std::bad_alloc(); }
        n = m;
    }
    void release() noexcept { for (double* b : {x, r, p, Ap}) cudaFree(b); x = r = p = Ap = nullptr; n = 0; }
    ~DeviceBuffers() { release(); }
};
// membre du solveur GPU ; ensure(n) au début de solve()
```
À combiner avec les critères C15/C19 du §4.

### 8-8. RCM par adjacence CSR — **mesuré**, compilé+testé
**Code actuel** : `memory_planner.cpp:336-460` parcourt toutes les faces pour chaque cellule visitée, soit O(Nc·Nf).

**Code proposé** :
```cpp
std::vector<size_t> adj_offset(num_cells + 1, 0);
for (size_t f = 0; f < n_faces; ++f)
    if (valid_edge(f)) { ++adj_offset[owner[f] + 1]; ++adj_offset[neighbour[f] + 1]; }
for (size_t c = 0; c < num_cells; ++c) adj_offset[c + 1] += adj_offset[c];
std::vector<uint32_t> adj(adj_offset[num_cells]);                     // remplissage par `fill`
std::vector<uint32_t> start_order(num_cells);                         // (degré, indice) : départs déterministes
std::iota(start_order.begin(), start_order.end(), 0u);
std::sort(start_order.begin(), start_order.end(), by_degree);
size_t head = ordering.size();                                         // la file BFS est la fin de `ordering`
ordering.push_back(component_start);
while (head < ordering.size()) {
    const uint32_t current = ordering[head++];
    for (size_t k = adj_offset[current]; k < adj_offset[current + 1]; ++k) { /* voisins non visités, triés par degré */ }
}
```
**Mesure** (permutation identique à l'ancien code, largeur de bande 26 853 → 690) :

| Cellules | Avant | Après |
|---|---|---|
| 1k | 2,76 ms | 0,35 ms |
| 8k | 166 ms | 3,2 ms |
| 27k | 1,958 s | 12,1 ms |

À 1M cellules, l'ancien code prendrait environ 45 min (estimé, par extrapolation quadratique).

### 8-9. Import OpenFOAM sans regex — **mesuré**, compilé+testé
Le code est donné en §4, C20. Les hashes du maillage importé sont identiques avant et après.

| Cas | Maillage | Avant | Après |
|---|---|---|---|
| c10 | — | 0,011 s | 0,002 s |
| c40 | 69k points, 197k faces, 64k cellules | 0,655 s | 0,147 s |
| c80 | 531k points, 1,56M faces, 512k cellules, 75 Mo | 5,98 s | 1,33 s |

Le gain sur c80 est de **×4,5**. Le « ×10 » annoncé précédemment est corrigé. Surtout, l'import ne crashe plus sur un long commentaire.

### 8-10. VTU binaire — **estimé**, non compilé
**Code actuel** : `vtu_writer.cpp:123-269`. Tous les DataArray sont en `format="ascii"`, avec `setprecision(12)` et le format scientifique, soit environ 20 octets par double.

**Code proposé** : format `appended` en `encoding="raw"`, sans dépendance à base64.
```cpp
// en-tête : <DataArray type="Float64" Name="T" format="appended" offset="OFF"/>
// ... puis, après </UnstructuredGrid> :
os << "  <AppendedData encoding=\"raw\">\n_";
auto put = [&](const void* p, std::uint64_t bytes) {
    os.write(reinterpret_cast<const char*>(&bytes), sizeof bytes);   // header_type="UInt64"
    os.write(static_cast<const char*>(p), static_cast<std::streamsize>(bytes));
};
for (const auto& blk : blocks) put(blk.data(), blk.size() * sizeof(double));
os << "\n  </AppendedData>\n";
// <VTKFile ... header_type="UInt64" byte_order="LittleEndian">, flux ouvert en std::ios::binary
```
Gain estimé : ÷2,5 en taille (8 octets contre environ 20 par double), et beaucoup plus en temps car il n'y a plus de formatage `%e`.

Deux précautions :
- les offsets doivent être calculés en première passe ;
- l'ASCII doit rester disponible en option pour le débogage.

---

## 9. Plan d'action priorisé

Chaque étape renvoie à des correctifs **déjà écrits** dans ce rapport, et le plus souvent déjà
compilés et testés en copie de travail (fichiers dans `correctifs/`, Annexe A).

| Étape | Contenu | Constats | Correctifs | Critère de sortie |
|---|---|---|---|---|
| **0. Rendre le CI fiable** (≈ 1 jour) | `.gitignore` + `git rm --cached build*` ; gardes MPI (`communication_avoiding.h`) et HDF5 ; label `phase13` sur phase9 ; pilote `run_validation.py` qui échoue sur rouge ; protection de `master` avec les vrais noms de check ; `-Werror=multichar` | D1, D2, D4, D5, D8, D10, D11, m1 | `s567/d4.diff`, `s567/run_validation.py`, `s567/workflows/`, §6 D1 | CI rouge sur `master` pour la bonne raison (phase9) ; build OK avec MPI/HDF5 OFF ; pilote rc = 1 sur `master` |
| **1. Justesse du couplage P-V** (≈ 2 à 3 jours) | Green-Gauss ; relaxation implicite ; ancrage par décalage ; flux de masse persistant ; sortie gradient nul ; statut du CG pression ; NaN Krylov | A-B1 à A-B4, A-M1, A-M2, A-M6/7, C1, C2 | `s3a/`, `s4_s8/s4_s8.patch` | Couette p0 ≠ 0 à 1e-8 ; Poiseuille u_max à 1e-3, ordre ≥ 1,8 ; phase9 verte avec les tolérances resserrées (§5.3) |
| **2. Justesse physique** (≈ 1 semaine) | Diffusion turbulente ; cp énergie ; signe et linéarisation radiatifs ; P1 π ; Rosseland ; SA fw ; lois de paroi ; k-ω/SST ; enthalpie ; Boussinesq ; BDF2 variable ; RK2 | B1 à B18, RK2 | `s3b/s3b.diff`, `s3b/b17.diff`, `s567/rk2.diff`, `s567/energy.diff` | Tests par bug 10/10 (§3.6.5) ; puis advection-diffusion Pe = 10, cavité Ra = 1e4, canal Re_τ = 395, facteur de forme décalé |
| **3. Géométrie** | Orientation contrôlée ; `compute_cell_geometry_oriented` partout ; FIXED_GRADIENT ; centres exacts ; non-orthogonalité | A-M3 à A-M5, m6 | `s3a/` | MMS ordre ≥ 1,8 sur maillage tétraédrique et déformé |
| **4. Robustesse solveurs et I/O** | **Restart GMRES non écrêté et politique adaptative remise à l'endroit** ; breakdown d'Arnoldi ; oracles GMRES ; CG vrai résidu ; tolérance relative ; doublons CSR ; LU grossier AMG ; restart atomique et borné ; HDF5 multi-champ ; import OpenFOAM sans regex ; ordre des halos | C3 à C13, C17 à C20, N-DUP, N-AMG, N1 à N4, §10.4.1 à 10.4.4 | `s4_s8/s4_s8.patch`, `pr415/gmres.diff` | Tests dédiés (NaN, ‖b‖ < 1, restart demandé respecté, système singulier, restart tronqué) ; `test_gmres_solver` 11/11 |
| **5. Validation réelle** | Tests tautologiques remplacés ; tests C++ et Python enregistrés ; `validation_report.py` corrigé ; matrice VMFL avec métriques. **Job pytest ajouté mais non requis** tant que les 6 échecs restants ne sont pas corrigés | D6 à D9 | `s567/d6/`, `s567/d9.diff`, `s567/validation_report.diff`, `s567/tpl/` | Au moins 1 cas VMFL PASS mesuré, publié en artefact CI ; pytest 0 échec avant de le rendre requis |
| **6. PRs** | #410 + correctifs d'abord (96/96) ; puis #415 rebasée avec ses correctifs §10.8 ; fermeture de #409, #411 et des doublons (§7.4) | D14, D15, §10 | `pr415/`, §7.4 | CI verte sur les deux workflows de #415 ; plus de doublons |
| **7. Performance** (priorités mesurées, §8) | 1) AMG en préconditionneur de pression (×2 à ×2,3), LU grossier factorisé une fois ; 2) RCM par CSR (×160) ; 3) matvec distribué sans try/catch (×14) ; 4) espaces de travail préalloués, BiCGStab à 2 matvecs ; 5) import OpenFOAM (×4,5). OpenMP ensuite ; bounds-check en dernier (0 à 8 %) | §8 | `s4_s8/s4_s8.patch` | Profil avant/après sur un cas de 1M de cellules |
| **8. Nettoyage** | Doublons, code mort P-V, `linear_algebra/`, traces de débogage | D12, m5, dbg | `s567/dbg.diff` | Aucun en-tête non inclus ; aucun `getenv`/`cerr` sous `core/linalg/` |

> **Règle pour la suite** : chaque correctif de §3 arrive avec **le test qui l'aurait détecté**
> (§5.2, liste « tests manquants »). Les tests actuels passent sur du code faux parce qu'ils
> n'exercent ni gradient de pression, ni flux de masse non nul, ni NaN, ni restart GMRES réel.

---

## 10. Revue détaillée de la PR #415

> Revue faite en lecture seule sur le dépôt : aucun commit, aucun push, aucun commentaire GitHub.
> Tout a été reproduit dans une copie de travail (`/tmp/cfdx_pr415n`, tête `f829e32f`) et dans
> une copie corrigée (`/tmp/cfdx_meta/scratch_415/fixc`, issue d'un `git archive`).
> Statut de chaque code proposé : **compilé + testé**, **compilé**, ou **non compilé**.
> Ce qui n'a pas été reproduit est marqué **PLAUSIBLE**.

### 10.1 Périmètre et relations avec #409, #410 et #411

| Élément | Valeur |
|---|---|
| Tête analysée | `f829e32f` (branche locale `pr/415-new`). La ref distante `pr/415` pointe encore l'ancienne tête `1443e313` |
| `master..415` | 53 commits, 15 fichiers, +2886 / −330 lignes |
| #409 (`0fa79e66`) | **entièrement contenue** dans #415 (0 commit hors #415), merge-tree propre |
| #411 (`a0fd6e66`) | **ancêtre** de #415 : 38 commits `master..411`, puis 11 commits `411..415` |
| #410 (`d89aa8dc`) | 24 commits absents de #415, **conflit sur 5 fichiers** : `cfdx-validation.yml`, `preconditioner.h`, `finite_volume_transport.h`, `steady_incompressible_solver.h`, `tests/unit/test_finite_volume_transport.cpp:72-76` |

Pour #410, le correctif de transport `d9d865a1` (`diag += F`) est **équivalent** à celui de #415 ;
seuls les commentaires diffèrent. Le seul conflit de fond est le test : #415 attend
`EXPECT_THROW` alors que #410 passe `extra_diagonal{10}` et vérifie `EXPECT_NEAR(diag, 10)`. Le
paramètre `extra_diagonal` existe déjà dans #415 (l. 170).

**Recommandation** (arbitrée avec le §7.2) :
- ne merger ni #409 ni #411 séparément : les fermer comme contenues dans #415 ;
- merger d'abord #410 + correctifs (seule base mesurée verte, 96/96) ;
- rebaser ensuite #415, en gardant le test de `b01f0498` comme **cas supplémentaire** à côté de `EXPECT_THROW`.

Une variante « #415 seule, #410 fermée » est techniquement possible. Elle n'est pas retenue : les 24 commits
de #410 absents de #415 devraient alors être repris un par un (§7.2).

### 10.2 Build, ctest et CI

| Contexte | Résultat |
|---|---|
| Build local de la PR (Release, sans MPI, HDF5 `/tmp/hdf5_inst`) | OK |
| ctest local de la PR | **90/93** : échecs `test_gmres_solver`, `test_phase9_acceptance`, et `test_performance_runtime` (*Not Run* : `mpi.h` absent, défaut D4 déjà présent sur master) |
| CI GitHub sur `f829e32f` | **Les 2 workflows échouent** (`test_gmres_solver` et `test_phase9_acceptance`, 2/103). Toutes les têtes précédentes ont échoué ou été annulées. État GitHub : `MERGEABLE` / `UNSTABLE` |
| ctest de la copie corrigée `fixc` | **92/93** : seul reste `test_performance_runtime` (*Not Run*, `mpi.h`, préexistant) |
| `phase9` sur `fixc` | **PASS**, `successful=6 failed=0` |

Commande de build utilisée :

```bash
PATH=/tmp/hdf5_inst/bin:$PATH cmake .. -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF -DHDF5_ROOT=/tmp/hdf5_inst
make -j && ctest --output-on-failure
```

La description de la PR affirme que le bug d'échappement `'\\n'` est corrigé. **C'est faux** :
cinq occurrences subsistent (voir 10.4.5).

### 10.3 Cause racine de l'échec COUPLED (phase9)

**Constat.** Un run instrumenté (`CFDX_DEBUG_COUPLED=1`) de phase9 sur la PR montre ce qui suit.

La matrice est saine :

```
COUPLED_MATRIX_SUMMARY n=512 nnz=9236 M_l2=11.54 G_l2=1.64 D_l2=1.56 C_l2=1.23
zero_rows=0 zero_cols=0 schur_diag_min=0.0498 supplied_schur_max_delta=0.965 gauge_row=384
```

Le restart GMRES, lui, décroît : 40 → 35 → 30 → … → 10. Le résidu vrai reste bloqué à
1,33e-6 (relatif 1,18e-6, pour une tolérance de 1e-10) jusqu'à l'itération 2000. On obtient
alors `MAX_ITER` ; la reprise en BiCGStab finit aussi en `MAX_ITER`, et une exception est levée
dès la **première** itération externe.

Deux défauts combinés en sont la cause.

1. **`gmres_solver.h:43`**, le restart demandé est écrêté en silence :
   ```cpp
   int current_restart = std::clamp(restart, controls.restart_min, controls.restart_max);
   ```
   Les valeurs par défaut de `KrylovControls` sont `restart_min=10`, `restart_max=40` et
   `adaptive_restart=true`. Le `GMRES(512)` que demande le commit `f829e32f` devient donc
   `GMRES(40)` : **ce commit n'a aucun effet**. C'est aussi le cas du `restart=64` du solveur de
   pression, et du `128` du repli identité.
2. **`krylov_controls.h:32-35`**, la politique adaptative est **inversée** :
   ```cpp
   if (cycle_reduction > c.target_reduction_per_cycle * 2.0)   // mauvaise réduction
       return std::max(c.restart_min, current - 5);             // -> on RÉTRÉCIT l'espace
   if (cycle_reduction < c.target_reduction_per_cycle * 0.5)   // bonne réduction
       return std::min(c.restart_max, current + 5);             // -> on l'AGRANDIT
   ```
   Sur un système de point-selle qui stagne, l'espace de Krylov est ainsi réduit à chaque cycle
   jusqu'à 10. C'est exactement la trajectoire 40 → 10 observée.

**Preuve par correction.** Avec les correctifs 10.4.1 et 10.4.2 :
- phase9 passe ;
- COUPLED converge en **48 itérations**, sans reprise BiCGStab ;
- l'erreur de profil vaut L2 = 2,0e-9 et Linf = 3,1e-9 ; la continuité vaut 1,4e-10 ;
- l'`ALGORITHM_INVARIANCE` avec SIMPLE donne `max_abs_dU = 2,73e-7`.

**Remarque.** Dans le harnais 10.6, COUPLED/Couette converge **même avec le binaire de la PR**,
alors que phase9 échoue sur ce binaire. Seule l'initialisation diffère : p = 0 dans le harnais,
p ≠ 0 dans phase9. Sans le correctif, la robustesse de COUPLED dépend donc de l'état initial.

### 10.4 Revue ligne à ligne et corrections proposées

#### 10.4.1 [BLOQUANT] Restart GMRES écrêté — **compilé + testé**

Code actuel, `pr/415:src/cfdx/core/linalg/gmres_solver.h:43-44` :

```cpp
int current_restart = std::clamp(restart, controls.restart_min, controls.restart_max);
current_restart = std::min<int>(current_restart, static_cast<int>(n));
```

Code proposé : le restart de l'appelant fait foi et n'est borné que par n. La fenêtre
adaptative est élargie pour le contenir.

```cpp
// The caller's restart is authoritative: it is only bounded by n. The
// adaptive window is widened to contain it instead of clamping it to the
// default [10, 40], which silently turned GMRES(512) into GMRES(40).
int current_restart = std::min<int>(restart, static_cast<int>(n));
controls.restart_max = std::max(controls.restart_max, current_restart);
controls.restart_min = std::min(controls.restart_min, current_restart);
```

(`controls` est déjà une copie locale dans `solve_gmres`, donc la modifier n'a aucun effet de
bord chez l'appelant.)

Test (à ajouter à `test_gmres_solver.cpp`, **non compilé** ; le correctif, lui, est compilé et testé via phase9) : sur un système non symétrique où `GMRES(40)`
stagne, `GMRES(n)` doit converger en au plus n itérations.

```cpp
run_case("gmres_requested_restart_is_not_clamped_to_default_window", []() {
    // Upwind 1D convection-diffusion, n = 64: full GMRES converges in <= n its.
    const std::size_t n = 64;
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        A.push_back(i, i, 2.5);
        if (i > 0) A.push_back(i, i - 1, -2.0);
        if (i + 1 < n) A.push_back(i, i + 1, -0.5);
    }
    A.finalize();
    Vector b(n, 1.0), x(n, 0.0);
    const auto r = solve_gmres(A, b, x, static_cast<int>(n), n, 1e-12);
    EXPECT_TRUE(r.status == SolverStatus::CONVERGED);
    EXPECT_TRUE(r.iterations <= n);   // impossible if restart were clamped to 40
});
```

#### 10.4.2 [BLOQUANT] Politique de restart adaptatif inversée — **compilé + testé**

Code actuel, `pr/415:src/cfdx/core/linalg/krylov_controls.h:32-35` (cité en 10.3).

Code proposé :

```cpp
if (!c.adaptive_restart) return std::clamp(current, c.restart_min, c.restart_max);
// Poor reduction per cycle = the Krylov space is too small: enlarge it.
if (cycle_reduction > c.target_reduction_per_cycle * 2.0)
    return std::min(c.restart_max, current + 5);
// Fast reduction: a smaller space is enough and cheaper (memory, orthogonalisation).
if (cycle_reduction < c.target_reduction_per_cycle * 0.5)
    return std::max(c.restart_min, current - 5);
return std::clamp(current, c.restart_min, c.restart_max);
```

Test à mettre à jour, `tests/unit/test_performance_runtime.cpp:37` :

```cpp
// poor reduction (0.5 per cycle) must enlarge the Krylov space
const int restart = choose_gmres_restart(30, 0.5, controls);
if (restart <= 30) { std::cerr << "check 7 failed\n"; return 1; }
```

Ce test n'est pas compilable ici (`mpi.h`, D4). La politique elle-même est validée
indirectement par phase9 et par le harnais.

#### 10.4.3 [BLOQUANT] Gestion du breakdown d'Arnoldi — **compilé + testé**

Code actuel : `pr/415:src/cfdx/core/linalg/gmres_solver.h:124-142` et `:204-253`.
Il présente cinq défauts.

- `H(j+1,j) = hnext` est posé **avant** le test de plancher (l. 125). Un `hnext` au niveau de
  l'arrondi entre donc dans la rotation de Givens.
- `vnext = w / hnext` est calculé dès que `hnext > 0` (l. 140-142). Cela injecte un vecteur de
  bruit O(1) dans la base.
- La branche `std::abs(diag) <= 1e-30` (l. 204) est morte en pratique, puisque `diag = rho`
  vaut au moins l'arrondi.
- Sur un breakdown non « heureux », la fonction rend un itéré inutilisable. Mesuré sur le cas
  `diag(1,0)` : `status=1`, `x = (1.5, 4.2e15)`.
- Des traces `getenv("CFDX_DEBUG_COUPLED")` + `std::cerr` (`GMRES_BREAKDOWN`,
  `GMRES_BREAKDOWN_TRUE_RESIDUAL`, `GMRES_CYCLE`, l. 214, 229, 252) se trouvent dans un en-tête
  de bibliothèque.

Code proposé : il remplace le corps de la boucle d'Arnoldi et la fin de cycle.
`<iostream>` et `<cstdlib>` sont remplacés par `<limits>`.

```cpp
            const double hnext = krylov_norm2(w.w, SolverPrecision::FP64, controls.reduction);
            // Breakdown test BEFORE normalisation: dividing a roundoff-level
            // w by hnext would inject an O(1) noise vector into the basis.
            double hcolumn_scale = 0.0;
            for (int i = 0; i <= j; ++i)
                hcolumn_scale = std::max(
                    hcolumn_scale,
                    std::abs(w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j))));
            const double breakdown_floor =
                128.0 * std::numeric_limits<double>::epsilon() *
                std::max(1.0, hcolumn_scale);
            const bool breakdown = !(hnext > breakdown_floor);
            w.H(static_cast<std::size_t>(j + 1), static_cast<std::size_t>(j)) =
                breakdown ? 0.0 : hnext;
            if (!breakdown) {
                double* vnext = w.v(static_cast<std::size_t>(j + 1));
                for (std::size_t k = 0; k < n; ++k) vnext[k] = w.w(k) / hnext;
            }

            for (int i = 0; i < j; ++i) {
                const double h0 = w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
                const double h1 = w.H(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j));
                w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) =
                    w.cs[static_cast<std::size_t>(i)] * h0 + w.sn[static_cast<std::size_t>(i)] * h1;
                w.H(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(j)) =
                    -w.sn[static_cast<std::size_t>(i)] * h0 + w.cs[static_cast<std::size_t>(i)] * h1;
            }

            const double a = w.H(static_cast<std::size_t>(j), static_cast<std::size_t>(j));
            const double b2 = w.H(static_cast<std::size_t>(j + 1), static_cast<std::size_t>(j));
            const double rho = std::hypot(a, b2);
            ++iterations;
            if (!(rho > breakdown_floor)) {
                // Column j is linearly dependent on the previous ones: the
                // operator is singular on the Krylov space (non-happy
                // breakdown). Drop the column; the least-squares solution
                // over the first `used` columns is still well defined.
                arnoldi_breakdown = true;
                break;
            }
            w.cs[static_cast<std::size_t>(j)] = a / rho;
            w.sn[static_cast<std::size_t>(j)] = b2 / rho;
            w.H(static_cast<std::size_t>(j), static_cast<std::size_t>(j)) = rho;
            w.H(static_cast<std::size_t>(j + 1), static_cast<std::size_t>(j)) = 0.0;
            const double gj = w.g[static_cast<std::size_t>(j)];
            w.g[static_cast<std::size_t>(j)] = w.cs[static_cast<std::size_t>(j)] * gj;
            w.g[static_cast<std::size_t>(j + 1)] = -w.sn[static_cast<std::size_t>(j)] * gj;
            ++used;
            estimated_residual = std::abs(w.g[static_cast<std::size_t>(j + 1)]);

            if (controls.residual_replacement &&
                controls.residual_replacement_period > 0 &&
                iterations % controls.residual_replacement_period == 0) {
                const double exact = true_residual();
                if (exact <= tol) {
                    result.status = SolverStatus::CONVERGED;
                    result.iterations = iterations;
                    result.residual = exact;
                    result.residual_relative = exact / std::max(b_norm, 1.0);
                    return result;
                }
            }

            if (estimated_residual <= tol || breakdown) {
                // breakdown with rho > floor is a happy breakdown: the
                // Krylov space is invariant and g[j+1] is the exact residual.
                arnoldi_breakdown = breakdown && estimated_residual > tol;
                break;
            }
        }

        if (used == 0) {
            // A*z0 is in span{}: no admissible direction at all.
            result.status = arnoldi_breakdown ? SolverStatus::DIVERGED
                                              : SolverStatus::MAX_ITER_REACHED;
            result.iterations = iterations;
            result.residual = beta;
            result.residual_relative = beta / std::max(b_norm, 1.0);
            return result;
        }

        std::fill(w.y.begin(), w.y.end(), 0.0);
        for (int i = used - 1; i >= 0; --i) {
            double sum = w.g[static_cast<std::size_t>(i)];
            for (int j = i + 1; j < used; ++j)
                sum -= w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) *
                       w.y[static_cast<std::size_t>(j)];
            const double diag = w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(i));
            w.y[static_cast<std::size_t>(i)] = sum / diag;   // diag = rho > floor by construction
        }

        for (int j = 0; j < used; ++j) {
            const double alpha = w.y[static_cast<std::size_t>(j)];
            const double* zj = w.zv(static_cast<std::size_t>(j));
            for (std::size_t i = 0; i < n; ++i) x(i) += alpha * zj[i];
        }

        beta = true_residual();
        if (beta <= tol) {
            result.status = SolverStatus::CONVERGED;
            result.iterations = iterations;
            result.residual = beta;
            result.residual_relative = beta / std::max(b_norm, 1.0);
            return result;
        }
        if (arnoldi_breakdown) {
            // Non-happy breakdown: restarting cannot enlarge the Krylov space.
            // Report the TRUE residual of the least-squares iterate.
            result.status = SolverStatus::DIVERGED;
            result.iterations = iterations;
            result.residual = beta;
            result.residual_relative = beta / std::max(b_norm, 1.0);
            return result;
        }

        const double reduction = beta / std::max(previous_cycle_residual, 1e-300);
        if (controls.adaptive_restart)
            current_restart = choose_gmres_restart(current_restart, reduction, controls);
        previous_cycle_residual = beta;
```

#### 10.4.4 [BLOQUANT] Oracles faux dans `test_gmres_solver.cpp` — **compilé + testé (11/11)**

Code actuel, `pr/415:tests/unit/test_gmres_solver.cpp:147-172` :
- La matrice `[[1,0,1],[0,1,1],[-1,-1,0]]` avec b = (1,1,3) est donnée comme singulière (le
  commentaire dit « u=v=1 ») et le test attend `DIVERGED`. Or **det = 2** : le système est
  régulier, et sa solution exacte est (−1,5 ; −1,5 ; 2,5). Mesuré avant correctif, sur les deux
  têtes de la PR : `status=0` (CONVERGED) avec cette solution. Le test est rouge **parce que le
  solveur a raison**.
- Pour le cas « non heureux » A = diag(1,0), b = (1,1), le test attend `x(1) = 0`. La solution
  des moindres carrés sur span{b} est pourtant **x = (1,1)**, avec un résidu vrai exactement
  égal à 1.

Code proposé : les trois cas complets sont dans `fixc/tests/unit/test_gmres_solver.cpp:147-208`.

```cpp
run_case("gmres_nonsingular_saddle_point_converges", []() {
    // [ 1  0  1 ] [u]   [1]      det = 2 : the system is regular.
    // [ 0  1  1 ] [v] = [1]      u + p = 1, v + p = 1, -u - v = 3
    // [-1 -1  0 ] [p]   [3]      => u = v = -1.5, p = 2.5
    SparseMatrix A(3, 3);
    A.push_back(0, 0, 1.0); A.push_back(0, 2, 1.0);
    A.push_back(1, 1, 1.0); A.push_back(1, 2, 1.0);
    A.push_back(2, 0, -1.0); A.push_back(2, 1, -1.0);
    A.finalize();
    Vector b(3); b(0) = 1.0; b(1) = 1.0; b(2) = 3.0;
    Vector x(3, 0.0);
    const auto result = solve_gmres(A, b, x, 3, 3, 1e-12);
    EXPECT_TRUE(result.status == SolverStatus::CONVERGED);
    EXPECT_NEAR(x(0), -1.5, 1e-10);
    EXPECT_NEAR(x(1), -1.5, 1e-10);
    EXPECT_NEAR(x(2), 2.5, 1e-10);
});

run_case("gmres_singular_saddle_point_reports_true_residual", []() {
    // [ 1  0  1 ]   det = 0, right null vector (1, 1, -1),
    // [ 0  1  1 ]   left null vector y = (1, -1, -1).
    // [ 1 -1  0 ]   b = (1, 1, 3): y.b = -3 != 0 -> inconsistent,
    //               global min ||b - A x|| = |y.b| / ||y|| = sqrt(3).
    SparseMatrix A(3, 3);
    A.push_back(0, 0, 1.0); A.push_back(0, 2, 1.0);
    A.push_back(1, 1, 1.0); A.push_back(1, 2, 1.0);
    A.push_back(2, 0, 1.0); A.push_back(2, 1, -1.0);
    A.finalize();
    Vector b(3); b(0) = 1.0; b(1) = 1.0; b(2) = 3.0;
    Vector x(3, 0.0);
    const auto result = solve_gmres(A, b, x, 3, 30, 1e-12);
    EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
    EXPECT_TRUE(std::isfinite(result.residual));
    EXPECT_TRUE(result.residual >= std::sqrt(3.0) - 1e-8);
    // GMRES minimises over the Krylov space only: >= sqrt(3), and
    // never above the initial residual ||b|| = sqrt(11).
    EXPECT_TRUE(result.residual <= std::sqrt(11.0) + 1e-12);
    for (std::size_t i = 0; i < 3; ++i) EXPECT_TRUE(std::isfinite(x(i)));
});

run_case("gmres_nonhappy_breakdown_reports_true_residual", []() {
    // A = diag(1, 0), b = (1, 1): K_2 = span{(1,1),(1,0)} is rank 1 under
    // A. The least-squares iterate over span{b} is x = (1, 1) and the
    // true residual is exactly 1 (the unreachable second component).
    SparseMatrix A(2, 2);
    A.push_back(0, 0, 1.0);
    A.finalize();
    Vector b(2); b(0) = 1.0; b(1) = 1.0;
    Vector x(2, 0.0);
    const auto result = solve_gmres(A, b, x, 2, 2, 1e-12);
    EXPECT_TRUE(result.status == SolverStatus::DIVERGED);
    EXPECT_NEAR(result.residual, 1.0, 1e-12);
    EXPECT_NEAR(x(0), 1.0, 1e-12);
    EXPECT_TRUE(std::isfinite(x(1)));
    EXPECT_TRUE(std::abs(x(1)) < 10.0);
});
```

#### 10.4.5 [MAJEUR] A-B1 non corrigé : gradient de Gauss faux côté voisin — **compilé + testé**

Code actuel, `pr/415:src/cfdx/physics/steady_incompressible_solver.h:166-169` :

```cpp
if (own.neighbour(f) >= 0) {
    const std::size_t n = static_cast<std::size_t>(own.neighbour(f));
    vf = 0.5 * (field(c) + field(n));
```

Quand c est le **voisin** de la face, `n = neighbour(f) = c` : la valeur de face vaut alors
`field(c)`, et la contribution de la face au gradient disparaît. Le gradient de pression est
donc faux sur la moitié des faces internes. Conséquence mesurée : Poiseuille SIMPLE donne
u_max = 0,593 au lieu de 1, et **aucun ordre de convergence** (−0,06).

Le bloc G du système couplé (l. ~647) utilise en outre `coeff = 0.5` en dur. Il faut prendre la
même interpolation que le gradient explicite retiré de b, sinon l'opérateur implicite et le
second membre ne sont pas cohérents.

Code proposé :

```cpp
// Linear interpolation weight of cell c on face f shared with `other`:
// phi_f = w*phi_c + (1-w)*phi_other. Exact for linear fields on non-uniform
// orthogonal meshes, and 0.5 on uniform ones.
inline double face_interpolation_weight(
    const FvGeometry& geometry, std::size_t f, std::size_t c, std::size_t other)
{
    const double dc = (geometry.face_centres[f] - geometry.cell_centres[c]).mag();
    const double dn = (geometry.cell_centres[other] - geometry.face_centres[f]).mag();
    if (!(dc + dn > 0.0) || !std::isfinite(dc + dn))
        throw std::runtime_error("face_interpolation_weight: degenerate face distance");
    return dn / (dc + dn);
}

// gauss_gradient_with_boundary, internal face:
if (own.neighbour(f) >= 0) {
    // A-B1: the opposite cell is neighbour(f) only on the owner side.
    const std::size_t other = owner
        ? static_cast<std::size_t>(own.neighbour(f))
        : static_cast<std::size_t>(own.owner(f));
    const double w = face_interpolation_weight(geometry, f, c, other);
    vf = w * field(c) + (1.0 - w) * field(other);
}

// solve_coupled_momentum_continuity, block G, internal face:
// Same face interpolation as gauss_gradient_with_boundary, so the implicit
// G block equals the explicit gradient removed from b above.
const double wc = face_interpolation_weight(geometry, f, c, ncell);
const double wn = 1.0 - wc;
A.push_back(c,        3*nc + c,     Sf.x * wc);
A.push_back(c,        3*nc + ncell, Sf.x * wn);
A.push_back(nc + c,   3*nc + c,     Sf.y * wc);
A.push_back(nc + c,   3*nc + ncell, Sf.y * wn);
A.push_back(2*nc + c, 3*nc + c,     Sf.z * wc);
A.push_back(2*nc + c, 3*nc + ncell, Sf.z * wn);
```

Test : sur un champ linéaire p = a·x + b·y + c, le gradient de chaque cellule intérieure doit
être exact à 1e-12 près, **y compris sur un maillage étiré**. Le résultat de Poiseuille est
donné en 10.6 : l'ordre passe à 1,87.

```cpp
run_case("gauss_gradient_exact_for_linear_field_on_stretched_mesh", []() {
    auto mesh = make_channel_mesh(8, 16, /*stretch_y=*/1.3);   // helper of channel_mesh.inc
    const auto geometry = FvGeometry(mesh);
    Field<double, Location::CELL> p(mesh.n_cells());
    for (std::size_t c = 0; c < p.size(); ++c)
        p(c) = 2.0 * geometry.cell_centres[c].x - 3.0 * geometry.cell_centres[c].y + 1.0;
    const auto g = gauss_gradient_with_boundary(p, mesh, geometry, linear_extrapolation_bcs(mesh));
    for (std::size_t c : interior_cells(mesh)) {
        EXPECT_NEAR(g(c).x, 2.0, 1e-12);
        EXPECT_NEAR(g(c).y, -3.0, 1e-12);
    }
});
```

(Ce test suppose des helpers `make_channel_mesh(…, stretch)`, `interior_cells` et
`linear_extrapolation_bcs` à écrire. Il est **non compilé** ; le correctif lui-même est
compilé et testé.)

#### 10.4.6 [MAJEUR] `'\\n'` non corrigé (contrairement à ce qu'annonce la PR) — **compilé + testé**

| Lieu (`pr/415:`) | Effet |
|---|---|
| `finite_volume_transport.h:420` et `:437` : `<< '\\n'` | littéral multi-caractères : le **nombre** `23662` est imprimé à la place d'un saut de ligne |
| `steady_incompressible_solver.h:1069` et `:1082` : `"...GMRES\\n"` | affiche `\n` littéralement, sans saut de ligne |
| `preconditioner.h:409` : `<< "\\n"` | idem |

Code proposé : `'\\n'` → `'\n'` et `"\\n"` → `"\n"` aux cinq endroits (compilé + testé). Garde-fou pour la CI (**non compilé**, échappement à valider dans le YAML) :

```bash
# .github/workflows/cfdx-validation.yml, step "lint"
! grep -rnE "'\\\\\\\\n'|\"[^\"]*\\\\\\\\n\"" src/ tests/ || { echo "escaped newline literal"; exit 1; }
```

Autre point : `finite_volume_transport.h:418` et `:435` écrivent **sans condition** sur `cerr`
une ligne `CFDX solver cascade` dès que n ≤ 256. C'est du bruit en production, à passer sous
`DiagnosticsControls` (10.4.9).

#### 10.4.7 [MAJEUR] C6 : inverse de bloc 4×4 mal dé-mise à l'échelle — **non compilé**

Code actuel, `pr/415:src/cfdx/core/linalg/preconditioner.h:152-154` puis `:192-197` :

```cpp
const double scale = 1.0 / row_scale[r];                     // B = S A, S = diag(1/row_scale)
...
// The inverse above is B^{-1} for B = S A. Hence A^{-1} = B^{-1} S.
inv_blocks_[c][4 * r + q] = aug[8 * r + 4 + q] * row_scale[q];   // multiplie par S^{-1} !
```

Le commentaire est juste (A⁻¹ = B⁻¹·S), mais le code multiplie par S⁻¹. Comme l'auto-contrôle
`inverse_residual > 1e-10` (l. 213) détecte l'erreur, `setup()` **renvoie `false` dès qu'une
ligne a une somme ≠ 1**. Le préconditionneur par blocs de cellule est alors abandonné en
silence au profit du repli, ce qui masque le bug. COUPLED ne l'utilise pas ici, donc phase9 n'en
dépend pas.

Code proposé :

```cpp
// The inverse above is B^{-1} for B = S A with S = diag(1/row_scale).
// Hence A^{-1} = B^{-1} S: scale column q by 1/row_scale[q].
for (std::size_t r = 0; r < 4; ++r)
    for (std::size_t q = 0; q < 4; ++q)
        inv_blocks_[c][4 * r + q] = aug[8 * r + 4 + q] / row_scale[q];
```

Test : un bloc diagonal `diag(2, 3, 5, 7)` sur une seule cellule. Aujourd'hui `setup()` doit
renvoyer `false` ; après correction, `setup()` renvoie `true` et `apply(e_k)` vaut `e_k / d_k`.

#### 10.4.8 [MINEUR] Diagonale de Schur écrasée même avec pression imposée — **non compilé**

Code actuel, `pr/415:steady_incompressible_solver.h:1052` :

```cpp
coupled_pressure_schur_diagonal[reference_cell] = 1.0;
```

La ligne de jauge n'existe que si aucune pression n'est imposée. Avec `has_fixed_pressure`,
cette ligne remplace une vraie diagonale de Schur par 1.

```cpp
if (!has_fixed_pressure)   // gauge row only exists when pressure is floating
    coupled_pressure_schur_diagonal[reference_cell] = 1.0;
```

#### 10.4.9 [MINEUR] `getenv` et `cerr` de diagnostic dans le code de production — **non compilé**

Occurrences (`pr/415:`) :
- `gmres_solver.h:214`, `:229` et `:252` (supprimées en 10.4.3) ;
- `preconditioner.h:406` ;
- `steady_incompressible_solver.h:764` et `:953` (`CFDX_DEBUG_COUPLED`) ;
- `steady_incompressible_solver.h:1213` (`CFDX_FREEZE_STATE`) ;
- `steady_incompressible_solver.h:1852` (`CFDX_DEBUG_CELL`).

`CFDX_FREEZE_STATE` (l. 1213-1219, 1412-1416, 1958-2003) sert au diagnostic seulement et n'a
aucun effet sur la solution, mais il **copie le champ de faces à chaque itération**.

Code proposé : un bloc de contrôles explicite, transmis au solveur et testable.

```cpp
struct DiagnosticsControls {
    bool coupled_matrix_summary = false;   // ex CFDX_DEBUG_COUPLED
    bool freeze_state_probe     = false;   // ex CFDX_FREEZE_STATE
    long debug_cell             = -1;      // ex CFDX_DEBUG_CELL (-1 = off)
    std::ostream* sink          = nullptr; // nullptr = silent
};
// SteadyIncompressibleControls { ...; DiagnosticsControls diagnostics; };
// usage : if (auto* os = controls.diagnostics.sink; os && controls.diagnostics.coupled_matrix_summary) *os << ...;
```

Les solveurs linéaires (`solve_gmres`) ne doivent **rien** écrire : ils renvoient déjà
`SolverResult`, qui suffit au diagnostic.

#### 10.4.10 [MINEUR] Messages d'erreur illisibles — **non compilé**

Code actuel, `pr/415:steady_incompressible_solver.h:1306-1307`, `:1361-1362`, `:1562-1563` :
`std::to_string(residual)`. Cela affiche `residual=0.000000` pour 1,3e-6, et masque exactement
la valeur qui a servi au diagnostic de 10.3.

```cpp
inline std::string format_linear_failure(const char* what, const SolverResult& r) {
    std::ostringstream os;
    os << what << " (status=" << static_cast<int>(r.status)
       << ", iterations=" << r.iterations
       << std::scientific << std::setprecision(3)
       << ", residual=" << r.residual
       << ", relative=" << r.residual_relative << ")";
    return os.str();
}
// throw std::runtime_error(format_linear_failure("coupled linear solve failed", coupled_result));
```

#### 10.4.11 [PLAUSIBLE] Rhie-Chow couplé : terme de compensation explicite absent — **non compilé**

Dans la continuité couplée (`pr/415:steady_incompressible_solver.h:~740-810`), la partie
implicite orthogonale −D·(p_N − p_P) est bien présente. En revanche, la compensation explicite
ρ·rfn·(∇p_old moyen)·S_orth, présente dans le Rhie-Chow ségrégué, ne l'est pas : seule la partie
non orthogonale va dans b. Les deux algorithmes ne résolvent donc pas exactement le même
problème discret. Cet écart n'a **pas été reproduit** : COUPLED retrouve Couette exactement, et
l'écart avec SIMPLE (`max_abs_dU = 2,73e-7`) reste du même ordre que la tolérance.

À vérifier avec un test d'invariance COUPLED/SIMPLE sur Poiseuille en maillage **étiré** avant
de corriger. Si l'écart est confirmé, ajouter à b, pour chaque face interne :

```cpp
// explicit RC compensation, consistent with the segregated interpolation
const Vec3 gp_f = wc * grad_p_old(c) + wn * grad_p_old(ncell);
b(3*nc + c) += rho * rfn_f * dot(gp_f, S_orth_f);   // sign to be aligned with the D block
```

#### 10.4.12 Boucle principale, branche COUPLED — remarque

`pr/415:steady_incompressible_solver.h:1292-1335` : `relax_momentum_equation` n'est pas
appliquée à COUPLED. Les champs sont relaxés explicitement **après** la résolution, par
`alpha_u` et `alpha_p`. Ce choix est défendable, mais un solveur couplé se relaxe
habituellement en implicite (diagonale/α), ou par un pas de pseudo-temps. En explicite, le
point fixe converge lentement, comme le montre le Poiseuille COUPLED qui n'atteint pas son
critère en 3000 itérations (10.6). À documenter, ou à passer en relaxation implicite.

#### 10.4.13 Optimisation : préconditionneur de Schur — **non compilé**

`CoupledBlockSchurPreconditioner` (`pr/415:preconditioner.h:259+`) approche S = C − D·diag(M)⁻¹·G
par sa **diagonale**. Or S est un Laplacien de pression (κ ~ h⁻²) : Jacobi sur S devient
inefficace quand on raffine. Le restart l'a déjà montré (10.3).

Proposition : assembler S creux (même motif que le Laplacien de pression) et appliquer S⁻¹ par
quelques balayages de CG préconditionné par Jacobi, ou par un V-cycle AMG quand il sera
disponible. `set_pressure_schur_diagonal` reste un simple repli.

```cpp
// apply(): y_u = diag(M)^{-1} r_u ; r_p' = r_p - D y_u ;
//          y_p ≈ S^{-1} r_p' by k PCG sweeps (k = 5..10, fixed -> linear preconditioner, GMRES ok) ;
//          y_u -= diag(M)^{-1} G y_p
```

**Attention** : avec un nombre de balayages fixe, le préconditionneur reste **linéaire**. Il est
alors compatible avec GMRES classique. Si la tolérance interne est adaptative, il faut passer à
FGMRES.

### 10.5 Bugs connus de master : que fait la PR ?

| Réf. audit | Sujet | Dans #415 |
|---|---|---|
| A-B1 | gradient de Gauss côté voisin | **toujours présent** (10.4.5) |
| A-B2 | relaxation de l'équation de quantité de mouvement | **corrigé** (`relax_momentum_equation`, l. 506-552) |
| A-B3 | transport `diag += F` | **corrigé** |
| A-B4 | jauge de pression | **corrigé** (décalage uniforme, l. 1569-1575) |
| C2 | statut du solveur de pression ignoré | **corrigé** (`rp.status` vérifié, exception l. 1557) |
| C1, B1–B6 | turbulence, énergie, rayonnement, SA, thermophysique, `krylov_reductions` | **toujours présents** : `git diff master pr/415` sur ces fichiers est vide |
| `'\\n'` | littéraux échappés | **toujours présent** (10.4.6) |

### 10.6 Harnais Couette / Poiseuille (PR contre correctifs)

Harnais : `/tmp/cfdx_meta/scratch_415/h415.cpp`.
- Canal nx = 8, ν = 0,1, α_u = 0,7, α_p = 0,3, au plus 3000 itérations, tolérance 1e-8.
- Couette : u_top = 1.
- Poiseuille : p_in = 0,8, p_out = 0, solution exacte u = 4y(1−y), donc u_max = 1.
- Couette sur ny = 16.

| Cas | PR (`f829e32f`) | Correctifs (`fixc`) |
|---|---|---|
| Couette SIMPLE, bornée | conv., 409 it, Linf 2,76e-7 | conv., 392 it, Linf 2,56e-7 |
| Couette SIMPLE, non bornée | conv., 409 it, Linf 2,75e-7 | conv., 392 it, Linf 2,55e-7 |
| Poiseuille SIMPLE ny = 16 | non conv. (3000), **u_max 0,593**, L2 0,303 | non conv. (3000), u_max 1,0127, L2 7,27e-3, Linf 1,66e-2 |
| Poiseuille SIMPLE ny = 32 | non conv., u_max 0,570, L2 0,316, **ordre −0,06** | non conv., u_max 1,0034, L2 1,99e-3, **ordre 1,87** |
| Couette COUPLED, bornée | conv., 16 it, Linf 4,2e-9 | conv., 16 it, Linf 4,26e-9 |
| Couette COUPLED, non bornée | conv., 16 it, Linf 4,2e-9 | conv., 16 it, Linf 4,24e-9 |
| Poiseuille COUPLED ny = 16 | non conv. (3000), u_max 0,99915, L2 3,08e-3 | identique |
| Poiseuille COUPLED ny = 32 | **exception `MAX_ITER` linéaire**, u_max 0,91 | non conv. (3000), u_max 0,99993, L2 8,97e-4, **ordre 1,78** |

Lecture du tableau :
- A-B1 explique à lui seul le Poiseuille SIMPLE faux de la PR.
- Le correctif restart explique l'exception COUPLED en ny = 32.
- Avec les correctifs, l'ordre 2 est retrouvé (1,87 / 1,78).
- **PLAUSIBLE** : la boucle non linéaire ne satisfait pas son critère en 3000 itérations sur
  Poiseuille. Le critère conjoint (`steady_incompressible_solver.h:1897-1902`) exige que
  six grandeurs soient toutes ≤ 1e-8, dont `velocity_change_inf` et `pressure_change_inf`, qui
  sont **non normalisées** et amorties par α. Ces deux-là stagnent autour de 1e-8 alors que
  l'erreur de discrétisation est à 1e-3. Il faut soit normaliser les incréments par
  α·|U_ref|, soit ne garder que les résidus normalisés et la continuité (défaut m4 de l'audit).

### 10.7 Verdict

**NON MERGEABLE en l'état.**
- La CI est rouge sur les deux workflows.
- Le commit qui devait réparer COUPLED (`f829e32f`, restart 512) n'a aucun effet.
- Un test verrouille un comportement faux (oracle `DIVERGED` sur un système régulier).
- A-B1 laisse Poiseuille SIMPLE faux d'un facteur 1,7.

Avec les correctifs 10.4.1 à 10.4.6, déjà compilés et testés, la PR devient verte localement
(92/93 ; le dernier échec est le D4 préexistant) et phase9 passe.

### 10.8 Changements requis, dans l'ordre

1. **Restart GMRES** : supprimer l'écrêtage `[10,40]` (10.4.1) **et** remettre la politique
   adaptative dans le bon sens (10.4.2). Mettre à jour `test_performance_runtime.cpp:37`.
2. **Breakdown d'Arnoldi** : test de plancher avant normalisation, rendu du résidu vrai, suppression
   des `getenv`/`cerr` dans `gmres_solver.h` (10.4.3).
3. **Oracles de `test_gmres_solver.cpp`** : système régulier → `CONVERGED` ; ajouter le vrai cas
   singulier et le cas `diag(1,0)` → x = (1,1), résidu 1 (10.4.4).
4. **A-B1** : interpolation pondérée côté propriétaire *et* voisin, et même pondération dans le
   bloc G couplé (10.4.5). Ajouter le test Poiseuille d'ordre ≥ 1,8.
5. **`'\\n'`** aux cinq endroits, plus le garde-fou grep en CI, et le silence du `cerr` « solver cascade » (10.4.6).
6. **C6** : `/ row_scale[q]`, avec un test de bloc diagonal (10.4.7).
7. Garder la diagonale de Schur quand la pression est imposée (10.4.8). Formater les messages
   d'échec en notation scientifique (10.4.10).
8. Remplacer les `getenv` par `DiagnosticsControls`, et retirer ou encadrer `CFDX_FREEZE_STATE` (10.4.9).
9. Rebaser sur `master` après le merge de #410 + correctifs (§7.2), en gardant le test `extra_diagonal` de `b01f0498` à côté du `EXPECT_THROW`.
10. Suivi dans des PR séparées : Rhie-Chow couplé (10.4.11, à confirmer d'abord), critère de
    convergence normalisé (m4), préconditionneur de Schur non diagonal (10.4.13), et les bugs de
    master B1–B6 et C1.

### 10.9 Checklist de merge

- [ ] `ctest` vert en local, en dehors de D4 (`test_performance_runtime`, `mpi.h`), qui doit être explicitement marqué `DISABLED`/`SKIP` si MPI est absent
- [ ] CI GitHub verte sur **les deux** workflows, sur la tête finale
- [ ] `test_phase9_acceptance` : `successful=6 failed=0`, COUPLED sans reprise BiCGStab
- [ ] Test « restart non écrêté » (10.4.1) présent et vert
- [ ] Oracles GMRES corrigés (11/11)
- [ ] Poiseuille SIMPLE et COUPLED : u_max à 1 % près, ordre observé ≥ 1,8 entre ny = 16 et 32
- [ ] `grep` du garde-fou `'\\n'` vide
- [ ] Aucun `std::getenv` ni `std::cerr` sous `src/cfdx/core/linalg/`
- [ ] Description de la PR corrigée : retirer « bug `\n` corrigé » tant que ce n'est pas vrai ; citer la cause racine (restart) et non le « restart 512 »
- [ ] #415 rebasée sur `master` après le merge de #410 + correctifs, test `extra_diagonal` conservé ; #409 et #411 fermées comme contenues
- [ ] Ref distante `pr/415` à jour (elle pointe encore `1443e313`)

**Artefacts de la revue** (dans `/tmp/cfdx_meta/scratch_415/`) :
- `p9dbg.out`, `p9dbg.err`, `p9dbg.clean` : phase9 de la PR, instrumentée ;
- `fixc/_b/p9.out` et `p9b.out` : phase9 corrigée ;
- `fixc/_b/ctest.log` ;
- `h_pr.out` et `h_fix.out` : harnais 10.6 ;
- `gmres.diff` et `other.diff` : diff complet des correctifs ;
- `ci_build_f829.log` et `ci_val_f829.log` : logs de la CI.

---

## Annexe A — Reproduction et correctifs

**Correctifs conservés à côté du rapport** (`~/audits/cfdx_dev/correctifs/`, copie des fichiers de
`/tmp/cfdx_meta/`, qui est volatil). Dans le corps du rapport, un chemin `/tmp/cfdx_meta/X`
correspond à `correctifs/<section>/X`.

| Dossier | Contenu | Sections |
|---|---|---|
| `s3a/` | Arborescence `cfdx/` corrigée (discrétisation, couplage P-V, géométrie) | §3.1, §3.2 |
| `s3b/` | `s3b.diff` (turbulence, énergie, rayonnement, temps), `b17.diff` (S2S de #398), `overlay_src/` (sources corrigées complètes), `tests_updated/` (oracles mis à jour) | §3.3 à §3.6 |
| `s4_s8/` | `s4_s8.patch` : solveurs linéaires, I/O, parallélisme et optimisations | §4, §8 |
| `s567/` | `d4.diff`, `d9.diff`, `energy.diff`, `format_residual.diff`, `dbg.diff`, `rk2.diff`, `validation_report.diff`, `case_io.diff`, `parallel_hdf5.diff`, `prod_e2e.diff` ; `d6/` (3 tests remplaçants) ; `workflows/` (`ci.yml`, `cfdx-validation.yml` et leurs diffs) ; `tpl/` ; `run_validation.py` (remplaçant complet) | §2, §5, §6, §7 |
| `pr415/` | `gmres.diff`, `other.diff` (correctifs complets de #415), `ctest_fixc.log`, logs CI `ci_build_f829.log` et `ci_val_f829.log` | §10 |
| `../AUDIT_v1.md` | Version 1 de ce rapport, pour mémoire | — |

**Environnements de reproduction** (volatils, sous `/tmp`) :
- Clone : `/tmp/cfdx_dev` (HEAD `27bae126`, URL de push désactivée, aucun fichier suivi modifié), branches des PRs : `pr/<N>`.
- Build : `/tmp/cfdx_dev/_audit_build` (MPI OFF, HDF5 1.14.4-3 compilé dans `/tmp/hdf5_inst`). Logs : `_audit_build2.log`, `_audit_ctest.log`.
- Copies de travail : `/tmp/cfdx_scratch`, `/tmp/cfdx_pr410` (#410 + correctifs, 96/96), `/tmp/cfdx_pytest_scratch`, `/tmp/cfdx_pr415n` (#415 tête `f829e32f`), `/tmp/cfdx_meta/scratch_415/fixc` (#415 corrigée).
- Programmes de reproduction : `/tmp/cfdx_meta/scratch_a/` (couplage P-V), `scratch_b/` (physique, facteurs de forme), `scratch_c/` (Krylov, NaN, VTU), `scratch_415/h415.cpp` (harnais Couette/Poiseuille).
- CI : `/tmp/cfdx_meta/runs.json`, `fail_36036549020.log`, `fail_36036548885.log`.
- Issues et PRs extraites : `/tmp/cfdx_meta/issue_*.md`, `pr_*.md`, `pr_*_files.txt`.

**Appliquer un correctif** (par le propriétaire, sur une branche de travail). Les diffs sont des
`diff -u` entre copies de travail : leurs préfixes de chemin varient (`/tmp/cfdx_dev/src/…`,
`d9_orig/…`, `orig_src/…`). Il faut donc choisir le `-p` fichier par fichier, et toujours
commencer par `--dry-run`. Exemple vérifié sur `master@27bae126` (dry-run seulement) :
```bash
git switch -c fix/audit-s4 master
patch -p3 --dry-run < ~/audits/cfdx_dev/correctifs/s4_s8/s4_s8.patch   # tous les hunks passent
patch -p3           < ~/audits/cfdx_dev/correctifs/s4_s8/s4_s8.patch
cmake --build build -j && ctest --test-dir build --output-on-failure
```
Les diffs faits contre #410 (`validation_report.diff`, `workflows/`) ou contre #415 (`pr415/`)
s'appliquent sur ces branches, pas sur `master`. `s3a/` et `s3b/overlay_src/` contiennent des
**fichiers complets** : il suffit de les comparer ou de les copier.

## Annexe B — Limites de l'audit

- **Statut par bloc.** Chaque code proposé porte son statut en titre : *compilé + testé*, *compilé*
  ou *non compilé*. Tous les blocs « Code actuel » sont copiés du source réel à `27bae126`
  (ou `pr/<N>` quand c'est indiqué) ; les identifiants des codes proposés sont ceux du dépôt.
- **Non exécuté** :
  - CUDA : pas de `nvcc`, donc C15, C19 et l'item 8-7 sont *non compilés* ;
  - HDF5 parallèle (MPI-IO) : non construit ;
  - MPI : seuls des tests à np = 2 et 4 ont tourné, sans mesure de scalabilité (item 8-5 estimé) ;
  - `actionlint` n'a pas été lancé sur les workflows proposés : les valider par un run sur une branche.
- **Gains de performance** : les items *mesurés* du §8 viennent de bancs en copie de travail
  (même compilateur, même options). Ils ne remplacent pas un profil sur un cas de 1M de cellules.
- **Validation Python** : 6 échecs pytest restent après correctifs (§5.4). Ils ne sont pas analysés jusqu'au bout.
- Les numéros de ligne correspondent à `27bae126` (et `f829e32f` pour #415). Ils bougeront avec les merges.
- `/tmp` est volatil : seuls `correctifs/` et ce rapport sont conservés durablement.
