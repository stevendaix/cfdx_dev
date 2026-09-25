
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
+    validate_implicit_source(source_i
