# Architecture CI CFDX — rationalisation du 8 novembre 2026

## 1. Trois questions, trois niveaux

| Question | Responsable | Preuve | Cadence |
|---|---|---|---|
| Le code compile-t-il et reste-t-il sûr ? | `cfdx-ci.yml` + `precommit.yml` | Release, DebugSanitizers, CTest court, Python, Ruff, suites CI | À chaque PR |
| Les contrats N1–N15 sont-ils respectés ? | `cfdx-numerical.yml` + `cfdx-maturity.yml` | Matrice N1–N15 et audits JSON | À chaque PR (sauf qualification N7 Re400 longue) |
| Les campagnes physiques, les grands maillages et la robustesse globale passent-ils ? | `cfdx-validation.yml` | CTest MPI/HDF5 complet, N9, Ghia, NACA, SST/thermique, distance, JUnit, rapport | Lundi, manuel, ou label `validation-total` |

L'étiquette **long** détermine les tests déplacés vers la validation complète ; elle ne
modifie pas leurs critères de convergence, leurs tolérances, leurs budgets internes,
leurs commandes ou leurs références physiques.

## 2. Graphe des validations

```text
 PR (sans validation-total)
 ├── Pre-commit: static policy + Python + change-aware smoke
 ├── CFDX CI: Build/Test Release ── CTest court + Python + N12
 │            Build/Test DebugSanitizers ── CTest non-slow
 ├── CFDX numerical maturity: N1 ... N15 (indépendants)
 └── Maturity records / audit / workflows ciblés

 PR avec label validation-total / schedule lundi / dispatch
 └── CFDX total validation
     └── Prepare: UNE compilation Release MPI/HDF5 + archive
         ├── Fast CTest (sur cette configuration MPI/HDF5) en parallèle
         ├── Long CTest shard 0 ──┐
         ├── Long CTest shard 1 ──┼── Rapport + vérification croisée des preuves
         ├── Long CTest shard 2 ──┤
         ├── Long CTest shard 3 ──┘
         ├── Wall-distance standards
         ├── Wall-distance N128 (9 méthodes, jobs séparés)
         ├── Thermique/rayonnement
         └── Application E2E
                 └── Total validation gate : PASS uniquement si tous PASS
```

Les jobs `long` restent indépendants (`fail-fast: false`) : un échec ne cache pas
les résultats des autres. `RUN_SERIAL` conserve son sens **au sein d'un runner**,
mais n'empêche pas les campagnes totalement indépendantes de progresser sur
des runners différents.

## 3. Contrat de couverture : aucune qualification oubliée

La liste CTest réelle (`ctest --show-only=json-v1`) est la source de vérité.
`scripts/ci/ctest_shard.py` :

1. Filtre exactement le label `long`.
2. Conserve ensemble les dépendances CTest déclarées.
3. Rejette les tests désactivés, dépendances incohérentes et partitions vides.
4. Affecte chaque test à exactement **une** partition (plan déterministe).
5. Génère la sélection CTest par regex échappée et ancrée (compatible avec CMake 3.28).
6. Utilise le temps observé de `CTestCostData.txt` lorsqu'il est fourni, sinon le `TIMEOUT` comme **poids de planification**, jamais comme durée mesurée.

Le mécanisme natif des **fixtures CTest** rajoute les tests de génération de
maillage requis lorsqu'une qualification NACA en a besoin. Ils peuvent donc
être réexécutés sur plusieurs runners indépendants ; cela **ne constitue pas un
doublon de test physique**. Cette exception est affichée dans la vérification.

`scripts/ci/ctest_timings.py` vérifie le XML JUnit de chaque shard : chaque test
attendu doit avoir été exécuté exactement une fois, sans failure/error/skip.
`scripts/ci/verify_shards.py` vérifie ensuite que les **quatre artefacts**
contiennent exactement le même plan, et que la couverture globale est complète.
Tout manque de preuve devient un **FAIL**, pas un PASS par défaut.

### Campagnes lourdes reclassées

- `test_ghia_cavity_order` (anciennement inclus dans le CTest PR : `TIMEOUT=900`).
- `test_n7_adaptive_relaxation_re400` (`TIMEOUT=1800`).
- `test_polyhedral_poisson_convergence` (`TIMEOUT=1800`).
- `test_n9_poiseuille_body_force_production` (`TIMEOUT=900`).
- `test_n9_negative_diagnostics` (`TIMEOUT=600`).
- `test_wall_distance_transport_stress65` (`TIMEOUT=600`).

Ce reclassement garantit que les workloads à gros budget ne se glissent plus
silencieusement dans les suites PR dites rapides. Un test toujours critique mais
court (par exemple SST/3D, Rhie–Chow MMS rapide, NACA `--quick`, transport 33³)
reste surveillé par les workflows PR selon ses labels existants. **Une PR sans
label n'exécute pas les six campagnes reclassées** : politique de promotion à
définir selon les exigences de la branche protégée.

### Observations locales sur le ZIP fourni

Après configuration CMake Release avec `CFDX_BUILD_TESTS=ON` et
`CFDX_ENABLE_LONG_VALIDATION=ON` : **374 tests**, dont **33 tests `long`**
(27 auparavant). La sélection des quatre regex est complète : **6/7/10/10**
tests longs planifiés, plus le test fixture de génération des maillages qui
est automatiquement inclus par CTest quand nécessaire. Les poids 8040/7920/
7980/8070 proviennent des `TIMEOUT` ou de valeurs de repli, et **ne sont pas
des secondes de calcul observées**.

## 4. Suppression des redondances N9

Avant, une PR portant `validation-total` pouvait déclencher **simultanément** :

- La matrice N9 dans le total-validation (`long`).
- `n9-campaign.yml`, nouvelle compilation et campagne N9.
- `n9-physical-validation.yml`, cinq compilations supplémentaires et cinq
  exécutions par modèle de Schur.

Après : le niveau 3 fait autorité. Les deux autres workflows sont conservés
pour **`workflow_dispatch` uniquement** : ils deviennent de vrais outils de
diagnostic, non des vérifications répétées sur chaque qualification totale.

## 5. Savoir où regarder en cas d'échec

| Symptôme | Où chercher | Artefact / diagnostic |
|---|---|---|
| Erreur compilation | `CFDX CI / Build and test (...)` ou `Prepare validation build` | `build.log`, cache compiler, `compile_commands.json` |
| Echec CTest rapide | `CFDX CI` ou `Fast regression` | `ctest-Release.xml`, `ctest-fast-results.xml`, `*-timings.json` |
| Echec Ghia/N9/NACA | `Long validation shard N/4` | `ctest-long-N.xml`, `long-validation-shard-N.log`, manifest de répartition |
| Test long manquant | `validation-report` ou `Final shard gate` | `long-shards-coverage.json` et message `FAIL-CLOSED` |
| Distance à la paroi | `wall-distance` / `wall-distance-fine` | CSV par méthode et N128 |
| Régression Python/PyFluent | `CFDX CI / Release` | Trace de `pytest` (la dépendance `ansys-units` est épinglée dans cette CI) |
| Résultat physique faux malgré convergence numérique | Campagne + rapport | Valeurs physiques et références indépendantes ; **pas** le seul statut CTest |

Un test au verdict PASS prouve seulement le contrat effectivement codé ; un
contrat de référence ou un prérequis ne prouve **pas** qu'une campagne physique
CFDX ait été qualifiée. Le registre de qualification doit continuer de faire
la distinction.

## 6. Évaluer le gain, sans confondre timeout et durée

Pour chaque nouvelle exécution GitHub Actions, relever :

- Durée réelle `prepare`, `fast-tests`, `long-tests[0..3]`, rapport, et `final-gate`.
- Durée totale en temps mural (du lancement au verdict final).
- Les 20 tests CTest les plus lents (`*-timings.json`).
- La part compilations/installations/cache versus calcul physique.
- Les régressions de durée entre deux HEADs à environnement équivalent.

**On ne peut pas conclure que les 2 heures seront devenues 20 minutes sans
exécution GitHub réelle.** Un test ayant `TIMEOUT=5400` peut être rapide ou
consommer réellement 90 minutes. Le parallélisme élimine l'attente d'un test
indépendant derrière un autre, pas le coût incompressible d'une simulation.
Aucun timeout CTest, aucun seuil de conservation ou de convergence physique
n'a été abaissé.

## 7. Lancer localement

```bash
cmake -S . -B build-validation -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
  -DCFDX_ENABLE_LONG_VALIDATION=ON
cmake --build build-validation --parallel 4

# Contrôles CI ultracourts, sans compilation CFD
python3 -m pip install pyyaml
python3 scripts/ci/ci_policy.py
python3 -m unittest discover -s tests/python -p 'test_ci_*.py' -v
python3 -m unittest discover -s tests/python -p 'test_ctest_*.py' -v

# Préparer les 4 shards / vue CTest sans exécuter une campagne
python3 scripts/ci/ctest_shard.py --build-dir build-validation \
  --label long --lanes 4 --output-dir build-validation/ci-shards
ctest --test-dir build-validation -LE long --parallel 4 --output-on-failure

# Exécution d'un seul shard *si un build exécutable existe*
ctest --test-dir build-validation \
  -R "$(cat build-validation/ci-shards/long-shard-0.regex)" \
  --parallel 3 --output-on-failure
```

**Attention :** `ctest -N`/`ctest --show-only=json-v1` valide le manifeste et la
sélection, **pas la compilation ni le succès des solveurs**.

## 8. Politique de fusion et périmètre

- Checks requis recommandés sur les PR : `Pre-commit`, Release et
  DebugSanitizers ; choisir explicitement les contrôles N-level à rendre
  bloquants suivant la politique de CFDX.
- `validation-total` est **opt-in** lors du travail courant. Avant une fusion
  critique ou une promotion numérique/physique, lancer la validation complète
  sur le **HEAD exact**, conserver la preuve `Total validation gate = PASS` et
  recommencer si le HEAD change.
- Si l'équipe veut **interdire techniquement** toute fusion sans qualification,
  il faut mettre en place une règle de fusion explicite et un check garanti sur
  toutes les PR : **un job facultatif sauté n'est pas une preuve PASS**.
- Ne pas rendre obligatoires les workflows filtrés par chemins ou manuels
  (TH, N9 diagnostics, CUDA) ; ils peuvent rester `pending` ou `skipped`.

### Périmètre de la modification

Change uniquement l'orchestration CI, la classification CTest des campagnes
coûteuses, les contrôles de non-omission et la documentation. N'introduit
**aucun changement** au solveur, aux résidus, aux tolérances, à la physique,
aux géométries ou aux critères de qualification. Les nouveaux tests Python
contrôlent le protocole, pas l'écoulement. Il faut toujours faire exécuter
la CI GitHub complète et les cas physiques pour valider ce changement sur le
HEAD de la PR.