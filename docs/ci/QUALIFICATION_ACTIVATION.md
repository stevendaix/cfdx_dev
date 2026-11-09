# CFDX - Activer une qualification en quelques clics

Pour **suivre tous les cas en parallèle**, consulter aussi [le portefeuille transversal](QUALIFICATION_PORTFOLIO.md). Le dashboard de chaque PR comprend les 25 cas (y compris ceux qui n'ont pas été exécutés), indépendamment du profil ciblé.

## 1. Depuis GitHub Actions (un seul menu)

1. Ouvrir **Actions > CFDX qualification - activate a profile**.
2. Cliquer **Run workflow** et choisir **la branche / le HEAD a qualifier**.
3. Choisir un **profil** (`naca`, `momentum`, `sst`, `thermal`, `wall-distance`,
   `linalg`, `n9`, `ghia`, `smoke` ou `total`).
4. Choisir la **portee** : `quick` ou `extended`.
5. Lancer. Les fichiers `plan.json`, `tests.txt`, `ctest.xml` et `timings.json`
   sont disponibles dans les **artifacts** du run cible.

`total` appelle **la meme** `cfdx-validation.yml` que `validation-total` et le
lancement hebdomadaire. Pour `total`, le champ `scope` est ignore : la validation
est toujours **complete**. Une campagne de qualification ne peut etre promue
qu'apres lecture des resultats sur le HEAD exact.

## 2. Depuis une pull request (une etiquette)

Pour creer d'un coup les etiquettes manquantes (une seule fois, avec `gh`
authentifie sur le depot) :

```bash
bash scripts/ci/setup_qualification_labels.sh
```

Le script ne supprime pas les etiquettes existantes et ne les modifie pas.
Ajouter ensuite **une seule** etiquette au choix :

| Etiquette | Domaine | Mode automatique |
|---|---|---|
| `qualification:smoke` | Sentinelles numeriques | `extended` (= court pour smoke) |
| `qualification:naca` | NACA 2D/3D, paroi, efforts/GCI | `extended` |
| `qualification:momentum` | Pression, Rhie-Chow, momentum | `extended` |
| `qualification:sst` | k-omega SST et turbulence | `extended` |
| `qualification:thermal` | Energie, CHT, P1/S2S, convection | `extended` |
| `qualification:wall-distance` | Distance a la paroi CTest | `extended` |
| `qualification:linalg` | AMG, Krylov, anisotropie, Schur | `extended` |
| `qualification:n9` | Matrices N9 et Coupled | `extended` |
| `qualification:ghia` | Cavite Ghia et raffinement | `extended` |
| `validation-total` | Qualification physique complete | Workflow canonique `cfdx-validation.yml` |

L'etiquette active reste en place lors d'un nouveau commit et relance la
qualification sur le nouveau HEAD (`synchronize`). Retirer l'ancienne etiquette
avant d'en activer une autre. Deux etiquettes `qualification:*` simultanees sont
un **FAIL de configuration**, jamais une selection silencieuse. Si
`validation-total` est present, il est prioritaire : le workflow cible se saute
pour eviter le doublon.

**Ne pas exiger le workflow cible comme check obligatoire sur toutes les PR** :
une PR sans etiquette n'a pas execute la qualification. Les verifications
systématiques sont celles de `cfdx-ci.yml`/`precommit.yml`. Pour imposer la
validation totale, utiliser la politique de fusion decrite dans
`CI_ARCHITECTURE.md`.

## 3. Profils courts vs approfondis

- **`quick`** : liste explicite de sentinelles CTest, interdiction des tests
  etiquetes `long`; compilation uniquement des executables C++ requis.
- **`extended`** : tests du domaine reconnus par nom, y compris ceux etiquetes
  `long`; compilation du build de validation complet (certains tests Python
  lancent des outils binaires indirectement).
- **`total`** : *pas seulement CTest*. Reutilise la validation complete, les
  4 partitions longues, les 9 methodes de distance N128, les tests thermiques,
  E2E, le rapport physique et la verification croisee des preuves.

Exemples de couverture **dans le manifeste Release de reference** (374 CTest,
33 `long`, MPI non active sur la machine de controle) :

| Profil | `quick` (CTest, prerequisites inclus) | `extended` (CTest, prerequisites inclus) |
|---|---:|---:|
| smoke | 4 | 4 |
| naca | 8 | 32 |
| momentum | 5 | 39 |
| sst | 5 | 21 |
| thermal | 4 | 38 |
| wall-distance | 3 | 17 |
| linalg | 4 | 48 |
| n9 | 3 | 14 |
| ghia | 2 | 7 |

Le nombre peut changer si les tests/CMake changent ; **le manifeste de chaque
run est souverain**. Si un test ancre obligatoire manque, le run echoue avant
la compilation. Les dependances CTest et fixtures de maillage sont ajoutees
au plan ; aucune selection vide ou desactivee n'est acceptee.

## 4. Execution en local

Une seule commande par phase, identique a la logique GitHub :

```bash
cmake -S . -B build-qualification -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_LONG_VALIDATION=ON \
  -DCFDX_BUILD_PYTHON=OFF

# Decouvrir les noms et les etiquettes disponibles
python3 scripts/ci/qualification.py list

# Voir les tests AVANT de payer leur compilation
python3 scripts/ci/qualification.py plan \
  --build-dir build-qualification --profile naca --scope quick

# Compiler seulement les cibles necessaires en mode quick
python3 scripts/ci/qualification.py build \
  --build-dir build-qualification --profile naca --scope quick --parallel 3

# Executer, puis verifier le JUnit et la couverture requise
python3 scripts/ci/qualification.py run \
  --build-dir build-qualification --profile naca --scope quick --parallel 3

# Passer a la campagne CTest NACA approfondie
python3 scripts/ci/qualification.py plan \
  --build-dir build-qualification --profile naca --scope extended
python3 scripts/ci/qualification.py build \
  --build-dir build-qualification --profile naca --scope extended --parallel 3
python3 scripts/ci/qualification.py run \
  --build-dir build-qualification --profile naca --scope extended --parallel 3
```

L'installation des dependances Python, MPI et HDF5 reste a effectuer selon la
configuration locale, comme dans `cfdx-validation.yml`. Le `total` local ne se
reduit pas a `ctest` : utiliser le workflow complet pour les preuves finales.

## 5. Resultats et limites (pas de faux verts)

- `plan.json` consigne les tests primaires, prerequis, tests longs et le perimetre.
- `tests.regex` contient des noms echappes et ancres (`^...$`).
- `ctest.xml` et `timings.json` recensent les tests reels ; un test attendu
  manque, passe en `skipped`, echoue ou est duplique => **FAIL**.
- La selection d'un test de qualification n'est **pas** sa validation physique.
  Par exemple, un test du *contrat de certificat GCI* peut passer pendant que
  les efforts NACA N256 ne sont toujours pas qualifies.
- Aucun seuil de convergence, tolerance, comparaison physique ou timeout
  n'est modifie par ce selecteur.
- L'etiquette de PR et le menu `Run workflow` ne donnent pas un statut PASS ;
  seuls les resultats et les preuves du HEAD effectivement teste le peuvent.

## 6. Ajouter un nouveau domaine

Le fichier `scripts/ci/qualification_profiles.json` est la source de verite :

1. Ajouter un objet avec `id`, `label=qualification:<id>`, `description`.
2. Definir `quick` (noms CTest exacts), `quick_targets` (cibles CMake necessaires),
   `extended_regex` (regex ancree) et `required_extended` (tests sentinelles longs).
3. Ajouter l'option `id` dans `.github/workflows/cfdx-qualification.yml`.
4. Lancer `python3 scripts/ci/ci_policy.py` et
   `python3 -m unittest discover -s tests/python -p 'test_ci_*.py'`.
5. Configurer CMake et verifier le plan reel pour `quick` et `extended`.

et reference independante correspondants.
