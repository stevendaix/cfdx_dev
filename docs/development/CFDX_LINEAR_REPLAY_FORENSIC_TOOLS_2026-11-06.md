# CFDX — captures CSR, rejeu de préconditionneurs et audit de convergence (2026-11-06)

**Statut : implémenté et testé localement, qualification CFD globale non terminée.**
**Base :** archive `CFDX_GMRES_TVD_CODE_COMPLET_2026-11-06.zip` (commit Git local de référence `5677e7f`).
Aucun push GitHub. Aucun choix automatique de nouveau préconditionneur ni changement de tolérance.

## Pourquoi ces outils

Les erreurs des AMG anisotropes sur les matrices de pression N128/N256 ne sont pas correctement caractérisées par des matrices fabriquées. Il faut capturer l'opérateur **réel après jauge**, le second membre et la solution initiale, puis comparer les méthodes sur un problème strictement identique. Le statut de solveur et le résidu vrai original doivent être certifiés indépendamment.

## Composants

| Composant | Fichier | Contrat |
| --- | --- | --- |
| Format CSR reproductible avec empreinte de contrôle | `src/cfdx/core/linalg/linear_system_capture.h` | FP64 bit-exact, empreinte FNV-1a anti-corruption, rejet des entrées non finies/incohérentes, écriture `.tmp` puis renommage ; série uniquement |
| Capture générique | `src/cfdx/core/linalg/linear_solver_dispatch.h` | `LinearSolverRequest.linear_capture_path` opt-in ; matrice originale avant scaling/AMG |
| Capture pression CFD | `src/cfdx/physics/steady_incompressible_solver.h` | `pressure_replay_prefix`, `pressure_replay_iterations` ; pression après élimination de jauge ; couvre aussi CG réutilisable |
| Commande NACA | `tests/validation/test_naca0012_geometry.cpp` | Variables d'environnement de capture pour les itérations choisies |
| Rejeu des préconditionneurs | `tools/diagnostics/cfdx_linear_replay.cpp` | CSV, temps total setup+solve, FGMRES et résidu vrai, échec explicite |
| Campagne batch | `tools/diagnostics/compare_linear_captures.py` | Comparaison de plusieurs captures, CSV avec `failure_reason`, arrêt en échec quand candidat non convergé |
| Auditeur de convergence | `tools/diagnostics/qualify_convergence_history.py` | Preuve à partir du JSON CFDX, résidus par composante, continuité, bilan masse, maintien de l'ordre 2 |

### Rejouer une matrice originale

Depuis la racine du projet, après compilation de la cible optionnelle :

```bash
cmake --build build --target cfdx_linear_replay -j 2
build/cfdx_linear_replay capture.cfdxlin native directional energy
```

Code retour 0 seulement si **toutes** les méthodes sélectionnées ont le statut CONVERGED et satisfont `||b-Ax||_2 / ||b||_2 <= tolerance` (ou la norme du résidu initial si `b=0`). Les valeurs proviennent du système original, non des résidus récursifs ou scaled. Le CSV contient le temps **setup+solve** mais n'est pas une mesure MPI.

```bash
python3 tools/diagnostics/compare_linear_captures.py ./captures \
  --binary build/cfdx_linear_replay --method native --method directional \
  --csv comparison.csv
```

Le programme en lot conserve les lignes des candidats défaillants et le diagnostic `failure_reason` de la hiérarchie AMG, puis quitte avec le code 1. Ne supprimer ni les erreurs, ni les méthodes en échec pour fabriquer une qualification verte.

### Capturer la pression pendant NACA

Créer au préalable le dossier de destination ; définir les deux variables :

```bash
export CFDX_NACA_PRESSURE_REPLAY_PREFIX=/tmp/cfdx/naca256_pressure
export CFDX_NACA_PRESSURE_REPLAY_ITERATIONS=1,20,21,240
# Lancer la cible NACA existante avec le maillage et le budget appropriés.
# Les fichiers obtenus sont /tmp/cfdx/naca256_pressure_i<iter>_c<corr>.cfdxlin
```

La capture crée seulement les fichiers des itérations réellement atteintes. La validation du sélecteur refuse les valeurs nulles, dupliquées et malformées. La capture n'est **pas autorisée sur MPI multirang** : un système distribué doit d'abord être assemblé de façon globale avec un protocole explicite. Les captures ne contiennent pas de coordonnées de cellules ni d'informations de forces ; ne pas les confondre avec une campagne aérodynamique complète.

### Auditer l'histoire CFD

```bash
python3 tools/diagnostics/qualify_convergence_history.py history.json \
  --momentum-relative 1e-8 --pressure-relative 1e-8 \
  --continuity-relative 1e-8 --mass-relative 1e-8 \
  --min-full-order 4 --output result.json
```

Les quatre seuils illustratifs ci-dessus **doivent être remplacés par les seuils existants du cas** ; ils ne constituent pas de nouvelles tolérances CFDX. Le succès suppose un verdict `converged` du solveur, des métriques finies, l'éligibilité explicite à l'ordre supérieur, une fraction d'ordre 2 égale à 1 et son maintien pendant les dernières itérations. Cet audit ne vérifie **ni les forces, ni GCI, ni les erreurs physiques SST**.

## Résultats locaux sur matrices de pression NACA0012 (itération 1)

Les maillages de qualification précédemment générés sont réutilisés ; l'itération de capture n'est PAS un état CFD convergé. FGMRES, même norme de résidu vrai originale cible `1e-10` :

| Maillage | Taille CSR | AMG natif | AMG directionnel | AMG énergétique |
|---|---:|---|---|---|
| N64 | 3 096 x 3 096 ; 15 222 nnz | 8 it., `1.31e-11` | 15 it., `3.52e-11` | Non mesuré dans ce replay |
| N128 | 12 336 x 12 336 ; 61 166 nnz | 9 it., `3.18e-11` | NOT_APPLICABLE : diagonale invalide, niveau grossier **6** | NOT_APPLICABLE : V-cycle échoué |
| N256 | 49 248 x 49 248 ; 245 214 nnz | 10 it., `2.82e-11` | NOT_APPLICABLE : diagonale invalide, niveau grossier **7** | NOT_APPLICABLE : V-cycle échoué |

Les temps ne sont pas une campagne statistique. Le benchmark batch retourne volontairement **1** lorsque directionnel échoue sur N128/N256. L'AMG natif fonctionne : ceci isole un défaut des hiérarchies expérimentales, pas un défaut linéaire universel.

## Qualification locale exécutée

- Release CTest ciblé : `test_linear_system_capture`, `test_incompressible`, `test_convergence_history` : **3/3**.
- Capture réelle de la pression après élimination de jauge : `test_steady_incompressible_solver` (compilation Clang O0 isolée) : **30/30** scénarios.
- Python : audit d'historique **9/9**, campagne batch **5/5** (total **14/14**).
- ASan/UBSan : capture C++ **7/7** scénarios, sans erreur dans les chemins testés.
- NACA : capture réelle et rejeu N64/N128/N256 ; PAS de convergence non linéaire démontrée par ces exécutions (budget 1–2 itérations choisi pour captures).

## Verrous et suite précise

1. Diagnostiquer les matrices grossières aux niveaux 6/7 avec leurs diagonales, connexion des agrégats et défauts de prolongation. Ne pas assouplir les critères ; ajouter un repli sûr, certifié sur le résidu original, si une hiérarchie invalide est rencontrée.
2. Rejouer ces mêmes matrices après correctif avec `native`, `directional`, `energy`, temps complet et complexité Galerkin ; comparer aussi les allocations et la mémoire.
3. Récolter les matrices aux itérations 20/21/240 sur N128/N256 (les captures d'itération 1 ne caractérisent pas l'instabilité tardive).
4. Poursuivre les simulations **complètes** N128/N256 en TVD van Leer + continuation et certifier Cp, Cf, Cd, Cl et GCI séparément ; vérifier SST et couplage COUPLED sur des campagnes physiques.
5. CI complète GitHub (Release, DebugSanitizers, Python, précommit et audit gate) sur le HEAD exact, puis MPI et maillages 3D. Ces validations ne sont PAS acquises.
