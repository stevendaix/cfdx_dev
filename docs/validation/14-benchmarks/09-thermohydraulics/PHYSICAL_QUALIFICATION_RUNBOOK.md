# TH-01 / TH-03 / TH-07 / TH-08 — runbook de qualification physique CFDX

## 1. Limite démontrée de la PR initiale

Les quatre GitHub Actions `th0*` ne lançaient que `validate_thermohydraulic_reference.py`, c'est-à-dire des contrôles de **métadonnées**. Aucune comparaison numérique CFDX/référence n'était exécutée. Un CI vert signifie **CONTRAT OK**, pas **MODÈLE VALIDÉ**.

Le présent patch :

- corrige plusieurs incohérences des données de références ;
- vérifie leur intégrité et les groupes sans dimension ;
- exécute en CI des prérequis réels de production (Ghia rapide, énergie, flottabilité, cavité de de Vahl Davis Ra=1000 à trois maillages) ;
- ajoute un post-traitement indépendant à quatre cas, qui rejette les valeurs physiquement incohérentes, les références non épinglées, l'extrapolation, les résidus non finis, les bilans mauvais et les séries de maillages incomplètes ;
- **ne crée pas** de mesures expérimentales, ne substitue pas un cas laminaire à un cas turbulent, et **n'affiche jamais un PASS de validation physique complète** sans données réelles et preuve du modèle correspondant.

## 2. Commandes reproductibles

```bash
python -m pytest -q tests/python/test_thermohydraulic_reference_contracts.py tests/python/test_thermohydraulic_physical_qualifier.py
python scripts/validate_thermohydraulic_reference.py --output build/thermohydraulic/reference-contracts.json
bash scripts/run_thermohydraulic_qualification.sh build/thermohydraulic/solver-exports
```

Si `build/thermohydraulic/solver-exports/TH-*.json` n'existent pas, la dernière commande **échoue par conception (code 2)** et écrit `build/thermohydraulic/physical-prequalification.json` avec les quatre raisons `INCOMPLETE_OR_FAIL`. Ne pas masquer ce code de retour dans une porte de qualification physique. Les quatre CI de contrat de références produisent séparément un rapport d'absence de preuve, sans prétendre avoir lancé ces simulations.

## 3. Conditions à implémenter/exécuter pour chaque cas

| Cas | Scénarios obligatoires (3 maillages minimum **par** scénario) | Référence indépendante et grandeurs de contrôle |
|---|---|---|
| TH-01 | `Re_tau_180`, `Re_tau_395`, `Re_tau_590` | Fichiers DNS MKM `.means`/`.reystress` ; profils `u_plus`, `reynolds_uv_plus`, `total_shear_plus`, friction `Re_tau`; U+, flux/pression, bilan cisaillement. |
| TH-03 | `Re_h_37500` | Données CMP32 Driver–Seegmiller ; `Cf`, `Cp`, `U/Uref` aux stations `x/h=-4,1,4,6,10` ; réattachement issu du **zéro du cisaillement** et non d'un contour. |
| TH-07 | `Gr_0`, `Gr_640000`, `Gr_950000`, `Gr_1600000` | Wagner/Wetzel : vitesses, températures, cisaillement des parois, flux de chaleur, bilan de flottabilité, `Re_b=4328` et `Pr=0,71`. |
| TH-08 | `isothermal`, `non_isothermal` | OCDE/AEN–Vattenfall : LDV/PIV, température normalisée, bilan massique et enthalpique, données temporelles pour LES/URANS. |

Avoir une simulation de **fluide laminaire** ou `test_natural_convection_de_vahl_davis` ne remplace aucun des scénarios turbulents. Une solution RANS stationnaire ne valide pas la RMS, les spectres ou des contraintes turbulentes résolues. Le résultat actuel du comparateur est au maximum `PREQUALIFICATION_PASS` : il reste à démontrer incertitude statistique, erreur de modélisation, GCI spatial complet et, pour TH-08, l'indépendance du pas de temps et le filtrage de l'instrumentation.

## 4. Export CFDX pour le post-traitement

Dans le répertoire de preuves fourni à `--evidence-dir`, créer `TH-01.json`, `TH-03.json`, `TH-07.json`, `TH-08.json`. Chaque fichier est un manifeste contenant les résultats **réellement calculés** du solveur avec la structure suivante (extrait illustratif, les points ne sont pas des données de CFDX) :

```json
{
  "schema_version": 1,
  "case_id": "TH-03",
  "solver": "CFDX",
  "production_run": true,
  "solver_revision": "<sha-git-exact-40-caracteres>",
  "model": "RANS_SST",
  "case_configuration_file": "TH-03/input.yaml",
  "case_configuration_sha256": "<sha256-de-64-chiffres-hex>",
  "solver_log_file": "TH-03/solver.log",
  "solver_log_sha256": "<sha256-de-64-chiffres-hex>",
  "scenarios": [
    {
      "id": "Re_h_37500",
      "runs": [
        {
          "mesh_id": "coarse",
          "cells": 50000,
          "characteristic_h": 0.004,
          "solver_executed": true,
          "converged": true,
          "true_residual_max": 1e-8,
          "mass_imbalance_relative": 1e-7,
          "qoi": {"reattachment_x_over_h": 6.0},
          "profiles": [
            {
              "name": "cf",
              "station": "lower_wall",
              "units": "Cf",
              "reference_source_url": "<url-de-la-source-primaire>",
              "reference_file": "TH-03/reference/cf.csv",
              "reference_sha256": "<sha256-de-64-chiffres-hex>",
              "simulation_file": "TH-03/coarse/cf.csv",
              "simulation_sha256": "<sha256-de-64-chiffres-hex>"
            }
          ]
        }
      ]
    }
  ]
}
```

C'est volontairement **incomplet** : le fichier réel exige au moins trois entrées `runs` ordonnées grossier → fin, des fichiers de profils pour toutes les grandeurs et stations obligatoires, ainsi que tous les scénarios. Les profils sont des CSV de colonnes `x,value,sigma` pour les références et `x,value` pour les valeurs CFDX. Ils doivent utiliser la même abscisse, le même repère, la même station, des normalisations fixées (`U/Uref`, `U/u_tau`, `Tstar`, etc.), des sources documentées et des hashes SHA-256 exacts. **Ne pas inventer d'incertitudes `sigma` :** les obtenir des expériences ou de blocs DNS/LES indépendants, puis en figer les règles avant de comparer aux résultats.

L'outil rejette par défaut les points hors domaine (pas d'extrapolation), un profil de moins de cinq mesures, un SHA incorrect, un flux non convergé, un résidu vrai supérieur à `1e-6`, un déséquilibre massique supérieur à `1e-4`, un déséquilibre thermique supérieur à `1e-3`, des profils avec RMS pondérée au-delà de `2 sigma`, un pic au-delà de `5 sigma`, une dérive fine/moyenne d'une quantité intégrale **ou de profil** supérieure à 5 %, un profil qui s'écarte davantage de sa limite lors du raffinement, des références qui changent entre maillages, un ordre observé hors `[0.5,5]` ou un GCI relatif fin supérieur à 5 % lorsqu'il est définissable. **Ces bornes sont des présélections d'audit et pas des seuils de validation physique issus des publications.** Toute incertitude statistique ou de modèle non couverte maintient le statut `NOT_FULLY_QUALIFIED`.

## 5. Erreurs et divergences des contrats précédents

- **TH-01 :** les fichiers officiels MKM sont sous `chan*/profiles/`, non directement à la racine du dossier MKM. Le fichier brut `chan180.means` indique `Re_tau=178.12` ; la PR avait `178.13`. Les étiquettes nominales `180/395/590` ne sont pas des mesures.
- **TH-03 :** `Re_h=37500` et `Uref=44.2 m/s`, `h=0.0127 m` impliquent `nu=1.49696e-5 m²/s` (contrôle automatique). La variante numérique NASA TMR à `Re_h=36000` ne doit pas être confondue avec le contrat expérimental.
- **TH-07 :** les pourcentages `[5.3,11.6,21]` et `[-8.3,-11.6,-19.9]` précédemment étiquetés « évolution de cisaillement pariétal » ne concordaient pas avec les tables. Le contrôle dérive maintenant explicitement, et **séparément**, les changements de `Re_tau` et de `tau_wall` ; le carré intervient dans leur relation.
- **TH-08 :** la valeur publiée `T*≈0.38` ne doit pas être identifiée à `6/(9+6)=0.40` (mélange volumique à propriétés constantes). Il faut conserver les deux et comparer les bilans enthalpiques avec les propriétés de l'eau et les débits retenus par la référence.

## 6. Verrou observé dans le code de production

Le fichier `apps/production_solver/main.cpp` contient un chemin SST **expérimental** (`--experimental-sst`) qui couple la viscosité turbulente au solveur de mouvement. Il ne s'agit pas d'un solveur turbulent déjà qualifié : il exige un cas HDF5 avec métadonnées de paroi et des champs initiaux `k`/`omega` (ou un redémarrage cohérent). Surtout, il **rejette explicitement les cas combinant SST avec énergie, flottabilité ou rayonnement** (guard dans la section `if (options.experimental_sst)`). Cela bloque notamment une qualification TH-07 par SST turbulent avec convection mixte. Les validations de production laminaire/thermique ne permettent pas de contourner ce refus. TH-01 et TH-03 demandent encore géométrie, maillages, conditions de turbulence et post-traitement des profils, absents de cette PR.

TH-08 demande en plus un protocole temporel documenté pour les fluctuations thermiques ; le comparateur actuel de **profils moyens** ne certifie pas ces quantités. La clôture de ces lacunes exige un développement de solveur et des campagnes indépendantes, pas un changement de seuil ni une simple modification du contrat JSON.

## 7. Méthode de diagnostic en cas d'échec

1. D'abord échouer sans diminuer les critères ; conserver l'archive exacte, le SHA du code, les fichiers de maillage, les coefficients physiques et les vrais résidus.
2. Vérifier masse, quantité de mouvement, énergie, les conditions d'entrée et les dimensions sans se limiter aux résidus internes.
3. Comparer l'évolution des intégrales sur au moins trois maillages ; vérifier gradient, cisaillement au mur, `y+`, reconstruction des flux, convection et correction de pression.
4. Comparer les données **brutes** point par point avec interpolation contrôlée et incertitude de mesure ; distinguer erreur numérique, erreur de maillage, erreur statistique et erreur de modèle.
5. Activer uniquement alors les gates DNS/LES RMS/spectres, stations d'écoulement et invariance temporelle correspondant réellement au modèle utilisé.

### Sources primaires

- MKM : <https://turbulence.oden.utexas.edu/data/MKM/>
- Driver & Seegmiller (CMP32) : <https://torroja.dmt.upm.es/turbdata/agard/docs/AGARD-AR-345.pdf>
- Wagner & Wetzel (2022) : <https://doi.org/10.1002/gamm.202200006>
- OCDE/AEN T-Junction (2011), page 48 : <https://www.oecd-nea.org/upload/docs/application/pdf/2021-03/csni-r2011-5.pdf>

## 8. Infrastructure native ajoutée au solveur CFDX

Le binaire `cfdx_production_solver` accepte désormais `--benchmark-export`. Il
écrit, à partir du **dernier état du solveur**, les trois fichiers suivants sous
`<output-dir>/benchmark_observables/` :

- `cell_fields.csv` : coordonnées/volumes des centres, `u_x/u_y/u_z`, pression
  du solveur, **pression statique** (corrigée de `2 rho k / 3` avec SST), `k`,
  `omega`, viscosité turbulente, température disponible, distance paroi,
  résidu de bilan massique de cellule ;
- `face_fluxes.csv` : géométrie orientée de chaque face, patch, flux de masse
  **autoritaire du solveur Rhie–Chow**, et cisaillement au mur ;
- `mass_conservation.json` : bilan entrée/sortie et bilan maximal par cellule.

Les résultats intermédiaires `result_*.vtu` contiennent maintenant `u_x`,
`u_y`, `u_z` et `p`. La pression statique de comparaison externe doit être
lue depuis `cell_fields.csv`, **pas** le champ `p` de calcul SST. Le cisaillement
issu de `face_fluxes.csv` est explicitement une **approximation unilatérale
d'ordre 1**, destinée au diagnostic du point de réattachement, et non un
`Cf` indépendamment qualifié. Sa convergence devra être comparée à la
reconstruction de face quadratique existante et aux références expérimentales.

### Lancer et auditer une campagne

1. Copier `benchmarks/thermohydraulics/campaign.example.json` et renseigner
   les vrais cas `.cfdx.h5` de trois grilles par scénario, sans modifier les
   identifiants de scénarios. Les chemins sont relatifs au fichier campagne.
2. `python scripts/run_thermohydraulic_cases.py --campaign <campagne.json> --check-only`
   vérifie les entrées et signale les capacités physiques non disponibles.
3. Après build : `python scripts/run_thermohydraulic_cases.py --campaign <campagne.json> --case TH-03 --solver build/th-qualification/cfdx_production_solver --output build/thermohydraulic/native-runs`
   ne lance que les grilles TH-03 pour lesquelles les entrées existent, avec
   `--experimental-sst`. Chaque calcul conserve `solver.log`,
   `run_evidence.json`, les SHA-256 et les fichiers CFDX d'origine.
4. Pour la marche descendante, faire le diagnostic de réattachement :
   `python scripts/extract_th03_wall_diagnostics.py --help`. Préparer
   séparément les profils bruts avec unités/normalisations et les incertitudes
   expérimentales épinglées ; aucun profil de référence n'est fabriqué.
5. Une fois les manifestes `TH-*.json` constitués, lancer
   `python scripts/qualify_thermohydraulic_benchmarks.py --evidence-dir <preuves> --require-pass`.
   **Ce mode exige obligatoirement** pour chaque grille une propriété
   `native_run_evidence_file` ainsi que `native_run_evidence_sha256`, pointant
   vers le `run_evidence.json` émis par le lanceur (chemins relatifs à
   `<preuves>`). Il vérifie les 5 artefacts natifs CFDX, `solver.log`, leur
   empreinte, la convergence, les métriques massiques globales et locales,
   et la cohérence avec le manifeste de profils. Les tests synthétiques du
   comparateur ne sont pas des preuves de simulation.

Les sorties `run_evidence.json` et les fichiers associés doivent être **copiés
sans modification** sous le répertoire de preuves, en conservant leurs chemins
relatifs ; les références indépendantes, les profils normalisés, les réglages
expérimentaux et la version exacte du solveur doivent y être archivés.

### Couverture effective et limitations non résolues

| Benchmark | Production CFDX dans ce correctif | Blocage de validation physique restant |
|---|---|---|
| TH-01 canal MKM | Export de flux et des champs ; comparateur | Entrée périodique 3D, pilotage du débit, DNS/LES temporel, contraintes résolues et statistiques indépendantes non disponibles |
| TH-03 marche descendante | **Lanceur SST expérimental** sur vrais HDF5, diagnostic signed `Cf`, bilans et preuves natives | Fournir géométrie et 3 maillages réellement convergés ; reconstruire précisément la traction, mesurer `Cp`, profils et références expérimentales |
| TH-07 convection mixte | Export des champs thermiques laminaires et prérequis énergie/flottabilité | **SST+énergie+flottabilité non couplés dans le pilote de production** ; `Pr_t`, bilans de chaleur pariétaux et référence indépendante non qualifiés |
| TH-08 jonction T | Export stationnaire et prérequis de conservation | **Transitoire URANS/LES, transport thermique turbulent, moyennes temporelles, variances, spectres et pas de temps** non qualifiés |

En particulier, ce correctif ne supprime pas le garde du solveur qui refuse
`--experimental-sst` avec énergie/flottabilité/rayonnement : supprimer ce garde
sans implémenter les termes turbulents cohérents dans le flux thermique et la
boucle couplée serait une régression physique, pas une amélioration.

### Tests de régression nouveaux

- `test_thermohydraulic_observables` (C++) : signe des flux, conservation,
  cisaillement diagnostique, correction de la pression SST, refus des NaN.
- `tests/python/test_thermohydraulic_native_cases.py` : campagnes trois
  maillages, capacités absentes, fichiers invalides, zéro signé du cisaillement.
- `tests/python/test_thermohydraulic_native_export_e2e.py` : petit cas HDF5
  réellement exécuté par `cfdx_production_solver` avec `--benchmark-export` ;
  flux nuls attendus, champs et VTU vérifiés. **Ce test laminaire n'est pas un
  des benchmarks TH physiques.**
- `tests/python/test_thermohydraulic_physical_qualifier.py` : mutations de
  mesures, empreintes, grilles, et attestation native stricte.
