# CFDX — pilotage transversal de toutes les qualifications

Ce tableau de bord est **un inventaire vérifiable et un journal de progression**, pas
une preuve scientifique automatique. Les 25 cas de
`docs/validation/CFDX_QUALIFICATION_REGISTRY.json` doivent pouvoir progresser
**en parallèle**. NACA n'est ni un verrou CI global ni le critère de réussite de
la turbulence, des thermiques, de Ghia ou de la distance à la paroi.

## Principe

Une ligne par cas suit quatre états **indépendants** :

1. `registry_status` : état historique (les anciens `PASS` ne certifient pas
   l'archive présente).
2. `contract_status` : CTest courts/contrats explicitement liés au cas.
3. `campaign_status` : CTest de solveur/campagne liés au cas, distincts de la
   preuve de qualification physique; `NO_DEDICATED_TEST` est visible.
4. `physical_evidence_status` : verdict de
   `scripts/qualification_evidence_gate.py` sur les **preuves réelles** (source,
   binaire, configurations, logs, maillages, oracles indépendants, QoI,
   convergence, bilans et ordre spatial). Les tests de métadonnées ne peuvent
   jamais le promouvoir.

Les rattachements explicites se trouvent dans
`docs/validation/QUALIFICATION_PORTFOLIO_LINKS.json` : 25 identifiants du
registre, ancrages CTest, workflows de référence, prochaine action par cas.
Le catalogue comporte également une liste *à inscrire au registre*, suivie
séparément : **de Vahl Davis, Rayleigh–Bénard, VOF avec transport, LES,
VMFL056 et NACA 3D**. Un cas de cette liste n'est pas automatiquement compté
comme qualifié.

Un nom présent dans le manifeste CTest ne prouve **aucune exécution**. Un
résultat `PASS` n'apparaît dans les colonnes CTest que si un JUnit réel fournit
le résultat correspondant. Manquant, désactivé, échoué, dupliqué ou
non exécuté restent distingués. Des workflows TH-01/03/07/08 de vérification
de référence ne sont pas des calculs CFDX de production.

## Vue continue GitHub / CI

- **Chaque PR :** `cfdx-qualification.yml` publie le portefeuille des 25 cas
  avec la totalité de l'inventaire CTest (compilation non nécessaire). Cela
  détecte automatiquement les tests référencés absents, et rend les `NOT_RUN`
  visibles, même sans étiquette PR.
- **Chaque PR Release/DebugSanitizers :** `cfdx-ci.yml` réutilise son JUnit
  CTest existant pour publier les vraies réussites/échecs **sans exécuter une
  seconde fois** les solveurs. Les tests longs ne sont pas enregistrés dans
  les builds rapides : leur absence dans ce manifeste est attendue; le
  contrôle strict des noms a lieu dans le job d'inventaire avec l'option
  `CFDX_ENABLE_LONG_VALIDATION=ON`.
- **Sur demande par domaine :** `qualification:<profil>` ou `workflow_dispatch`
  continue à produire le JUnit et, désormais, la projection de ce résultat
  sur les 25 lignes. Le résultat du domaine ne vaut pas validation globale.
- **Campagne totale :** `cfdx-validation.yml` agrège les JUnit de la suite
  rapide et des 4 partitions longues, puis publie la même matrice avec les
  `NOT_RUN`/`FAIL` éventuels. La CI échoue toujours sur ses contrôles existants,
  sans changement des critères physiques.

Artefacts typiques : `qualification_portfolio.md`, `.json`, `.csv`,
`physical_evidence_gate.json`. Le CSV peut être filtré par domaine ou
comparé entre deux exécutions. L'absence du dossier `qualification_evidence/`
est présentée comme **0 cas munis de preuves vérifiables**, jamais comme
réussite implicite. Une archive source peut être identifiée par SHA-256,
mais cette empreinte seule **n'atteste pas le binaire utilisé**.

## Reproduction locale : tous les cas sans coût de calcul CFD

```bash
cmake -S . -B build-qualif -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
  -DCFDX_ENABLE_LONG_VALIDATION=ON -DCFDX_ENABLE_MPI=OFF
ctest --test-dir build-qualif --show-only=json-v1 > manifest.json

python3 scripts/qualification_portfolio.py \
  --ctest-manifest manifest.json \
  --evidence-dir qualification_evidence \
  --output-dir out/qualification-portfolio \
  --strict-links
```

Après des calculs CTest, incorporer **leurs** fichiers JUnit (un ou plusieurs,
chacun unique), sans modifier leurs seuils :

```bash
ctest --test-dir build-qualif --output-on-failure --no-tests=error \
  -R '^(test_mms_scalar_diffusion|test_radiation_vv|test_s2s_radiation_vv)$' \
  --output-junit "$PWD/tests-domaines.xml"
python3 scripts/qualification_portfolio.py \
  --ctest-manifest manifest.json --junit tests-domaines.xml \
  --evidence-dir qualification_evidence \
  --previous-report out/qualification-portfolio/qualification_portfolio.json \
  --output-dir out/qualification-portfolio-suivant
```

L'option `--previous-report` affiche chaque changement de statut entre deux
captures. `--strict-links` sert à contrôler les ancrages de test sur une
configuration CMake **avec** les campagnes longues; `--require-evidence` est
réservée à une porte de promotion de qualification physique, et doit échouer
si un seul cas n'a pas ses preuves complètes.

## Plan de travail parallèle (sans attendre NACA)

| Chantier | Court (PR / quotidien) | Prolongé (indépendant) | Preuve finale |
|---|---|---|---|
| Laminaire et MMS | Couette, Poiseuille, champs manufacturés | Séries 3 grilles, flux et ordre | Solution analytique, bilans, erreur |
| Écoulements internes | Ghia, channel, BFS quick | Ghia Re 100/400/1000, canal et BFS raffinés | Profils, rattachement, conservation |
| Turbulence | Positivité et contrats SST | Canal turbulent, plaque plate, BFS Driver, NACA RANS | NASA/DNS, Cf, U+, y+, contraintes |
| Énergie / convection | T1, T2, CHT, conduction | de Vahl Davis Ra 1e3..1e6 et CHT multi-région | Nu, profils T, bilans énergie/interface |
| Rayonnement | P1, DOM, S2S contrats | VMFL056 et configurations optiques réalistes | Flux, rayonnement incident, réciprocité |
| Parois et maillages | Distance analytique et géométrie | N128, stress 65³, méthodes approximatives | Erreurs L2/L∞, résidus vrais, temps |
| NACA / pression-vitesse | N9, flux Rhie–Chow, traction | N64/N128/N256 même TVD, GCI, 3D | Cp/Cf/forces, momentum, similitude |
| Multiphase / LES | Contrat de courbure VOF | Advection conservatrice, cas LES DNS | Bilan de phase, énergie/statistiques |

Les lignes `READY`, `BLOCKED` et `PLANNED` n'autorisent pas la suppression
des tests existants. Les nouveaux cas physiques exigent référence indépendante,
contrat d'acceptation, paramètres et maillages figés, solveur réel, preuve
numérique attachée au SHA d'exécution et revue indépendante.

## Premier constat sur le ZIP du 8 novembre 2026

- 391 CTest enregistrés en Release avec validations longues, dont 33 `long`.
- 25 cas au registre : 3 `PASS` historiques, 10 `READY`, 6 `DIAGNOSTIC`,
  4 `BLOCKED` et 2 `PLANNED`.
- 17 sur 25 n'ont pas encore d'ancrage de campagne de solveur **distinct**
  dans cette cartographie (ce qui n'interdit pas l'existence de tests de
  composants, de scripts ou de campagnes de référence séparées).
- 0/25 dossiers de preuves physiques fournis dans ce ZIP. Un contrôle
  d'absence de preuves doit retourner `NOT_QUALIFIED` pour chacun.
- Dix tests CTest de composants sur plusieurs domaines ont été compilés
  et exécutés localement : 10/10 PASS (dont Poiseuille, diffusion MMS,
  turbulence, transport scalaire, contraintes visqueuses, thermique,
  rayonnement, distance et courbure VOF). La CI GitHub complète **n'a pas** été
  exécutée et aucun `PASS` physique n'est inféré de ce contrôle.

Protocole détaillé :
`docs/validation/CFDX_QUALIFICATION_EVIDENCE_PROTOCOL_2026-11-08.md`.
