# CFDX — Protocole strict de qualification physique et numerique (2026-11-08)

## Portee, limite et regle de promotion

Cette branche renforce les **preuves necessaires**. Aucune campagne physique non executee ne devient PASS. Le registre historique comporte 25 cas; un statut historique PASS ne certifie pas le HEAD actuel. Un dossier synthetique qui passe les auto-tests prouve le verificateur, pas le solveur CFD.

L architecture exige trois niveaux independants : **V** (verification des equations/flux par MMS analytique), **S** (verification de convergence du solveur, bilans et etats acceptes), **P** (validation physique sur references preexistantes). Les verifications de traceabilite ne suffisent pas a etablir V ni P; elles bloquent les declarations sans preuve.

Pour une release, aucune relaxation de seuil ni suppression d un test lent ne transforme un diagnostic en PASS. Les cas longs sont une porte de qualification distincte : un smoke CI rapide ne remplace pas la campagne physique.

## Controles automatiques ajoutes

1. `reference_artifact` — reference.json chargee et empreinte compatible.
2. `provenance` — SHA-256 binaires/configuration/log/code oracle/reference/ZIP source; commit Git propre si Git.
3. `raw_evidence` — maillage, historique CSV et oracle volumes finis distincts, SHA et verification du dernier iterate.
4. `case_contract` — equations/modeles, Re, Mach, incidence et normalisations conformes aux references.
5. `solver` — 3 maillages converges, meme schema final, iterations/temps reels.
6. `residuals` — residus non lineaires et residus FV independants; TVD final obligatoire pour NACA.
7. `qoi_stability` — 3 derniers etats acceptes: variations Cd/Nu/profils inferieures a 10 % du budget de reference.
8. `conservation` — bilans locaux/integrales : masse, quantite de mouvement, chaleur, rayonnement selon cas.
9. `force_decomposition` — Cd pression + visqueux = Cd total; integration surfaces independante et controle volume.
10. `reference` — toutes les grandeurs physiques imposees, valeurs, incertitudes et unites independantes.
11. `mesh_series` — qualite geometrique positive; serie de maillages et raffinement compatibles; NACA couche limite.
12. `spatial_order` — 3 points, tendance convergente de chaque QoI, ordre observe, GCI physique, extrapolation.

Si un seul controle manque ou echoue, le statut est `NOT_QUALIFIED`. Le seul succes autorise est `EVIDENCE_GATES_PASS`, volontairement distinct d une certification complete des equations.

## Obligations par cas (25/25)

| Identifiant | Statut registre (historique) | QoI physiques exigees | Bilans indispensables |
|---|---|---|---|
| `LAM-COUETTE` | PASS | u_mean, velocity_profile_l2, wall_shear | masse |
| `LAM-POISEUILLE` | PASS | volume_flow, velocity_profile_l2, pressure_gradient | masse |
| `INC-GHIA` | PASS | u_centerline_rms, v_centerline_rms, vortex_position_error | masse |
| `VER-MMS` | READY | field_l2, field_linf, face_flux_error | oracle MMS/geometrique |
| `FORCE-VMFL036` | READY | Cd, Cd_pressure, Cd_viscous | masse |
| `EXT-NACA0012` | DIAGNOSTIC | Cd, Cl, Cm, Cp_profile_l2, Cf_profile_l2 | masse, momentum x/y |
| `INT-CHANNEL` | DIAGNOSTIC | bulk_velocity, friction_factor, pressure_gradient | masse |
| `SEP-BFS` | DIAGNOSTIC | reattachment_length, velocity_profile_l2 | masse |
| `TUR-FLATPLATE` | PLANNED | Cf_profile_l2, Uplus_profile_l2 | masse |
| `TUR-CHANNEL` | BLOCKED | velocity_profile_l2, reynolds_stress_l2 | masse |
| `SEP-BFS-DRIVER` | BLOCKED | reattachment_length, wall_shear_l2 | masse |
| `TH-MIXED-CONVECTION` | BLOCKED | Nu_mean, temperature_l2, pressure_drop | masse, energie |
| `TH-TJUNCTION` | BLOCKED | mixing_temperature_l2, pressure_drop, energy_balance | masse, energie |
| `TUR-NACA0012` | PLANNED | Cd, Cl, Cm, Cp_profile_l2, Cf_profile_l2 | masse, momentum x/y |
| `TH-CONV` | READY | Nu_mean, temperature_l2 | masse, energie |
| `CHT-INTERFACE` | READY | interface_heat_flux, interface_temperature_jump | energie, interface |
| `RAD-P1` | READY | heat_flux, incident_radiation | energie, rayonnement |
| `RAD-DOM` | READY | heat_flux, angular_error | energie, rayonnement |
| `RAD-S2S` | READY | heat_flux, view_factor_reciprocity_error | energie, rayonnement |
| `N9-S8-NACA0012` | READY | Cd, Cl, Cm, Cp_profile_l2, Cf_profile_l2 | masse, momentum x/y |
| `T1-SOLID-COND` | READY | heat_flux, temperature_l2 | energie |
| `T2-FLUID-ENERGY` | READY | energy_flux, temperature_l2 | masse, energie |
| `INC-TRANSPOSE-STRESS-MMS` | DIAGNOSTIC | momentum_error_l2, shear_traction_error_l2 | masse |
| `SCALAR-TRANSPORT-NONORTH-AFFINE` | DIAGNOSTIC | scalar_error_l2, flux_error_l2 | scalaire |
| `WALL-DISTANCE-GEOMETRY` | DIAGNOSTIC | distance_l2, distance_linf, near_wall_linf | oracle MMS/geometrique |

## Structures des fichiers de preuve

Une preuve valide pour `LAM-COUETTE` se trouve dans `qualification_evidence/LAM-COUETTE/evidence.json`, avec `reference.json` dans le meme dossier et des fichiers attaches haches. Chaque cas a son propre dossier et trois maillages.

`evidence.json` (structure, valeurs indicatives **pas une campagne CFD**) :

```json
{
  "schema": "cfdx-physical-evidence/v1",
  "case_id": "LAM-COUETTE",
  "setup": {"model":"laminar","reynolds":100,"mach":0,"alpha_degrees":0},
  "provenance": {
    "source_kind": "archive",
    "source_revision": "<SHA-256_REEL_DE_L_ARCHIVE>",
    "artifacts": {
      "source_archive":{"path":"source.zip","sha256":"<SHA256>"},
      "binary":{"path":"solver.bin","sha256":"<SHA256>"},
      "configuration":{"path":"configuration.json","sha256":"<SHA256>"},
      "solver_log":{"path":"execution.log","sha256":"<SHA256>"},
      "reference":{"path":"reference.json","sha256":"<SHA256>"},
      "oracle_source":{"path":"oracle_fv.py","sha256":"<SHA256>"}
    }
  },
  "meshes": ["<N16 complet>", "<N32 complet>", "<N64 complet>"]
}
```

Chaque `meshes[]` doit inclure `id`, `status=CONVERGED`, `iterations`, `runtime_seconds`, `final_scheme`, `h`, `mesh_family`, `mesh_quality`, `residuals`, `residual_limits`, `independent_residuals`, `balances`, `quantities` et 3 artifacts `mesh/history/independent_oracle` avec leur SHA-256. Les CSV de l historique doivent porter la derniere iteration et **toutes les colonnes de QoI exigees** ; le PDF/une capture d ecran n est pas une preuve numerique suffisante.

`reference.json` contient `case_id`, `setup` strictement identique, `citation` physiquement traçable, `independent_of_solver=true` et un objet `quantities` pour chaque QoI, avec `value`, `unit`, `uncertainty`, `absolute_tolerance`, `relative_tolerance`. Aucun etalon n est choisi a partir de la solution CFDX a qualifier.

L oracle sous `cfdx-independent-fv-oracle/v1` doit contenir `case_id`, `mesh_id`, `residuals`, `balances` (et `Cd_control_volume` pour cas force). Le hash verifie les octets ; **la veritable independance doit etre etablie par revue du code oracle**, une autre formulation/discretisation ou des cas MMS. Une valeur stockee deux fois ne constitue pas un controle independant.

## Commandes de qualification

```bash
# Smoke tests rapides : auto-tests des detecteurs et contre-exemples synthetiques.
python -m unittest -q tests.python.test_naca_mesh_similarity_gate tests.python.test_naca_force_mesh_qualification tests.python.test_naca_surface_force_oracle tests.python.test_qualification_evidence_gate

# Gate de release, qui refuse les dossiers manquants.
python scripts/qualification_evidence_gate.py --registry docs/validation/CFDX_QUALIFICATION_REGISTRY.json --evidence-dir qualification_evidence --output coverage.json --require-all
```

### Campagne NACA, export des faces et verificateur GCI

**N64, N128 et N256 doivent etre resolus separement** avec le meme schema TVD final. Les paths ci-dessous sont exemplaires. Le demarrage UPWIND/TVD continuation est permis, mais le dernier etat accepte doit etre entierement TVD.

```bash
# Pour chaque N choisi : exemple N64.
CFDX_NACA_CONVECTION_TVD=1 \
CFDX_NACA_WALL_TRACTION_FACE=1 \
CFDX_NACA_PHYSICAL_HISTORY_CSV=out/N64.history.csv \
CFDX_NACA_FACE_TRACTION_RAW_CSV=out/N64.faces.csv \
CFDX_NACA_CONVERGED_FORCE_CSV=out/N64.force.csv \
./build/test_naca0012_qualification --single meshes/naca0012_laminar_n64.h5

python scripts/naca_surface_force_oracle.py --raw-faces out/N64.faces.csv --certificate out/N64.force.csv --output out/N64.surface-check.json

python scripts/naca_force_mesh_qualification.py \
  --mesh-quality out/mesh_quality.json \
  --cert-64 out/N64.force.csv --history-64 out/N64.history.csv --raw-faces-64 out/N64.faces.csv \
  --cert-128 out/N128.force.csv --history-128 out/N128.history.csv --raw-faces-128 out/N128.faces.csv \
  --cert-256 out/N256.force.csv --history-256 out/N256.history.csv --raw-faces-256 out/N256.faces.csv \
  --output out/NACA.grid-gci.json
```

Attention : `GCI_CANDIDATE` n est pas une validation experimentale. Les profils `Cp`, `Cf`, la decomposition pression/cisaillement, le bilan de quantite de mouvement autour du volume de controle et la concordance Reynolds/Mach/alpha sont des portes independantes du calcul des forces total.

### Diagnostic A/B de la pression a la paroi

`reconstruct_wall_face_traction()` accepte maintenant un `dp_dn` explicite ; en l absence de ce parametre, l hypothese homogene precedente subsiste. `CFDX_NACA_WALL_PRESSURE_NORMAL_MODE=cell_gradient` avec `CFDX_NACA_WALL_TRACTION_FACE=1` active un essai A/B. Il NE CORRIGE PAS la condition limite pression du solveur et ne doit pas etre promu sans identite de la paroi, MMS et convergence spatiale independante.

### Benchmark N15 : provenance a l execution

Le benchmark N15 ne declare plus `PASS` uniquement lorsque le sous-processus retourne 0. Le dossier doit provenir d un commit Git **propre** de la racine demandee, ou fournir `--source-archive /chemin/vers/source.zip`. Le binaire effectivement lance et la fiche de cas sont identifies par SHA-256. Sans source explicite, le rapport peut contenir les mesures mais la sortie est `UNVERIFIED_SOURCE`, code retour 2, et non un PASS. Un ZIP n est pas une preuve de correspondance avec le binaire : l evidence release exige une reconstruction/compilation sur la source identifiee.

```bash
python scripts/n15_benchmark.py run --case cas_n15.json --repo . --source-archive /chemin/vers/source.zip --output out/benchmark.json -- /chemin/vers/solveur
```

### Mode de panne attendue / provenance

Les anciennes campagnes PR #420/#428 peuvent etre historiquement PASS. Elles restent `NOT_QUALIFIED` sur cette copie tant que leurs etats, maillages, references et logs n ont pas ete reattaches et controles sur le **meme SHA source**. Le mode `archive` exige le hash du vrai fichier source, le mode `git` un commit complet propre.

### Suites necessaires avant toute qualification physique generale

- Achever N64/N128/N256 en TVD avec residus non lineaires et FV independants, bilans, profils et forces corrects ; traiter les familles isotropes et anisotropes comme deux jeux de preuves, pas une seule courbe GCI si non similaires.
- Comparer Cp et Cf et Cd pression/cisaillement a des references laminaire Re=1000 d une part, NASA RANS Re=6 millions d autre part : ne jamais melanger ces references.
- Verification MMS sur Rhie-Chow, pression, nonorthogonalite, traction paroi et esquisse de couche limite ; chaque expression non couverte reste explicitement non verifiee.
- Refaire Ghia, De Vahl Davis, SST, CHT, P1, DOM, S2S, conduction, distance a la paroi et autres cas du registre sur plusieurs resolutions reelles avec dossiers reproductibles.
- Garder des campagnes longues dans les gates de release en plus des smoke tests CI rapides. Les 59 tests Python actuels n executent pas 25 physiques.
