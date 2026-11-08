# CFDX CI — point d'entrée

**Objectif : une PR obtient un diagnostic rapide, la qualification longue reste
exhaustive, chaque validation est identifiable par son nom, son déclencheur et
ses preuves.** Les critères numériques/physiques ne sont pas relâchés.

Pour les responsabilités, les déclencheurs, le détail des partitions et les
commandes locales : [Guide de la CI](../../docs/ci/CI_ARCHITECTURE.md).

| Niveau | Workflow(s) | Déclenchement | But et verdict |
|---|---|---|---|
| **1 — PR logicielle** | `precommit.yml`, `cfdx-ci.yml` | Chaque PR | Contrôles de politique CI, Release, DebugSanitizers, CTest court, Python, N12 rapide |
| **2 — Maturité numérique** | `cfdx-numerical.yml`, `cfdx-maturity.yml` | PR/push | N1–N15 et registres de maturité; N7 Re400 long dans le niveau 3 |
| **3 — Qualification complète** | `cfdx-validation.yml` | `validation-total`, lundi, manuel | Build MPI/HDF5, CTest complet, 4 shards `long`, wall-distance, thermique, E2E et rapport |
| **Diagnostic N9** | `n9-campaign.yml`, `n9-physical-validation.yml` | Manuel **uniquement** | Reproduction ciblée; N9 est déjà qualifié au niveau 3 |
| **Fonctions spécialisées** | GUI, docs, CUDA, TH, contrats, autres | Selon leur workflow | Ne sont pas toutes des qualifications physiques CFDX |

**Une seule autorité pour la campagne longue :** `cfdx-validation.yml`.
Les deux workflows N9 manuels ne tournent plus en parallèle lorsque le label
`validation-total` est ajouté à une PR.

**Attention aux checks exigés dans GitHub :** les jobs facultatifs ou filtrés ne
doivent pas être imposés comme checks requis sans une politique de merge qui
garantit leur déclenchement. Le label n'est pas une preuve qu'une campagne est
passée : seuls les résultats du **HEAD testé**, les artefacts et le `final-gate`
le sont. Pour exiger systématiquement la qualification avant fusion, configurer
explicitement le ruleset/processus de fusion de l'équipe (voir guide).

**Aucune promesse de durée avant mesure :** le planificateur équilibre des
`TIMEOUT` (bornes supérieures), pas des temps réels. Les artefacts JUnit et
`*-timings.json` permettent enfin de mesurer la vraie durée par test.

## Activer une qualification sur demande

**Un point d'entree unique :** Actions > **CFDX qualification - activate a profile**,
choisir le profil et `quick`/`extended`, puis `Run workflow`. Sur PR, une seule
etiquette `qualification:naca`, `qualification:sst`, `qualification:thermal`,
`qualification:momentum`, `qualification:linalg`, `qualification:wall-distance`,
`qualification:n9`, `qualification:ghia` ou `qualification:smoke` active le
profil approfondi et le relance apres chaque commit. `validation-total`
conserve la validation physique canonique et a priorite sur ces etiquettes.

Guide et commandes locales : [Activation des qualifications](../../docs/ci/QUALIFICATION_ACTIVATION.md).
**Attention :** un profil CTest au vert n'est pas une preuve physique de
qualification complete (NACA N256, GCI, SST, thermique, etc.).
