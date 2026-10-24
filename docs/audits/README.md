# ÉSPACE D'AUDIT MÉTHODE — docs/audits/ — Fiches Méthode

> **Constitution technique de cfdx_dev.** Ce dossier applique les **4 Règles d'Or** décidées dans le Plan Directeur de Qualification et d'Audit (2026-10-12). Il est la porte d'entrée obligatoire de tout développement numérique.

## Scope et autorité

Chaque fiche est un **audit mathématique par méthode** : forme continue, discrétisation Volume Fini exacte, traitement des frontières, pièges connus, ordre de précision attendu, références (DOI/ISBN).

Ce dossier **ne duplique aucun registre existant** (cf. `docs/validation/CFDX_VV_GOVERNANCE.md §16.3–§16.5— pas de matrice en doublon, un seul propriétaire par fait) :

| Dimension | Propriétaire autoritaire | Rôle de la fiche |
|---|---|---|
| Statut des cas de qualification | `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` | cible les campagnes exigées par la méthode |
| Traçabilité exigence → preuve (#461) | `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json` | identifie les exigences N x |
| Terminologie / gates V&V | `docs/validation/CFDX_VV_GOVERNANCE.md` | s'y réfère, ne redéfinit rien |
| Contrat mathématique d'une méthode | **la fiche de ce dossier** | source unique du « quoi vérifier » avant d'écrire le code |

## Les 4 Règles d'Or (non négociables)

1. **Zéro implémentation sans Fiche Méthode.** Aucune équation ou algorithme n'est codé sans fiche approuvée dans `docs/audits/`, référence DOI/ISBN à l'appui.
2. **Zéro fusion sans benchmark et GCI.** Toute méthode est comparée quantitativement (L2/L∞) à une solution analytique, un benchmark établi (Ghia) ou un code tiers (OpenFOAM/SU2/Fluent) **sur le même maillage exact**, avec étude de convergence (GCI).
3. **Zéro régression silencieuse.** Bilans (masse, QDM, énergie) à la tolérance machine ; bornitude des champs (VOF ∈ [0,1], k ≥ 0…) vérifiée **avant** clipping, avec télémetrie.
4. **Zéro divergence parallèle.** CPU vs MPI vs GPU reproductibles à la tolérance machine stricte (1e-12) sur le même maillage ; aucun fallback silencieux.

## Processus

1. **Avant de coder** → copier `_FICHE_METHODE_TEMPLATE.md’ en `NN-nom-de-la-methode.md’, remplir toutes les sections, citer les DOI/ISBN.
2. La fiche passe en revue dans la même PR que l'implémentation (ou une PR préalable séparée pour les méthodes lourdes).
3. La PR doit cocher la checklist de fusion (`.github/pull_request_template.md’) et référencer la(les) fiche(s) concernée(s).
4. Le job CI `audit-gate’ vérifie mécaniquement : structure des fiches, présence d'une fiche référencée pour toute PR touchant le noyau numérique.
5. Après fusion → mettre à jour `QUALIFICATION_REGISTRY.md’, le registre JSON autoritaire et l'issue #34.

## Index des fiches

| Fiche | Méthode | État |
|---|---|---|
| [01-diffusion-laplacien.md](01-diffusion-laplacien.md) | Opérateur diffusion / Laplacien FV | ✓ rédigée — code existant à auditer contre elle |
| [02-schemas-convection.md](02-schemas-convection.md) | Schémas de convection (UPWIND…QUICK, correction différée) | ✓ rédigée — code existant à auditer contre elle |
| [03-interpolation-rhie-chow.md](03-interpolation-rhie-chow.md) | Interpolation de Rhie–Chow / colocalisé | ✓ rédigée — code existant à auditer contre elle |
| [04-couplage-pression-vitesse.md](04-couplage-pression-vitesse.md) | SIMPLE / SIMPLEC / PISO / PIMPLE / COUPLED | ✓ rédigée — code existant à auditer contre elle |
| [05-schur-preconditionnement.md](05-schur-preconditionnement.md) | Complément de Schur : PCD, LSC, BFBT, MGR, AMG | ✓ rédigée — code existant à auditer contre elle |
| [06-gradients-reconstruction.md](06-gradients-reconstruction.md) | Gradients Green–Gauss / moindres carrés (N2) | ✓ rédigée — PR #511 en attente |
| [07-benchmarks-references.md](07-benchmarks-references.md) | Oracles de validation : Ghia, Poiseuille, Couette, MMS, GCI | ✓ rédigée |
| [08-codes-de-reference.md](08-codes-de-reference.md) | SU2, OpenFOAM, Fluent, Basilisk — comparaisons croisisées | ✓ rédigée |
| [09-turbulence-rans.md](09-turbulence-rans.md) | Turbulence RANS : k-ε, k-ω, SST, SA + lois de paroi | APPROUVÉE (2026-10-12, PR #861) ; reste requis (Phase 3, audit du code Phase 10) |
| [10-thermique-cht.md](10-thermique-cht.md) | Énergie ρ·cp, CHT multi-régions, Boussinesq | APPROUVÉE (2026-10-12, PR #861) ; reste requis (Phase 3, audit du code Phase 11) |
| [11-rayonnement.md](11-rayonnement.md) | Rayonnement : Rosseland, P1, S2S, DOM | APPROUVÉE (2026-10-12, PR #861) ; reste requis (Phase 3, audit du code Phase 12) |
| [12-parallelisme-gpu-ooc.md](12-parallelisme-gpu-ooc.md) | Parité CPU/MPI/GPU/OOC (Règle d’Or 4) | APPROUVÉE (2026-10-12, PR #861) ; reste requis (Phases 6–7, N13) |
| [23-variable-viscosity-stress.md](23-variable-viscosity-stress.md) | Contrainte visqueuse transposée à viscosité variable (F11) | Implémentation + oracle affine ajoutés ; CI exacte et MMS/GCI polyédrique en attente |

Les fiches à venir (VOF/PLIC, maillage dynamique, FSI, GUI) seront créées **avant** toute implémentation correspondante (Phase 4 de la roadmap). Les fiches 09–12 couvrant les Phases 1–3 du Plan Directeur sont approuvées (2026-10-12, PR #861) : l'audit du code contre chaque fiche reste requis avant toute promotion multiphysique (statuts BLOCKED du registre de qualification).

## Règle d'exception

Aucune. Une fiche incomplète bloque la PR. Une fiche peut être marquée « PARTIAL » uniquement sur les sections « Statut CFDX », jamais sur les sections mathématiques (1–6).
## Note d'infrastructure CI (2026-10-12)

Pendant la qualification N9, un défaut **préexistant sur master** a été identifié et documenté (issue #858, correctif PR #857) : les workflows `cfdx-numerical.yml` (contenu dupliqué , YAML invalide) et `cfdx-validation.yml` (scalaire non quoté dans `application-e2e`) sont incapables de démarrer depuis les fusions du 5 octobre. Conséquence directe : chaque push produit deux « startup failures » et la « CFDX total validation » (schedule + label `validation-total`) est silencieusement inopérante. La chaîne d'évidence de la campagne N9 n'est **pas** affectée : elle repose sur `n9-campaign.yml` et `n9-physical-validation.yml`, workflows valides et indépendants. Aucune tolérance , aucun test, aucune politique n'est modifiée par le correctif #857.


**Suivi (2026-10-12 : 17 h 00Z).** Depuis la rédaction de cette note : `#854` (constitution technique) a été fusionnée dans `master` (squash `c1a5b6e`, 16 h29Z) ; la preuve de guérison de `cfdx-numerical.yml` est acquise — le workflow a exécuté le correctif `#857` de bout en bout (run `#1826`, 7 min 45 s) ; et un conflit add/add sur `.github/workflows/precommit.yml` (step « Ruff diagnostics » ajouté indépendamment des deux côtés de la pile) a été résolu côté base (commit `4fe52f3a`, sur-ensemble fonctionnel couvrant aussi `scripts/n9_physical_campaign.py`), ce qui rétablit la fusionabilité de la PR empilée `#855` sans modifier le HEAD à qualifier.
