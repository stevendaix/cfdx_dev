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

Les fiches à venir (VOF/PLIC, maillage dynamique, FSI, turbulence, rayonnement) seront créées **avant** toute implémentation correspondante (Phases 3–4 de la roadmap).

## Règle d'exception

Aucune. Une fiche incomplète bloque la PR. Une fiche peut être marquée « PARTIAL » uniquement sur les sections « Statut CFDX », jamais sur les sections mathématiques (1–6).
