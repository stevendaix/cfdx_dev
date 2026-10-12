# QUALIFICATION_REGISTRY — Tableau de suivi exigences → preuves

> **Tableau de bord du Plan Directeur de Qualification et d'Audit (constitution technique, 2026-10-12).**
> Ce document est **synthétique et non autoritaire** (gouvernance §16.2–§16.3) : les sources de vérité machine-readable restent `docs/validation/CFDX_QUALIFICATION_REGISTRY.json’ (statut des cas) et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json’ (traçabilité exigences #461). Toute promotion doit modifier le registre JSON **et** ce tableau dans la même PR (gouvernance §16.7).

Légende — Implémenté / Test / Oracle / SHA : ✓ = prouvé, ~ = partiel, ✗ = absent. Un SHA validé n'est compté que s'il désigne un commit **mergé** sur lequel la campagne a réellement été exécutée.

## Phase 1 — Noyau incompressible (verrouillage avant toute multiphysique)

| Exigence | Fiche méthode | Implémenté | Test exécutable | Oracle quantitatif | SHA validé | Gap principal |
|---|---|---|---|---|---|---|
| N9 Couplage P-V (campagne physique complète : Couette, Poiseuille, Ghia Re=100, skewed ; 6 modèles Schur ; vrais résidus ; chemin `.cfdx.h5’ complet) | [04](docs/audits/04-couplage-pression-vitesse.md), [03](docs/audits/03-interpolation-rhie-chow.md) | ~ | ✓ (négatifs N9.7 #851 ; **pilote de campagne scripts/n9_physical_campaign.py** + workflow n9-campaign, label validation-total) | ~ (systèmes synthétiques ; Ghia non démontré sur HEAD) | ✗ | Pilote de campagne N9 fourni (9 tests requis + équivalence MPI requise quand le build l’expose + pytest body-force opt-in ; rapport JSON exact-HEAD, gates inchangés) : **campagne 6×4 à exécuter sur HEAD exact** ; équivalence COUPLED/ségrégé à documenter — **#601 reste ouvert jusqu’à campagne verte** |
| N8 Schur/AMG/PCD (politique d'acceptation quantitative + MPI) | [05](docs/audits/05-schur-preconditionnement.md) | ✓ (6 modèles) | ~ | ~ | ✗ | Seuils mesurés sur population CFD réelle manquants ; AMG/Schur/PCD/MGR **non qualifiés MPI** — matrice #481 à remplir |
| N11 Conservation/bornitude (frontières + MPI + télémetrie clipping) | [01](docs/audits/01-diffusion-laplacien.md) | ~ | ✓ (Poisson MPI #814) | ~ | ✓ (Poisson) | Bilans QDM/énergie production ; signes/orientations ; bornitude par variable avant clipping — **#480 reste ouvert** |
| N10 Maillages difficiles (non ortho, skew, aspect ratio) | [01](docs/audits/01-diffusion-laplacien.md), [06](docs/audits/06-gradients-reconstruction.md) | ~ | ✓ (Poisson polyédrique #852) | ✓ (Kuhn ; L2 pondérée + vrai résidu ≤ 1e-8) | ✓ (Kuhn) | Familles skew/non-ortho/high-aspect non couvertes — **#639 reste ouvert** |
| M0 `.case_hash’/`.geometry_hash’ dans le schéma HDF5 | [08](docs/audits/08-codes-de-reference.md) | ~ (écrits : Python SHA-256 dans case_io.py, C++ FNV-1a dans case_hdf5_io.cpp, rejet fail-closed à la lecture) | ~ (unitaires ; round-trip end-to-end manquant) | ✗ | ✗ | Évidence round-trip/end-to-end + cohérence des deux systèmes de hash (SHA-256 Python vs FNV-1a C++) à documenter : le lecteur C++ accepte à ce jour les fichiers Python sans hashes FNV — décision à prendre |

## Phase 2 — Dette et PR en attente (audits du 2026-10-12)

| PR | Constat d'audit | Décision exigée |
|---|---|---|
| #828 (probes CSV) | ✓ Durcit les gardes (rejet itérations dupliquées, métadonnées partielles, monotonie stricte) ; tolérances intactes (1e-12/1e-8) ; tests de régression ajoutés. Scope honnête (pas de qualification solveur revendiquée). | **Fusionner si CI vert** sur le HEAD exact |
| #571 (MMS transient) | ⚠️ Change T=0.05→0.02 « pour réduire le runtime CI » — en tension directe avec la règle §5 de `docs/development-and-validation-rules.md’ (« ne pas réduire l'horizon physique uniquement pour économiser du temps CI »). Aucune tolérance changée. | **Ne pas fusionner en l'état** : exiger la démonstration que seuls les subtests d'ordre spatial sont affectés et que les mêmes gates passent à T=0.05 ET T=0.02 dans un run documenté ; sinon fermer |
| #511 (gradient N2) | ✓ Renforce réellement les gates (ordre ≥ 1.8, vrais ratios 6/4 et 8/6, rejet de stencil déficient, diagnostics) ; ⚠️ branche en retard sur master ; seuillage pivot vs estimateur de rang SVD/QR assumé différé. | **Rebaser sur master, re-auditer contre la fiche 06, relancer la campagne exact-HEAD** ; ne pas fermer |

## Phase 3 — Multiphysique (aucune entrée promue avant verrouillage Phase 1)

| Domaine | Fiche méthode | Benchmark exigé | Statut |
|---|---|---|---|
| Turbulence RANS + lois de paroi | [09](docs/audits/09-turbulence-rans.md) | canal turbulent, couche limite plate (`TUR-CHANNEL’, `TUR-FLATPLATE’) | BLOCKED (Phase 1) — fiche rédigée, à approuver |
| Thermique/CHT multi-régions | [10](docs/audits/10-thermique-cht.md) | bilans d'énergie à l'interface fluide–solide | BLOCKED (Phase 1) — fiche rédigée, à approuver |
| Rayonnement couplé | [11](docs/audits/11-rayonnement.md) | cavité chauffée, bilans S2S/P1/Rosseland | BLOCKED (Phase 1) — fiche rédigée, à approuver |
| GPU/OOC | [12](docs/audits/12-parallelisme-gpu-ooc.md) | équivalence CPU/MPI/GPU 1e-12, OOC au-delà de la VRAM sans fallback | BLOCKED (Phase 1) — fiche rédigée, à approuver |

## Phase 4 — Avancé (VOF/PLIC, maillage dynamique, FSI, GUI) — BLOCKED (Phase 3)

| Domaine | Benchmark exigé | Référence bibliographique |
|---|---|---|
| VOF PLIC multidimensionnel + tension superficielle | conservation stricte de la fraction volumique, batteries d'interface | Rider & Kothe 1998, DOI 10.1006/jcph.1998.5906 ; Basilisk/Popinet 2003, DOI 10.1016/S0021-9991(03)00298-5 |
| Maillage dynamique | transfert conservatif des champs, topologie | — |
| FSI partitionné (Aitken) | benchmarks FSI couplés | Küttler & Wall 2008, DOI 10.1007/s00466-008-0255-5 |
| GUI/CLI/API | cycle de vie utilisateur complet | — |

## Phase 5 — Certification

- [ ] Comparaison ligne à ligne `spec.md’ ↔ ce registre.
- [ ] Campagne de reproductibilité cross-platform (Linux CPU, Linux GPU, WSL2).
- [ ] Rapport final, fermeture #34, tag `v1.0.0’ (ou `v0.9.0-beta’).

---
*Mis à jour : 2026-10-12 — création (Phase 0 du Plan Directeur) ; révision Phases 1–3 : pilote de campagne N9 (scripts/n9_physical_campaign.py + workflow n9-campaign), fiches méthode 09–12, état M0 documenté (PR empilée sur governance/constitution-technique). Prochaine revue : à chaque fusion sur le noyau.*
