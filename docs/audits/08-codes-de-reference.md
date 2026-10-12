# FICHE MÉTHODE 08 — Codes de référence et positionnement produit (SU2, OpenFOAM, Fluent, Basilisk)

**En-tête d'évidence**
- Scope : politique de comparaison croisée avec les codes tiers et positionnement produit de CFDX. S'étend `docs/validation/CFDX_VV_GOVERNANCE.md’ §7 (classification des comparaisons) — n'en duplique pas la terminologie.
- Statut : APPROVED (2026-10-12)
- Issues liées : #34, #461 ; `docs/validation/N14_REFERENCE_CODE_MATRIX.json’, `FLUENT_VMFL_MATRIX.md’, `SU2_NACA0012_READER.md’
- Dernier commit vérifiéà: 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : aucune comparaison croisée n'est un oracle d'acceptation en soi ; c'est un élément de preuve.

## 1. Références sources
- T. D. Economon, F. Palacios, S. R. Copeland, T. W. Lukaczyk, J. J. Alonso, « SU2: An Open-Source Suite for Multiphysics Simulation and Design », *AIAA Journal* 54(3):828–846, 2016 — DOI : 10.2514/1.J053813.
- H. G. Weller, G. Tabor, H. Jasak, C. Fureby, « A tensorial approach to computational continuum mechanics using object-oriented techniques », *Computers in Physics* 12(6):620–631, 1998. (OpenFOAM)
- S. Popinet, « Gerris: a tree-based adaptive solver for the incompressible Euler equations in complex geometries », *J. Comput. Phys.* 190(2):572–600, 2003 — DOI : 10.1016/S0021-9991(03)00298-5 ; Basilisk (évolution de Gerris, basilisk.fr, S. Popinet).
- W. J. Rider, D. B. Kothe, « Reconstructing volume tracking », *J. Comput. Phys.* 141(2):112–152, 1998 — DOI : 10.1006/jcph.1998.5906. (PLIC, référence VOF)

## 2. Rôle de chaque code de référence
| Code | Usage CFDX | Classement gouvernance §7 |
|---|---|---|
| **OpenFOAM** | référence FVM colocalisé incompressible ; comparaisons Poiseuille/Ghia/canal sur le **même maillage exact** | même méthode |
| **SU2** | référence compressible/aéro (NACA0012) ; lecture de solutions existantes | méthode proche, but différent |
| **Fluent/CFX/STAR-CCM+** | oracles VMFL (`FLUENT_VMFL_MATRIX.md’), cible UX | même méthode |
| **Basilisk** | référence VOF/PLIC + adaptativité quadtree (batterie de benchmarks d'interface) | méthode différente |

Règle d'or #2 : comparaison = **même maillage exact**, même problème physique, réglages documentés (schémas, solveurs, critères d'arrêt, QoI). Divergence documentée = donnée utile, pas échec silencieux.

## 3. Positionnement produit (confort d'utilisation)
CFDX vise un **code open source moderne** avec un confort d'utilisation comparable aux codes commerciaux (Fluent) :
- **Workbench** GUI/CLI/API cohérents : setup déclaratif, arborescence physique/solveurs type Fluent, présets (laminaire/RANS/k-ε/SST…) ; erreurs **fail-closed en langage clair** (jamais de NaN silencieux, jamais de valeur par défaut cachée).
- **Fichier de cas unique** `.cfdx.h5’ vérsionné (schema + hashes `case_hash’/`geometry_hash’ — à compléter, M0) ; restart N→M ; probes CSV sûrs (#828).
- **Python API** (`scripts/cfdx.py’) + visualisation PyVista ; rapports de sélection numérique lisibles ; télémetrie conservation/bornitude dans le rapport de cas.
- L'ergonomie ne dégrade jamais la rigueur : chaque preset déclare les tolérances et oracles qu'il revendique (registre de qualification).

## 4. Comparaison croisée — contrat d'exécution
1. Geler version/commit du code tiers, maillage, CL/CI, propriétés, QoI.
2. Exécuter les deux codes avec critères de convergence documentés.
3. Comparer les QoI appariés ; différences absolue/relative ; distinguer l'incertitude de référence de l'erreur CFDX.
4. Classer (même méthode / proche / différente) ; publier la superposition graphique.
5. Toute divergence > seuil déclaré → investigation documentée, pas d'ajustement de seuil a posteriori.

## 5. Pièges connus
- « Même nom d'option ≠ même numérique » (gouvernance §11) : vérifier la sémantique des schémas/solveurs entre codes.
- Basilisk n'est pas un FVM cell-centré classique : comparer des QoI intégraux, pas des champs point à point.
- Les références Fluent VMFL sont héritées de publications : geler la provenance (`FLUENT_VMFL_MATRIX.md’).

## 6. Ordre de précision attendu et validation exigée
Aucun ordre propre — la validation exigée est la campagne croisée elle-même : ≥ 1 cas par physique majeure (incompressible, turbulence, thermique, VOF) contre le code de référence le plus approprié, avec graphique et métriques dans le registre.

## 7. Statut d'implémentation CFDX
- Existant : `docs/validation/N14_REFERENCE_CODE_MATRIX.json’, `SU2_NACA0012_READER.md’, `FLUENT_VMFL_MATRIX.md’, `EXT-NACA0012_QUALIFICATION_CASE.md’.
- Gaps : comparaisons OpenFOAM sur même maillage pas encore déployées comme campagne CI standard ; `case_hash’/`geometry_hash’ manquants (M0) ; script de superposition standardisé à fournir (Partie 5 du Plan Directeur).

## 8. Journal
- 2026-10-12 — création.
