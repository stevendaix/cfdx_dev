# FICHE MÉTHODE 12 — Parallélisme (MPI, GPU) et Out-of-Core — contrat de parité (Règle d'Or 4)

**En-tête d'évidence**
- Scope : contrat de reproductibilité CPU séquentiel vs MPI (décomposition de domaine, échange de halos) vs GPU (CUDA), redémarrage N->M, et mode Out-of-Core (OOC) au-delà de la VRAM. C'est la traduction exécutable de la **Règle d'Or 4** : zéro divergence parallèle, aucun fallback silencieux.
- Statut : DRAFT (2026-10-12) — contrat à faire approuver dans la PR ; le code N13/Phase 6-7 existant doit être audité contre cette fiche avant toute promotion.
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json ; politique : docs/validation/N13_EQUIVALENCE_POLICY.json ; preuves : docs/validation/N13_CROSS_BACKEND_EQUIVALENCE.md
- Issues liées : #34 (roadmap, Phases 6-7)
- Dernier commit vérifié : f5a728f205cb45ad6186c9fab92554ae4af06bdd (master, 2026-10-12 — présence des sources et tests confirmée)
- Limitations connues : équivalence non exécutée depuis l'environnement d'audit ; statut déclaratif jusqu'à campagne sur HEAD exact.

## 1. Références sources
- U. Trottenberg, C. W. Oosterlee, A. Schuller, *Multigrid*, Academic Press, 2001 — ISBN : 978-0-12-701070-0. (analyse de la convergence indépendante du nombre de processus, halo/garde-fou des schémas)
- MPI Forum, *MPI: A Message-Passing Interface Standard* (versions 3.1 et suivantes) — standard cité par désignation (pas de DOI : document normatif).
- NVIDIA, *CUDA C++ Programming Guide* (documentation de conception kernels/warps — citée par désignation).
- Contrat interne : docs/validation/N13_EQUIVALENCE_POLICY.json (tolérances et politiques) + docs/validation/N13_CROSS_BACKEND_EQUIVALENCE.md (preuves existantes).

## 2. Formulation continue
Aucune équation nouvelle : la physique est celle des fiches 01 à 11. Cette fiche définit le **contrat d'invariance** : pour un cas donné (même maillage exact, même configuration numérique déclarée, mêmes initialisation et critères d'arrêt), les champs produits par les exécutions séquentielle CPU, MPI (P processus, décomposition de domaine quelconque) et GPU doivent satisfaire :
$$\|\phi_{backend}-\phi_{ref}\|_{\infty}\le 10^{-12}\times\text{échelle du champ},\qquad \text{ou reproductibilité bit-à-bit documentée}$$
Le choix (bit-à-bit vs tolérance stricte 1e-12) doit être **déclaré par module** dans la politique N13, jamais implicite.

## 3. Discrétisation Volume Fini
Contraintes sur l'implémentation (pas de reformulation, mais des exigences) :
- **Sommes et réductions** : l'ordre de sommation (produits scalaires, résidus, normes, bilans de conservation) doit être déterministe — réduction reproductible (ordre fixe ou arbre fixe), car une réduction flottante non déterministe entre backends détruit la parité 1e-12.
- **Halos MPI** : échange de toutes les variables dépendantes avant chaque assemblage ; les coefficients de face aux interfaces de décomposition doivent être calculés de façon identique des deux côtés (une seule formule, aucune variante par côté).
- **GPU** : kernels sans atomiques non déterministes (ou avec atomiques documentés comme bit-à-bit stables) ; layouts mémoire identiques aux champs CPU pour la comparaison exacte ; pas de fusion de passes qui changerait l'ordre des opérations flottantes.
- **OOC** : le débordement VRAM -> RAM/disque doit préserver l'ordre des opérations exactement (streaming sans réordonnancement flottant) ; chaque page/tile spillée est télémetrée.

## 4. Traitement des frontières
- Les conditions aux limites (fiches 01-11 §4) s'appliquent à l'identique sur tout backend ; les patches partitionnés par la décomposition doivent retrouver exactement les mêmes coefficients que l'exécution mono-domaine.
- Communication : frontières de décomposition != frontières physiques ; aucune CL physique ne doit être appliquée par erreur sur une interface MPI (et réciproquement).
- Redémarrage N->M : le champ relu sur un maillage/redécoupage différent doit reproduire le champ de référence à la tolérance déclarée (test test_n13_restart_equivalence).

## 5. Pièges connus
- **Réduction non déterministe** : premier suspect de toute divergence 1e-12 ; exiger des réductions reproductibles, documenter l'ordre.
- **Fallback silencieux interdit** (Règle d'Or 4) : si MPI ne s'initialise pas, si le GPU est absent, ou si la VRAM est insuffisante sans OOC activé -> **échec fort avec message clair**, jamais de repli discret sur un autre backend ; l'activation OOC doit être explicite et télémetrée (chaque spill visible dans le rapport).
- Différence de padding/alignement GPU qui change l'ordre des sommes (warps) — geler les layouts.
- Précision mixte (FP32 vs FP64) interdite sans déclaration explicite et preuve que la tolérance 1e-12 tient ; par défaut FP64 partout.
- Halos incomplets (oubli d'une variable auxiliaire : gradients, flux de Rhie-Chow) — parfaite illustration du checkerboarding parallèle ; liste des variables échangées à documenter par module.
- L'arithmétique FMA/contractions du compilateur GPU peut différer du CPU : geler les drapeaux de contraction (ou documenter bit-à-bit impossible et se replier sur la tolérance 1e-12 déclarée).

## 6. Ordre de précision attendu et validation exigée
- Aucun ordre propre : l'exigence est l'**équivalence** sur le même maillage exact, à la tolérance N13 (1e-12) ou bit-à-bit déclarée.
- Oracles : tests d'équivalence dédiés sur cas représentatifs (Poisson, Poiseuille, cavité) : test_n13_mpi_equivalence (décomposition P>1), test_n13_cuda_equivalence (GPU), test_n13_restart_equivalence (N->M) ; OOC : test_phase7_ooc_validation (au-delà de la VRAM, sans fallback, avec télémetrie de spill) ; test_phase8_cuda_poisson (backend CUDA sur l'opérateur pivot).
- Conservation : les bilans (masse, énergie) doivent tenir à la tolérance machine **sur chaque backend** (fiche N11 + Règle d'Or 3) ; un bilan qui ne se ferme qu'en séquentiel est un défaut bloquant.
- Campagne : cross-backend sur un même SHA exact, rapport JSON machine-readable, exécutée avant toute promotion des modules distribués (politique N13).

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : sources MPI/GPU/OOC du noyau (src/cfdx/...) — inventaire exact à produire lors de l'audit Phase 3 ; politique : docs/validation/N13_EQUIVALENCE_POLICY.json ; preuves : docs/validation/N13_CROSS_BACKEND_EQUIVALENCE.md.
- Tests exécutables : tests/validation/test_n13_mpi_equivalence.cpp, test_n13_cuda_equivalence.cpp, test_n13_restart_equivalence.cpp, test_phase7_ooc_validation.cpp, test_phase8_cuda_poisson.cpp.
- Écarts connus code <-> fiche : à auditer lors de la Phase 3 — réductions déterministes et fallbacks silencieux sont les deux premiers points de contrôle.
- Gaps et prochaines actions : équivalence GPU sur les solveurs de production (au-delà de Poisson) ; campagne OOC au-delà VRAM avec télémetrie ; reproductibilité cross-platform (Linux CPU, Linux GPU, WSL2) réservée à la Phase 5 ; promotion BLOCKED derrière Phase 1.

## 8. Journal des modifications de la fiche
- 2026-10-12 — création (Phase 3 du Plan Directeur, PR empilée sur governance/constitution-technique).
