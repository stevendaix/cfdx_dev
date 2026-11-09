# CFDX V13 — BFS : noyaux bilatéraux Fp, état global transactionnel et triplet de paroi SST

**Date :** 2026-11-09
**Base reproductible :** sources complètes CFDX V12, elles-mêmes dérivées de `cfdx_dev-master (2)(1).zip` (SHA-256 `c3a0d40206829596dca61871c10f6079724c6d0e11bbf63a47d4e894e914e62c`).
**Statut :** PREUVES ALGÉBRIQUES CIBLÉES + CONSTRUCTION DU SOLVEUR DE PRODUCTION ; BFS Driver–Seegmiller NON QUALIFIÉ.

## 1. Priorités et hypothèses testées

1. Différencier le noyau droit `ker(Fp)` du noyau gauche `ker(Fp^T)` pour un auxiliaire PCD non symétrique.
2. Quantifier la part incompatible du RHS via `b_incompatible = W (W^T b)` sans la cacher derrière la convergence projetée.
3. Raccorder un garde optionnel au résidu physique indépendant `U,p,k,omega,nu_t` réassemblé APRÈS le transport SST.
4. Séparer les résidus de candidat et les états réellement commis par une transaction complète.
5. Mesurer conjointement les références de face `omega,nut,G` et les valeurs des cellules de paroi, y compris les cellules partageant plusieurs faces.

Références : PETSc KSP <https://petsc.org/main/manual/ksp/> (sous-espaces gauche/droit d'un système singulier) ; OpenFOAM backwardFacingStep2D <https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/> (conditions aux limites SST et mesures de Driver–Seegmiller).

## 2. Certificat algébrique des captures *réelles* Fp

Deux captures `CFDXFP01` de 712 cellules, 3 442 coefficients CSR, ont été mesurées par une SVD dense LAPACK bornée aux petits systèmes ; les bases sont stockées dans `CFDXNB13` avec un FNV64 du CSR complet. **Aucune projection n'est installée dans le FGMRES COUPLED de production.**

| Capture | `dim ker(Fp)` | `dim ker(Fp^T)` | `||W W^T b||/||b||` | `||Fp x - b_compatible||/||b||` | `||Fp x - b||/||b||` |
|---|---:|---:|---:|---:|---:|
| Fp initial | 2 | 2 | 0.0627946816675 | 5.45e-13 | 0.0627946816675 |
| Fp avec Dirichlet d'entrée expérimental | 1 | 1 | 0.0627946816675 | environ 5e-13 | 0.0627946816675 |

Le nombre de modes est défini numériquement au seuil SVD relatif `1e-12` et certifié par les défauts `||Fp Z||/||Fp||_F` et `||Fp^T W||/||Fp||_F` de l'ordre de `1e-16`. Il ne s'agit pas d'une preuve d'une singularité *exacte*, indépendante de la précision flottante.

**Interprétation mathématique :** le second membre capturé n'appartient pas à l'image résolue de Fp au seuil choisi. La meilleure solution projetée n'annule pas le résidu original, dont le plancher est ~6,28 %. Déclarer Fp convergé à `1e-8` à partir du seul résidu projeté serait faux. Cette incompatibilité ne prouve pas encore que la cause racine est exclusivement la condition limite ; il faut auditer flux Fp/Kp, forme de Schur, interpolation et connectivité des blocs flottants.

Le lecteur C++ `pcd_fp_bisided_basis_v13.h` refuse un changement de coefficient même quand le graphe CSR est identique, et certifie séparément les deux bases. Son calcul de projection est un **composant de recherche réutilisable**, non branché aux équations physiques. `--require-compatible` retourne **3** sur les captures incompatibles (échec attendu).

## 3. Acceptation indépendante de l'itération externe

Le résidu `GLOBAL_NONLINEAR_V8` est réévalué sur l'état complet après SST. L'option *désactivée par défaut* `v13_joint_globalization=true` exige : solveur stationnaire, SST, audit global V8, réassemblage momentum après SST et budget de retries non nul. Pour une méthode et un opérateur de convection inchangés, l'état candidat n'est accepté que si aucun des six canaux (`U` maximal, `U` global, continuité, `k`, `omega`, fermeture `nu_t`) ne dépasse l'enveloppe d'accroissement `1.20*max(r_precedent,10*tolerance_physique)`. Les **critères finaux** ne sont ni modifiés ni contournés ; l'enveloppe de garde n'est PAS un nouveau seuil de convergence.

En cas de refus, le mécanisme existant `NonlinearRetryableFailure` restaure le candidat, y compris U, p, flux et SST. Ce n'est ni un Jacobien monolithique Newton–Krylov, ni une preuve de contraction globale.

**Correction de preuve de V13 :** `GLOBAL_NONLINEAR_V8` est un état *provisoire*. Une nouvelle ligne `GLOBAL_NONLINEAR_COMMITTED_V13` est émise seulement APRÈS `transaction.commit()`. Le validateur strict `scripts/qualify_bfs_committed_v13.py` distingue ces deux événements, vérifie la cohérence de leurs six champs et refuse toute attribution de qualification physique à un simple commit.

**Rejeu production BFS 712 cellules, Block-local, ordre deux, SST 800 :** une première itération a calculé FGMRES extérieur 139 itérations, SST convergé après 213, puis a commis un état NON convergé : `U=0.5275372423`, `p=0.001760025754`, `k=1.9203796e-13`, `omega=4.3727723e-14`, `nut_closure=0`; convergence globale `false`, statut du solveur `exit=1` après une itération demandée. Les tours suivants n'ont **pas** fermé les équations dans les budgets précédents ; aucune amélioration A/B globale n'est revendiquée.

## 4. Parois : audit collectif omega, nut, production G

Les nouveaux diagnostics dans `sst_openfoam_wall_coupling_v8.h` regroupent les contributions face par propriétaire et calculent la référence pondérée par aire pour le triplet. Une cellule reliée à deux faces n'est comptée qu'une fois dans l'inventaire propriétaire. Sur la première itération du BFS : **70 faces, 69 cellules, 1 cellule multi-face**, écart maximum omega propriétaire ~`0.99767`, écart maximum nut proxy propriétaire ~`0.99983`, écart G proxy initial `0`. Le proxy cellulaire n'est pas une condition limite de face équivalente ; ne pas attribuer directement ces écarts à une erreur physique.

Le traitement simultané de `omegaWallFunction`, `nutUSpaldingWallFunction` et la production `G` **n'a pas été implanté** dans le SST C++ ; ni `kqRWallFunction`, ni les règles de pondération de coin, ni les contraintes patch ne sont encore reproduites au solveur. Aucune concordance OpenFOAM n'est revendiquée.

## 5. Preuves et couverture

- **Compilation et édition de liens :** `cfdx_production_solver` Release diagnostique avec `-O0` uniquement sur `main.cpp` (non un benchmark optimisé).
- **5/5 CTest ciblés :** Schur PCD, ancrage V12, SVD/C++ bilatéral V13, garde global V13, audit paroi triplet V13.
- **99 Python réussis, 6 ignorés** : BTS BFS/PCD, SVD, validation committed, certificats existants.
- Les C++ ont rejoué les **deux captures Fp d'origine**, contrôlé bases gauche/droite et refusé une base périmée après changement des coefficients.
- Le script strict refuse les états provisoires, les données non finies, les champs incohérents et les duplications de commit.
- La campagne globale BFS de 20 540 cellules et les six Schur sur temps longs n'a **pas** été achevée en V13. `x_r/h`, `Cf`, profils de vitesse et contraintes de Reynolds n'ont pas de valeurs CFDX certifiées.
- **Non exécutés :** CI complète, DebugSanitizers, MPI, GPU, trois maillages OpenFOAM natifs, comparaison expérimentale NASA, Newton–Krylov monolithique de production, Jacobien exact de SST.

## 6. Méthode de reprise

1. Vérifier tous les résidus originaux Fp et le rang après chaque reconstruction de l'opérateur ; comparer la condition gauche/droite en fonction des patches et des composantes flottantes.
2. Corriger l'incompatibilité `b` / espace colonne de Fp à la **source de l'assemblage auxiliaire**, en démontrant l'identité des flux et des conditions limites, pas en cachant `W^T b` par projection.
3. Raccorder un Schur robuste de production (PCD/MGR/AMG et variantes) et tester le vrai résidu du COUPLED avec norme indépendante à chaque pas.
4. Ajouter une globalisation fondée sur les résidus complets qui démontre la diminution de `U,p,k,omega` sur les états **commis**, et non seulement un rejet sous une enveloppe heuristique.
5. Implémenter, valider séparément puis conjointement les conditions aux limites SST (`omega`, `nut`, production) avant l'étude expérimentale sur trois maillages natifs.

Aucun seuil de convergence n'a été assoupli et aucun push GitHub n'a été effectué.
