# FICHE MÉTHODE 10 — Benchmark A/B de l'action de Schur PCD (système couplé Ghia)

**En-tête d'évidence**
- Scope : mesure isolée, reproductible et peu coûteuse du coût et du comportement linéaire d'actions de Schur PCD alternatives sur les systèmes couplés pression-vitesse **réellement assemblés** par la qualification Ghia (cavité entraînée 32×32, Re=100). Benchmark de diagnostic : ni changement du préconditionneur de production, ni qualification physique.
- Statut : PROVISIONAL (2026-10-19) — collecte d'évidence benchmark uniquement ; n'autorise aucune modification de production (`pcd_schur.h`, défauts du solveur) sans décision explicite post-comparaison quantitative sans régression physique.
- Issues liées : #914 (audit Fp vs Kp), #390 (robustesse convergence) ; PR #966 ; consolidation #929 ; axe mesuré #965.
- Invariants : aucun budget, oracle, tolérance ou gate modifié ; `test_ghia_cavity_n8_pcd` reste l'unique juge du cas physique.

## 1. Objectif et références

Départager, sur les VRAIS systèmes couplés de la production, le **nombre d'applications de l'action de Schur** de leur **coût unitaire** (construction/setup vs application), comme demandé sur #914, sans nouvelle campagne longue.

- Kay, D., Loghin, D., and Wathen, A. J. (2002). « A preconditioner for the steady-state Navier–Stokes equations. » *SIAM Journal on Scientific Computing*, 24(1), 237–256. DOI : 10.1137/S1064827500384937. Référence de fond sur le préconditionnement pressure-convection-diffusion ; elle ne valide ni l'action CFDX particulière ni les résultats mesurés ici.

## 2. Formulation continue

Système d'Oseen couplé (linéarisation de Navier–Stokes incompressible stationnaire) :

$$ \begin{bmatrix} F & B^T \\ B & 0 \end{bmatrix} \begin{bmatrix} u \\ p \end{bmatrix} = \begin{bmatrix} f \\ g \end{bmatrix} $$

Le complément de Schur en pression $S = -B\,F^{-1}B^T$ n'est jamais formé : le préconditionneur applique une approximation de $S^{-1}$. PCD approxime $S \approx -F_p$ où $F_p$ est un opérateur de pression (masse, laplacien et convection-diffusion de pression) construit sur la vitesse de linéarisation ; l'action de production applique $-F_p^{-1}$ (choix **PROVISIONAL** documenté dans `pcd_schur.h`, jamais qualifié sur un flux massique non nul — Ghia est précisément ce cas). Alternatives mesurées : $-\mathrm{diag}(RC)\,F_p^{-1}$ (mise à l'échelle par la diagonale Rhie–Chow), action diagonale $-\mathrm{diag}(RC)^{-1}$ (référence SIMPLE-diag), et substitution de l'opérateur alternatif $K_p$ (variante PROVISIONAL documentée).

## 3. Discrétisation Volume Fini

Maillage cartésien 32×32, volumes finis colocalisés, couplage Rhie–Chow. Le solveur de production assemble le système monolithique et le résout par FGMRES couplé, préconditionné par `CoupledBlockSchurAMGPreconditioner` (factorisation Full, approximation de vitesse par blocs, approximation de Schur PCD via `PcdSchurApproximation` : apply = un solve $F_p$ imbriqué puis multiplication par $-1$). Chaque itération Krylov couplée coûte donc un solve imbriqué $F_p$ (production : FGMRES(512)+ILU0, ~53–56 itérations, budget 200, tol 1e-8) — chaîne de coût établie par le diagnostic : 874,7 s vs budget ctest 720 s, contre 5,6 s / 497 itérations pour SIMPLE sur le même cas.

## 4. Traitement des frontières

Cavité de Ghia : lid Dirichlet $U=(1,0)$, parois sans glissement, pas de pression imposée. Jauge de pression : la ligne de la cellule de référence est imposée à l'identité dans la matrice couplée par l'assembleur (avant finalize), la colonne correspondante est éliminée du second membre, et la valeur de référence est injectée au RHS et à l'état initial. Si une pression est fixée (`has_fixed_pressure`), la jauge ne s'applique pas ; le benchmark distingue explicitement les deux cas et vérifie structurellement la convention du système capturé.

## 5. Pièges connus

- Une action plus rapide sur les étages capturés est une évidence sur le coût du solveur linéaire, pas sur la convergence non-linéaire ni la justesse physique.
- Une action non-baseline qui ne converge pas est enregistrée comme mesure échouée, jamais substituée ni traitée comme PASS.
- L'état initial de production n'est **pas nul** (x0 = [U_old ; p_old], composante de jauge = valeur de référence) : le rejeu doit utiliser le x0 capturé, sinon la comparaison n'est ni fidèle ni équitable.
- Le point de capture doit être l'appel Krylov exact (après traitement de jauge et setup du préconditionneur) : capturer en amont ne garantit la fidélité que par lecture de code, pas par construction.
- Une reconstruction indépendante des opérateurs PCD ne suffit pas : le rejeu utilise les opérateurs de production capturés ; la reconstruction n'est qu'une garde (cross-check ≤ 1e-12), jamais une source.
- Les conventions de jauge doivent être identiques dans la matrice, le second membre, l'état initial, le préconditionneur et l'action de Schur.
- Aucun résultat de ce benchmark n'autorise une modification de l'action PCD de production : un suivi explicite, une qualification couplée complète et l'absence de régression d'oracle physique sont requis.
- Le workflow branche-only est une collecte d'évidence, pas destiné au merge en l'état ; l'enregistrement CTest et le hook read-only d'assemblage exigent une revue séparée avant toute intégration.

## 6. Protocole de mesure

1. Capturer, au point exact d'appel Krylov (hook additif nul par défaut, zéro overhead), la matrice couplée, le second membre, l'état initial de production, la diagonale Rhie–Chow de Schur pression (épinglée et brute, copiée avant l'épinglage) et les opérateurs PCD de production (masse, laplacien, convection-diffusion), sur 5 itérations non-linéaires du VRAI solveur (config miroir exacte de `--pcd`, budget non-linéaire réduit à 5 car la capture est le seul objectif).
2. Rejouer les étages #2 et #5 à travers la baseline et les actions alternatives avec la même matrice, le même second membre, le même état initial, les mêmes contrôles Krylov externes (FGMRES restart min(4096, n), budget 5000, tol 1e-10, `KrylovControls` fixe) et la même tolérance.
3. Recalculer indépendamment le résidu vrai ||b − A x|| / ||b||. Rapporter séparément : itérations couplées, applications de l'action de Schur, itérations et temps des solves imbriqués, temps d'action vs temps non-action (produits matrice-vecteur et orthogonalisation Krylov), temps setup imbriqué vs setup Schur, temps total, et temps de reconstruction des opérateurs.
4. Gardes de validité : la baseline doit reproduire la convergence linéaire de production (résidu vrai ≤ 1e-10) sur chaque étage ; les variantes convergées doivent s'accorder sur la solution à 1e-6 relatif ; les opérateurs reconstruits doivent matcher les capturés à 1e-12 ; la ligne de jauge doit être structurellement l'identité (1 nnz, valeur 1.0) avec RHS et x0 cohérents ; le slot `laplacian_solve` de l'action single-inverse reçoit un callback qui throw (toute invocation du mauvais opérateur échoue bruyamment).
5. Séparer le verdict de **validité** (les gardes ci-dessus, évaluées après le bilan complet) du rapport de **performance** (mesures agrégées par variante, jamais gating).

## 7. Statut d'implémentation

- Sources : `src/cfdx/physics/steady_incompressible_solver.h` (hook `CoupledAssemblySnapshot` + `controls.coupled_assembly_callback`, additif, nul par défaut) ; `src/cfdx/core/linalg/pcd_schur.h` **NON MODIFIÉ**.
- Test : `tests/validation/test_pcd_schur_action_benchmark.cpp`, enregistré uniquement sous `CFDX_ENABLE_LONG_VALIDATION`, TIMEOUT ctest 1200 s.
- Workflow : `.github/workflows/pcd-action-benchmark-validation.yml` (branche-only, pas destiné au merge en l'état).
- Variantes mesurées : `fp_ilu0` (baseline production), `fp_amg` (axe #965), `kp_amg` (PROVISIONAL), `fp_rc_scaled`, `rc_diag` (référence SIMPLE-diag).
- Traçabilité : PR #966, #965, #914, #390 ; consolidation #929.
- Journal : 2026-10-19 — création (f2207d12/7c2465dd), restructurée au format des fiches méthode lors de la réponse à l'audit de la PR #966.
