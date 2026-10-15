# FICHE MÉTHODE — Restart/budget du solveur PCD imbriqué (chemin COUPLED+PCD, Ghia N8-PCD)

> Fiche 18. Créée à l'occasion de l'audit du blocage Ghia 32×32 COUPLED+PCD **après** le correctif #898 (fiche 16) : avec le restart couplé complet, le plateau a persisté, et le défaut se situe dans les solveurs imbriqués — exactement le cas annoncé par la fiche 16 §7 (« si le plateau persiste avec le restart complet, le défaut est alors dans la qualité de l'approximation PCD (solveurs imbriqués, tolérances internes) — audit suivant »).

**En-tête d'évidence (gouvernance §16.5)**
- Scope : solveur linéaire **imbriqué** Fp (convection-diffusion pression) de l'action PCD `-Fp⁻¹` (`pcd_schur.h`), invoqué par `CoupledBlockSchurAMGPreconditioner::apply` via `solve_fp` avec les contrôles `pcd_convection_diffusion_*` ; qualification N8-PCD (Ghia, Re=100, 32×32 et préflight 16×16). Non couvert : l'ordre Fp vs Kp (note PROVISIONAL de `pcd_schur.h` — la qualification sur cas convectif non nul reste ouverte), les chemins SIMPLE/SIMPLEC, les tolérances et gates de qualification (inchangés), Kp (`pcd_laplacian_*`, jamais invoqué par l'action single-inverse).
- Statut : DRAFT (2026-10-15)
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (INC-GHIA) ; NUMERICAL_METHOD_CAPABILITY_MATRIX.json (pcd : partial — non avancé par cette fiche)
- Issue / PR liées : #870, #908, #898 (fiche 16), run de diagnostic 38743667702 (preuves)
- Dernier commit vérifié : correctif porté par la PR #902 (côté test uniquement)
- Limitations connues : verdict quantitatif 32×32 en attente du run de diagnostic / validation-total au HEAD portant le correctif (test derrière CFDX_ENABLE_LONG_VALIDATION).

## 1. Références sources (obligatoire : DOI ou ISBN)
- Y. Saad, M. H. Schultz, « GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems », SIAM J. Sci. Stat. Comput. 7(3), 1986, pp. 856–869 — DOI: 10.1137/0907058 (restart et troncature de l'espace de Krylov ; le taux réductible par cycle est borné par un polynôme de degré m sur le spectre préconditionné).
- Y. Saad, « A Flexible Inner-Outer Preconditioned GMRES Algorithm », SIAM J. Sci. Comput. 14(2), 1993, pp. 461–469 — DOI: 10.1137/0914028 (FGMRES tolère un préconditionneur variable **entre** itérations, mais chaque application reste une application d'opérateur : une inverse interne résolue à précision fixe best-effort fait de l'action un bruit non linéaire pour l'externe, qui planche au-dessus de ce bruit).
- D. Silvester, H. Elman, D. Kay, A. Wathen, « Efficient preconditioning of the linearized Navier–Stokes equations for incompressible flow », J. Comput. Appl. Math. 128(1–2), 2001, pp. 261–274 — DOI: 10.1016/S0377-0427(00)00515-X (approximation PCD du complément de Schur).

## 2. Formulation continue
Système couplé point de selle `A = [Auu G ; D C]`, action PCD `S⁻¹ ≈ -Fp⁻¹` avec Fp l'opérateur pression convection-diffusion (une seule inverse ; contrat `pcd_schur.h`). Le solveur externe (FGMRES, restart complet 4096, contrôles dédiés `coupled_gmres_controls` — #898) applique ce préconditionneur à chaque itération : chaque application doit être une application **d'opérateur linéaire déterministe** pour que l'externe puisse converger à sa gate. Une inverse Fp best-effort (jamais CONVERGED) viole ce contrat et impose un plancher de résidu externe ≈ au bruit de l'interne.

## 3. Discrétisation Volume Fini
Maillage cavité 32×32×1 (1024 cellules ; une seule couche en z → pas de faces z internes ; stencil effectif 5 points). Fp = convection upwind **bornée** (flux de masse ρU·S, Pe_cell ≈ 3 à Re=100) + diffusion ρν = 0.01, ligne unité de jauge à la cellule de référence 528 (`impose_pcd_reference_row` ; les opérateurs tout-Neumann sont singuliers en constant, la jauge les rend non singuliers).

Chaîne de solveur interne (statu quo) : FGMRES + ILU0, requête `pcd_convection_diffusion_linear_solver` (restart hérité 40), budget `pcd_convection_diffusion_max_iterations = 100`, tolérance relative 1e-8, contrôles `KrylovControls{}` **par défaut** (clamp [10, 512], adaptatif actif, remplacement de résidu période 50). Le budget de 100 itérations est donc consommé en 2–3 cycles tronqués, et la politique adaptative **rétrécit** le restart (40 → 35 → 30) quand un cycle réduit bien — complice de la troncature.

## 4. Traitement des frontières
Cas tout-parois (lid entraîné) : toutes les frontières de pression sont Neumann (zéro-gradient) ; la jauge est la ligne unité à la cellule de référence 528, imposée aux opérateurs de pression PCD eux-mêmes (pas d'épinglage d'entrée/sortie de l'action — défaut déjà documenté et corrigé dans `pcd_schur.h`). Les flux de masse de bord sont tous nuls (lid parallèle à ses faces normales, parois fixes). La preuve du jumeau montre que le plancher n'est **pas** un défaut de jauge ni de frontière : Fp converge en ~55 its non tronqué, la même matrice cale tronquée.

## 5. Pièges connus (Gotchas)
- Signature du défaut : la **première** résolution couplée passe (flux de masse nul → Fp ∝ Kp, bénin, `COUPLED_LINEAR_TRUE_RESIDUAL relative=9.70666e-11` au run 38743667702) ; le défaut n'apparaît qu'à la deuxième itération non linéaire, quand le flux de masse devient non nul. Un quick/preflight à flux nul ne peut pas le détecter.
- Chaque apply best-effort porte une erreur relative ~1e-6..1e-8 dans l'action : l'externe (même FGMRES flexible) planche au-dessus du bruit de son préconditionneur ; les retries de relaxation reproduisent le même comportement (les alpha ne changent pas le système linéaire).
- `KrylovControls{}` par défaut : `restart_max = 512` borne tout restart demandé ; l'adaptatif rétrécit le restart après un bon cycle. Toute requête imbriquée doit être consciente de ces contrôles.
- Ne pas conclure « PCD inadéquat » avant d'avoir éliminé la troncature (fiche 16 §5 : éliminer le restart d'abord ; c'est maintenant fait pour l'externe **et** l'interne).
- Diagnostic de production : `COUPLED_NESTED_PCD_BEST_EFFORT operator=Fp` n'est imprimé que si `coupled_matrix_summary` est actif ET si le status ≠ CONVERGED — son absence massive est la preuve de solves internes propres.

## 6. Ordre de précision attendu et validation exigée
- Sans objet pour l'ordre spatial : correctif de configuration Krylov, aucune discrétisation touchée.
- Validation exigée : contrat N8-PCD inchangé — gate couplée `coupled_linear_tolerance=1e-10`, gates profil Ghia (max 0.15 / RMS 0.075 pour 32×32), continuité et résidu momentum 1e-7, préflight 16×16. Budget interne porté à 200 (besoin mesuré ~55, marge ×3.6 ; les solves convergés sortent tôt — le budget n'est dépensé que dégradé, jamais rallongé sur un solve sain).
- Risque résidual explicite : la qualité PCD elle-même (ordre Fp vs Kp) sur convection non nulle n'est **pas** qualifiée par cette fiche (note PROVISIONAL de `pcd_schur.h`). Si l'externe cale encore avec l'interne propre, l'audit suivant est cette qualification — la courbe externe du run de diagnostic tranchera.

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : `tests/validation/test_ghia_cavity.cpp` (politique restart/budget du solveur imbriqué Fp, miroir du bloc #898 pour les modes `--pcd`/`--pcd-preflight`). `steady_incompressible_solver.h` NON modifié (un plancher serveur écraserait un choix explicite — raison déjà actée par #898) ; `pcd_schur.h` non modifié ; Kp non modifié (jamais invoqué par l'action single-inverse).
- Preuve : jumeau TypeScript FP64 (pattern #872) — réplique exacte de l'assemblage Fp (`make_mass_flux` + `assemble_scalar_equation` upwind borné + jauge), de l'ILU0 (`advanced_preconditioners.h`) et de `solve_gmres` (contrôles par défaut). 3 champs de vitesse proxy du U¹ réel × 5 secondes de RHS : statu quo MAX_ITER_REACHED à 3e-7..2e-5 relatif (la bande de production 1.4e-8..2.9e-6 est couverte) ; restart ≥ budget → CONVERGED en 53–56 its, 15/15.
- Écarts connus code ↔ fiche : le jumeau utilise des champs de vitesse proxy (le U¹ réel n'est pas disponible hors run) ; verdict quantitatif attendu du run de diagnostic au HEAD portant le correctif. Si le 32×32 passe, lever DRAFT.
- Gaps et prochaines actions : qualification Fp vs Kp sur convection non nulle (PROVISIONAL) ; statut pcd du NUMERICAL_METHOD_CAPABILITY_MATRIX.json inchangé (partial) jusqu'au verdict.

## 8. Journal des modifications de la fiche
- 2026-10-15 — issue #908 / PR #902 — création, à partir des preuves du run de diagnostic 38743667702 et du jumeau FP64.
