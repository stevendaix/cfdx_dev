# FICHE MÉTHODE 22 — Couplage pression–vitesse et flux Rhie–Chow en formulation COUPLED

**En-tête d'évidence**
- Scope : contrat discret du bloc monolithique incompressible `[M G; D C]`, flux massique aux faces, réponse pression orthogonale et correction différée non orthogonale. Ne couvre pas la qualification physique complète des solveurs SIMPLE/SIMPLEC/COUPLED.
- Statut : DRAFT — correction F03 en cours, PR #1078 ; non approuvée pour qualification.
- Registre : audit des équations #1069, constats F03 / E029 / E033 ; suivi #1076.
- Source principale : `src/cfdx/physics/steady_incompressible_solver.h`.
- État de preuve : revue statique et régression d'opérateur ajoutée ; CI exact-HEAD en cours. Aucune conclusion de convergence physique n'est revendiquée.

## 1. Références sources
- P. M. Gresho, R. L. Sani, *Incompressible Flow and the Finite Element Method*, Wiley, 1998. Référence générale sur le couplage pression–vitesse et les modes de pression.
- C. M. Rhie, W. L. Chow, “Numerical Study of the Turbulent Flow Past an Airfoil with Trailing Edge Separation”, *AIAA Journal*, 21(11), 1525–1532, 1983. DOI: 10.2514/3.8284.
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980. Référence de base pour l'algorithme SIMPLE et le flux aux faces.

## 2. Formulation continue
Pour un fluide incompressible de masse volumique constante :
\[
\nabla\cdot\mathbf U=0,\qquad
\rho\,\nabla\cdot(\mathbf U\otimes\mathbf U)
=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
\]
Le bloc linéarisé monolithique est \([M\;G;D\;C][U,p]^T=[b_U,b_p]^T\). Le flux de continuité doit utiliser la même réponse pression–vitesse que celle représentée par les blocs de couplage. La formulation du flux dépend du sens de la variable passée au reconstructeur : le prédicteur ségrégé \(H/A\) n'est pas la vitesse physique \(U\) déjà résolue par le bloc COUPLED.

## 3. Discrétisation Volume Fini
Pour une face interne orientée de \(P\) vers \(N\), la réponse pression discrète est décomposée en saut orthogonal implicite et correction différée non orthogonale :
\[
R_{p,f}=\rho\left[D_f(p_N-p_P)+R_{f,\mathrm{nonorth}}\right].
\]
Avec le prédicteur ségrégé, le flux standard utilise \(\phi_f=\rho(H/A)_f\cdot S_f-R_{p,f}\). Pour le bloc monolithique, \(U\) est la vitesse physique ; le contrat doit inclure la compensation interpolée \(\rho(dAU\,\nabla p)_f\cdot S_f\) :
\[
\phi_f=\rho U_f\cdot S_f+\rho(dAU\,\nabla p)_f\cdot S_f-R_{p,f}.
\]
Le bloc de continuité doit assembler la même compensation, ici différée au niveau de l'itération de Picard courante, en plus de la réponse orthogonale implicite et du terme non orthogonal différé. Les lignes voisines doivent utiliser le même flux avec des signes opposés. Cette identité doit être vérifiée par comparaison entre résidu matriciel et bilan des flux reconstruits, pas par réévaluation deux fois du même helper.

## 4. Traitement des frontières
- Face interne : interpoler les composantes de \(dAU=V/A_P\) et du gradient selon les règles explicites du solveur ; appliquer un signe opposé dans la ligne du voisin.
- Pression imposée : employer la valeur de pression de frontière pour le saut orthogonal et la correction non orthogonale ; appliquer également la compensation \(dAU_P\nabla p_P\cdot S_f\) si le contrat du flux COUPLED l'active.
- Pression à gradient nul : ne pas inventer de saut de pression imposé ; conserver exactement la même convention de frontière entre le flux et la continuité.
- Patch EMPTY : ne contribue ni au flux ni au bilan de continuité.
- Toute géométrie dégénérée, diagonale non positive ou coefficient non fini doit échouer explicitement.

## 5. Pièges connus
- Confondre la vitesse physique résolue par COUPLED et le prédicteur \(H/A\) du chemin ségrégé.
- Ajouter la réponse pression au flux final sans la compensation correspondante dans le bloc de continuité, ou l'inverse.
- Corriger uniquement le flux post-solve : cela rend le diagnostic différent de la matrice qui a été résolue.
- Utiliser une identité testée par deux appels au même helper comme oracle indépendant.
- Omettre le signe du flux sur la ligne voisine ou le traitement des frontières à pression imposée.
- La compensation est actuellement différée ; la stabilité et la convergence doivent être mesurées, sans ajuster les tolérances ou budgets.

## 6. Ordre de précision attendu et validation exigée
- Oracle orthogonal indépendant : expression directe du flux attendu sur un canal à deux cellules, avec pression constante et affine, et vérification du terme de compensation séparément.
- Oracle affine sur maillage adapté avec conditions de frontière compatibles avec le champ affine ; contrôler que la réponse pression et la compensation s'annulent à l'ordre attendu.
- Mode checkerboard : vérifier que la pression alternée reste couplée au flux ; comparer avec le mode cellule-centre sans correction.
- Conservation : antisymétrie interne, flux de frontière, résidu local par cellule et accord entre résidu du bloc et divergence du flux autoritaire.
- Tester Release et DebugSanitizers sur le HEAD exact ; exécuter le test ciblé N9/S9 et les gates audit/maturité. Toute qualification physique nécessite en plus les cas COUPLED pertinents et leurs métriques indépendantes.
- Aucun relâchement de tolérance, gate, budget d'itérations ou politique de fallback n'est acceptable.

## 7. Statut d'implémentation
- Le chemin COUPLED est modifié dans la PR #1078 pour distinguer le contrat du flux à partir de la vitesse physique de celui du prédicteur ségrégé.
- La compensation \(dAU\nabla p\) est ajoutée au flux autoritaire et au second membre de continuité, y compris les faces à pression imposée.
- Une régression d'opérateur indépendante du helper de flux est ajoutée dans `tests/validation/test_n9_s9_rhie_chow_consistency.cpp`.
- La PR #1078 a depuis été fusionnée. Cela ne ferme pas la qualification F03 : les résultats de validation totale comportent encore des échecs sans lien causal établi avec ce changement.
- Le complément ajoute un indicateur de divergence L1 sur toute la grille pour le mode de pression alternée, comparant le flux corrigé au flux interpolé sans correction. Il mesure la réponse de l'opérateur, pas l'amortissement d'une solution convergée ; il doit passer en CI exact-HEAD avant d'être considéré comme preuve.

## 8. Journal des modifications
- 2026-10-24 — création pour l'audit #1069 et la PR #1078 ; contrat et critères de preuve documentés, qualification en attente.

## 9. Evidence update — 2026-10-24

Added regression coverage on the COUPLED physical-velocity flux contract:
- affine pressure / constant velocity on a nonuniform four-cell mesh, with compatible fixed-pressure inlet/outlet and zero-gradient walls;
- independent face-level pressure-jump/gradient-compensation oracle on each internal face;
- cellwise divergence of the authoritative reconstructed flux;
- pressure-gauge invariance on the COUPLED path, for both zero-gradient walls and fixed-pressure boundaries with consistently shifted boundary values;
- explicit H/A and physical-velocity APIs, with a direct comparison showing their difference equals the specified compensation;
- an actual coupled block solve on the nonuniform mesh, capturing the production matrix/RHS and comparing every continuity-row residual against the independently reconstructed authoritative flux divergence;
- a manufactured recovery gate starting from exact constant velocity and incorrect pressure, checking recovery of \(\mathbf U=(2,0,0)\) and \(p=1+2x\);
- a quantitative checkerboard central-face test comparing naive interpolated flux (0) with pressure-jump-corrected flux (-2).

Remaining before promotion:
- run the added grid-level checkerboard flux-response metric on exact HEAD;
- add or identify an explicit end-to-end regression proving the production SIMPLE/SIMPLEC solver continues to call the H/A predictor API and COUPLED calls the physical-velocity API;
- execute and review the applicable physical N9 matrix and total validation.

These additions strengthen the operator-level and block-level evidence, but the exact-HEAD Release/DebugSanitizers and numerical maturity campaigns must pass before promotion. Keep the PR draft until those results are green.

## 10. CI audit follow-up — 2026-10-26

The first exact-HEAD CI run for follow-up PR #1133 failed, so the PR remains draft. The failure review found a production-selection blocker in the PCD setup: its matrix-only pressure operators called `assemble_scalar_equation` with non-orthogonal correction enabled but without a current scalar field. The PCD operators now explicitly disable that deferred correction, which is consistent with their documented matrix-only contract; this does not alter the physical transport operator.

The same run exposed two independent regression-test/diagnostic defects:
- the skewed one-cell energy audit did not pass its current temperature to the energy assembly, despite requiring it for the non-orthogonal correction;
- scalar solve results left the physical residual-before/after-relaxation fields as NaN. The scalar solver now records the physical infinity-norm residual for the linear candidate and returned field, with relative values normalized by the RHS infinity norm; the 1x1 path uses the same relative-residual contract.
- the skewed variable-viscosity transpose-stress regression exposed normal-distance interpolation of viscosity on a skewed mesh. Internal-face viscosity is now reconstructed from least-squares cell gradients at the face centre, preserving affine-field exactness when the stencil supports it.

These are code/test corrections, not relaxed acceptance criteria. On the exact-HEAD run for the pre-rebase commit, Release CTest failed only `test_variable_viscosity_transpose_stress`; DebugSanitizers failed that test plus `test_mpi_deterministic_fused`. The CHT-time audit, energy-boundary audit, and finite-volume transport unit test passed on that run. The viscosity reconstruction has since been changed to use least-squares gradients with boundary ghosts excluded, and the PR has been rebased onto current master while preserving its H04 relative-tolerance fix. The new exact-HEAD CI has not yet completed. The 19-byte LeakSanitizer allocation from an unknown external module remains unresolved; it is neither suppressed nor waived. Release, DebugSanitizers, and full validation must pass on one exact HEAD before promotion.


## 11. Exact-HEAD validation follow-up — 2026-10-26, HEAD `3e09c881`

The latest N9 physical-validation run failed on the PCD Schur model's skewed Couette case. The test reached the configured nonlinear stagnation detector at iteration 58:
- momentum relative residual: `1.1e-5`;
- continuity (L_\infty): `1.63e-3`;
- pressure relative residual: `2.27e-4`;
- velocity change: `3.9e-5`;
- pressure change: `8.0e-5`.

The nested PCD (F_p) solves repeatedly report `status=1` at 1000 iterations despite relative residuals around (10^{-14}), while the outer coupled solve's independently recomputed linear true residual is around (10^{-12}). This is not sufficient evidence of nonlinear physical convergence: the continuity residual remains materially above the (10^{-7}) physical gate. Do not raise the stagnation window, relax the physical gates, or reinterpret the nested residual as proof of convergence.

Other exact-HEAD evidence observed so far:
- Release build and test job: passed.
- Pre-commit: passed.
- GUI CI, documentation, audit gate, numerical maturity audit, and N10 public fixture verification: passed.
- N9 physical validation: failed as described above.
- DebugSanitizers and the total-validation workflow: still running at the time of this update.
- The earlier compile blocker from the unqualified `Vec3` name was corrected; Release now builds.

Next investigation must isolate why the PCD coupled solve stalls specifically on skewed Couette, compare the reconstructed authoritative face flux and cellwise continuity against the assembled continuity block, and inspect the nested (F_p) termination/status semantics. Any fix must preserve the strict physical residual gate and be requalified on exact HEAD. PR #1133 remains draft and must not be merged until the physical failure is resolved and required CI completes.


## 12. Résultats CI et instrumentation de diagnostic — 2026-10-26

HEAD évalué par les campagnes complètes : `33318dae25eb71ff8a860ac5fc0893d26eda1427`.

- `CFDX CI / Release` : succès.
- `Pre-commit Fast Guards`, `CFDX documentation`, `cfdx-audit-gate`, `CFDX numerical maturity audit` et `N10 Public Fixture Verification` : succès.
- `N9 physical validation` et `N9 qualification campaign` : échec identique sur `COUPLED skew Couette`, à l'itération 58.
- `CFDX GUI CI`, `DebugSanitizers`, `CFDX numerical maturity` et `CFDX total validation` : encore en cours lors de cette lecture.

Mesures de l'échec N9 : `momentum_rel=1.1e-5`, `continuity_linf=1.63e-3`, `pressure_rel=2.27e-4`, `velocity_change=3.9e-5`, `pressure_change=8.0e-5`. Le résidu vrai du système linéaire couplé est de l'ordre de `1e-12`. Les appels imbriqués `Fp` affichent `status=MAX_ITER_REACHED` après 1000 itérations, mais avec un résidu relatif mesuré proche de `1e-14`; le contrat actuel les accepte comme applications de préconditionneur si l'itéré reste fini et exploitable. Ce statut seul ne prouve donc pas la cause de l'échec physique et ne justifie pas de modifier la tolérance.

### Suite de diagnostic

1. Ne pas relâcher les seuils de continuité/momentum ni le budget de stagnation.
2. Le bloc couplé assemble la correction non orthogonale explicitement à partir de (p_{old}), tandis que le flux physique d'acceptation est reconstruit avec le gradient de la pression mise à jour. C'est un point à vérifier quantitativement, pas encore une cause démontrée.
3. Le commit `7c72fb51e510ec4aff63b5234b39e2b3f27065e6` enrichit le message d'échec de `test_n9_physical_matrix` avec les canaux indépendants déjà calculés : continuité reconstruite depuis la vitesse, écart flux/vitesse et identité Rhie–Chow. Ces données doivent permettre de distinguer une incohérence de flux d'une stagnation physique réelle au prochain run.
4. Après la collecte de ces mesures, comparer le résidu de continuité de la ligne algébrique capturée au résidu des flux reconstruits sur le même état ((U,p)). Toute correction doit préserver les gates actuelles et être accompagnée d'un test de non-régression ciblé.


## 13. Requalification après instrumentation — 2026-10-26, HEAD `618951c1339a3aeda2cc24b560157b8e313127a4`

Le rerun sur le HEAD courant confirme que l'instrumentation est active mais que le défaut physique n'est pas corrigé.

Campagne [N9 physical validation, run 1158](https://github.com/stevendaix/cfdx_dev/actions/runs/40380695906) :
- `COUPLED skew Couette` stagne à l'itération 58 ; `momentum_rel=1.1e-5`, `continuity_linf=1.630e-3`, `pressure_rel=2.27e-4`.
- `corrected_flux_continuity_linf=1.630e-3`, `reconstructed_velocity_continuity_linf=1.027e-3`.
- `flux_velocity_mismatch_linf=3.51e-4`, `flux_velocity_mismatch_normalized=4.293e-3`, `rhie_chow_identity_linf` arrondi à zéro.
- Le résidu vrai linéaire couplé est de l'ordre de `1e-12`. Les applications imbriquées Fp atteignent encore 1000 itérations avec un résidu relatif affiché autour de `1e-14`.

Interprétation limitée aux preuves disponibles : l'identité Rhie–Chow vérifie l'identité locale instrumentée, mais elle ne prouve ni l'égalité entre la ligne de continuité assemblée et la divergence du flux autoritaire, ni la convergence non linéaire. La continuité reconstruite depuis la vitesse est inférieure à celle du flux corrigé d'environ (6.03\times10^{-4}) en norme (L_\infty). Cet écart doit être localisé cellule par cellule avant de modifier l'assemblage. La valeur normalisée de l'écart flux/vitesse reste une métrique diagnostique, pas un seuil d'acceptation.

État CI exact-HEAD :
- Release : succès.
- Pre-commit, GUI CI, documentation, audit gate, numerical maturity audit, numerical maturity et N10 : succès.
- N9 physical validation et N9 qualification campaign : échec physique ci-dessus.
- DebugSanitizers : échec indépendant de `test_mpi_deterministic_fused` ; les 242 autres tests de ce job passent. Le défaut MPI doit rester séparé de F03, sans suppression de sanitizer ni changement de test.
- CFDX total validation : annulé après l'échec d'une étape requise ; ce n'est pas une validation réussie.

### Actions exigées avant correction

1. Capturer sur le même état ((U,p)) la ligne algébrique de continuité ((A x-b)_p), la divergence du flux corrigé et la divergence reconstruite depuis la vitesse ; rapporter les maxima, cellules d'argmax et valeurs signées. Ne pas comparer des champs issus de deux itérations différentes.
2. Vérifier les signes et la compensation de non-orthogonalité entre l'assemblage couplé (correction explicite fondée sur (p_{old})) et la reconstruction du flux d'acceptation (pression mise à jour). La divergence des métriques est un indice à tester, pas encore une preuve de cause.
3. Pour Fp, distinguer le résidu déclaré par le solveur de la terminaison au budget et vérifier le résidu vrai recalculé avec exactement la matrice/RHS internes. Ne pas accepter le seul résidu affiché comme preuve de convergence physique ; ne pas durcir ou relâcher le contrat best-effort sans preuve.
4. Garder les gates physiques et le détecteur de stagnation inchangés. Relancer N9, Release, DebugSanitizers, les audits et la validation totale après une correction causale. La PR reste draft.

## 13. Diagnostic pre-relaxation du bloc COUPLED — 2026-10-26

La campagne exacte sur HEAD `618951c1` confirme que le système linéaire couplé atteint un résidu vrai relatif inférieur à `1e-12`, alors que le flux reconstruit après mise à jour non linéaire conserve une divergence de `1.63e-3`. Les canaux indépendants mesurés au point d'échec sont : divergence issue de la vitesse `1.027e-3`, écart de flux interne `3.51e-4`, écart normalisé `4.293e-3`, identité de l'opérateur Rhie–Chow `0`.

Une instrumentation ciblée a été ajoutée dans `solve_steady_incompressible` pour mesurer la divergence du flux calculé depuis la solution monolithique **avant** les relaxations séparées de U et p. Elle journalise aussi `alpha_u` et `alpha_p`. Cela doit trancher entre :
- une incohérence du bloc de continuité et du flux déjà présente sur la solution linéaire ;
- une perte de continuité introduite par la relaxation non linéaire indépendante ;
- un effet du terme non orthogonal différé évalué avec l'ancien gradient de pression.

Le diagnostic est conditionné par `coupled_matrix_summary`, n'affecte pas les champs, les critères d'acceptation ni les budgets. Les gates physiques restent inchangées. Le prochain run doit être lu avant de modifier les signes ou la formulation du bloc : l'inspection algébrique seule ne justifie pas un changement de signe qui pourrait rompre les cas orthogonaux déjà qualifiés.

Sur la dernière campagne disponible :
- Release : succès ;
- pre-commit, documentation, GUI, audit gate, maturité numérique et fixtures N10 : succès ;
- N9 physique et campagne N9 : échec sur Couette non orthogonal ;
- DebugSanitizers : échec uniquement sur la fuite de 19 octets de `test_mpi_deterministic_fused` ;
- validation totale : échec du rapport de validation (sept suites Python en échec) et campagne longue annulée, en plus du blocage N9.

La PR reste draft. Ne pas fusionner tant que le diagnostic exact-HEAD n'a pas identifié la cause et que les contrôles requis ne sont pas verts.


## 14. Décomposition cellulaire de la continuité — 2026-10-26

Le diagnostic `COUPLED_CELL_OPERATOR_AUDIT` est ajouté sur le candidat monolithique exact, avant toute relaxation séparée de U et p. Pour la cellule argmax du résidu de ligne et celle argmax de la divergence du flux (si différente), il publie séparément la somme des contributions matricielles vitesse/pression, le second membre et le résidu algébrique, puis les composantes du flux conservatif : interpolation de vitesse, correction de pression orthogonale, correction non orthogonale différée et compensation Rhie–Chow.

Le diagnostic recalcule le flux à partir du même état U,p et des mêmes diagonales de quantité de mouvement. Il est conditionné par `coupled_matrix_summary`, n'altère ni le solveur ni les critères d'acceptation. Les résultats doivent être lus sur le prochain run exact-HEAD avant toute modification de signe, coefficient, tolérance ou budget. Le statut `MAX_ITER_REACHED` des sous-résolutions Fp reste une piste secondaire séparée ; aucun changement de statut n'est inclus dans ce commit.
