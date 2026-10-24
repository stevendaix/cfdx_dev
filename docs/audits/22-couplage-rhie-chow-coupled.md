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
- À la rédaction de cette fiche, la CI Release/DebugSanitizers et les campagnes de cohérence exact-HEAD ne sont pas encore toutes terminées. PR #1078 reste draft ; ne pas merger avant inspection des résultats et revue du test affine/checkerboard.

## 8. Journal des modifications
- 2026-10-24 — création pour l'audit #1069 et la PR #1078 ; contrat et critères de preuve documentés, qualification en attente.

## 9. Evidence update — 2026-10-24

Added regression coverage on the COUPLED physical-velocity flux contract:
- affine pressure with constant velocity on a rectilinear but nonuniform four-cell mesh, with compatible fixed-pressure inlet/outlet and zero-gradient walls;
- independent face-level pressure-jump/gradient-compensation oracle on every internal face;
- cellwise divergence of the authoritative reconstructed flux for the constant-velocity manufactured state;
- pressure-gauge invariance on the COUPLED path for zero-gradient walls and a separate fixed-pressure-boundary test where the boundary pressure is shifted consistently.

Still required before promotion:
- the exact production matrix/RHS evidence hook now supports a cellwise comparison of the assembled continuity block residual against independently reconstructed authoritative flux divergence after an actual block solve;
- a manufactured coupled velocity-pressure solution on a nonuniform mesh;
- the checkerboard central-face test now quantifies the naive interpolated flux (0) versus the pressure-jump-corrected flux (-2) and checks both named contracts; a grid-level pressure-mode amplitude/suppression metric remains to be added;
- explicit regression that SIMPLE/SIMPLEC retain the H/A predictor contract. The call sites now use named APIs (make_rhie_chow_flux_from_hbya / make_rhie_chow_flux_from_physical_velocity), and the nonuniform affine test compares the two contracts directly.

These additions improve the operator-level evidence but do not yet establish the complete coupled block contract or physical qualification. Keep the PR draft until the remaining independent checks and exact-HEAD CI pass.
