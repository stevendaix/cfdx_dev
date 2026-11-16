# CFDX V26 - Jv multi-amplitude, objectif L2 SST et verrous de convergence

**Date : 10 novembre 2026. Source : CFDX V25 exacte.** Version experimentee en opt-in, **desactivee par defaut**. Aucun seuil de convergence, schema FV, condition limite ou champ de production n'a ete assoupli. Aucune qualification physique BFS ou NASA revendiquee.

## 1. Motivation et reference V25

Sur le BFS de diagnostic (712 cellules, memes HDF5 et parois V25), a l'iteration interne SST 10, le preconditionneur spatial fournit une direction Newton dont le residu directionnel GMRES vaut 0.16613194317, inferieur au critere **lineaire seulement** 0.25. Sa recherche de pas V25 trouve un alpha de 0.00087718878836 qui reduit legerement le merite local Linf de 0.92596107454 a 0.92595972137. Le Newton-Jacobi V22, conserve comme solution de repli, propose un meilleur merite : 0.88990745136.

Cette observation ne prouve ni un defaut de Jacobien, ni une mauvaise condition aux limites, ni une convergence globale. V26 distingue la stabilite de Jv, le merite lisse et le **certificat physique Linf**.

## 2. Code de production V26

### 2.1 Audit directionnel sur l'operateur non lineaire reassemble

Nouveau `sst_audit_log_jvp_v26()` dans `src/cfdx/physics/sst_globalization_v26.h` : sur l'etat SST de reference `x=[k,omega]`, pour la direction de Newton logarithmique `d`, l'audit mesure

```
 x_plus_i  = x_i * exp(+h * d_i / max_j(abs(d_j)))
 x_minus_i = x_i * exp(-h * d_i / max_j(abs(d_j)))
 Jv_h      = (F(x_plus)-F(x_minus)) / (2*h)
 h = 1e-3, 1e-4, 1e-5
```

La fonction F est **le vecteur signe de residus volumes finis reassembles** deja utilise par JFNK (y compris contraintes proprietaires non projetees en tant qu'equations du Jacobien). Les essais utilisent aussi des differences unilaterales et deux comparaisons entre amplitudes. Les sept evaluations (base + six perturbations) sont uniquement des diagnostics. Elles ne sont ni une estimation rigoureuse de l'erreur du Jacobien, ni un certificat SST. Un point de non-regularite reste possible hors de ces trois echelles.

### 2.2 Objectif lisse L2 issu des memes lignes FV

Pour chaque equation phi (k ou omega), les residus de ligne locaux restent definis par le schema de production :

```
 rho_phi_i = abs(A_i*phi-b_i) / (sum_j abs(A_ij) *
                      max(phi_floor, max_{j dans stencil i} abs(phi_j)))
 M_inf = max_i(rho_k_i, rho_omega_i)     [porte physique EXISTANTE]
 M_2 = hypot(sqrt(sum_i rho_k_i^2/n),
             sqrt(sum_i rho_omega_i^2/n)) / sqrt(2)
```

`M_2` suit la meme politique native/compensee que `M_inf` (`sst_compensated_residual_gate`). Les valeurs RMS sont calculees dans le controle deja present `sst_local_residual_precision_audit`, a partir des memes matrices et etats, sans discretisation parallele.

### 2.3 Recherche double avec aucune degradation du critere final

Nouveau `sst_globalize_log_direction_v26` : les 13 candidats logarithmiques (au plus) sont reassembles **une seule fois** chacun. La V26 produit deux propositions depuis ces memes essais :

* V25 : candidat minimisant `M_inf` sous sa decroissance stricte precedente.
* V26 : candidat minimisant `M_2` **a condition** de diminuer simultanement `M_2` et le `M_inf` physique, chacun avec la garde de flottants `(1-1e-9*alpha)`.

Ensuite le solveur compare les propositions Picard, Newton-Jacobi, blocs locaux, spatial, V25 et V26 avec le **meme** test `sst_nonlinear_step_reduces_merit` sur les vrais residus reassembles et les parois actualisees. Ni `M_2` seul, ni un residu Krylov ne peuvent marquer un pas comme converge. Echec de tous les pas => rollback. Un candidat V26 rejete ne supprime jamais le recours V22.

Le drapeau `sst_smooth_globalization_v26_experimental=true`, dans `numerics.raw_settings`, exige egalement l'activation des modes V22, V24 et V25. Sans ce drapeau, le chemin V25 est preserve.

## 3. Resultats reproductibles BFS sur le MEME HDF5 a 712 cellules

Condition de comparaison : un seul champ numerique differe entre les configurations JSON HDF5 (le nouveau drapeau V26). Toutes les donnees HDF5 de maillage, BC et tolerances sont identiques. Budget SST diagnostique **10** iterations internes, tolerance physique **1e-8** non modifiee. Comparaison au premier essai exterieur avant retries.

| Indicateur a l'iteration SST 10 | V25 (binaire V26, option desactivee) | V26 activee |
|:--|--:|--:|
| Residuel GMRES direction spatiale | 0.16613194317 | 0.16613194317 |
| Direction spatiale avec reduction Linf | oui, alpha = 0.00087718879 | oui, meme alpha V25 |
| Linf SST avant selection | 0.92596107454 | 0.92596107454 |
| Meilleur Linf spatial V25 | 0.92595972137 | 0.92595972137 |
| L2 equilibre initial | non journalise | 0.17800728499 |
| Meilleur L2 spatial parmi les essais | non journalise | **0.17800747609** (hausse) |
| Proposition V26 a double decroissance | sans objet | **refusee** |
| Merit Linf apres selection finale | **0.88990745136** | **0.88990745136** |
| Methode selectionnee | Newton-Jacobi V22 | Newton-Jacobi V22 |
| Etat exterieur accepte/converge | non | non |

La variation du L2 pour la direction spatiale est positive (+1.9110e-7, soit environ +1.07e-6 en relatif). La recherche par objectif lisse **n'apporte donc pas de nouveau pas** sur ce cas. La selection conserve le pas V22 deja le meilleur en Linf. Il serait incorrect d'annoncer un gain de convergence grace a V26.

### 3.1 Audit des trois pas logarithmiques de Jv

| h logarithmique | Norme 2 de Jv_h | Desaccord unilaterales vs centree |
|--:|--:|--:|
| 1e-3 | 0.017251845567 | 6.8516556682e-5 |
| 1e-4 | 0.017251847453 | 8.1698807518e-6 |
| 1e-5 | 0.017251865832 | 2.9850167301e-5 |

Desaccord relatif entre Jv a 1e-3 / 1e-4 : `8.1049739281e-7`. Entre 1e-4 / 1e-5 : `6.0224799895e-6`. Le maximum des diagnostics vaut `6.8516556682e-5`. Les trois amplitudes ne montrent pas une divergence majeure de Jv, ce qui **oriente** vers une mauvaise adequation de la direction au merite plutot que vers une erreur evidente de differenciation. Ceci n'est pas une verification analytique de tout le Jacobien, ni de l'equivalence aux fonctions de paroi d'OpenFOAM.

## 4. Tests executes et qualification

* Compilation et lien de `cfdx_production_solver` : reussis, `Clang -O0` diagnostique, sans benchmark Release optimise.
* CTest cibles V17/V18/V19/V21/V22/V23/V24/V25/V26 et protocoles Python V24/V25/V26 : **12/12 reussis**.
* Nouveau test C++ V26 : **7/7** sous-cas, dont calcul de RMS a partir de la matrice FV, Jv multiscale, source non reguliere, selection Linf/L2 differente, refus des pas L2-only, rejet des projections invalides et dependances de drapeaux.
* Nouveau protocole Python V26 : **4/4** sous-cas, detection fail-closed d'une fausse qualification, trace manquante ou Linf degrada.
* BFS 712 cellules : les deux calculs comparatifs ont **rejete leur premiere tentative exterieure** ; les runs de preuve sont volontairement limites au premier rejet. Aucune demonstration a 800 iterations, aucun etat exterieur pleinement converge.
* CI complete, NACA, BFS 20 540 cellules, Schur/PCD/LSC/BFBt/AMG/MGR, alignement complet omegaWallFunction et comparaison Driver-Seegmiller NASA : **non testes dans V26**.
* Option experimentale desactivee par defaut, pas de changement de BC, pas de push GitHub.

## 5. Reproduction

```sh
cmake -S . -B build-v26 -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_TESTS=ON \
    -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v26 --target cfdx_production_solver test_sst_globalization_v26
ctest --test-dir build-v26 -R 'test_sst_globalization_v26|test_bfs_sst_v26_protocol' --output-on-failure
# Dans numerics.raw_settings du meme HDF5 de 712 cellules :
# sst_smooth_globalization_v26_experimental=true
# (et V22, V24, V25 deja actives ; sst_transport_max_iterations=10, tol=1e-8)
CFDX_SST_DIAGNOSTICS=1 build-v26/cfdx_production_solver \
    --mesh bfs_v26_on_10.cfdx.h5 --output-dir bfs_out \
    --iterations 1 --experimental-sst
python3 scripts/diagnose_bfs_sst_v26.py v26_on_final.log
```

Les journaux contiennent les tentatives exterieures suivantes et des retries : les comparaisons retiennent uniquement l'**essai initial**. Un timeout de campagne n'est ni une convergence ni une preuve d'echec a 800 iterations.

## 6. Decision et suite

Le diagnostic multiscale n'indique pas de plancher de differenciation numerique dominant dans cette direction a l'iteration 10. Le L2 dynamique indique au contraire que la direction Newton spatiale peut etre une **direction de montee** pour la norme lisse, malgre sa bonne resolution du systeme gele et sa tres faible diminution du Linf. Il faut maintenant proposer des **directions** plus adaptees au merite dynamique (eventuellement restauration du sous-probleme actif L_inf, amortissement par region de confiance ou direction melangee certifiee), et non simplement prolonger la recherche de pas sur cette direction. Conserver la surveillance de `k`, `omega`, `G`, `nu_t` et des contraintes proprietaires, puis valider le BFS a 800 iterations, trois etats exterieurs, grand BFS et variantes Schur avant tout certificat NASA.

### Bibliographie methodologique

- PETSc SNES 3.26.0, `https://petsc.org/release/manual/snes/` et `https://petsc.org/release/manualpages/SNES/SNESLINESEARCHBT/` : methodes de line search, objectif L2 et separation du solveur non lineaire; V26 n'implemente pas PETSc.
- OpenFOAM Foundation v11 `https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html` : `updateCoeffs()` et `manipulateMatrix()` appliquent des poids de faces distincts; equivalent CFDX non certifie.
