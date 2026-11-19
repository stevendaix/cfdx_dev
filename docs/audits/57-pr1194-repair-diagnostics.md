# PR #1194 — repair and temporary diagnostic campaign

## Purpose

This follow-up PR is based on PR #1194 and has two deliberately separate goals:

1. Repair integration defects that prevent the existing regression and SST BFS case from running.
2. Run a bounded GitHub Actions campaign and retain the complete logs, generated cases, solver outputs, residual histories, JSON evidence, build metadata and CTest output as downloadable artifacts.

The workstation result mentioned in the parent PR discussion (converged k-omega SST solution, reported 6–7% physical discrepancy and long runtime) is not available in this repository snapshot. This PR does not claim to reproduce it. The CI campaign uses the repository's generated TH03 debug case and is explicitly diagnostic, not experimental validation.

## Repairs

- Update ScalarEquation aggregate initializers in the TVD deferred-defect regression test to account for the implicit_relaxation boolean added to the struct.
- Serialize sst_adaptive_relaxation as the boolean string true, because the production solver's raw-setting contract accepts only true or false. The previous generator emitted enabled, causing the SST case to fail before entering the solver.

No physical convergence tolerance, residual threshold, iteration acceptance criterion or physical reference value is relaxed.

## CI campaign

The temporary workflow .github/workflows/bfs-pr1194-repair-diagnostics.yml:

- builds the production solver and the repaired regression target in Release with HDF5;
- runs the focused TVD regression and BFS smoke test;
- runs a same-mesh comparison of SIMPLE, SIMPLEC, PCD-auto, BlockLocal-AMG, LSC-AMG and BFBT-AMG;
- runs the PR-specific --bfs-converged, --bfs-fieldsplit and --bfs-pimple modes on a common generated case;
- applies a 600-second wall-clock cap per solver execution, retaining timeout exit codes and partial output;
- uploads the complete diagnostic directory as a 30-day artifact, even when solver runs fail or time out.

The workflow is intentionally diagnostic: solver non-convergence and timeouts are recorded as negative evidence and do not get silently converted into success. Build/test failures and solver outcomes are distinguishable in the uploaded evidence.

## Interpretation rules

- EXIT_ZERO is process status, not a physical convergence certificate.
- The generated debug mesh is not a Driver–Seegmiller validation mesh.
- A small mesh can diagnose configuration, dispatch, stability and cost, but cannot establish grid-independent separation or reattachment.
- No comparison to the unavailable workstation result is valid until the original case, mesh, inlet profile, initial fields, SST settings and solver controls are recovered or checked into an immutable reproducibility package.
- Keep all failed variants and timeouts in the report. Do not select only the best run.

## Artifact contents

The artifact contains the exact commit and runner identity, CMake configuration, build log, focused CTest log and exit code, each campaign log and output directory, generated HDF5 cases, solver run directories, direct-mode logs/exit codes, CMake cache, compile commands and CTest temporary output.

## Removal

This workflow is temporary. After the diagnostic campaign has been reviewed and the required evidence has been transferred into a permanent audit report, remove the temporary workflow or replace it with a smaller permanent regression.

## Driver–Seegmiller experimental reference audit

The diagnostic workflow now retrieves the three *raw* NASA Turbulence Modeling Resource files, records their SHA-256 hashes, and runs the existing independent parser/contract tests:

- `cf.exp.dat`: lower-wall skin-friction coefficient and its tabulated experimental uncertainty;
- `cp.expnew.dat`: top- and bottom-wall pressure coefficients, preserving both original and uniformly shifted pressure origins;
- `profiles.exp.dat`: mean velocity and Reynolds-stress profiles, including the stations (x/H=-4,1,4,6,10).

Reference definition: NASA TMR 2DBFS, Driver & Seegmiller (1985), DOI 10.2514/3.8890. The reference Reynolds number is approximately (Re_H=36,000), based on step height; the upstream momentum-thickness Reynolds number is 5,000. The measured reattachment position is (x_r/H=6.26\pm0.10). The reference velocity (U_{ref}) is taken near the channel centre at (x/H=-4), not assumed equal to a nominal inlet velocity. NASA also notes that the plotted (C_p) data are shifted to make (C_p\simeq0) near (x/H=40), while the downloadable file retains both origins.

The comparison contract is deliberately multi-observable:

1. (x_r/H): obtain from the first downstream zero crossing of a qualified streamwise wall-shear distribution; compare the signed error with the experimental interval.
2. (C_f(x/H)): compare at the experimental abscissae using the supplied pointwise experimental uncertainty; report maximum and RMS normalized errors.
3. (C_p(x/H)): compare both top and bottom walls and explicitly select the original or shifted NASA pressure reference. A pressure-origin mismatch is a failed comparison, not a tunable offset.
4. (U/U_{ref}), (V/U_{ref}), and (-\overline{u'v'}/U_{ref}^2): compare profiles at (x/H=1,4,6,10) using interpolation in (y/H). The parser also reads the normal Reynolds stresses; these must not be silently omitted from a future complete qualification.
5. Grid convergence: keep separate from experimental agreement. Require at least three systematically refined, geometrically similar grids and report observed order/GCI for key quantities; a single mesh cannot establish discretization uncertainty.

The workflow currently performs a **reference-data integrity and parser-contract audit**, not a fabricated CFD-to-experiment comparison. The 192-cell CI case has a synthetic inlet profile and is a solver diagnostic only. Its SIMPLE result is not globally converged (momentum residual remains (O(10^{-2})), eddy-viscosity closure changes materially); its COUPLED/PCD result fails the first coupled linear solve. Therefore no (C_f), (C_p), velocity-profile or Reynolds-stress agreement is claimed for those runs. The fail-closed tool `qualify_bfs_nasa_fields_v16.py` emits `NO_CONVERGED_CFDX_RESULT_SUPPLIED` until a genuine prediction dataset is provided.

### Experimental setup mismatch to resolve before physical comparison

Do not identify the OpenFOAM tutorial's nominal (U=44.2\,m/s), (H=0.0127\,m) and (Re_H\simeq36,000) as a complete reproduction of the NASA validation boundary conditions. NASA TMR describes a naturally developing turbulent inlet boundary layer, approximately (1.5H) thick, with momentum-thickness Reynolds number 5,000 and back pressure adjusted to (M\simeq0.128). The inlet development length is deliberately selected to recover the measured upstream profile and skin friction. The current compact CI generator's 1/7-power-law profile and 5% turbulence intensity are not the measured Driver–Seegmiller inlet data; they must remain labelled surrogate/debug inputs, not validation inputs.

### Corrected CI contract

The focused build now explicitly builds `test_bfs_qualification` as well as the production solver and TVD regression executable. This addresses the review finding that CTest's `test_bfs_quick` wrapper launches a separate executable which was not built by the original workflow. The workflow also runs the independent NASA reference-parser tests and archives the raw experimental files, checksums, parser output and exit status.

References:
- NASA TMR 2DBFS: https://tmbwg.github.io/turbmodels/backstep_val.html
- Driver & Seegmiller (1985), AIAA Journal 23(2), 163–171, DOI: https://doi.org/10.2514/3.8890
- NASA NTRS record: https://ntrs.nasa.gov/search.jsp?R=19850039701
- OpenFOAM tutorial case: https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/

## Analyse des résultats CI initiaux (commit 7380f1e)

La première campagne complète a bien compilé le solveur Release et le test TVD, mais elle a révélé une erreur de construction du smoke test : `test_bfs_quick` lance `test_bfs_qualification`, qui n'avait pas été demandé dans la liste des cibles à construire. CTest a donc rapporté `Not Run` (exécutable absent), et non un défaut du solveur. La liste des cibles du workflow a été corrigée pour construire explicitement cet exécutable.

Résultats de la matrice numérique sur le maillage debug commun (192 cellules, 1500 itérations maximum) :

| Variante | Temps observé | Observation finale | Statut |
|---|---:|---|---|
| SIMPLE | 30,09 s | résidu momentum global 2,23 %, composante active max 6,33 %, continuité 8,53e-4 | non convergé |
| SIMPLEC | 2,72 s | divergence numérique vers l'itération 119, champs et résidus de l'ordre de 1e19 | rejeté |
| COUPLED/PCD-auto | 4,38 s | FGMRES 2000 itérations, résidu relatif vrai 6,27e-3 | solveur linéaire non convergé |
| COUPLED/BlockLocal-AMG | 48,98 s | momentum global 2,28 %, continuité 3,88e-4 | non convergé |
| COUPLED/LSC-AMG | 149,77 s | momentum global 5,26 %, continuité 4,81e-4 | non convergé |
| COUPLED/BFBT-AMG | 189,06 s | momentum global 4,40 %, continuité 3,61e-4 | non convergé |

Ces chiffres ne classent pas universellement les algorithmes : ils décrivent uniquement ce maillage, cette initialisation uniforme, ces contrôles et cette exécution. En particulier, les temps des variantes couplées incluent des résolutions internes différentes et ne doivent pas être interprétés sans nombre d'itérations externes, coût par action de préconditionneur et preuve de convergence.

Le mode FieldSplit Full Schur échoue au premier système couplé : FGMRES 445 itérations, résidu relatif vrai 3,31e-3. Le mode COUPLED/PCD échoue après 2000 itérations avec 6,27e-3. Ces deux résultats localisent un défaut de résolution linéaire initial, avant toute comparaison expérimentale.

Le mode PIMPLE progresse nettement plus loin. À 1500 itérations :
- momentum global relatif : 8,13e-4 ;
- pire composante momentum active : 6,64e-3 ;
- continuité normalisée : 1,54e-5 ;
- résidu global k : 4,63e-6 ;
- résidu global omega : 8,08e-7 ;
- transport SST : candidat non certifié / `converged=false` ;
- gate global U-p-k-omega : `converged=false`.

Le mode PIMPLE est donc le meilleur candidat de progression sur ce petit test, mais il n'a pas convergé. Le faible résidu de continuité et les faibles résidus linéaires de pression ne compensent ni le défaut momentum par composante ni le défaut de fermeture SST.

### Conséquence physique

Aucune des variantes ne fournit un état global convergé : le post-traitement fail-closed doit donc s'arrêter avant le calcul des quantités d'intérêt. Il serait incorrect de comparer les résidus de ces runs à l'expérience comme s'il s'agissait d'une solution RANS établie, et incorrect de déduire (x_r/H), (C_f) ou (C_p) à partir d'un état rejeté ou incomplet.

Le point de recollement expérimental (6.26\pm0.10) est un critère de validation de la solution, pas une condition numérique à imposer au solveur. La campagne doit d'abord produire un état U-p-k-omega certifié, puis des contraintes pariétales et profils cohérents, puis seulement effectuer les comparaisons multi-observables et l'étude de maillage.


The repository audit gate now classifies this file explicitly as a campaign evidence dossier (rather than a numerical-method fiche). Its exemption from the seven-section method-fiche template is covered by a regression test; the dossier still carries the experimental DOI and explicit diagnostic/qualification status.
