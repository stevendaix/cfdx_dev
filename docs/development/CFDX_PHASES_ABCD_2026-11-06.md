# CFDX — Lots A/B/C/D : diagnostic, AMG anisotrope, flux et couplage

**État au 6 novembre 2026.** Base vérifiée : ZIP utilisateur `cfdx_dev-master(2).zip` et lot précédent de corrections directionnelles. Ce document ne prétend pas qualifier le HEAD distant GitHub.

## Principe de qualification

Le résidu physique **vrai** de l'opérateur d'origine, l'erreur par rapport à une solution manufacturée, la conservation et le temps **total** de résolution sont des critères distincts. Aucun critère de convergence n'a été assoupli. Les nouvelles variantes restent opt-in et les paires CG / préconditionnement non symétrique sont refusées.

## Phase A : diagnostics reproductibles — implémenté, raccordement automatique CFD encore absent

- `src/cfdx/core/linalg/anisotropy_analysis.h` : analyse CSR des diagonales, contrastes hors-diagonale, connexions fortes, graphe non réciproque, ramifications, modes dominants géométriques (directions principales), alignement voisin, plus mauvaises lignes.
- Les indicateurs géométriques optionnels (centres cellules, rapport d'aspect, non-orthogonalité, skewness) restent séparés des indicateurs algébriques : **chaque équation** doit être analysée sur sa propre matrice.
- `tools/diagnostics/cfdx_anisotropy_report.cpp` : CLI MatrixMarket, import optionnel de centres/qualité CSV, sortie JSON déterministe, aucune dépendance PETSc/HYPRE externe.

Exemple :

```bash
cmake --build build --target cfdx_anisotropy_report
build/cfdx_anisotropy_report pressure.mtx --equation pressure \
  --centers cell_centers.csv --aspect cell_aspect_ratio.csv \
  --faces face_quality.csv
```

`centers.csv` : trois composantes par cellule. `cell_aspect_ratio.csv` : un scalaire par cellule. `face_quality.csv` : non-orthogonalité en degrés, puis skewness par face. Les sorties géométriques ne sont pas fabriquées quand les fichiers manquent. Il reste à raccorder cette instrumentation à l'export automatique des matrices de NACA/SST.

## Phase B : AMG directionnel — expérimentations disponibles, pas de promotion automatique

Le **lot précédent** conserve son Line-Gauss-Seidel multiplicatif, son semi-coarsening algébrique et son repli documenté pour directions obliques complexes. Le présent lot ajoute :

- `AMGInterpolationPolicy::EnergyMinimizingAggregation` avec candidats de prolongation à **portée deux**, corrections de type descente d'énergie projetées, contrainte `P*1=1` et rejet des pas qui augmentent l'énergie `trace(P^T A P)` ;
- `PreconditionerModel::EnergyMinimizingAMG`, accessible uniquement en FGMRES sur les opérateurs elliptiques à diagonale positive et symétriques mesurés ;
- test Neumann anisotrope avec projection explicite de la constante, tests d'orientations x, y et obliques et résidus vrais.

Sur la matrice manufacturée 16×16 retenue, l'énergie du prolongateur est passée de **10 733,7** à **4 903,27** ; le mode constant est préservé à 1e-12. **Ce résultat de qualité de prolongation ne signifie pas que le solveur est plus rapide.** Dans nos six cas 256/1024 inconnues, l'AMG energy-minimizing est systématiquement plus lent que l'AMG natif (voir CSV). Il ne doit donc pas être sélectionné automatiquement. Une qualification future doit couvrir les matrices où les autres méthodes stagnent, ainsi que le coût de construction des niveaux grossiers.

Usage C++ :

```cpp
cfdx::core::LinearSolverRequest request;
request.preconditioner = cfdx::core::PreconditionerModel::EnergyMinimizingAMG;
request.gmres_restart = 50;
// En Neumann pur, choisir aussi NullSpaceModel::Constant.
```

## Phase C : discrétisation — certification supplémentaire, reste incomplet en COUPLED

La base existante contient la reconstruction WLS tenant compte des CL, le flux diffusif décomposé avec correction non orthogonale, Rhie–Chow et une reconstruction corrigée en skewness. Dans le présent lot :

- `tests/validation/test_pressure_anisotropic_affine_flux.cpp` teste l'identité entre la partie orthogonale implicite, le flux différé et la réponse de pression tensorielle analytique sur des faces orientées/étirées ;
- les campagnes de skewness et du Laplacien non orthogonal sont rejouées : dix cas de skewness passés et campagne Laplacien marquée PASS ;
- les contrôles de géométrie dégénérée et d'orientation de face restent fail-closed.

**Verrou non résolu :** activer le gradient de pression corrigé en skewness dans COUPLED exigerait son Jacobien dans la matrice implicite, cohérent avec les flux des équations de quantité de mouvement. Le garde de validation existant continue de refuser cette combinaison. Le correctif de l'opérateur COUPLED complet et sa qualification sur maillages fortement non orthogonaux ne sont **pas** terminés.

## Phase D : couplage — Vanka et MGR corrigés, variantes restantes à qualifier

- Nouveau `VankaPressurePatchPreconditioner` : blocs locaux pression–vitesses depuis les deux graphes D/G, LU pivoté par patch, balayages multiplicatifs aller/retour, borne de taille, repli pour variables de vitesse sans couplage pression et échec explicite si pivot singulier.
- `PreconditionerModel::VankaPressurePatch` est opt-in sur la matrice couplée 4N et FGMRES **uniquement**. Aucun changement silencieux des choix du solveur de production.
- `NativeMGRPreconditioner::setup/update_values` : reconstruction **transactionnelle** ; une hiérarchie/factorisation rejetée ne modifie plus l'ancienne. L'application MGR est maintenant sûre en cas d'aliasing entrée/sortie.
- En contrepartie, `update_values()` MGR reconstruit aujourd'hui toute la hiérarchie : **coût potentiellement plus élevé**, à optimiser après ajout d'un véritable clonage d'état AMG numérique.

Sur une matrice couplée manufacturée à 120 inconnues, tolérance de résidu vrai 1e-12 :

| Cas | Sans préconditionneur | Vanka | Temps total comparatif |
|---|---:|---:|---|
| Normal | 50 itérations | 15 itérations | 3,47 ms → 0,26 ms |
| Étirement coefficients ×1000 | 24 itérations | 6 itérations | 0,98 ms → 0,082 ms |

Temps mesurés **une fois** localement sur système synthétique ; la préparation Vanka est incluse dans son temps, pas de promesse de gain sur le CFD industriel. Erreur relative de la solution du cas étiré : environ 3,0e-9 (bien qu'avec résidu vrai ≈1,4e-13) ; c'est pourquoi le benchmark vérifie séparément ces deux grandeurs.

`PCD`, `LSC`, `BFBt`, `Block-Schur` existent à des niveaux de maturité différents dans la base de code. Le présent lot n'en démontre pas la qualification industrielle. `Braess-Sarazin` et l'interpolation adaptative complète selon la distance algébrique ne sont pas intégrés ici.

Usage C++ :

```cpp
cfdx::core::LinearSolverRequest request;
request.preconditioner = cfdx::core::PreconditionerModel::VankaPressurePatch;
request.krylov = cfdx::core::KrylovModel::FGMRES;
// solve_linear_system(A_4N, rhs, x, LinearProblemKind::CoupledPressureVelocity, ...)
```

## Preuves exécutées sur cette base locale

- CMake Release / CTest : **5/5 exécutables PASS** (`test_anisotropy_energy`, `test_vanka_pressure_patch`, `test_mgr_transactional`, `test_pressure_anisotropic_affine_flux`, `test_mgr_preconditioner`). Dans les nouveaux exécutables : 8+5+4+3 = 20 cas.
- ASan + UBSan : 7 anciens cas A/B avant l'ajout du test de qualité de face, 4 cas MGR et 5 cas Vanka : **16 cas instrumentés PASS**. La dernière extension du test A/B a passé Release, mais n'a pas été rejouée sous sanitizers.
- Non-orthogonal skew : **10/10**. Campagne Laplacien non orthogonal : PASS.
- Benchmark contrôlé : **18/18** résolutions elliptiques et **4/4** couplées satisfont les limites de résidu vrai et d'erreur de solution retenues. Le CSV est livré avec les journaux. Mesures à faible répétition, non statistiques.
- La CLI MatrixMarket a produit un JSON valide et reproductible sur la matrice de démonstration (`anisotropic_rows=2`, `worst_offdiagonal_contrast=1000`).
- `git apply --check` vérifié séparément sur la base corrigée précédente et la base ZIP originale (voir livrables).

## Verrous encore ouverts pour annoncer « convergence finale CFDX »

1. Récupérer/générer les maillages NACA N64/N128/N256 et les jeux de conditions physiques associés ; les `.h5` N128/N256 ne sont pas dans l'extraction de travail ni dans le build testé.
2. Exécuter SIMPLE, SIMPLEC et COUPLED jusqu'à convergence simultanée des **résidus non linéaires physiques**, de la masse globale et locale, des efforts et du critère de stationnarité, sans relâcher les seuils.
3. Comparer `Cd`, `Cl`, `Cp`, `Cf`, `Cm`, bilans, champ proche paroi et GCI trois maillages ; ne pas inférer cette convergence des certificats linéaires.
4. Construire le Jacobien cohérent de la correction de pression skewness pour COUPLED, le tester par différentiation finie et par MMS.
5. Comparer PCD/LSC/BFBt/MGR/Vanka puis développer Braess–Sarazin sur de vraies matrices CFD ; mesurer le temps total incluant préparation/actualisations.
6. Qualifier les profils industriels anisotropes (rotation des directions, couches prismatiques, rapport d'aspect 10¹–10⁴, SST, thermique) et exécuter les CI GitHub Release/DebugSanitizers exactes.

**Conclusion :** les ajouts A/B/D et les tests C sont intégrés et vérifiés localement, mais la phase C couplée et l'objectif de convergence physique finale demeurent ouverts. Ne pas merger en prétextant une qualification complète.