# CFDX — Suite des phases A–D : instrumentation, AMG, skewness et Schur

**Base de référence :** sources A–D du lot précédent (`2188f2f` dans la reconstruction locale), elles-mêmes dérivées du ZIP utilisateur. Le HEAD GitHub distant **n'est pas** modifié ou qualifié par ce document.

## Contrats invariants

- Le résultat d'un solveur doit satisfaire le **résidu physique vrai** `||b-Ax||/||b||` de la matrice d'origine, et la précision de la solution manufacturée doit être contrôlée indépendamment.
- Aucun seuil de convergence, aucune limite d'itérations et aucune protection de conservation non linéaire ne sont relâchés.
- Aucun nouveau préconditionneur n'est sélectionné automatiquement : les nouvelles méthodes nécessitent un choix explicite et FGMRES quand elles ne sont pas certifiées SPD.
- Aucun résultat sur un système manufacturé ne prouve la convergence physique complète du NACA0012, de SST ou de COUPLED sur maillage industriel.

## Phase A — Diagnostic raccordé aux équations assemblées

`LinearSolverRequest` accepte désormais `diagnose_anisotropy`, `anisotropy_equation` et un pointeur géométrique temporairement emprunté. `solve_linear_system` mesure l'anisotropie sur **l'opérateur original**, avant mise à l'échelle et avant le préconditionnement, et renvoie `LinearSolveReport::anisotropy_json`.

Exemple :

```cpp
cfdx::core::LinearSolverRequest request;
request.diagnose_anisotropy = true;
request.anisotropy_equation = "pressure";
// request.anisotropy_geometry = &geometry; // facultatif, doit vivre pendant l'appel
const auto result = cfdx::core::solve_linear_system(
    matrix, rhs, x, cfdx::core::LinearProblemKind::PressurePoisson,
    request, 500, 1e-10);
std::cout << result.anisotropy_json << '\n';
```

`IncompressibleSolverControls::anisotropy_report_callback` expose en option les rapports `momentum_x`, `momentum_y`, `momentum_z`, puis `coupled_full_4N`, pour **la matrice couplée réellement assemblée**. Les trois blocs scalaires utilisent les centres de cellules et la géométrie réelle des faces internes (non-orthogonalité, skewness). Le système 4N complet n'hérite **pas** fictivement d'une géométrie à N cellules. L'instrumentation ne se déclenche pas sans callback. Pour une campagne automatisée, l'appelant écrit les résultats JSON où il le souhaite.

**Reste** : intégrer un export de telles preuves dans les campagnes physiques NACA/SST ; l'analyse opt-in n'est pas une nouvelle correction géométrique du maillage.

## Phase B — Interpolation de distance deux à coût plafonné

`AMGInterpolationPolicy::EnergyMinimizingAggregation` : support maximal du prolongateur ramené de 12 à 8 coefficients par ligne, trois pas projetés au lieu de cinq, choix initial de pas depuis le produit projeté direction/gradient et une approximation diagonale de la courbure. Le **vrai** `trace(P^T A P)` est ensuite toujours évalué ; seules les corrections énergétiquement décroissantes sont acceptées, `P 1 = 1` et la sauvegarde du prolongateur original restant obligatoires. Cette procédure est une **minimisation approximative**, et non un solveur d'optimisation énergétique exact.

Test : neuf scénarios `test_anisotropy_energy` (géométrie, CSR, interpolation, Neumann projeté et FGMRES oblique). Le benchmark `cfdx_benchmark_energy_vanka` certifie les résidus vrais et les erreurs de solution, y compris 1 024 inconnues.

Exemple sur 1 024 inconnues, anisotropie oblique de 1 000 (exécutions indépendantes, non stabilisées thermiquement) : ancienne expérimentation `113,26 ms` / `13` itérations, nouvelle expérimentation `59,97 ms` / `8` itérations. Le natif reste bien plus rapide dans le même benchmark : `11,89 ms` / `6` itérations. **Cette variante n'est donc pas promue automatiquement.** La baisse de temps ne constitue pas encore une mesure statistiquement robuste.

Les idées de limitation du prolongateur et d'interpolation à portée deux sont cohérentes avec les recommandations HYPRE BoomerAMG, mais CFDX utilise ici une implémentation native différente, sans appeler HYPRE.

## Phase C — Jacobien-direction du gradient corrigé en skewness

`gauss_pressure_gradient_jacobian_action` différencie exactement la reconstruction Gauss de production **pour géométrie fixe et conditions aux limites homogénéisées**. Elle est confrontée à la différence finie du même opérateur pour les deux cellules skewées du test `test_pressure_skew_gradient_affine`, sans paramètre de tolérance de convergence modifié.

**Verrou toujours ouvert :** cette dérivée du gradient n'est **pas** le Jacobien complet des flux continuité/momentum/Rhie–Chow. L'option `pressure_skew_corrected_gradient` dans le chemin COUPLED doit rester protégée tant que le flux implicite complet, la conservation, la convergence et les tests N9/NACA ne sont pas qualifiés. Il serait incorrect de l'activer uniquement parce que la dérivée d'un sous-opérateur passe les tests.

## Phase D — Correction Block-Schur et Braess–Sarazin opt-in

Un défaut d'algèbre dans le préconditionneur Schur existant était présent : lorsque `CoupledSchurVelocityApproximation::Diagonal` était demandé, le Schur `C-DM^-1G` utilisait l'inverse diagonal, tandis que le correcteur de vitesse continuait parfois à consommer l'inverse du bloc 3x3 complet. Les deux actions utilisent maintenant le **même inverse approché**. Les diagonales de quantité de mouvement non positives ou non finies sont rejetées.

`PreconditionerModel::BraessSarazin` est un préconditionnement **de type** Braess–Sarazin à inverse diagonal de vitesse `M=alpha diag(Auu)` et Schur pression AMG, `alpha` fini dans `[0.25,4]`, compatible uniquement avec FGMRES couplé. Le modèle est accepté dans le chemin COUPLED de production seulement si `CoupledSchurModel::BlockLocal` est demandé. Il n'est pas activé par défaut.

```cpp
controls.coupled_linear_solver.krylov = cfdx::core::KrylovModel::FGMRES;
controls.coupled_linear_solver.preconditioner = cfdx::core::PreconditionerModel::BraessSarazin;
controls.coupled_linear_solver.braess_sarazin_alpha = 1.0;
controls.coupling.schur_model = cfdx::physics::CoupledSchurModel::BlockLocal;
```

Le test `test_braess_sarazin` compare l'assemblage de `C-D(alpha diagAuu)^-1G` à un calcul indépendant, vérifie six résolutions manufacturées et la précision de la solution et refuse les combinaisons non qualifiées.

Les tests historiques `test_pcd_schur`, `test_lsc_bfbt_schur`, `test_coupled_block_schur_amg` et `test_schur_quantitative_qualification` sont exécutés en Release. **Limite critique :** le test de qualification quantitative LSC/BFBt passe son assertion mais affiche encore des erreurs relatives à l'oracle Schur exact importantes (jusqu'à ~39). La précision de ces approximations est donc **non démontrée** sur ces cas. Ne pas déduire leur qualification de la seule couleur verte du test. Les comparaisons avec matrices réellement assemblées et le couplage C skew restent ouverts.

## Validation réellement exécutée sur la copie locale

- CMake/Ninja en configuration Release : nouveaux exécutables compilés.
- `ctest` ciblé : **9/9** (diagnostic/AMG, Braess, Vanka, PCD, LSC/BFBt, Schur, flux skew).
- `test_anisotropy_energy` : **9/9 scénarios** ; `test_braess_sarazin` : **3/3 scénarios** ; `test_pressure_skew_gradient_affine` : PASS pour affine, différence finie et Rhie–Chow.
- `cfdx_benchmark_energy_vanka` : **22/22 résolutions** avec vérification du résidu vrai et de l'erreur analytique ; fichier CSV joint au lot.
- **Pas de CI GitHub complète**, pas de compilation Sanitizers validée dans cette passe, pas de validation finale NACA0012 N128/N256/512 ni de vérification des efforts, des bilans et du GCI.

## Suite impérative pour qualification physique

1. Rejouer Release/DebugSanitizers et tous les tests Python sur le HEAD exact qui recevra le patch ; réparer les échecs reproductibles.
2. Assembler le Jacobien conservatif complet de Rhie–Chow et de la continuité, incluant la correction skew ; comparer produit `Jv` et différences finies de **tout le résidu**.
3. PCD/LSC/BFBt et Vanka/Braess : mesurer `||r||` de la matrice **CFD réelle**, erreur de Schur indépendamment de son algèbre interne, conservation, temps setup/apply, mémoire, et absence de stagnation.
4. Sur NACA N64/N128/N256 (puis N512), converger avec les critères inchangés, mesurer Cd, Cl, Cm, Cp, Cf, tractions pression/cisaillement, bilans et GCI ; conserver N9, Ghia, SST et thermique dans les gates physiques.

Références méthodologiques : documentation officielle hypre BoomerAMG (coarsening/interpolation tronquée et MGR) et PETSc PCFIELDSPLIT/Schur. **Aucun code de PETSc ou HYPRE n'est lié dans ce lot.**
