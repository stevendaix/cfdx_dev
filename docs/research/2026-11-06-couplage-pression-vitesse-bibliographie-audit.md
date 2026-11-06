# CFDX — Audit bibliographique et expérimentation sur le couplage pression–vitesse

**Date :** 2026-11-06. **Base inspectée :** `cfdx_dev-master(3).zip` transmis dans la conversation. **Statut :** recherche et prototype expérimental ; **aucune qualification N9/NACA/SST nouvelle n'est déclarée**. Les commits Git et le HEAD de la PR #1182 ne sont pas identifiables à partir de ce ZIP sans métadonnées `.git`.

## 1. Distinction fondamentale

- **Couplage physique et temporel** : SIMPLE, SIMPLEC, SIMPLER/SIMPLE-R, PISO, PIMPLE, fractionnement, COUPLED.
- **Opérateurs algébriques** : bloc de quantité de mouvement `M`, gradient `G`, divergence `D`, couplage pression `C` (dont la stabilisation colocalisée).
- **Préconditionnement** : Schur, PCD, LSC, BFBT, bloc local, AMG sélectif, Vanka/Braess–Sarazin, Schwarz.

L'identité des noms ne garantit ni la même discrétisation ni les mêmes critères. Pour comparer les méthodes, **geler maillage, propriétés, BC, convection, non-orthogonalité, jauge et vraies tolérances physiques**, et isoler ensuite la seule variable d'algorithme.

## 2. Constats vérifiables dans le code fourni

| État | Emplacement | Lecture technique |
|---|---|---|
| SIX stratégies | `src/cfdx/physics/pressure_velocity_algorithms.h`, enum `PressureVelocityAlgorithm` | SIMPLE, SIMPLEC, PISO, PIMPLE, FRACTIONAL_STEP, COUPLED. |
| SIX modèles Schur | même en-tête, enum `CoupledSchurModel` | BlockLocal, PCD, LSC, BFBT, SIMPLE, SIMPLEC ; ces choix sont des modèles de préconditionnement, **pas** six couplages physiques. |
| PISO/PIMPLE : correcteurs | `steady_incompressible_solver.h`, boucle `for (corr...)` | Le prédicteur `HbyA` est initialement calculé une fois par itération externe ; la boucle de corrections garde `phiHbyA` figé, ce qui ne représente pas toutes les corrections voisines décrites dans les variantes PISO itérées. |
| Contrat symbolique | `pressure_velocity_system.h` | `PressureVelocitySystemContract` est une description et validation de métadonnées, **pas** un oracle numérique des blocs produits. |
| Rhie–Chow | `steady_incompressible_solver.h`, `make_rhie_chow_mass_flux` et `make_full_rhie_chow_pressure_rows_*` | Plusieurs chemins opérateurs existent ; comparer la dérivée réelle du flux aux blocs du Jacobien est obligatoire, notamment pour COUPLED et skewness. |
| PCD atypique | `src/cfdx/core/linalg/pcd_schur.h`, commentaire et `apply` | Action **provisoire** `-Fp^{-1} r` ; elle ne correspond pas à la factorisation PCD habituelle `Kp^{-1} Fp Mp^{-1}` (à signe près). L'appellation PCD actuelle ne suffit donc pas pour la qualification bibliographique. |
| Qualification incomplète | `docs/audits/04-couplage-pression-vitesse.md`, §7 | La fiche « APPROVED » est datée du 2026-10-12 mais marque toujours l'équivalence COUPLED/ségrégé à faire ; « APPROVED » signifie ici fiche, pas résolution des réserves numériques. |
| Oracles existants | `tests/validation/test_n9_physical_matrix.cpp` et `test_steady_incompressible_solver.cpp` | Tests Couette, Poiseuille, Ghia et chemins couplés existants. Ils doivent être rejoués sur l'état exact du code proposé avant promotion. |

### Détail PISO : manque détecté

Après l'étape de correction, le code calcule `U = HbyA - rAU * V * grad(p)` mais, par défaut, il **n'actualise pas** `HbyA` à partir des nouvelles vitesses avant le correcteur suivant :

\[
H_i/A_{ii}=A_{ii}^{-1}\left(b_i+V_i\nabla p_{\mathrm{ancien}}-\sum_{j\ne i}A_{ij}u_j\right).
\]

Lorsque `u_j` change, l'approximation de `H/A` doit changer si l'on veut prendre en compte les termes voisins. Sur certains cas orthogonaux, une nouvelle correction de pression avec prédicteur figé peut principalement recalculer un résidu déjà annulé par la première correction : le nombre de correcteurs ne prouve pas une meilleure approximation. Il faut vérifier cet effet sur un flux **non nul** et un maillage **non orthogonal**.

## 3. Prototype strictement optionnel intégré dans ce lot

`CouplingControls::piso_neighbor_correction = false` par défaut.

Si **et seulement si** l'algorithme vaut PISO ou PIMPLE, et qu'un correcteur supplémentaire est prévu, le prototype :

1. récupère la vitesse corrigée du correcteur courant ;
2. recalcule `HbyA` à partir de la **matrice de quantité de mouvement relaxée et gelée** `ex/ey/ez`, du champ corrigé et du gradient de pression source de cette matrice ;
3. reconstruit `phiHbyA` à partir du nouveau `HbyA` avant la prochaine équation de pression ;
4. conserve le schéma de reconstruction Rhie–Chow et les BC existants ;
5. expose `piso_neighbor_predictor_rebuilds` et `piso_neighbor_predictor_flux_delta_linf` dans `IncompressibleIteration` ;
6. rejette explicitement l'option avec SIMPLE, SIMPLEC, FRACTIONAL_STEP ou COUPLED.

**Ce prototype n'est pas la démonstration d'un PISO complet** : pas de preuve d'ordre temporel, de réduction d'erreur de splitting, ni de gain sur Ghia/NACA. Aucun critère physique, budget ou tolérance n'est abaissé. Le comportement existant reste inchangé lorsque l'option est désactivée.

**Exemple d'activation directe depuis un appel C++ :**

```cpp
cfdx::physics::IncompressibleSolverControls controls;
controls.algorithm = cfdx::physics::PressureVelocityAlgorithm::PISO;
controls.coupling.n_pressure_correctors = 3;
controls.coupling.piso_neighbor_correction = true; // EXPÉRIMENTAL UNIQUEMENT
```

## 4. Relecture critique de la bibliographie

| Référence | Résultat et portée | Décision CFDX |
|---|---|---|
| Étude OpenFOAM de 2017, DOI `10.1016/j.cam.2016.03.037` | Rhie–Chow original vs correction de dépendance à la sous-relaxation (Majumdar) ; SIMPLE, SIMPLE-C, SIMPLE-R, PISO ; cavité et NACA0012. **Le papier appelle SIMPLE-R « SIMPLE-Revised », qui recouvre la famille SIMPLER.** | Ne pas compter SIMPLE-R et SIMPLER comme deux inventions sans analyser la formulation exacte ; priorité à l'invariance en `alpha_u`. |
| Majumdar (1988), DOI `10.1080/10407788808913607` | Rôle de la sous-relaxation dans l'interpolation de quantité de mouvement. | Ajouter un oracle d'invariance de l'état convergé vis-à-vis de `alpha_u`, en imposant **les mêmes vraies portes**. |
| Aguerre et al. (2020), COMPLEX, DOI `10.1016/j.compfluid.2019.104396` | Correction voisine approchée par Taylor au premier ordre, gradient de correction contraint par masse ; Fourier et essais laminaires/turbulents. | Priorité secondaire : dériver la matrice locale et vérifier la conservation avant prototype. Ne pas confondre avec le simple recalcul de `HbyA` ci-dessus. |
| ESIMPLE (2021), DOI `10.1016/j.buildenv.2021.108170` | Jusqu'à 3× dans le meilleur cas étudié, bâtiments/HVAC ; pas de garantie pour NACA/SST. | Recherche P3 ; ne pas transférer un facteur de vitesse à CFDX. |
| Issa (1986), DOI `10.1016/0021-9991(86)90099-9` | Couplage PISO par correction d'opérateur. | Prototype de recalcul voisin, puis démonstration de cohérence de l'itération interne. |
| Guermond, Minev, Shen (2006), DOI `10.1016/j.cma.2005.10.010` | Distinction pressure-correction, velocity-correction, consistent splitting ; analyse d'erreur aux frontières. | Préparer rotational pressure correction **séparée**, avec BC et ordre temporel clairement spécifiés. |
| Aoussou, Lin, Lermusiaux (2018), DOI `10.1016/j.jcp.2018.06.062` | Pressure-correction itérée : tolérance de splitting relativement à l'erreur temporelle, IMEX BDF/RK et accélération Aitken. Cas étudiés sur maillages décalés structurés : la transposition sur volumes finis colocalisés doit être validée. | Évaluer une porte de splitting indépendante du résidu linéaire, pas un `nCorrectors` arbitraire. |
| Block-SAMG (2018), DOI `10.1016/j.compfluid.2018.02.034` | Coarsening *sélectif* et interpolation bloc, équation de pression prise comme matrice primaire, smoother ILUC0 ; implémenté dans foam-extend. | Prototype bloc sur le **vrai** système `M,G,D,C` après la correction du Jacobien RC. Ne pas confondre avec un AMG scalaire appliqué à `Fp`. |
| Heinlein, Klawonn, Knepper, Saßmannshausen (2026), DOI `10.1553/etna_vol65s156` | Schwarz monolithique à deux niveaux/GDSW, étude Reynolds/CFL, parallélisme ; **éléments finis + Trilinos**. | Projet moyen terme si MPI est qualifié. Référence 2026 publiée (et non seulement prépublication arXiv 2025). |
| PETSc, `PCFIELDSPLIT` et `PCLSC` | Composition Schur explicite et préconditionnement least-squares commutator, avec contrats de signe et opérateurs. | Utiliser comme **référentiel mathématique**, pas prétendre utiliser PETSc/HYPRE si les bibliothèques ne sont pas installées. |
| OpenFOAM User Guide 6.3 | PISO/SIMPLE, correcteurs non orthogonaux, jauge de pression et relaxation. | Vérifier que les correcteurs de non-orthogonalité ne sont pas confondus avec les correcteurs PISO. |

### Point de vocabulaire à corriger

**Le « PCD » du code n'est pas encore un PCD complet.** Dans `pcd_schur.h`, il est expressément défini comme `-Fp^{-1}` avec justification empirique sur une matrice de test. La formule bibliographique la plus classique est une composition `Kp^{-1} Fp Mp^{-1}` (ou une autre forme explicitement équivalente après conventions). Les frontières de `Fp`/`Kp`, la jauge et les échelles doivent être auditées indépendamment. Il faut :

- mesurer `||S^{-1}r - \widetilde S^{-1}r||` sur des petits systèmes où `S` peut être calculé exactement ;
- faire varier le flux convectif, le Reynolds, l'anisotropie et la jauge ;
- garder `Fp^{-1}` sous son nom mathématique tant qu'un PCD classique n'est pas reproduit ;
- ne pas remplacer l'actuel opérateur dans la résolution de production sans preuve de convergence et de coût total.

## 5. Matrice de qualification nécessaire

| Porte | Cas | Données à enregistrer |
|---|---|---|
| Identité opérateurs | petits volumes finis, orthogonal/skewed, dirichlet/Neumann | `||J_RC v - (phi(p+eps v)-phi(p))/eps||`, `||A x-b||` par bloc, jauge |
| SIMPLE vs SIMPLEC vs PISO/PIMPLE | Couette, Poiseuille, Ghia Re 100 et Re 1000 | vrais résidus QDM, continuité **locale**, flux autoritaire, temps total ; champs finaux à critères communs |
| Sensibilité à `alpha_u`, `alpha_p` | cavité orthogonale et skewed | invariance des solutions à convergence complète, pas seulement nombre d'itérations |
| PISO transitoire / splitting | solution manufacturée non stationnaire, Δt/2, Δt/4 | erreur U/p en norme L2 et L∞, ordres expérimentaux, effet du nombre de correcteurs |
| COUPLED / Block-AMG/Schwarz | matrice capturée de CFDX, NACA skewed et canal anisotrope | vraie action Schur, vraies itérations FGMRES, mémoire, setup vs apply, temps jusqu'aux portes physiques |
| SST + NACA N64/N128/N256 | cas Re et BC déclarés | résidus masse/QDM/SST indépendants, Cd/Cl/Cp/Cf/Cm, forces pression/cisaillement, GCI, wall-distance |

**Refus de promotion :** `MAX_ITER`, résidu linéaire excellent mais masse/SST mauvais, pression checkerboard, différence d'opérateur RC/COUPLED, évolution non monotone inexpliquée des erreurs temporelles, ou dépendance non physique à `alpha_u` doivent rester des échecs explicites. Les prototypes bibliographiques ne constituent pas des qualifications physiques.

## 6. Commandes d'exécution locale

```bash
cmake -S . -B build -G Ninja -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_RUNTIME_TESTS=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_pressure_velocity_system test_simplerc_schur test_pcd_schur test_steady_incompressible_solver -j2
ctest --test-dir build -R '^(test_pressure_velocity_system|test_simplerc_schur|test_pcd_schur|test_steady_incompressible_solver)$' --output-on-failure
```

La campagne complète N9 est plus coûteuse et doit se faire sur le HEAD exact. Ce document n'en présume pas les résultats.

## 7. Sources primaires

- https://doi.org/10.1016/j.cam.2016.03.037
- https://doi.org/10.1080/10407788808913607
- https://doi.org/10.1016/j.compfluid.2019.104396
- https://doi.org/10.1016/j.buildenv.2021.108170
- https://doi.org/10.1016/0021-9991(86)90099-9
- https://doi.org/10.1016/j.cma.2005.10.010
- https://doi.org/10.1016/j.jcp.2018.06.062
- https://doi.org/10.1016/j.compfluid.2018.02.034
- https://doi.org/10.1553/etna_vol65s156
- https://petsc.org/release/manualpages/PC/PCLSC/
- https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/
- https://www.openfoam.com/documentation/user-guide/6-solving/6.3-solution-and-algorithm-control

## 8. Résultats vérifiés de cette copie locale

Compilation **Release GCC 14 / C++17**, MPI désactivé. Quatre cibles compilées et exécutées avec `ctest --output-on-failure` :

| Cible | Résultat | Portée |
|---|---:|---|
| `test_pressure_velocity_system` | PASS | Contrat déclaratif et contrôles associés. |
| `test_steady_incompressible_solver` | PASS, **31/31 cas internes** | Contient les deux nouveaux cas PISO (restriction de l'option et réactualisation du prédicteur). Les assertions vérifient aussi que `rhie_chow_identity_linf < 1e-10` pour le cas opt-in. |
| `test_simplerc_schur` | PASS, **10/10 cas** | Opérateurs Schur SIMPLE/SIMPLEC sur petits systèmes. |
| `test_pcd_schur` | PASS, **11/11 cas** | Confirme notamment la conformité de **l'action inverse simple actuellement définie** ; ne qualifie pas un PCD complet. |

`git apply --check` sur les **préimages extraites du ZIP d'origine** passe ; l'application produit des fichiers de code identiques octet par octet à ceux inclus dans le paquet.

**Important :** le scénario à deux cellules n'est pas une validation de convergence physique. Cette qualification n'autorise ni un gain annoncé pour PISO/PIMPLE, ni un remplacement du PCD actuel. NACA N128/N256, Ghia complet, SST et convergence temporelle ne sont pas qualifiés par ces quatre cibles. L'option PISO/PIMPLE reste désactivée par défaut.

### Complément important extrait de l'article de 2017

Le papier *Influence of momentum interpolation methods…* (§4 et §6) confirme que **SIMPLE-R désigne SIMPLER/SIMPLE-Revised**, et rapporte que, sur son cas NACA avec turbulence, SIMPLE-R ne parvient pas au résidu demandé (résidus de vitesse stagnants vers `1e-4`). Les auteurs constatent généralement de meilleures performances avec SIMPLEC dans **leurs** essais, sans effet universel sur tous les maillages. La correction de Majumdar réduit notamment la dépendance de la traînée à la sous-relaxation, mais ne rend pas toute la chaîne RANS indépendante de celle-ci. **Conséquence CFDX :** la priorité de SIMPLE-R passe derrière la reconstruction flux/pression, l'invariance de solution à `alpha_u`, la convergence SST et la qualification de SIMPLEC.

### 8.1. Exécution physique supplémentaire : Ghia rapide

La cible `test_ghia_cavity_quick` compile et passe en **5,38 s**. L'exécutable `test_ghia_cavity --quick` résout en SIMPLE la cavité **Re=100, 32×32**, avec les portes physiques et la comparaison de profil codées dans le test :

| Grandeur mesurée | Valeur |
|---|---:|
| Itérations externes | 497 |
| Défaut de continuité maximal | `3.11284e-12` |
| Résidu physique de quantité de mouvement | `7.92388e-11` |
| Norme de continuité normalisée | `1.55177e-10` |
| Écart interne `phi_RC` – flux interpolé depuis U | `1.96259e-3` |
| Statut CTest | PASS |

Le dernier écart compare **deux reconstructions différentes** ; il ne prouve pas à lui seul une violation de la conservation (le défaut de continuité du flux autoritaire est petit), mais son évolution avec le raffinement spatial doit être étudiée avant toute conclusion sur la précision du Rhie–Chow colocalisé. Ce test est celui de **SIMPLE inchangé par défaut**, non un benchmark PISO expérimental. La campagne complète Ghia Re=1000, NACA et SST reste à exécuter.
