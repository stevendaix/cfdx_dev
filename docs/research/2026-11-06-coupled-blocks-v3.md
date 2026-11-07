# CFDX — Couplage pression–vitesse : lot expérimental V3 et preuves

Date : 6 novembre 2026. Base : archive **CFDX_couplage_PCD_Jacobien_sources_completes_2026-11-06.zip** (SHA-256 `81d1c565c62ef40976c25935b431a9a8fc4640883f8ad3bd869cee8fa1e1b634`).

## Statut et règle de promotion

**Prototype expérimentaux compilés et testés** ; **aucun changement du solveur/préconditionneur sélectionné par défaut**. Les sous-domaines Schwarz ne sont pas raccordés au dispatch de production. L'algorithme PCD complet reste un choix explicite, non qualifié physiquement. Les seuils de convergence et les références physiques existants n'ont pas été abaissés.

Les tests de cette livraison portent sur les **sources reconstruites depuis l'archive** et non sur le HEAD distant de la PR #1182. La CI globale n'a pas été exécutée. Les campagnes NACA N128/N256, SST/CHT, Ghia complète, MPI et performances réelles restent ouvertes.

## Briques intégrées

1. **Commutateur PCD à masse pression consistante**, fichier `src/cfdx/core/linalg/pcd_schur.h`. L'action expérimentale est `-Mp^{-1} Fp Kp^{-1} r`. Une masse complète nécessite `pressure_mass_solve` ; le résidu réel de la résolution de `Mp` est vérifié après chaque application. Mises à jour des matrices ET des solveurs groupées de façon transactionnelle, avec refus des callbacks obsolètes. La variante historique `-Fp^{-1}` reste par défaut. Une diagonale strictement positive ET une symétrie numérique sont nécessaires à la masse ; cela n'est pas une preuve générale de SPD. Le caractère physiquement acceptable du commutateur pour les opérateurs FV/limites/pressions imposées demeure à établir.
2. **Certificat directionnel du Jacobien**, fichier `src/cfdx/physics/pressure_velocity_jacobian_certificate.h`. Produit `Jv` par les blocs `(Auu,G,D,C)` comparé aux différences centrales d'un **résidu indépendant fourni par callback**. Directions déterministes, contrôle des entrées et mesures par ligne. Le coût est de deux évaluations du résidu par direction. Il s'agit d'un contrôle échantillonné, et le callback du résidu CFD de production doit encore être raccordé pour qualifier toutes les branches.
3. **Certificat d'inverse Schur**, fichier `src/cfdx/core/linalg/schur_inverse_quality_certificate.h`. Mesure `||S P_S^{-1}r - r||`, avec `S = C - D Auu^{-1}G`. Le callback de résolution de `Auu` est certifié par le résidu réel avant usage. Cet oracle distingue une mauvaise action Schur de l'échec d'un solveur de référence.
4. **Certificat structurel COUPLED CSR**, fichier `src/cfdx/physics/coupled_assembly_certificate.h`. Contrôle dimensions, valeurs finies, non-vacuité des lignes, contenus des blocs `Auu,G,D,C`, jauge pression et accord des diagonales Schur, puis rapporte `Ax-b` **au point de départ**. Les défauts linéaires rapportés ici ne constituent pas des mesures de conservation physique.
5. **Certificat flux COUPLED/Rhie–Chow**, fichier `src/cfdx/physics/coupled_flux_equivalence_certificate.h`. Confronte les lignes de continuité réellement assemblées à la divergence des flux conservatifs **reconstruits indépendamment**, en respectant owner/neighbour, conditions limites et cellule de jauge. Le contrôle de conservation est indépendant et optionnel. Cas `absolute_tolerance=0` et résidus exactement nuls correctement géré.
6. **Schwarz à recouvrement par cellules**, fichier `src/cfdx/core/linalg/cell_overlapping_schwarz.h`. Implémentation **sérielle de référence** du préconditionneur, pour l'ordre global `[Ux(0..Nc-1), Uy, Uz, p]`, sans PETSc ni HYPRE installés. Décomposition en blocs de cellules adjacentes, extension sur le graphe de couplage CSR (jusqu'à 8 sauts), LU locaux pivotés et normalisés par ligne, mode RAS et additive pondérée par multiplicité. `setup` et `update_values` sont transactionnels. Le maximum par défaut est 64 inconnues par sous-domaine et une limite dure de 256 ; un sous-système singulier échoue explicitement. Ce préconditionneur n'est **pas un Schwarz multigrille/à espace grossier**, n'a pas de halo MPI, et n'est **pas activé dans les campagnes physiques**.

## Certificats sur l'assemblage de production

`tests/validation/test_n9_s9_rhie_chow_consistency.cpp` utilise la capture originale de l'assemblage du solveur COUPLED, sur un écoulement non nul et des pressions imposées. Il ajoute :

- Contrôle structurel de la matrice passée au Krylov : `gradient_nnz=84`, `divergence_nnz=36` ; contrôle structurel accepté.
- Reconstruction indépendante de la divergence des flux de face et comparaison aux lignes de pression évaluées à l'état final : `matrix_flux_linf = 8.88178e-16`, `continuity_linf = 8.88178e-16` ; certificat accepté sous la tolérance originelle de `1e-8`.
- Le résidu de continuité **au point initial** du certificat structurel vaut `2.14159` ; c'est une valeur diagnostique et **n'invalide pas** la comparaison à l'état final.

Ces résultats n'extrapolent pas l'équivalence à d'autres maillages, frontières ou NACA.

## Campagne exécutée — machine locale

| Exécution | Résultat | Limite |
|---|---:|---|
| CTest ciblé : PCD, Jacobien, Schur, CSR, flux, Schwarz, solveur stationnaire et N9/S9 | **8/8 réussies** | Pas une CI complète |
| Tests unitaires six composants | **50/50 réussis** | Matrices tests et contrôles négatifs |
| `test_steady_incompressible_solver` | **32/32 réussis** | Compilation du binaire en `-O0`, modèle Release avec `-DNDEBUG` |
| `test_n9_s9_rhie_chow_consistency` | **18/18 réussis** | Même stratégie de compilation `-O0`, binaire Release `-DNDEBUG` |
| Les six tests unitaires (avant dernier renforcement des cas limites) | **6/6 binaires sous ASan+UBSan** | Voir logs actualisés du livrable final |
| Ghia `--quick` 32x32 Re=100 | **Non concluante : limite d'exécution du binaire `-O0`** | Compilation réussie et progression à 200 itérations, mais pas de résultat convergé ; ne pas compter comme PASS |

## Test algébrique Schwarz avec couplage intercellulaire non nul

Sur une chaîne synthétique de trois cellules et douze inconnues, avec couplage spatial `-1.8` et correction locale calculée depuis le même RHS :

- une application sans recouvrement : `||A z - b||_2 = 10.1161` ;
- une application RAS avec recouvrement d'un saut : `||A z - b||_2 = 2.54127`.

Il ne s'agit **ni** d'un gain en temps CPU global, **ni** d'un benchmark NACA, **ni** d'une preuve d'accélération de convergence Krylov. Pour les sous-domaines locaux singuliers — probables dans des formulations pression non stabilisées — l'algorithme refuse l'initialisation ; aucun défaut physique n'est dissimulé par régularisation arbitraire.

## Vérifications supplémentaires avant le merge

1. Brancher les certificats directionnels sur le **véritable résidu de production** : `Jv` contre différences finies aux faces skewées, limites pression et maillages anisotropes ; test de précision et d'échelle du pas fini.
2. Comparer les valeurs et l'application Schur sur matrices capturées d'un COUPLED convergé avec débit non nul ; calculer le coût de setup, application, itérations GMRES/FGMRES, et le résidu vrai.
3. NACA N64/N128/N256 avec mêmes tolérances physiques (masse locale, forces pression/cisaillement, Cp/Cf, GCI) ; refuser les résultats non convergés.
4. Ghia Re=100/1000, canal 3D, couche limite anisotrope, transitoire PISO/PIMPLE et bilans SST ; sans promotion sur un cas synthétique.
5. CI Release complète, DebugSanitizers complet, PyFluent et essais MPI. En particulier vérifier les régressions connues hors couplage.
6. Une extension Schwarz sérieuse devra ajouter partitionnement de graphe, solveurs locaux plus économiques, coarse level, interface MPI et budget mémoire mesuré, puis démontrer son gain contre Block Jacobi et AMG.

## Bibliographie, éléments réutilisés

- Kay, Loghin, Wathen, *A preconditioner for the steady-state Navier–Stokes equations*, SIAM Journal on Scientific Computing 24(1), 237–256, 2002 : justification du commutateur pression et limites théoriques.
- Elman, Silvester, Wathen, *Finite Elements and Fast Iterative Solvers*, 2e éd., Oxford, 2014 : complément de Schur / PCD / LSC et influence des opérateurs.
- Cai & Sarkis, *A restricted additive Schwarz preconditioner for general sparse linear systems*, SIAM Journal on Scientific Computing, 21, 792–797, 1999.
- PETSc, `PCASM` : https://petsc.org/release/manualpages/PC/PCASM/ ; `PCASMType` : https://petsc.org/main/manualpages/PC/PCASMType/ ; `PCASMSetOverlap` : https://petsc.org/release/manualpages/PC/PCASMSetOverlap/ . Les références décrivent des options PETSc ; le prototype CFDX est indépendant et limité à une implémentation série pédagogique.
- PETSc, `KSP` et préconditionnement Schwarz : https://petsc.org/main/manual/ksp/ .
- PETSc, `SNES` : documentation de la vérification des Jacobien et produits directionnels : https://petsc.org/release/manualpages/SNES/SNESTestJacobian/ .
- Guermond, Minev et Shen (2006), *An overview of projection methods for incompressible flows* : reste pertinent pour la qualification future des corrections temporelles, non implémentée dans ce lot.

## Reproduire depuis la source

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF
cmake --build build --target test_pcd_schur test_pressure_velocity_jacobian_certificate test_schur_inverse_quality_certificate test_coupled_assembly_certificate test_coupled_flux_equivalence_certificate test_cell_overlapping_schwarz test_n9_s9_rhie_chow_consistency test_steady_incompressible_solver
ctest --test-dir build -R '^test_(pcd_schur|pressure_velocity_jacobian_certificate|schur_inverse_quality_certificate|coupled_assembly_certificate|coupled_flux_equivalence_certificate|cell_overlapping_schwarz|n9_s9_rhie_chow_consistency|steady_incompressible_solver)$' --output-on-failure
```

Le test N9/S9 et celui du solveur stationnaire ont été compilés séparément en `-O0` avec les autres macros et options de Release inchangées, afin d'éviter un surcoût de compilation `-O2` non nécessaire à l'exécution des tests. Cela doit être revalidé avec la compilation Release complète de la CI.
