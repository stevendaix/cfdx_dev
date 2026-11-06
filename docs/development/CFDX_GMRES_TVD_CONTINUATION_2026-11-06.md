# CFDX — accélération GMRES et convergence TVD par continuation (2026-11-06)

**Statut : expérimental, opt-in, qualification partielle.** Base : archive `cfdx_source_optimise.zip` du dialogue, pas le HEAD GitHub. Aucun push/merge. Aucune tolérance physique ou linéaire changée.

## Problème numérique, correction et contrat

1. **GMRES/Arnoldi :** le MGS en deux passes recopiant à chaque étape une colonne de l'espace de Krylov en `Vector` temporaire a été remplacé par un produit scalaire sur vue contiguë. La politique de réduction MPI, sa plage d'appartenance et l'arithmétique FP32/FP64 sont conservées. La solution reste certifiée avec le vrai résidu de l'opérateur original.
2. **Convection :** le schéma implicite conservatif UPWIND reste le même, tandis que la correction explicite d'ordre supérieur est multipliée par une fraction `lambda` entre 0 et 1. Le flux reconstruit aux faces internes applique la même fraction. Pour chaque face `o/n`, les contributions différées `(-lambda F (phi_high - phi_up), +lambda F (phi_high - phi_up))` s'annulent en somme. Aucun schéma partiel ne peut être déclaré convergé.
3. **Plan de continuation :** le nouveau contrôleur impose un préchauffage UPWIND, une rampe puis un minimum de plusieurs itérations à 100 % d'ordre supérieur. Il refuse les budgets insuffisants, les paramètres invalides ou un solveur COUPLED/transitoire/SST non qualifié pour ce chemin. Défaut `enabled=false`.
4. **Résidu de référence :** le ratio entre résidu courant et résidu de la première itération n'est pas valide lorsque ces deux itérations résolvent des opérateurs convectifs différents. Le chemin opt-in ne l'utilise plus comme porte finale; il conserve la certification indépendante du résidu physique RHS-normalisé par composante, de la pression, de la continuité et des variations de champ. Le moniteur reçoit alors comme limite du maximum synthétique dimensionnellement normalisé la tolérance **relative déjà demandée**, sans modifier les limites physiques.
5. **Traçabilité :** `high_order_fraction` et `full_order_eligible` figurent dans les traces incompressibles. Le CSV de qualification NACA ajoute `high_order_fraction` et ne remplace pas les contrôles de force ni de GCI.

## Tests exécutés sur la copie locale

- CTest Release ciblé : `test_gmres_solver`, `test_incompressible`, `test_finite_volume_transport`, `test_convection_continuation` : **4/4 executables**.
- Tests internes : GMRES **29/29**, convection continuation **5/5**, transport FV **25/25**, incluant les fractions 0, 0.25, 0.5, 1 et les flux assemblage/reconstruction; validation des contrôles SIMPLE/COUPLED.
- ASan + UBSan en Clang sur les exécutables autonomes : GMRES **29/29**, continuation **5/5**, sans signalement dans les chemins exécutés. FV ASan n'a pas été terminé : la tentative de compilation autonome a échoué à l'édition des liens (AMG hors unité), *pas* pendant le test numérique. La CI sanitizers complète n'a pas été exécutée.
- Benchmark GMRES (Clang `-O3`, 9 répétitions dont une chauffe, médianes de 8, même machine) :

| Inconnues | Itérations | Résidu physique vrai | Baseline | Optimisé | Accélération |
|---:|---:|---:|---:|---:|---:|
| 1 024 | 111 | 8.40823e-13 | 7.22338 ms | 6.46950 ms | 1.12× |
| 2 304 | 111 | 8.62683e-13 | 18.729 ms | 16.9106 ms | 1.11× |
| 4 096 | 103 | 9.64677e-13 | 32.591 ms | 26.7554 ms | 1.22× |

Ces microbenchmarks synthétiques ne prédisent pas le temps total NACA/N256 ni les performances MPI/GPU.

## NACA0012, Re = 1000, N64, TVD van Leer, SIMPLE

Commande de replay :

```bash
CFDX_NACA_CONVECTION_VANLEER=1 \
CFDX_NACA_TVD_CONTINUATION=1 \
CFDX_NACA_IMPORT_INITIAL_STATE_CSV=/path/to/upwind_checkpoint.csv \
CFDX_NACA_PHYSICAL_HISTORY_CSV=history.csv \
./build/test_naca0012_qualification --single naca0012_laminar_n64.h5 0.45 0.2 1000
```

Le CSV de redémarrage doit correspondre au maillage exact; le redémarrage est une **initialisation**, pas une certification. Conditions : 0.45/0.2, résidu physique et continuité relative à 1e-8, contrôle par composante, pression linéaire tolérance 1e-10, continuation 40+80 itérations et maintien 4 itérations. La fraction est 1 à partir de l'itération 120. Sans correction du résidu de référence, les 1000 itérations épuisent le budget; au dernier pas, momentum RHS relatif = 2.75002e-9, composante max = 7.66654e-9, continuité = 1.36246e-12, variations dU = 5.907e-9 / dP = 1.80929e-9, mais le ratio du résidu initial (un autre opérateur) vaut 1.80807e-5.

**Avec le correctif**, `converged=1` à l'itération **726** et `high_order_fraction=1` : temps solveur 35.0856 s, Cd=0.0792732, Cl=-1.81662e-12, Cd_p=0.0327556, Cd_v=0.0465177, continuité rapportée par NACA=5.34566e-10. Le mode `--single` prouve cette convergence **itérative** sur un maillage; il ne valide ni la valeur de Cd ni la convergence spatiale.

## Ce qui reste obligatoire, sans assouplir les critères

- Reproduire la même qualification TVD ordre complet sur N128/N256 et calculer Cp/Cf/Cd, bilan de masse cellulaire, temps complet et GCI (trois grilles convergées). Réexaminer les gradients/tractions de paroi indépendamment.
- Quantifier la stagnation N256 (couplage non linéaire et qualité du maillage), tester pseudo-transitoire CFL adaptatif et continuation non orthogonale avec restauration complète, puis comparer SIMPLE/SIMPLEC/COUPLED.
- Phase B : qualification AMG oblique/variable 3D sur matrices de CFD réelles, avec setup + solve + mémoire + vraie erreur et comparaison à PETSc/HYPRE officiels quand disponibles.
- Phases C/D : assemblage Jacobien conservatif complet de skewness en COUPLED, Schur approché sur matrices réelles, SST/thermique, tests Python et CI GitHub sur le HEAD exact.

## Reproduire les validations locales

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON
cmake --build build --target test_gmres_solver test_incompressible test_finite_volume_transport test_convection_continuation -j 2
ctest --test-dir build -R '^(test_gmres_solver|test_incompressible|test_finite_volume_transport|test_convection_continuation)$' --output-on-failure
cmake --build build --target cfdx_benchmark_gmres_arnoldi
./build/cfdx_benchmark_gmres_arnoldi 4096
```

**Attention :** le benchmark compare un exécutable ancien et un nouveau compilés avec les mêmes drapeaux; le seul exécutable courant ne mesure pas automatiquement le gain par rapport à l'ancienne version.
