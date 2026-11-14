# CFDX V28 — sous-espace Newton/Cauchy, rayon adaptatif, relance en stagnation

Date : 2026-11-10. Base de travail : sources V27. Statut : **expérimental / non qualifié physiquement**. Aucune tolérance de convergence finale ne change ; pas de push GitHub.

## Diagnostic qui motive la modification

En V27, sur le BFS 712 cellules, le candidat Cauchy à la dixième itération SST contracte le vrai résidu dynamique (`Linf=0.921785`, `L2=0.175168`) mais perd devant Newton–Jacobi (`Linf=0.889907`). En V27, le rayon `0.5` est réinitialisé à chaque essai et la recherche spatiale n'est exécutée qu'aux itérations 1, 10, 20, ... . En campagne prolongée, un pas Cauchy à l'itération 20 fait descendre `Linf` de `0.827229` à `0.824035`, puis les trois itérations 21–23 sont refusées : le déclenchement sur stagnation est absent.

## Développement V28

- `sst_adaptive_trust_v28.h` construit trois mélanges **logarithmiques** des états candidats Cauchy et Newton–Jacobi (`25/75`, `50/50`, `75/25`). Il tente au plus deux rayons par mélange et mesure la pente de l'objectif lisse **sur l'opérateur FV projeté et réassemblé**, pas sur le préconditionneur. La taille du pas est bornée également **après** la projection des fonctions de paroi.
- `Phi = 0.5 * L2_dynamic^2`, `rho=(Phi(x)-Phi(x+s)) / (m(0)-m(s))`. Le modèle `m` est un modèle directionnel linéaire sondé, pas le modèle quadratique complet PETSc. Un candidat exige simultanément positivité, convergence de la projection de paroi, `rho > 1e-3`, décroissance de `L2` et du **Linf physique existant**.
- `SSTV28AdaptiveTrustState` conserve un rayon par invocation du transport SST (`0.5` initial, min `1/64`, max `1.5`). **Seul un état effectivement choisi par la recherche physique finale avec alpha=1** peut augmenter le rayon ; cela inclut un candidat V27 qui gagne réellement. Si `rho < 0.25`, contraction `x0.5`; si `rho >= 0.75` et le pas **effectivement projeté** sature le rayon, expansion `x2`. Pas d'expansion sur une proposition perdante.
- **Correction de cohérence de la politique adaptative :** un mélange V28 rejeté lors des essais ne diminue **pas** le rayon avant la sélection finale : la proposition Cauchy V27 ou Newton–Jacobi peut encore gagner. La diminution du rayon intervient uniquement après le rejet de **toute** la mise à jour SST physique. Cela évite d'invalider une éventuelle correction V27 mesurée dans l'ancien rayon.
- Déclenchement de récupération ciblée : lorsqu'un pas a été refusé, V28 autorise au plus deux recalculs extraordinaires du préconditionnement spatial et du sous-espace avant la protection de stagnation historique. Le déclenchement n'est pas global et n'est **pas** actif sans l'option V28.
- Le chemin historique Picard, la direction Newton–Jacobi V22, les variantes V23/V24 et le test final des résidus physiques sont conservés. Le candidat V28 n'est jamais imposé s'il n'est pas le meilleur de la comparaison finale. Les états sont construits dans des copies locales et la remise en état du solveur reste active.
- Option HDF5 : `numerics.raw_settings.sst_adaptive_trust_v28_experimental=true`. Elle requiert V27 et toutes ses dépendances V26/V25/V24/V22. Désactivée par défaut.

## Résultats — statut de la validation

**BFS 712 cellules, dix itérations SST, première tentative extérieure, mêmes jeux de données HDF5 :** candidat Cauchy V27 `Linf=0.92178503383`, `L2=0.17516771348`; mélange V28 à 50 % Newton `Linf=0.91267743812`, `L2=0.17254039398`, `rho=0.95837839491`. Newton–Jacobi reste sélectionné (`Linf=0.88990745136`) : **pas de convergence SST ni BFS**.

**BFS 712 cellules, budget diagnostique porté à 40 (non atteint en pratique à cause de la stagnation) :** témoin sans V28 : dernier état intérieur `Linf=0.82403522517` sélectionné à l'itération 20 ; refus 21–23, rollback du calcul extérieur. V28 avec relance de stagnation : des pas supplémentaires aux itérations 22 et 24 contractent le défaut jusqu'à environ `0.8232783`, puis de nouvelles itérations sont refusées et le calcul extérieur est rejeté. Après un contrôle supplémentaire du calendrier des mises à jour du rayon et un nouveau rejeu de production, les mêmes valeurs sont retrouvées : rayon 0,5 à l'itération 20 ; réduction à 0,25 au rejet total de l'itération 21 ; proposition Cauchy V27 effectivement retenue à l'itération 22. Les nouvelles réductions arrivent de V27 à l'itération 22 et de V25 à l'itération 24, **pas d'un mélange V28 finalement retenu**. Ces résultats sont des **diagnostics sur un transport intérieur non convergé**, pas une preuve de contraction extérieure, ni une qualification du BFS 20 540 cellules. Conserver les journaux complets en preuves.

Tests ciblés : `test_sst_adaptive_trust_v28`, `test_bfs_sst_v28_protocol`, tests de compatibilité V24–V27 ; compilation et liaison de production diagnostique GCC 14.2.0 `-O0`. Les premiers résultats ont été obtenus avec tous les seuils originaux ; aucune méthode d'acceptation du résidu linéaire n'est substituée à la norme FV finale.

## Comparaison bibliographique (faits vs transposition)

| Référence | Mécanisme publié / visible | Usage retenu dans CFDX |
| --- | --- | --- |
| OpenFOAM Foundation v11 `omegaWallFunctionFvPatchScalarField.C` | `updateCoeffs()` actualise `G` et `omega`, `manipulateMatrix()` utilise `matrix.setValues(..., weights)` | Garder les dérivées / projections de paroi dans les évaluations physiques, mais **aucune équivalence OpenFOAM revendiquée** |
| OpenFOAM user guide `relaxationFactors` | Sous-relaxation des champs / équations pour les itérations stationnaires | Voie de comparaison future distincte de la région de confiance |
| SU2 Setup Guide RANS | Pseudo-temps implicite, `CFL_ADAPT`, `CFL_REDUCTION_TURB`, ILU/FGMRES, critères non linéaires | Inspiration pour un futur essai de continuation pseudo-temporelle **séparé**, sans copier les réglages d'un écoulement compressible |
| Fluent Theory Guide (Ansys) | Sous-relaxation explicite/implicite, pseudo-temps local pour ségrégé et global pour certains solveurs couplés | Principes documentés ; Fluent est propriétaire, son code source n'a pas été audité |
| PETSc `SNESNEWTONTR` | Rapport réel/prédit, rayon adaptatif, variantes Cauchy/dogleg et paramètres `eta1/eta2/eta3` | Politique inspirée des principes ; modèle V28 linéaire local, pondérations et contrôles propres à CFDX ; **non équivalent à SNES** |

Sources primaires :
- OpenFOAM source : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html
- OpenFOAM relaxation : https://doc.openfoam.com/2606/tools/processing/numerics/solvers/under-relaxation/
- SU2 : https://su2code.github.io/docs_v7/Setup-Guide-RANS-External-Aero/ ; https://su2code.github.io/docs_v7/Linear-Solvers-and-Preconditioners/
- Ansys Fluent Theory : https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_time_theory_seg.html ; https://ansyshelp.ansys.com/public/views/secured/corp/v252/en/flu_th/flu_th_sec_pseudo_auto.html
- PETSc : https://petsc.org/release/manualpages/SNES/SNESNEWTONTR/ ; https://petsc.org/release/manualpages/SNES/SNESNewtonTRSetUpdateParameters/

## Qualifs restant ouvertes

1. 800 itérations SST sur BFS réduit, puis au moins trois états extérieurs effectivement acceptés avec **résidu global indépendant** et vrais bilans ; le budget seul ne constitue pas une solution.
2. Comparer coûts CPU, réassemblages et scénarios AMG/PCD/LSC/BFBt/MGR sur BFS de 20 540 cellules ; chaque variante doit fournir sa preuve propre. La première pression compensée V19 ne certifie pas toutes ces variantes.
3. Équivalence `omegaWallFunction` : `updateCoeffs` / `manipulateMatrix`, coins, production `G`, viscosité `nu_t`, balises de mise à jour ; confrontation OpenFOAM contrôlée.
4. Validations NASA/Driver–Seegmiller sur solutions **globalement convergées**, avec bon champ de pression statique et comparaison de `Cf`, `Cp`, profils de vitesse et contraintes de Reynolds ; tests NACA et thermiques de non-régression.
5. Vérification multi-cœur/GPU/MPI à faire : ce lot reste CPU séquentiel diagnostique, il ne qualifie pas les performances.
