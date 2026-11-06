# CFDX V6 — Galerkin numérique, stagnation NACA et forces aux parois

**Statut : vérification numérique ciblée PASS ; NACA N256 et GCI FAIL-CLOSED.**
Base locale V5 : `493f50f`. Aucune qualification physique complète NACA 3 grilles ni merge distant n'est revendiqué.

## Changements dans le code

1. `GalerkinCSRProduct` réutilise ses buffers de lookup, AP et PtAP entre les actualisations. La table d'indices est remise à zéro uniquement après un calcul interrompu ; tous les marqueurs posés pendant un calcul réussi sont effacés avant son retour. Les échecs ne modifient jamais la matrice CSR cible. 40 actualisations comparées à l'implémentation Galerkin de référence et un échec injecté sont testés.
2. Le détecteur commun `ConvergenceMonitor` n'annonce plus une stagnation non linéaire si un défaut de conservation hors tolérance continue à diminuer au-dessus du seuil d'amélioration. Un vrai plateau des deux métriques reste un arrêt `STAGNATED`. Les tolérances numériques ne changent pas.
3. Le pilote NACA exporte séparément les historiques non linéaires (même sur échec) et un certificat de forces (uniquement après convergence). Les anciens certificats, champs convergés, profils et historiques sont effacés avant chargement d'un nouveau maillage. Un précontrôle refuse toute destination égale à une entrée ou déjà utilisée comme autre sortie.
4. Un outil `naca_mesh_similarity_gate.py` vérifie les facteurs de raffinement tangentiel/radial et de première couche normale ; `naca_force_mesh_qualification.py` exige certificats convergés, historiques cohérents et similarité avant tout GCI classique. Une exécution non qualifiée retourne le code 2 et écrit un diagnostic `BLOCKED`.

## Vérifications de convergence CFD

| Paramètre | N64 | N128 | N256 |
|---|---:|---:|---:|
| Cellules | 3 096 | 12 336 | 49 248 |
| Itérations non linéaires de la campagne de base | 335 | 808 (continuation N64) | 100 au maximum testé (continuation N128) |
| Statut | CONVERGED | CONVERGED | **NOT CONVERGED** |
| Cd, gradients aux cellules | 0.0864586 | 0.127905 | absent |
| Cd, traction reconstruite à la face (A/B) | 0.1341017143 | 0.1472343475 | absent |
| Cd pression, traction de face | 0.0280127683 | 0.0442498614 | absent |
| Cd visqueux, traction de face | 0.1060889460 | 0.1029844861 | absent |
| Cl, traction de face | -3.45e-10 | -4.76e-10 | absent |

Les deux reconstructeurs de force utilisent la **même solution CFD convergée** par niveau. La variation N64->N128 est +47,9 % pour la formule historique, contre +9,8 % pour la reconstruction de face. Cette comparaison ne valide PAS la valeur physique de Cd, car il n'existe ni N256 convergé ni validation indépendante du cisaillement et de la pression.

La comparaison des profils à la face par régions de corde révèle que `x/c` dans `[0, 0.05]` explique +0.0129757 des +0.0162371 de différence de Cd pression N128-N64 (**environ 80 %**). La contribution visqueuse totale diffère de -0.00310446 ; elle présente des compensations spatiales importantes.

### N256 : faux plateau corrigé, verrou physique toujours ouvert

À paramètres constants, l'ancienne logique annonçait `STAGNATED` à l'itération 58, pendant que la continuité s'améliorait encore. Après correctif, 100 itérations sont acceptées et l'arrêt est `MAX_ITERATIONS` ; les historiques sont **identiques jusqu'à l'itération 58 incluse**. À l'itération 100 : continuité normalisée 6.80743e-08 ; résidu de momentum relatif à la référence 0.0291828 ; résidu linéaire de pression ~3.97e-11. **Aucune force finale ne doit être publiée.** Les essais ne qualifient pas une convergence N256 en temps fini.

### Maillage et GCI

| Raffinement | Surface | Première couche normale |
|---|---:|---:|
| N64 → N128 | 2 | 3.2343977 |
| N128 → N256 | 2 | 3.2656101 |

Le modèle radial à puissance 1.8 ne raffine pas la première couche comme la surface. Même avec trois calculs convergés, un GCI classique à un unique facteur r=2 ne peut être certifié sur cette famille. Le script de qualification retourne `BLOCKED` pour les raisons indépendantes **N256 sans certificat de convergence** et **non-similarité du raffinement de paroi**. Un plan séparant sensibilité tangentielle et normale ou une famille réellement similaire est nécessaire.

## Tests reproductibles

```sh
cmake -S . -B build -DCFDX_BUILD_TESTS=ON
cmake --build build --target test_galerkin_csr_symbolic_numeric test_naca0012_qualification -j 2
ctest --test-dir build -R 'test_galerkin_csr_symbolic_numeric|test_convergence_monitor|test_naca_force_mesh_qualification|test_naca_force_certificate_fail_closed|test_naca_mesh_similarity_gate' --output-on-failure
python tests/python/test_naca_force_mesh_qualification.py
python scripts/naca_mesh_similarity_gate.py /path/to/naca0012_mesh_quality.json
```

`naca_mesh_similarity_gate.py` et `naca_force_mesh_qualification.py` retournent un code non nul sur les fichiers NACA existants : c'est leur **résultat attendu**, pas une régression du test unitaire.

Résultats exécutés : **18/18 CTest Release**, **11/11 scénarios unitaires Galerkin sous ASan/UBSan**, **tests monitor sous ASan/UBSan réussis**, **6/6 tests Python de la porte forces/maillage**. La suite complète CFDX sous sanitizers et les longs cas NACA ne sont pas qualifiés.

## Benchmark apparié de rafraîchissement numérique Galerkin

Les 11 paires de mesure (C++ optimisé, un cœur CPU fixé, 25 mises à jour par passage) donnent :

| Matrice fine | V5, médiane ms/actualisation | V6, médiane ms/actualisation |
|---|---:|---:|
| 16 384 inconnues | 1.18472 | 1.20186 |
| 65 536 inconnues | 4.98134 | 4.88790 |

Les checksums des matrices sont identiques. V6 ne fait plus d'allocation des trois buffers après le premier appel ; la vitesse n'a **pas** de gain significatif démontré au niveau de ce microbenchmark. Ces temps ne représentent pas une simulation NACA complète.

## Méthodes et références

- PETSc `MatProductSymbolic` / `MatProductNumeric` et `MatPtAP` : séparation de la structure et des calculs numériques, utile pour la mise à jour Galerkin répétée. L'implémentation reste du C++ CFDX original, sans dépendance PETSc.
- HYPRE BoomerAMG : stratégies de grossissement HMIS/PMIS, interpolation à distance deux, plafonnement des coefficients et lisseurs anisotropes — à qualifier **séparément** sur les captures de pression avant activation globale.
- Pour les efforts de paroi, les reconstructions de face ont déjà une vérification affine MMS et la reconstruction quadratique un test d'ordre deux ; les résultats physiques NACA doivent compléter ces preuves avant de promouvoir une formule comme calcul de référence.

Documentation primaire : https://petsc.org/main/manualpages/Mat/MatPtAP/ ; https://petsc.org/release/manualpages/Mat/MatProductNumeric/ ; https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html .

## Verrous de qualification

1. Revoir le bord d'attaque et le traitement des faces / normales / pression : 80 % de l'écart de pression N64-N128 est dans les premiers 5 % de corde.
2. Qualifier la reconstruction de la traction visqueuse sur une famille de maillages avec raffinement normal contrôlé et contre un cas analytique/expérimental indépendant.
3. Stabiliser et converger N256 sans abaisser les critères de momentum ni de conservation. Comparer les schémas de convection et relaxations seulement sur des campagnes A/B contrôlées ; enregistrer toute modification de formulation.
4. Exécuter une étude de convergence spatiale cohérente multi-paramètres, puis qualifier Cd/Cp/Cf/Cl/Cm et un GCI justifié.
5. Qualifier les temps jusqu'à convergence complète : la vitesse du produit Galerkin seul n'est pas un gain CFD prouvé.
