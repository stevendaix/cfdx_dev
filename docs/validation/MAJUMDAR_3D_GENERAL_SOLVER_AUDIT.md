# Audit Majumdar / Rhie–Chow — solveur 3D général

**Statut :** audit statique du chemin de production effectué ; qualification numérique 3D non acquise.  
**Périmètre :** `steady_incompressible_solver.h`, interpolation Rhie–Chow générique, tests N9/S9.  
**Conclusion courte :** la pondération par la diagonale de quantité de mouvement relaxée est déjà utilisée dans le chemin ségrégé général 3D. Il ne faut donc pas recopier aveuglément le correctif spécifique de `axisymmetric_vmfl036_solver.h` dans ce solveur. En revanche, la preuve actuelle ne qualifie pas toute la chaîne 3D sur maillage non orthogonal.

## 1. Faits vérifiés dans le code

Dans `solve_steady_incompressible` :

1. Le chemin ségrégé appelle `relax_momentum_equation(ex/ey/ez, ..., alpha_u)` avant de résoudre les équations de quantité de mouvement.
2. Après cette relaxation, le code construit les trois tableaux `rAU[d][c] = 1.0 / eqs[d]->diagonal[c]`.
3. Ces tableaux sont transmis à `make_rhie_chow_mass_flux(...)`.
4. La routine générique projette les réponses directionnelles sur la normale de face :
   `rfn = cx e_x² + cy e_y² + cz e_z²`, avec `cx/cy/cz = 0.5 (V rAU)_P + 0.5 (V rAU)_N`.
5. Le terme de pression combine le saut orthogonal `(p_N-p_P)/d` et la contribution différée non orthogonale `grad(p)·S_non`.

Comme la diagonale est déjà relaxée au moment où `rAU` est calculé, on a, pour une relaxation diagonale standard,
`a_P^r = a_P/alpha_u` et `rAU = 1/a_P^r = alpha_u/a_P`.
Le coefficient de Rhie–Chow du chemin ségrégé reçoit donc bien la réponse inverse de la diagonale relaxée. Il ne faut **pas** ajouter un facteur `alpha_u` supplémentaire dans `make_rhie_chow_mass_flux` : cela appliquerait la relaxation deux fois.

Sources :
- `src/cfdx/physics/steady_incompressible_solver.h` — relaxation des équations, construction de `rAU`, appel à `make_rhie_chow_mass_flux`.
- `tests/validation/test_n9_s9_rhie_chow_consistency.cpp` — test de mise à l'échelle de la correction lorsque `rAU` est multiplié par `alpha_u`.
- `src/cfdx/physics/axisymmetric_vmfl036_solver.h` — implémentation distincte du solveur axisymétrique ; ne pas supposer que ses opérateurs géométriques sont interchangeables avec le chemin 3D.

## 2. Limite importante : le chemin COUPLED

Le chemin `COUPLED` ne passe pas par `relax_momentum_equation` dans la même séquence. Il construit `coupled_rAU` à partir des diagonales du bloc de quantité de mouvement, puis appelle le même constructeur de flux. La relaxation du solveur couplé est gérée différemment (mise à jour explicite des champs / contrôle de couplage).

Il serait donc incorrect de conclure que la même relation de Majumdar, interprétée pour une équation ségrégée sous-relaxée, est déjà démontrée pour `COUPLED`. Ce chemin doit être analysé selon la formulation du système bloc ; il ne faut pas lui appliquer une correction `alpha_u` sans dérivation et test indépendant.

## 3. Ce que les tests actuels prouvent — et ne prouvent pas

Le test `n9_s9_rhie_chow_correction_scales_linearly_with_rAU_majumdar` vérifie la linéarité sur une face orthogonale. Un nouveau test, `n9_s9_rhie_chow_oblique_3d_face_operator_and_relaxation_scaling`, exerce directement l'opérateur de flux de production sur un hexaèdre 3D cisaillé : il vérifie la contribution non orthogonale avec un champ de pression linéaire exact, des réponses directionnelles anisotropes et la mise à l'échelle par `alpha_u`. Cela améliore la preuve opérateur, mais ne vérifie pas à lui seul :

- le passage effectif de `alpha_u` à la diagonale relaxée dans une itération complète du solveur 3D ;
- la cohérence entre le flux de Rhie–Chow, la matrice de correction de pression et la correction de vitesse sur un maillage 3D oblique/non orthogonal ;
- l'invariance de la solution convergée vis-à-vis de `alpha_u` à tolérances numériques comparables ;
- le comportement du chemin `COUPLED`, distinct du chemin ségrégé.

Les tests d'invariance de décalage de pression et de mode de pression alterné sont utiles, mais ne remplacent pas ces vérifications d'intégration.

## 4. Actions de qualification requises

- [x] Ajouter un test opérateur 3D sur maillage non orthogonal, avec réponses directionnelles anisotropes et vérification analytique du flux de pression.
- [ ] Ajouter une vérification d'intégration démontrant que le `rAU` passé au flux ségrégé provient bien de la diagonale relaxée pour plusieurs valeurs de `alpha_u`.
- [ ] Comparer les bilans de continuité du flux conservatif et du flux reconstruit à partir de la vitesse corrigée sur un cas 3D oblique.
- [ ] Faire un sweep `alpha_u` sur un cas 3D validé, en exigeant la convergence de toutes les métriques physiques avant de comparer les quantités d'intérêt.
- [ ] Auditer séparément `COUPLED` : dérivation documentée, puis test de cohérence du bloc et de sa réponse pression–vitesse. Ne pas transposer automatiquement la formule ségrégée.

## 5. Impact vérifié et limites

Le test ajouté ne change aucune formule de production. Il est conçu pour échouer si le terme différé non orthogonal est omis ou si la réponse directionnelle est projetée de manière incohérente ; il vérifie aussi que multiplier les trois `rAU` par `alpha_u` multiplie le flux de pression par le même facteur. La validation effective dépend du résultat CI de la PR. Aucun résultat de compilation ou d'exécution ne doit être annoncé tant que ce job n'a pas terminé.

## 6. Décision

**Pas de modification de formule dans le solveur 3D ségrégé à ce stade** : l'inspection du code indique que la diagonale relaxée est déjà transmise à Rhie–Chow. Une modification supplémentaire risquerait de compter deux fois `alpha_u`.

La décision reste toutefois **« implémentation statiquement présente, qualification 3D incomplète »**, et non « Majumdar entièrement validé en 3D ». La prochaine modification de code doit être guidée par les tests d'intégration ci-dessus, notamment sur géométrie non orthogonale, plutôt que par une copie du patch axisymétrique.

## 7. Note d'intégration (Lot A, branche `integration/numerics-2026-10`)

Cette fiche est intégrée dans le candidat d'intégration #929 avec trois adaptations explicites :

1. **Oracle du test oblique adapté au contrat #923.** L'opérateur de production de la branche d'intégration (`rhie_chow_pressure_flux_internal`) retourne le flux tensoriel complet `rho * Σ_d dAU_d · grad_d(p) · Sf_d` (coefficient orthogonal implicite `pressure_face_orthogonal_coefficient` + flux différé `pressure_face_nonorthogonal_flux`). L'oracle initial de #924 (`rho * rfn * grad(p)·Sf` avec `rfn = cx·e_x² + cy·e_y² + cz·e_z²`) validait l'opérateur pré-#923, qui projetait la réponse directionnelle sur la normale de face. L'oracle attendu est donc `rho * (cx·g_x·Sf_x + cy·g_y·Sf_y + cz·g_z·Sf_z)`. L'intention du test est préservée : l'oracle tensoriel ne peut pas être satisfait sans le terme différé non orthogonal sur la face cisaillée, et l'assertion de mise à l'échelle `alpha_u` (linéarité en `rAU`) est inchangée.
2. **Séquences `\n` littérales de la section 5 réparées** (même famille de corruption que celle corrigée par le commit `ef2ef62` du candidat) : les faux `\n` sont convertis en vrais sauts de ligne, sans autre modification du contenu.
3. **Retrait du label `long` de `test_n9_s9_rhie_chow_consistency`** (labels restants : `n9;qualification;validation;validation-total`). Le test reste exécuté par la CI standard et par la partition `validation-total` ; seul le doublon de partition longue est retiré (optimisation de budget CI). Aucun contenu de test, oracle ou tolérance n'est modifié.
