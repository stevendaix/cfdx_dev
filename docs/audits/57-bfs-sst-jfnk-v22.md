# CFDX V22 — opérateur SST vectoriel et Newton–Krylov globalisé sur le BFS

**Date :** 10 novembre 2026. **Base :** copie V21 exacte, reconstruite par patch depuis l'archive originale CFDX. **Statut :** implémentation expérimentale, sans qualification physique et sans push. Le chemin historique est conservé. La tolérance non linéaire `1e-8` et les seuils de conservation restent inchangés.

## 1. Contexte scientifique et problème confirmé

Sur le BFS réduit à 712 cellules, la V21 synchronise les 69 contraintes propriétaires `omegaWallFunction` avec l'état courant de `k`, `omega` et de la production de paroi `G`. Le défaut des lignes propriétaires atteint zéro mais la deuxième itération extérieure stagne encore, avec `k_returned_nonlinear_raw=0.008735979683`, `omega_returned_nonlinear_raw=24.41970936` et `omega_returned_row_scaled=1451.800745`. Le défaut dominant se déplace de la contrainte sur la cellule 357 vers les équations intérieures (`k` en 357, `omega` en 363). La convergence d'une équation gelée n'est donc pas suffisante.

**Point ouvert distinct :** la fonction de paroi V17–V21 n'est pas numériquement équivalente à OpenFOAM Foundation v11. Dans `omegaWallFunctionFvPatchScalarField.C`, `updateCoeffs()` met à jour `G` et `omega` avec des poids; `manipulateMatrix()` appelle `matrix.setValues(...)` avec des poids, tandis que la formulation propriétaire CFDX V18–V21 remplace entièrement la ligne. Une égalité de la cible locale ne prouve pas une égalité de l'opérateur FV. Sources : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html.

## 2. Modifications du solveur *réel* de production

### 2.1 Un seul assemblage SST par état

`sst_reassembled_equations_v22(mesh,...,k,omega,...)` extrait l'assemblage FV existant et retourne exactement les équations de transport de `k` et de `omega` (mêmes flux de convection upwind du transport SST, diffusion `gamma_k/gamma_omega`, limiteurs de production, destruction implicite, diffusion croisée, `F1/F2`, conditions limites, production de paroi et ligne propriétaire `omega`). L'évaluateur scalaire antérieur appelle cette même fonction : pas de nouveau schéma de physique introduit.

Le nouvel évaluateur vectoriel `sst_coupled_residual_v22` retourne chaque ligne signée de l'opérateur *réassemblé* :

`Rk_i(k,omega) = sum_j A_k,ij(k,omega) k_j - b_k,i(k,omega)`

`Rw_i(k,omega) = sum_j A_w,ij(k,omega) omega_j - b_w,i(k,omega)`

Les normes de ligne, échelles locales et diagonales sont disponibles pour les contrôles et le préconditionneur; l'accumulation des défauts signés est faite en `long double` avant conversion. L'ancienne certification du résidu scalaire reste fondée sur ses propres critères stricts (elle n'est pas remplacée par une norme facile à satisfaire).

### 2.2 Newton–Krylov sur le vecteur complet, avec positivité

L'option opt-in `sst_jfnk_v22_experimental` utilise les inconnues transformées `z_k=log(k/k_ref)` et `z_omega=log(omega/omega_ref)` ; un état inadmissible ou au plancher est refusé. Pour un produit Jacobien-vecteur, le code appelle **le résidu FV complet à l'état perturbé**, de sorte que les dérivées directionnelles incluent les dépendances en `k`, `omega`, `F1/F2`, les gradients/cross diffusion et les cibles de paroi recalculées, sous `U` et flux de masse figés au niveau de cette itération SST.

Le solveur réutilise le GMRES CFDX existant de `coupled_newton_krylov.h` ; PETSc et HYPRE ne sont ni installés ni simulés. Un préconditionneur Jacobi **figé** est dérivé des diagonales FV transformées en variables log ; il conditionne l'itération GMRES, mais ne modifie ni `F` ni le jugement final de convergence. Le critère linéaire interne `0.25` est une *forcing tolerance* Newton approximative et ne remplace **jamais** le seuil SST physique `1e-8` (voir PETSc `SNES`: https://petsc.org/release/manual/snes/).

Le mécanisme d'acceptation compare la proposition Picard historique et la proposition Newton (si Krylov a passé son contrôle) avec huit longueurs de pas. Un essai est accepté seulement si les champs restent admissibles et si le **mérite SST physique réassemblé** décroît ; le dernier état admissible est restauré en cas de rejet. Chaque diagnostic `SST_JFNK_V22` fournit `selected=picard|jfnk|none` : `proposed=1` seul **ne garantit pas** qu'un pas Newton a été choisi.

La nouvelle option ne peut pas être combinée à `sst_anderson_m2_v20_experimental`, `sst_nonlinear_line_search`, ou `sst_picard_damping != 1`. Elle est **désactivée par défaut** et requiert le chemin explicite expérimental de la simulation SST.

## 3. Rejeu sur le BFS 712 cellules : ce qui est réellement mesuré

Comparaison A/B sur une copie de `matched_b.cfdx.h5`, 712 cellules, même équation finale, même tolérance `1e-8`, même budget de correction pression/momentum. **Seul le budget d'observation SST a été limité à 10 itérations internes**, pour mesurer les premiers pas sans extrapoler à 800. Les deux modes sont non convergés et correctement rejetés par l'itération extérieure. Le journal V22 définit le mérite comme le maximum des défauts locaux du résidu reassemblé : il est distinct de la borne absolue brute de chaque équation.

| Mesure à l'itération SST 10 | Picard de référence (V22 option off) | V22 option on, Jacobi figé |
|---|---:|---:|
| Défaut brut `k` réassemblé | 1.877483714e-3 | 3.960597672e-6 |
| Défaut brut `omega` réassemblé | 1.701857945 | 1.769146644 |
| Défaut `k` normalisé par ligne | 0.794290718 | 0.004930225 |
| Défaut `omega` normalisé par ligne | 1546.506591 | 84.677114 |
| Défaut `k` relatif local | 0.869484293 | 0.889907451 |
| Défaut `omega` relatif local | 0.028364589 | 0.200784731 |
| Tolérance physique | 1e-8 | 1e-8 |
| État extérieur convergé | **NON** | **NON** |

**Lecture nécessaire :** certaines grandeurs (`k` brut, défauts normalisés par ligne) s'améliorent fortement, mais les défauts relatifs locaux et le `omega` brut ne montrent **pas** une meilleure solution générale à l'itération 10. Ces métriques ont des normalisations distinctes; aucune comparaison du seul chiffre le plus favorable ne peut qualifier le solveur.

Le journal de production V22 avec `selected` confirme :

- SST 1 : `proposed=0`, donc **pas Picard sélectionné**, `alpha=0.125`, mérite `28.902052503 -> 5.8448837064` ; l'essai Newton a été refusé car le vrai résidu directionnel GMRES vaut `0.2653795 > 0.25`.
- SST 10 : **pas Newton sélectionné**, `alpha=1`, mérite `0.92596107454 -> 0.88990745136` ; résidu directionnel GMRES `0.24793702144 <= 0.25`.
- Les essais complets avec rollback/retry ont été interrompus pour borner le temps. Le résultat est un **diagnostic incomplet de début d'itération**, sans convergence à 800 itérations, sans deuxième état extérieur accepté et sans rapport physique certifiable.

## 4. Contrôles exécutés

- Compilation et édition des liens du binaire C++ de production : **réussies**, configuration *diagnostique* (`-O0` pour `main.cpp`), pas de benchmark Release optimisé.
- CTest ciblés : **5/5** (SST V19 équilibré, Anderson V20, synchronisation V21, nouveau JFNK V22, blending historique).
- Pytest ciblés : **16/16** (diagnostics BFS V17/V19/V20/V21/V22); test de distinction explicite Picard/JFNK.
- Nouveau test FV C++ : vecteur `[Rk,Romega]` contre le résidu scalaire du même assemblage sur une cellule hexaédrique avec parois.
- Tests JFNK C++ : résolution d'un système couplé synthétique, préconditionneur diagonal (racine inchangée), refus des planchers actifs, rejet d'un résidu incomplet, absence de mutations.
- Lint `ruff` : **non exécuté**, module absent de l'environnement.
- CI générale, grand BFS à 20 540 cellules, convergence 800 itérations, six variantes Schur/PCD/AMG/MGR, NASA : **non exécutées / non acquises**.

Le diagnostic JSON `scripts/diagnose_bfs_sst_v22.py` refuse explicitement toute qualification NASA. Le journal `selected` n'était pas présent dans le premier diagnostic; un rejeu a confirmé la famille effectivement choisie.

## 5. Suite technique prioritaire

1. **Préconditionneur SST structuré** : exploiter des blocs `k-omega` 2x2 par cellule ou Schur sur les équations de transport, avec dérivées des sources/conditions de paroi, et un critère de convergence GMRES du résidu directionnel vrai. Ne pas confondre `0.25` (linéaire) avec `1e-8` (physique).
2. **Performance** : éviter de réassembler intégralement FV à chaque `Jv` quand l'opérateur analytique/frozen sparse plus correction de sources peut reproduire les différences directionnelles. Mesurer nombre d'assemblages, allocations, temps CPU et GPU; ne pas promettre une version GPU tant qu'elle n'est pas implémentée.
3. **BFS 712 long** : refaire la deuxième itération externe (V21 acceptait le premier état), examiner cellules 357, 363, 426, imposer trois états externes **réellement acceptés**, sans changer les seuils.
4. **BFS 20 540 et six COUPLED** : tester les six variantes Schur indépendamment avec mêmes données, mêmes tolérances, même critères; aucune extrapolation d'ILUT petite taille.
5. **Physique OpenFOAM/NASA** : implémenter et tester séparément l'interpolation pondérée des équations et coins `updateCoeffs` / `manipulateMatrix`; comparer `x_r/h`, `Cf`, `Cp` reconstruit en pression statique, profils et tensions de Reynolds seulement après convergence itérative + spatiale.

**Conclusion :** V22 crée le contrat vectoriel qui manquait pour attaquer la stagnation intérieure SST. Une direction Newton est effectivement admissible sur un pas réel du BFS, mais V22 n'obtient pas la convergence globale. Les critères physiques restent inchangés ; aucun push GitHub n'a été effectué.
