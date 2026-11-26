# Audit mathématique complémentaire — SST et qualification SPD

Date : 2026-11-26

Branche : `feat/bfs-convergence-petsc-openfoam` (PR #1194)

## Périmètre

Cet audit complémentaire cible trois opérateurs introduits ou activés dans le chemin BFS/SST de la PR #1194 :

1. le pseudo-temps dans le transport SST stationnaire ;
2. la correction matricielle couplée (k)-(omega) ;
3. le contrat SPD utilisé par la sélection CG/AMG.

Il s'agit d'un audit statique des équations et des assemblages. Les campagnes physiques BFS/SST et l'exécution locale des tests restent nécessaires.

## 1. Pseudo-temps SST : séparation candidat / résidu stationnaire

L'équation stationnaire physique s'écrit schématiquement :

\[
\mathcal{R}_{\mathrm{SST}}(\phi)=0,
\qquad \phi\in\{k,\omega\}.
\]

Le pseudo-temps ajoute un terme artificiel :

\[
\frac{\phi-\phi^{old}}{\Delta\tau}
\]

et le système candidat devient :

\[
\left(A+M_{\tau}\right)\phi^{new}
=b+M_{\tau}\phi^{old}.
\]

Ce terme peut stabiliser la génération d'un candidat, mais il ne fait pas partie de l'équation stationnaire. Le résidu de certification doit donc être calculé sur l'opérateur non décalé.

**Défaut trouvé dans le chemin `sst_pseudo_transient_dt` :** les sources étaient modifiées par le pseudo-temps avant la construction des matrices `eqk_phys` et `eqw_phys`. Les résidus présentés comme physiques incluaient donc le terme artificiel.

**Correction sur la branche :** les quatre champs source/destruction physiques sont maintenant copiés avant l'ajout du pseudo-temps. Les systèmes candidats utilisent les termes décalés ; les systèmes de certification utilisent les sources non décalées. Le chemin séquentiel omega-first réassemble également l'équation physique de k après mise à jour de omega.

Test de régression ajouté : `sst_pseudo_time_does_not_enter_steady_residual_certificate`.

## 2. Correction couplée (k)-(omega) : dimension, RHS et point fixe

Les équations scalaires assemblées sont :

\[
A_k k=b_k,\qquad A_\omega\omega=b_\omega.
\]

Les vecteurs (b_k) et (b_\omega) doivent être les RHS complets des objets `ScalarEquation`, car ils contiennent notamment les contributions de frontière et de relaxation.

Pour ajouter une approximation du couplage des sources explicites autour de l'état courant ((k_0,\omega_0)), posons :

\[
C_{k\omega}=-\frac{\partial S_k}{\partial\omega},
\qquad
C_{\omega k}=-\frac{\partial S_\omega}{\partial k}.
\]

Le système tangent affine doit être :

\[
A_k k+C_{k\omega}\omega
=b_k+C_{k\omega}\omega_0,
\]
\[
C_{\omega k}k+A_\omega\omega
=b_\omega+C_{\omega k}k_0.
\]

Le décalage du RHS garantit que le résidu du système couplé, évalué à l'état de linéarisation, est exactement le résidu des deux équations scalaires. Sans ce décalage, l'ajout des blocs croisés déplace artificiellement le point fixe.

**Trois défauts trouvés dans l'ancienne implémentation :**

- les différences finies étaient calculées par rapport à (log\omega) et (log k), puis utilisées comme coefficients multipliant les champs dimensionnels ; il manquait respectivement la division par $\\omega$ et par $k$, ainsi que la multiplication par le volume de cellule pour obtenir un coefficient de matrice intégré ;
- les côtés droits étaient reconstruits à partir de `sk` et `sw`, au lieu des RHS assemblés `eqk.rhs` et `eqw.rhs`, omettant ainsi les contributions de frontière et de relaxation ;
- les termes croisés étaient ajoutés sans correction affine du RHS, ce qui changeait le point fixe.

**Correction sur la branche :** les dérivées sont converties en dérivées par rapport aux variables dimensionnelles, avec le signe correspondant au transfert de la dérivée de source à gauche. Le solveur construit son RHS tangent depuis les RHS assemblés et ajoute le décalage affine évalué à l'état courant. Le plafonnement des corrections est désormais relatif à chaque variable, et non un même incrément absolu pour (k) et (omega).

Test de régression ajouté : `sst_coupled_2x2_tangent_preserves_assembled_scalar_fixed_point`.

Cette correction rend le tangent cohérent au point de linéarisation. Elle ne transforme pas pour autant cette approximation partielle des dérivées de sources en Jacobien exact du système SST complet.

## 3. Qualification SPD : symétrie et diagonale positive insuffisantes

Une matrice symétrique à diagonale positive n'est pas nécessairement définie positive. Contre-exemple :

\[
A=\begin{bmatrix}1&2\\2&1\end{bmatrix},
\qquad \lambda(A)=\{3,-1\}.
\]

Elle est symétrique et ses deux coefficients diagonaux sont positifs, mais elle est indéfinie. Ces deux propriétés ne suffisent donc pas à autoriser CG ni un cycle AMG dont le contrat est SPD.

Pour les problèmes génériques, la branche exige maintenant en plus un certificat suffisant fondé sur une borne de Gershgorin strictement positive :

\[
\min_i\left(a_{ii}-\sum_{j\ne i}|a_{ij}|\right)>\varepsilon.
\]

Pour une matrice symétrique, cette condition garantit que toutes les valeurs propres sont strictement positives. C'est un critère suffisant mais non nécessaire : certains systèmes SPD non strictement dominants seront volontairement refusés par ce chemin conservateur.

Le contrat historique des opérateurs canoniques PressurePoisson/Diffusion n'est pas remplacé par cette preuve générique ; les cas à noyau constant nécessitent toujours un traitement explicite de l'espace nul.

Tests de régression ajoutés :

- `generic_cg_and_amg_reject_symmetric_positive_diagonal_indefinite_matrix` ;
- `strict_gershgorin_margin_certifies_generic_symmetric_positive_definite_matrix`.

## 4. État des preuves

| Élément | Statut |
|---|---|
| Séparation résidu stationnaire / pseudo-temps | Corrigée dans le code ; test ajouté |
| Tangente couplée (k)-(omega) dimensionnelle | Corrigée dans le code ; test ajouté |
| Préservation du point fixe de la tangente | RHS affine corrigé ; test ajouté |
| RHS complet incluant frontières/relaxation | Le solveur utilise maintenant `ScalarEquation::rhs` |
| Qualification SPD générique | Critère Gershgorin strict ajouté ; tests ajoutés |
| Exécution C++ / CI de ces nouveaux tests | À confirmer sur le commit final |
| Qualification physique BFS/SST | Non acquise ; campagne convergée et étude de maillage requises |

## 5. Suite de l'audit

Avant toute promotion des options SST expérimentales :

1. exécuter les tests unitaires ciblés et la suite Release ;
2. comparer, sur un même état, le résidu scalaire et le résidu du tangent couplé au point de linéarisation ;
3. vérifier par différences finies centrées indépendantes les dérivées des sources (S_k) et (S_\omega), y compris dans les branches de limiteur ;
4. mesurer le résidu SST entièrement réassemblé après chaque candidat, sans pseudo-temps ;
5. vérifier les bilans locaux, la positivité et l'indépendance au maillage avant de comparer (x_r/h) à Driver–Seegmiller.

Aucune tolérance physique n'a été relâchée et aucun résultat de convergence physique n'est revendiqué par ces corrections.
