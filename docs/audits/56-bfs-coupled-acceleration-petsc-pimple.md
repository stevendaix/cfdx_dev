# Accélération du couplage pression–vitesse BFS : PCD, FieldSplit et PIMPLE

## 1. Objet et domaine d'application

Cette fiche décrit les réglages expérimentaux de résolution du problème incompressible stationnaire à marche descendante (BFS) introduits par `bfs_convergence_accelerator.h`.

Le périmètre couvre uniquement la stratégie de couplage pression–vitesse, la sélection du préconditionneur, les correcteurs de pression et la relaxation. Il ne constitue ni une validation du modèle SST, ni une qualification du réattachement Driver–Seegmiller. Les options restent opt-in ; aucun statut physique n'est promu.

## 2. Formulation continue

Le système incompressible stationnaire s'écrit :

\[
\nabla\cdot\mathbf{u}=0,
\qquad
\rho(\mathbf{u}\cdot\nabla)\mathbf{u}
=-\nabla p+\nabla\cdot\boldsymbol{\tau}+\mathbf{f}.
\]

Après linéarisation et discrétisation, le système couplé vitesse–pression prend la forme selle :

\[
\begin{bmatrix}
F & B^T\\
B & 0
\end{bmatrix}
\begin{bmatrix}
\delta u\\ \delta p
\end{bmatrix}
=
\begin{bmatrix}
r_u\\ r_p
\end{bmatrix},
\qquad
S=-BF^{-1}B^T.
\]

Le complément de Schur exact est une référence algébrique. Un préconditionneur PCD ou FieldSplit n'est pas l'opérateur exact et doit être évalué sur le système réellement assemblé.

## 3. Discrétisation Volume Fini

Les équations de quantité de mouvement et de continuité sont assemblées sur le maillage de volumes finis existant. La correction de pression doit rester cohérente avec le flux de masse Rhie–Chow utilisé pour le bilan de continuité.

Les stratégies PCD et NativeFieldSplit agissent sur le système linéaire couplé. FGMRES est retenu car les préconditionneurs et leurs applications peuvent être non symétriques ou variables. Le mode sans fallback refuse explicitement une substitution silencieuse du préconditionneur demandé.

Les correcteurs PIMPLE répètent les corrections pression–vitesse dans une itération externe ; ils ne changent ni le modèle continu ni la discrétisation spatiale.

## 4. Traitement des frontières

Les opérateurs de frontière restent ceux du solveur de production. Une stratégie de préconditionnement ne crée pas de nouvel opérateur de frontière et ne doit pas convertir implicitement une entrée débit/pression en entrée vitesse.

Les patches explicitement classifiés par le maillage sont conservés. Seuls les patches de type `UNKNOWN` peuvent être classifiés depuis les conditions physiques du cas. Les contradictions entre classification explicite et condition physique restent rejetées.

Une entrée `ZERO_GRADIENT` n'est admise par l'adaptateur courant que si aucune vitesse non nulle n'est prescrite. Les entrées pression/débit restent refusées tant qu'un opérateur dédié n'est pas implémenté.

## 5. Pièges connus

- Un préconditionneur plus robuste peut réduire le nombre d'itérations linéaires sans améliorer le résidu non linéaire physique.
- Un résidu Krylov faible ne prouve ni la convergence du système couplé réassemblé, ni la conservation des flux sur l'état final.
- Une relaxation adaptative peut améliorer la stabilité, mais ne doit jamais modifier les tolérances de convergence physiques fournies par l'appelant.
- Le schéma UPWIND est robuste mais diffusif ; il ne constitue pas une discrétisation d'ordre deux et ne peut pas servir de preuve de précision spatiale.
- Un nombre accru de correcteurs PIMPLE augmente le travail par itération externe ; le coût total doit être mesuré à critères physiques identiques.
- Aucun résultat de cette fiche ne démontre à lui seul l'accord expérimental du BFS turbulent ou une valeur qualifiée de \(x_r/h\).

## 6. Vérification et éléments de preuve

Les contrôles logiciels doivent vérifier séparément :
1. le plan linéaire effectivement résolu (Krylov, préconditionneur, stratégie de Schur) ;
2. l'absence de fallback non demandé ;
3. la conservation des tolérances physiques et des conditions aux limites ;
4. la convergence non linéaire réassemblée sur l'état final ;
5. les bilans de masse et les observables physiques indépendantes.

Les tests de configuration et de dispatch sont des tests de contrat logiciel. Une campagne BFS/SST ne peut être déclarée qualifiée qu'après convergence des équations, conservation, indépendance au maillage et comparaison aux références expérimentales.

## 7. Statut d'implémentation

**Statut : EXPÉRIMENTAL, opt-in, non qualifié physiquement.**

PCD, NativeFieldSplit/Full Schur, réglages PIMPLE et relaxation adaptative sont des choix numériques de résolution. Ils ne constituent pas une preuve de convergence ou de précision. Les critères physiques restent ceux du cas appelant ; les modes non disponibles doivent échouer explicitement plutôt que basculer silencieusement vers une autre méthode.

## Références

- M. Benzi, G. H. Golub, J. Liesen, “Numerical solution of saddle point problems”, *Acta Numerica*, 14, 2005, pp. 1–137. DOI: 10.1017/S0962492904000212.
- H. C. Elman, D. J. Silvester, A. J. Wathen, *Finite Elements and Fast Iterative Solvers*, 2nd ed., Oxford University Press, 2014. ISBN: 9780199678792.
