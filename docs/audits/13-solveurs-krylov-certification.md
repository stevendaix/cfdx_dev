# 13 — Solveurs Krylov, certification du résidu et précision mixte

**Statut :** fiche d'audit des implémentations existantes — correction G02–G06/G09 en cours.  
**Périmètre :** CG préconditionné, GMRES/FGMRES, BiCGStab, réduction MPI, arithmétique mixte et oracle ExactSchur. Cette fiche ne qualifie pas à elle seule un solveur de production.

## 1. Problème algébrique

On résout le système linéaire
\[
A x=b
\]
avec une méthode itérative. L'itéré (x_k) n'est convergé que si le résidu du système demandé
\[
r_k=b-Ax_k
\]
satisfait le contrat de tolérance. Une norme de résidu préconditionnée, une estimation de Hessenberg ou un résidu récursif ne constitue pas à elle seule une certification de (r_k).

## 2. Critère d'arrêt

Le contrat de référence est
\[
\|b-Ax_k\|_2 \le \max(\mathrm{atol},\mathrm{rtol}\,\|b\|_2).
\]
Quand (b=0), un itéré initial non nul n'implique pas une divergence : il faut réduire le résidu ( -Ax_k) jusqu'au seuil absolu, ou jusqu'à un seuil explicitement défini par le contrat de l'API. Les valeurs non finies doivent être rejetées avant les branches « second membre nul ».

Le résidu final doit être recalculé avec l'opérateur original en FP64 lorsque cet opérateur est disponible. Si une API matrix-free ne fournit pas d'application haute précision, le résultat doit indiquer que cette certification est indisponible ; il ne faut pas annoncer une certification FP64 non réalisée.

## 3. Préconditionnement et invariance d'échelle

Pour CG avec préconditionneur (M), le produit (r^T M^{-1}r) sert à la récurrence, mais sa racine n'est pas interchangeable avec \(\|r\|_2\). Les tests de convergence doivent employer le résidu du système original.

Pour toute échelle commune (s>0), le système (sAx=sb) possède la même solution. La classification de convergence et l'erreur de solution doivent donc rester invariantes à l'échelle, dans les limites de la précision flottante et du conditionnement. Les tests de rupture de BiCGStab doivent comparer les produits scalaires à des échelles relatives construites à partir des normes des vecteurs concernés, non à une constante absolue indépendante des unités.

## 4. Arnoldi et réduction distribuée

Les normes, projections d'Arnoldi et passes de réorthogonalisation doivent partager le même domaine owned et la même politique de réduction globale. Les valeurs ghost peuvent être mises à jour dans les vecteurs, mais ne doivent pas être comptées plusieurs fois dans les produits scalaires. Une validation CPU monoprocessus ne suffit pas à certifier MPI : la campagne distribuée doit inclure partitions différentes, rafraîchissement des halos et collectives.

## 5. Précision mixte

L'application de l'opérateur d'itération peut être FP32, mais le calcul du résidu de certification doit employer les coefficients FP64 d'origine. La norme calculée en double d'un produit matrice-vecteur FP32 ne restaure pas les bits perdus. Si le résidu certifié dépasse le seuil, le solveur doit poursuivre, raffiner ou retourner un statut non convergé ; il ne doit pas déclarer CONVERGED sur le seul résidu de travail.

## 6. Contrats de panne et oracle

Les statuts de rupture doivent conserver un diagnostic de résidu réel ou une valeur non finie explicite. Un zéro par défaut ne doit jamais être interprété comme une solution. Les oracles comme ExactSchur doivent rejeter NaN/Inf avant toute reconnaissance d'un second membre nul, puis vérifier le résidu final avant de retourner le succès.

## 7. Cas de régression obligatoires

- Résolution exacte connue avec (A,b) multipliés par plusieurs puissances de dix.
- (b=0) et (x_0\ne0), avec solution atteignable.
- Matrice FP64 dont un coefficient diffère de moins d'un ulp FP32 de sa valeur arrondie ; certification contre l'opérateur FP64.
- Ruptures algébriques réelles et systèmes bien conditionnés mais de petite/grande échelle.
- Plage owned distincte des ghosts, projections Arnoldi et réorthogonalisation.
- Entrées NaN/Inf dans les oracles et diagnostic de résidu renseigné sur chaque sortie.
- Validation Release et DebugSanitizers sur le même HEAD ; test MPI multi-rang lorsque le chemin distribué est activé.

Aucune tolérance ne doit être élargie pour faire passer ces tests. Les campagnes doivent distinguer les tests effectivement exécutés de ceux qui sont indisponibles ou ignorés.

## 8. Références

- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, 2003. DOI: [10.1137/1.9780898718003](https://doi.org/10.1137/1.9780898718003).
- H. A. van der Vorst, *Iterative Krylov Methods for Large Linear Systems*, Cambridge University Press, 2003. ISBN: 978-0-521-81828-5.
- R. Barrett et al., *Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods*, SIAM, 1994. ISBN: 978-0-89871-328-2.

## 9. Statut CFDX / G01–G09

Cette fiche est un contrat mathématique et une liste d'essais. La présence du code ou d'un test unitaire ne vaut pas qualification complète. Les constats G01–G09 restent à confirmer sur le HEAD exact par les tests de régression et CI ; toute portée non exécutée (notamment MPI multi-rang ou opérateur matrix-free sans FP64) doit rester explicitement limitée.
