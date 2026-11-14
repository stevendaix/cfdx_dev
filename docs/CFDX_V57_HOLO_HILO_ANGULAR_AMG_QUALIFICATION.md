# CFDX V57 — HO angulaire intégralement hi+lo avant la projection; HO–SMM/AMG globalisé

12 novembre 2026 — **branche CPU de recherche**, aucun changement du solveur DOM de production.

## Résumé exécutable et limite de qualification

La V56 stockait deux doubles dans le point fixe HO–LO mais transmettait **seulement `x.hi`** à `DomMomentSchurSystemV44::project()` ; les intensités, les courants de face et la fermeture de réflexion HO étaient calculés en double. V57 ajoute un transport HO directement en deux composantes avant la projection, et propage ensuite ces composantes jusqu'au correcteur SMM/AMG, y compris les puissances aux faces. L'auditeur reconstruit indépendamment les équations angulaires depuis la connectivité du maillage, et non depuis le CSR ou le bilan LO construit.

Dans la famille de maillages 8×8×10 / 64 directions testée, les cinq profondeurs 0,060 / 0,030 / 0,012 / 0,006 / 0,003 m passent les résidus HO, LO et énergétiques `1e-9` après 17–18 états extérieurs acceptés, avec **résidu linéaire hi+lo strict <= 1e-13**. Ces résultats sont propres à la branche expérimentale CPU et à ces cinq géométries. Ni V48/V50 en double, ni le DOM, ni le CHT en production ne sont certifiés. Les tests longs doivent être explicitement sélectionnés.

## 1. Opérateur HO réellement utilisé

Pour chaque direction m et cellule i, la matrice FV UPWIND historique A_m est prélevée **en lecture seule** par `angular_matrix_research_v57()`, sans recréer une diffusion artificielle ni changer son graphe. Les seconds membres sont reconstruits en arithmétique à deux doubles :

```
 b_(m,i) = V_i [ k_s,i * (J_hi+J_lo)_i + k_a,i * I_b,i ]
         + sum_(f incoming at i) |Omega_m.S_f| [ eps_f I_b,f
                        + (1-eps_f) * (Hplus/Hminus)_f * (R_hi+R_lo)_f ]
 A_m I_m = b_m
```

Les entrées du CSR d'origine restent `double` pour refléter la discrétisation effectivement assemblée, mais **les produits, la somme de la ligne et les inconnues `I_m` sont en `DomHiLoV55`**. La solution est calculée sur l'ordre topologique déduit des non-zéros du CSR. Une dépendance circulaire provoque un **refus explicite** : ce balayage est validé sur la famille testée, pas sur tout maillage polyédrique. Un fallback GMRES/IR pour ces cas reste à développer.

Le véritable défaut `b_m-A_m I_m` est recalculé sur les deux composantes, sans se fier à la précision apparente d'un solveur. Le moment J, le moment de réflexion R, le tenseur d'Eddington et le flux HO à chaque face utilisent ensuite les intensités *hi+lo* ; la projection vers `double` n'a lieu **qu'à la frontière du tenseur E du correcteur WLS V54**, ce qui est une limite encore ouverte. Les opacités, géométries et températures d'entrée sont elles aussi stockées en double.

## 2. Réassemblage indépendant et falsification

`audit_dom_angular_hilo_v57()` refait direction par direction : le flux de face selon le signe de `Omega.S`, l'absorption, la diffusion, l'émission, l'énergie locale, les moments et l'état de réflexion. Il n'invoque ni `lo.consistent_source`, ni la matrice `A_m`, ni la correction HO−LO. Le diagnostic compare séparément le courant obtenu aux faces à celui de l'évaluateur ; ce **n'est pas une preuve MPI**.

Contrôles de réfutabilité inclus :

- injection d'une perturbation de 2 % **uniquement** dans `lo` : la solution et les moments angulaires doivent changer ;
- intensité corrompue après résolution : l'auditeur indépendant doit la refuser ;
- courant de face remplacé après résolution : le comparateur face par face doit le refuser ;
- équilibre homogène analytique à 420 K, milieux optiques variables, parois grises réfléchissantes : fermeture sur 27 et 640 cellules ;
- seuil linéaire volontairement assoupli et seuil physique assoupli : erreurs explicites ;
- arrêt après un seul pas extérieur : aucun certificat physique.

La conservation de la face HO−LO est vérifiée, mais elle est une **identité construite**, pas un oracle indépendant du transport. Le défaut du véritable point fixe HO est distinct et doit tendre vers zéro.

## 3. Couplage HO–SMM/AMG globalisé

1. `DomAngularHiloTransportV57` résout toutes les ordonnées en deux composantes, obtient les courants et les moments, puis l'auditeur FV indépendant est exécuté.
2. Les tenseurs `E` actualisés, les conditions grises et les opacités fabriquent un nouvel opérateur SMM/WLS. **Il est reconstruit pour chaque état et chaque essai d'acceptation.**
3. Un **seul** courant HO−LO orienté est calculé par face, injecté avec signes opposés dans les deux cellules. Le second membre corrigé reste `hi+lo`, sans repli implicite.
4. GMRES préconditionné par AMG directionnel donne une direction. Le raffinement V56 la corrige et recalcule `||b-Ax||_2/||b||_2` dans l'arithmétique des inconnues stockées. Un échec du seuil strict `1e-13` refuse cette direction.
5. Une recherche de pas éprouve à nouveau **les équations HO et LO**. Le mérite `max(R_HO,R_LO)` doit contracter ; une contraction individuelle des deux équations ou le succès complet des seuils peut également accepter un candidat. Aucun seuil physique n'est augmenté.
6. Avant la publication expérimentale, `audit_dom_angular_v54()` est également réexécuté sur une projection double pour contrôle croisé du modèle historique. Puis l'opérateur V57 et l'auditeur hi+lo sont réassemblés sur **les deux composantes réellement conservées**. Le certificat ne contient pas de champ double publiable en production.

Les témoins `max_tensor_refresh`, `max_face_refresh`, `angular_sweeps`, `refreshed` et les histoires HO/LO/énergie permettent de détecter un coefficient non actualisé ou un faux certificat.

## 4. Résultats de qualification CPU V57 (Release diagnostique, -O1 sur les deux nouvelles unités)

| Profondeur (m) | Cellules | Acceptés | Itérations Krylov cumulées | HO max final | LO L2 final | énergie max finale | linéaire hi+lo final |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.060 | 640 | 17 | 856 | 1.81915e-10 | 2.85174e-10 | 1.81915e-10 | 9.73054e-15 |
| 0.030 | 640 | 17 | 1005 | 2.38761e-10 | 4.80093e-10 | 2.38761e-10 | 2.00799e-14 |
| 0.012 | 640 | 18 | 1166 | 1.66979e-10 | 5.22588e-10 | 8.97905e-11 | 8.40612e-15 |
| 0.006 | 640 | 18 | 1518 | 1.44992e-10 | 5.75007e-10 | 5.76041e-11 | 1.07653e-14 |
| 0.003 | 640 | 18 | 2724 | 1.02962e-10 | 5.26603e-10 | 3.03821e-11 | 1.13113e-14 |

Les résultats de ce tableau devront être remplacés par ceux du dernier lancement si les sources sont modifiées. Une première campagne de cinq géométries a réussi en 46,94 s hors compilation ; ces mesures ne sont pas des benchmarks CPU/GPU industriels. En 27 cellules à 0,03 m : 17 acceptations, 150 itérations Krylov, HO 2.53514e-10, LO 2.51799e-10.

### Oracle angulaire figé — précision

Sur 640 cellules, sans itérations extérieures, le défaut angulaire géométrique indépendant reste dans `[3.27e-16, 4.06e-16]` sur les cinq profondeurs. La réévaluation CSR du transport hi+lo est de l'ordre de `1e-31`, **mais cela ne suffit pas à valider la convergence extérieure**. L'initialisation du cas extrême présente encore un défaut de point fixe de `8.26e-3`, qui ne disparaît qu'avec la boucle globalisée.

## 5. Comparabilité, risques et absence de qualification industrielle

La famille reste un petit bloc 3D régulier transformé/rotaté (8×8×10), 64 ordonnées, avec absorption/diffusion faibles et des températures de paroi de 408/420/440 K. Ce n'est ni une géométrie industrielle complexe, ni un maillage non conforme, ni une campagne d'ordre spatial, ni une comparaison expérimentale externe.

**Limites qui restent ouvertes :**

- Le tenseur E transmis au correcteur est encore `double` ; le coefficient de face SMM/WLS n'est donc pas intégralement stocké en double-double. Tester ensuite E hi+lo et la sensibilité du préconditionneur AMG au coefficient de face.
- Les matrices FV angular CSR conservent des coefficients assemblés en double. Les cycles d'advection sur des maillages arbitraires ne sont pas traités par le balayage et sont refusés explicitement.
- Le stockage des intensités de toutes les ordonnées exige une empreinte proportionnelle à Ncell × Nangles, sans optimisation bloc/sweep GPU ; viser un ordonnancement compact et des buffers réutilisables.
- Les campagnes MPI/GPU natifs, une vraie CHT–DOM non conforme à émission radiative résolue **face par face**, les calculs transitoires longs et les vérifications physiques indépendantes demeurent à réaliser.
- ASAN/UBSAN ciblé (unités V57 instrumentées) ne représente **pas** une instrumentation de toute la bibliothèque `cfdx_core` ; la CI complète et le solveur de production n'ont pas été qualifiés.

## 6. Reproduction et options

```sh
# Sans longs travaux (ports V57 + oracles négatifs) :
sh scripts/run_radiation_v57_qualification.sh /tmp/cfdx-v57 --quick
# Ajout V53/V54/V55/V56 :
sh scripts/run_radiation_v57_qualification.sh /tmp/cfdx-v57 --extended
# La campagne V57 de 5 profondeurs doit être activée volontairement :
sh scripts/run_radiation_v57_qualification.sh /tmp/cfdx-v57 --long
# Refus explicite code 3 (pas une fausse simulation parallèle):
sh scripts/run_radiation_v57_qualification.sh /tmp/cfdx-v57 --mpi
sh scripts/run_radiation_v57_qualification.sh /tmp/cfdx-v57 --gpu
```

Une compilation à `-O1` des deux unités V57 est un compromis de diagnostic reproductible ; elle **ne modifie pas** l'optimisation habituelle du noyau de production.

## Références

- PETSc KSP, convergence, normes préconditionnées et `KSPMonitorTrueResidual`, https://petsc.org/main/manual/ksp/ ; un GMRES estimé n'est pas une preuve de résidu vrai sur l'état final.
- Carson & Higham, *A New Analysis of Iterative Refinement and Its Application to Accurate Solution of Ill-Conditioned Sparse Linear Systems*, SIAM J. Sci. Comput. 39 (2017), DOI: https://doi.org/10.1137/17M1122918. Les garanties de cet article dépendent de ses hypothèses ; elles ne qualifient pas une implémentation CFDX.

**Décision :** qualification numérique **expérimentale CPU sur les cinq cas testés**, pas promotion en solveur industriel. Aucun seuil assoupli. Aucun push GitHub.
