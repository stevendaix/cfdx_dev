@@ -0,0 +1,137 @@
# CFDX V65 — BFS152, globalisation V47 et comparaison stricte de residus physiques

**Date :** 13 novembre 2026. **Base :** CFDX V64, reconstruite depuis l'archive BFS V6 d'origine.
**Statut :** diagnostic C++ experimental, non active par defaut. **Qualification physique : NON.**

## 1. Question et methode

Le calcul V64 a montre une contraction du momentum de 1,898668 (iteration 1) a
0,00225539 (iteration 40) avec COUPLED/BlockLocal historique, malgre des
oscillations tardives. Le V47 historique s'arretait a l'iteration 8, car il
comparait un etat precedent et un candidat avec des coefficients/rAU/flux de
reference differents. Le but V65 n'etait pas de creer un autre Newton ni un
Schur, mais de comparer explicitement le MEME evaluateur etat-seul V51 a
l'ancien et au nouvel etat pour une tentative.

Le mode opt-in exige simultanement :
`CFDX_FPHYS_CANDIDATE_V51=1`, `CFDX_FPHYS_V65_SAME_OPERATOR=1`,
ainsi que les controles de globalisation experimentale V13, V14 et V47 dans
le HDF5 (sans V48). Il refuse les maillages de plus de 152 cellules.

A chaque essai d'une taille de pas beta :

1. Sauvegarder l'etat accepte X_old (U,p,k,omega) et l'etat candidat X_trial.
2. Recalculer integralement F_phys(X_old), en demandant deux evaluations
   bit-a-bit reproductibles, puis F_phys(X_trial), avec la MEME fonction V51.
3. Calculer les six normes brutes par bloc. Utiliser les echelles du seul
   etat X_old pour les deux vecteurs, sans permettre au candidat d'ameliorer
   artificiellement son merite par une nouvelle normalisation.
4. Reconstruire les flux de face Rhie-Chow a partir du rAU propre a chaque
   evaluation, et verifier les identites de flux et de divergence (1e-10
   pour le controle numerique independant ; garde local, pas tolerance finale).
5. Conserver les gardes V13/V14, sans relachement : croissance de bloc <=1,20
   et merite candidat / merite ancien <=0,995 ; la premiere iteration est la
   reference de comparaison historique. Un essai refuse restaure l'etat,
   puis reduit le pas par le V47 original. Aucun etat Newton n'est active.
6. N'enregistrer `GLOBAL_SAME_OPERATOR_COMMITTED_V65` qu'apres la transaction
   complete U,p,flux,SST et son certificat de production.

Le flux V51 **reconstruit de l'etat seul** n'est pas identique, en general,
au flux historique de la transaction COUPLED (`old_authoritative_phi_gap`,
`trial_authoritative_phi_gap`). Ces ecarts sont publies et non dissimules.
L'opt-in ne remplace donc pas la correction globale du porteur de flux de
production ; il teste la comparaison du merite, pas l'equivalence des
operateurs conserves entre deux transactions.

## 2. Trois campagnes de meme geometrie et physique

Le pilote a 152 cellules et 634 faces est un cas `TH03_debug_only` avec
traitement parietal SST historique. **Ce n'est pas encore le BFS
Driver-Seegmiller physiquement valide avec fonctions de paroi actives.**
Les conditions limites, SST, discretisation, geometrie et HDF5 sont
identiques sauf les drapeaux de globalisation. Le binaire V65, avec ses
options desactivees, reproduit textuellement les quarante certificats de
production V64.

| Methode | Etats acceptes | Dernier momentum apres SST | Continuite | Fin |
|---|---:|---:|---:|---|
| COUPLED/BlockLocal historique | 40 | 0,00225539 | 4,12874e-5 | Non converge (budget 40) |
| V47 historique V64 | 7 | 0,04063251 | 1,349925e-3 | Echec a 8 |
| V47 + merite V65 | 22 | 0,17435534 | 4,891534e-3 | Echec a 23 |

**Comparaison isochrone iteration 22 :** momentum de base 0,00352414 et
continuite 1,52253e-4, contre 0,17435534 et 4,89153e-3 pour V65. Ne pas
confondre les 22 transactions V65 avec un gain de vitesse de convergence.

### 2.1 Cause precise du rejet V65 a l'iteration 23

Sur les tentatives evaluees, le pas complet beta=1 donne un merite relatif
0,640693 (ameliore) mais une croissance maximale du bloc `omega` de
126,95 (superieure a 1,20). Le pas beta=1/1024 satisfait le garde par bloc,
mais son ratio de merite est 0,998296 (superieur au plafond 0,995).
Le solveur refuse donc les directions avec les criteres existants : les
petits pas ne contractent pas assez, les grands pas destabilisent SST.
Une tentative intermediaire ne produit pas de trace V65 (rejet au stade
anterieur) ; l'auditeur enregistre exactement 10 sondes V65 rejetees a
l'iteration 23, et non 11 faux certificats complets.

Les ecarts de reconstruction Rhie-Chow sur les sondes tracees sont typiquement
~1e-13, sans prouver que le flux retenu en production est identique au flux
etat-seul. L'ecart relatif avec le flux porteur de production atteint
4,88e-3 parmi les sondes V65. Cela appelle une investigation du contrat des
flux et des coefficients dans la transaction et dans le merite.

## 3. Qualification et protections

- 11/11 CTest selectionnes : V65 merit C++, V65 auditeur Python (10 tests,
  dont faux certificats), V64 read-only, V47, V49, V50, V51 et SST V41.
- Compilation du solveur `cfdx_production_solver` reussie en configuration
  diagnostique `CFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON`, pas de benchmark
  Release optimise.
- Parite binaire V65 options OFF : les 40 lignes de production V13, V8,
  post-SST et `INCOMPRESSIBLE_ITER` sont identiques au temoin V64.
- L'auditeur sans `--allow-partial` echoue volontairement sur la campagne
  V65 interrompue (0/1 qualification) ; avec cette option, il publie un JSON
  diagnostic dont `complete_40_iteration_run=false` et
  `physical_qualification=false`.
- Le pilote complet BFS712/BFS20540, les tests de parois SST actives, la CI
  complete, NACA et les campagnes longues rayonnement n'ont pas ete executes.
- Aucun seuil de convergence modifie ; aucune option experimentale activee
  par defaut ; aucun push GitHub.

## 4. Conclusion et chantier recommande apres V65

La correction de comparabilite est utile comme **instrument de diagnostic**, mais
la politique de merite V13/V14 couplee au V47 actuel est trop restrictive et
beaucoup plus lente sur ce pilote. **Conserver COUPLED historique par defaut.**

Le prochain test causal doit confronter, sur le meme etat de l'iteration 22 :
(1) les coefficients rAU et viscosite avant/apres SST ; (2) les flux porteurs
historiques aux flux etat-seul du candidat ; (3) les contributions omega au
rejet du pas complet, cellule par cellule ; et (4) la pente de merite physique
pour plusieurs directions deja disponibles, en se gardant de changer les
seuils pour obtenir un faux PASS. Tester ensuite une politique de pas coherente
sur un cas dont la comparaison des operateurs est fermee, avant BFS712.

## 5. Reproduction

```bash
cmake -S . -B build-v65 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_ENABLE_MPI=OFF -DCFDX_BUILD_PYTHON=OFF \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v65 --target cfdx_production_solver \
  test_comparable_physical_merit_v65 -j2
ctest --test-dir build-v65 -R 'v65' --output-on-failure

# Preparation des cas identiques ; campagne longue explicitement opt-in :
python3 scripts/run_bfs_globalization_v65.py \
  --base-case path/to/bfs152_v64_budget40.cfdx.h5 \
  --out /tmp/cfdx-bfs-v65-prepared --iterations 40
python3 scripts/run_bfs_globalization_v65.py \
  --base-case path/to/bfs152_v64_budget40.cfdx.h5 \
  --out /tmp/cfdx-bfs-v65-comparison --iterations 40 \
  --solver ./build-v65/cfdx_production_solver --execute --allow-partial
```

Le fichier source BFS152 HDF5 et les journaux des campagnes executees sont
fournis dans l'archive de livraison V65, sous `preuves/`.
