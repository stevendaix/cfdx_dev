@@ -0,0 +1,144 @@
# CFDX V71 — audit de co-itération momentum–pression–flux–SST

## Contrat du lot

La V70 a projeté pression et vitesse *après* les deux transports SST privés. Le flux
se refermait à environ `1e-13` sur les coefficients figés, mais la divergence
`F_phys` réassemblée restait `6.84e-5` après trois cycles et le mérite V65
se dégradait. La V71 teste un ordre d'actualisation distinct sur BFS152 :

1. Évaluer les six équations sur un seul état `X=(Ux,Uy,Uz,p*,k,omega)`.
2. Appliquer le prédicteur momentum existant en diagonale,
   `U <- U - gamma*diag(A_momentum)^{-1}*(A_momentum*U-b)`, sans modifier `p,k,omega`.
   La direction emploie **exactement** les diagonales `rAU` et les résidus signés
   de ce même assemblage. `Uz=0` est préservé pour une grille extrudée 2D.
3. Recalculer entièrement le résidu physique, viscosité SST et flux RC inclus.
4. Faire la correction couplée `U,p` avec la réponse de pression/flux Rhie–Chow
   du V70 et les `rAU` de **l'état prédit**, et non ceux du carrier historique.
5. Réévaluer `F_phys`, puis résoudre le **véritable transport SST** sur son flux
   reconstruit et reconstituer de nouveau `F_phys`.
6. Comparer les six résidus au **même ancien état**, avec les échelles V65 figées,
   vérifier explicitement les flux de face et le bilan de masse réassemblé.
7. Sur la face interne `3` (cellules 0/12), journaliser `phi`, `rAUx`, `nu_t`,
   les résidus `omega` et les différences de convection, diffusion, production,
   destruction et diffusion croisée via l'oracle discret indépendant V16/V66.

La V71 est activée seulement par `CFDX_FPHYS_V71_JOINT_PICARD=1`, en combinaison
avec V51, V65, V68 et V69. La variante V70 doit rester désactivée lors de la
comparaison. Uniquement pour le pilote `BFS152`, tentative 23/retry 8,
`gamma` dans `{0,1/32,1/8}` et fraction de pression dans `{1/2,1/4}`.

**Important :** ces six co-itérations sont **privées** et sans voie de commit.
Même un candidat satisfaisant le mérite V65 n'est pas qualifié tant que le flux
conservatif, les certificats momentum et SST, le replay et le rollback du
nouvel état ne sont pas également validés. Aucun changement de seuil ou de
modèle physique n'est autorisé.

## Tests et falsifications

- `test_bfs_joint_momentum_sst_v71` : identités de pression/flux,
  prédicteur diagonal, retransport SST avec le vrai flux, traces de l'état final,
  refus de `rAU` incohérent, dépassement de fraction et manipulation cachée
  de `omega` par la réponse de pression.
- `test_bfs_v71_evidence_falsification` : huit journaux falsifiés rejetés,
  aucune preuve de convergence à partir d'un essai privé.
- Régressions V69 et V70 : `test_bfs_sst_flux_picard_v69`,
  `test_bfs_sst_rc_conservative_v70`, `test_bfs_v70_full_rollback`.
- Audit `scripts/audit_bfs_v71_joint.py` : identité stricte des traces V65,
  des états acceptés, du transport privé V69 et des traces d'itération.

## Distinction physique / algébrique

La projection avec `rAU` figé ne fournit **pas** la conservation sur le flux
Physique Rhie–Chow après réassemblage du momentum/SST. Une réduction du résidu
`omega` à elle seule ne prouve **pas** la convergence physique COUPLED. Une
convergence de transport SST pour un flux donné n'implique pas la convergence du
transport pour le flux réévalué après changement de `nu_t`.

## Qualification restante

L'essai s'effectue sur le pilote `TH03/BFS152` et sa fonction de paroi
historique, sans les fonctions de paroi OpenFOAM de production. La qualification
physique Driver–Seegmiller, les BFS712/BFS20540, les études de maillage, NACA,
les campagnes radiatives longues et la CI complète restent distinctes.

## Rejeu BFS152 — résultats obtenus

Build CMake diagnostique : `Release` avec optimisations compilateur réduites
`-O0 -DNDEBUG` (y compris sur le coeur), et
`CFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON`. Il ne s'agit **pas** de mesures
de vitesse industrielle. Le binaire du solveur de production est compilé.

Le pilote réel a reproduit **149/149 décisions V65**, **22/22 états acceptés**,
**6/6 traces V69** et **164/164 traces INCOMPRESSIBLE_ITER** exactement à partir
de la référence fournie. Code de sortie final du solveur : `2`, échec historique
à l'itération 23. **Aucun nouvel état V71 accepté ou engagé.**

Six co-itérations indépendantes ont été effectuées à l'itération 23,
`retry=8` (`beta=1/256`), toujours à partir du même candidat :

| Gamma momentum | Fraction pression | Merit V65 | Croissance maximale | Masse F_phys L_inf | Omega F_phys L_inf | Garde V65 |
|---:|---:|---:|---:|---:|---:|---|
| 0 | 1/2 | 1.086784 | 12.8694 | 3.87401e-4 | 3.15529e-2 | NON |
| 0 | 1/4 | 0.952308 | 5.90678 | 6.10722e-4 | 1.46679e-2 | NON |
| 1/32 | 1/2 | 1.082802 | 4.56199 | 3.95882e-4 | 1.13604e-2 | NON |
| **1/32** | **1/4** | **0.928348** | **0.928348** | **6.24857e-4** | **7.07376e-4** | **OUI** |
| 1/8 | 1/2 | 1.070647 | 8.96282 | 4.21982e-4 | 2.26949e-2 | NON |
| 1/8 | 1/4 | 0.866005 | 11.1907 | 6.67313e-4 | 2.79319e-2 | NON |

**Attention :** la meilleure contraction du mérite (0.8660) n'est *pas*
admissible, car le bloc SST croît 11.19 fois. Le seul candidat qui satisfait
**les deux gardes V65 inchangés** est gamma=1/32, pression=1/4.
Il fait passer le résidu brut omega de 4.72060e-3 à 7.07376e-4, mais sa
**divergence physique de 6.24857e-4 excède de plusieurs ordres de grandeur**
la condition conservative absolue de 1e-10. Le candidat ne peut pas être
engagé. Le succès V65 n'est pas un certificat de convergence globale.

### Audit local de la face 3 pour l'essai le moins mauvais admissible V65

La face 3 est bien connectée entre les cellules 0 et 12 :

- Flux physique face 3 : 4.647820e-4 vers 3.428729e-4 kg/s.
- Cellule 0 : `rAUx` 1291.714 -> 1296.439 ; `nu_t` 0.00914235 -> 0.00902796.
- Cellule 12 : `rAUx` 516.474 -> 524.387 ; `nu_t` 0.00621995 -> 0.00623654.
- Cellule 12 : résidu omega signé 4.720596e-3 vers 6.354730e-4.
- Cellule 12, différences des termes du bilan omega : convection -0.536974,
  diffusion +0.089586, production -0.186802, destruction +0.633541,
  diffusion croisée -0.003437 (unités brutes de l'équation volumes finis).
  Leurs compensations expliquent la variation du résidu signé.

Le défaut maximal de fermeture des **12 lignes de bilans signés** omega
(cellules 0 et 12, six essais) vaut **4.05e-13**. L'outil d'audit indépendant
vérifie cette identité et rejette les falsifications. Cette précision valide
la décomposition, pas le modèle de turbulence vis-à-vis des expériences.

### Conséquence numérique

La première mise à jour coordonnée momentum/pression/SST peut faire baisser
simultanément le mérite et omega, **sans fermer la masse**. Le prochain
verrou ne doit plus être présenté comme un problème purement SST ou de
préconditionneur Schur : il faut faire converger **le point fixe pression–flux–rAU**
sur le même état que le transport SST. Le remplacement direct du flux porteur
ou une seule correction de pression après SST ne satisfont pas ce contrat.

Une éventuelle continuation doit répéter, sur une copie privée,
`momentum -> pression/flux -> SST -> F_phys`, **avec une acceptation uniquement
sur conservation physique après réassemblage, résidus des six blocs et
certificats linéaires/rollback entièrement réémis**. Un projeté conservatif
à coefficients gelés ne suffit pas. Aucun nouveau Schur ou Jacobien n'est
justifié par ces résultats.

## Vérifications de livraison

- 5/5 CTest ciblés réussis : V71 C++, V71 Python, SST V69 et projection/rollback V70.
- Test C++ V71 ASAN/UBSAN exécuté avec succès.
- 10 falsifications Python rejetées (journal modifié, publication illégale,
  `rAU` invalide, `NaN`, bilan omega incorrect et paire de cellules manquante).
- Patch cumulatif BFS V6 -> V71 et incrémental V70 -> V71 appliqués sur des
  extractions neuves : **2190/2190 fichiers identiques en contenu et droits
  d'exécution**.

**Limites :** comparaison sur un pilote `TH03_debug_only/BFS152` avec SST
historique de paroi, pas sur le Driver–Seegmiller physique validé. Aucune
qualification BFS712/BFS20540, fonctions de paroi OpenFOAM actives,
longue campagne NACA/rayonnement, ou CI complète dans ce lot.
