@@ -0,0 +1,78 @@
# CFDX V70 - Conservation du flux COUPLED/SST sur BFS152

## Portee et regle de certification

Cette livraison repond au diagnostic V68/V69 : la solution SST est precise
avec le flux porteur COUPLED, mais le residu de omega est tres different sur
le flux Rhie-Chow reassemble a partir de l'etat candidat. V69 avait trouve
un candidat *prive* beta=1/256 admissible aux gardes V65, sans certifier
sa fermeture de masse ni une transaction de production.

V70 teste une projection candidate-propre sous la forme

    delta U_P = -V_P rAU_P grad(delta p)_P
    p' = p + delta p
    phi' = RC(U',p';rAU_P)

Le coefficient rAU est fige UNIQUEMENT pendant la construction de la
reponse de pression. Pour chaque etat projete, l'evaluateur physique V51
reconstruit tous les coefficients, la viscosite turbulente et les flux.
Le SST est alors transporte de nouveau avec ce flux, sur une copie privee.

La conservation provisoire de la matrice de pression n'est pas un
certificat de conservation final. Le bilan par cellule, la somme des
flux de bord, le residu physique complet, la reproductibilite de
Rhie-Chow et les gardes de croissance/merite V65 doivent etre reverifies.
La tolerance absolue de divergence de la certification experimentale
reste 1e-10 (aucun assouplissement).

Une correction complete de pression peut sortir du bassin de fermeture de
l'operateur RC reassemble. V70 peut alors examiner une continuation
fractionnaire 1,1/2,...,1/128, au plus trois cycles, a seule fin de
localiser une direction qui reduit la vraie divergence. Une fraction
partielle n'est jamais presentee comme projection conservative totale.

## Non-commit et reprise en cas d'echec

Le noyau existant conserve ses six preconditionneurs Schur (ainsi que
l'option selfp existante), ses equations SST, ses tolerances et ses
conditions de paroi. Le mode V70 est experimental et desactive par defaut.
Il s'active uniquement avec les gardes et les diagnostics V51, V65, V68,
V69 deja disponibles : `CFDX_FPHYS_V70_CONSERVATIVE_SST_RC=1`.

Les sondes V70 ne modifient jamais le vrai U,p,phi,k,omega,nu_t. Une
verification explicite compare les champs de production a leur snapshot.
En cas de fuite, le chemin de rejet natif restaure tous les champs. Les
tests V70 couvrent aussi un objet SSTCoupledStateTransaction ainsi que le
schema de rollback reel du COUPLED (U,p et snapshots phi/SST).

**Important :** le code V70 n'engage PAS un candidat de production.
Un vrai commit de ce nouveau candidat demanderait en plus la recreation
du certificat momentum/continuite V8 associe a U',p',phi',k',omega',nu_t'
ainsi qu'une preuve de convergence reassemblee et un rollback complet.
Il serait errone d'utiliser le certificat obtenu avant projection.

La projection dense du prototype est volontairement limitee a 152
cellules. Elle n'est ni un solveur pression industriel, ni un benchmark
de performance, ni une methode a extrapoler au BFS712/20540.

## Reproduction courte

```
cmake -S . -B build-v70 -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON
cmake --build build-v70 --target test_bfs_sst_rc_conservative_v70 \
  test_bfs_v70_full_rollback test_bfs_v70_evidence_falsification
ctest --test-dir build-v70 --output-on-failure -R 'bfs_sst_rc_conservative_v70|bfs_v70_full_rollback|bfs_v70_evidence_falsification'
```

Les journaux BFS152 et l'audit de non-regression sont livres dans l'archive
complete. La commande `scripts/audit_bfs_v70_conservative.py` compare les
traces V65/V69 avec celles du rejeu V70 et rejette toute falsification d'un
commit ou d'une conservation physique non demontree.

## Statut

Aucune modification des seuils, aucun push GitHub, aucune validation
industrielle de BFS152, BFS712 ou BFS20540. Les cas de qualification
thermique/rayonnement V57 sont conserves mais non relances longuement.
