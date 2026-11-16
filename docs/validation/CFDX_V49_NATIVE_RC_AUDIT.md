# CFDX V49 - Audit de l'operateur Rhie-Chow natif et des residus SST sur BFS152

**Date : 12 novembre 2026.** Lot experimental non intrusif, issu de la vraie V48. **Aucune convergence physique BFS n'est revendiquee.** Pas de push GitHub, pas de modification des seuils de production, et nouvelle instrumentation OFF par defaut.

## 1. Contexte et probleme

V47 avait obtenu sept etats exterieurs puis s'etait arrete ; V48 utilisait une homotopie des coefficients inverses de momentum et avait obtenu onze etats avant l'arret au douzieme. V48 a conserve un verrou de non-qualification a beta < 1, mais l'etat (U,p,k,omega) reste susceptible d'avoir plusieurs evaluations du flux en fonction de l'ancienne diagonale momentum. Les bilans physiques V48 portent sur un flux de continuation : ces bilans ne demontrent pas automatiquement les residus de l'operateur natif.

Objectif V49 : mesurer, sans changer la solution de production, l'ecart entre le flux conserve par V48 et un flux Rhie-Chow dont les coefficients momentum sont reconstruits par point fixe **a partir de l'etat courant seulement**. V49 est une etape vers F_phys(X), **pas encore un oracle complet ni un Newton-Krylov global**.

## 2. Operateur de diagnostic

Pour un etat fixe U,p,k,omega (et sa fermeture nu_t, ses conditions limites et coefficients de paroi courants), V49 initialise un flux **deterministe** par interpolation des vitesses de cellule et des conditions limites, sans utiliser phi_committed ou rAU_committed. Puis il recherche un point fixe :

```
phi_0 = flux_interpole(U, conditions_limites)
for j = 0..95:
    A_U(phi_j) = assemblage_momentum(U,p,nu_t,phi_j,conditions_limites)
    rAU_j = diag(A_U(phi_j))^{-1}
    phi_RC = RhieChow_production(U,p,rAU_j,conditions_limites)
    verifier la difference RELATIVE en norme infinie
    si la difference est suffisamment petite:
        reassembler A_U(phi_RC), reconstruire de nouveau phi_RC
        n'accepter la fermeture que si le CONTROLE FINAL est aussi satisfait
    sinon phi_{j+1} = 0.5*phi_j + 0.5*phi_RC
```

Critere du point fixe de diagnostic `||phi_RC-phi_j||inf/max(||phi_RC||inf,||phi_j||inf,1e-300) <= 1e-10` ; budget 96 et amortissement interne 0,5. **Ce critere ne remplace aucune tolerance CFDX et ne garantit pas l'unicite du point fixe.** Une premiere implementation utilisait `max(1,||phi||)` : elle a ete corrigee avant la livraison, car les flux BFS152 sont de l'ordre de 0,0076 kg/s et une normalisation a 1 aurait surestime la precision relative.

Si la fermeture echoue, `closed=0` et aucun champ flux/rAU candidat n'est promu. Si la fermeture reussit, V49 reassemble ensuite le momentum natif et re-evalue le transport de k/omega avec ce **meme** flux, puis mesure la divergence conservative. Les fonctions de paroi et termes sources sont ceux de l'etat produit par l'iteration courante, sans re-resolution couplée de ces champs. Le resultat est un diagnostic a etat gele, pas encore un F(X) mathematiquement complet et qualifie.

Le mode `CFDX_NATIVE_RC_ORACLE_V49_REPEAT=1` lance deux evaluations independantes du meme etat : les nombres de pas, le flux, les trois champs rAU et le defaut doivent etre **strictement egaux**. Aucune memoire d'iteration anterieure n'est fournie aux callbacks de l'oracle.

## 3. Mesures de la production sur le BFS152

**Rejeu avec la version finale et l'echelle relative stricte.** Le binaire a compile en Release de diagnostic, avec `CFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON`. Les lignes de commit V13 pour les **trois premiers etats** sont textuellement identiques a celles du journal V48 original.

| Tentative | Iterations point fixe | Defaut relatif point fixe | Ecart L1 relatif phi_V48 / phi_natif | U natif | Continuité native | k natif | omega natif | Etat |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| Ext. 1, retry 0 | 31 | 2,3611e-12 | 4,86895e-3 | 1,24797 | 8,56622e-3 | 6,34297e-2 | 5,20243e-2 | Commit V48 |
| Ext. 2, retry 0 | 32 | 3,1437e-12 | 1,51744e-3 | 1,02257 | 9,91167e-3 | 1,17431e-2 | 1,45320e-2 | Rejet V13 |
| Ext. 2, retry 1 | 32 | 2,8155e-12 | 6,19269e-3 | 0,85837 | 5,16033e-3 | 1,06363e-1 | 1,41501e-1 | Commit V48 |
| Ext. 3, retry 0 | 40 | 7,2584e-13 | 1,41131e-3 | 0,236045 | 6,31831e-3 | 1,21915e-2 | 4,25473e-2 | Commit V48 |

**Comparaison pour les etats effectivement committed** :

| Etat | U V48 -> natif | Continuite V48 -> native | k V48 -> natif | omega V48 -> natif |
|---|---|---|---|---|
| 1 | 1,898668 -> 1,247969 | 1,266899e-2 -> 8,566224e-3 | 1,877770e-11 -> 6,342968e-2 | 1,691289e-13 -> 5,202434e-2 |
| 2 | 0,504888 -> 0,858369 | 9,040057e-3 -> 5,160332e-3 | 2,623040e-11 -> 1,063631e-1 | 1,223308e-13 -> 1,415012e-1 |
| 3 | 0,236013 -> 0,236045 | 6,326746e-3 -> 6,318315e-3 | 1,757976e-11 -> 1,219147e-2 | 2,703967e-13 -> 4,254726e-2 |

Ces differences ne sont pas un simple bruit FP64 : elles concernent des equations de turbulence evaluees avec deux flux differents. Cela **ne prouve pas a lui seul** que la physique native est correcte ni que V48 est incorrect ; cela montre la non-equivalence entre les operateurs utilises pour certifier les residus. Une convergence physique de CFDX ne peut pas etre revendiquee avant d'identifier l'operateur discret unique a qualifier.

**Reproductibilite :** le rejeu de l'etat 1 avec `CFDX_NATIVE_RC_ORACLE_V49_REPEAT=1` a retrouve bit a bit toutes les sorties internes de l'oracle ; le commit V13 est inchange. Les deux autres executions utilisent l'oracle non intrusif. Le code retour du solveur avec `--iterations 3` est `1`, car les criteres de convergence finale ne sont pas atteints : cela n'est ni une erreur de l'instrument ni une reussite de qualification.

## 4. Tests et controle logiciel

- Build CMake Release diagnostic avec MPI desactive ; `cfdx_production_solver` compile et lie (le `main.cpp` est compile a `-O0` par l'option diagnostique : pas de benchmark O2).
- **8/8 CTest cibles PASS** : nouvel evaluateur C++ V49, verificateur Python V49, et regressions V48, V47 (2), V14, V13, V41.
- **13/13 cas Python PASS** : preuve absente, double essai, residu non fini, echec explicite, fausse qualification, fausse fermeture, reference modifiee, echelle manquante, drapeau de rejeu invalide, etc.
- **4/4 oracles natifs fermes** dans le rejeu 3 etats, avec tolérance relative stricte `1e-10`.
- **3/3 commits V13 identiques** au journal V48 d'origine (lignes textuelles).
- **Rejeu strict** d'un etat : flux et rAU identiques au bit pres entre deux evaluations successives.

Aucun `git push` n'a ete effectue. Pas de CI totale, pas de sanitizers complets, pas de qualification GPU/MPI, pas de run BFS712 ni BFS20 540, pas de convergence finale.

## 5. Sources bibliographiques et interpretation

- OpenFOAM `simpleFoam/pEqn.H` (code historique accessible) : `rAU`, `HbyA`, `phiHbyA`, correction de pression et flux conservatif. https://github.com/OpenFOAM/OpenFOAM-2.3.x/blob/master/applications/solvers/incompressible/simpleFoam/pEqn.H . V49 **ne reproduit pas** cet algorithme SIMPLE, mais s'en sert pour exiger une coherence de flux/momentum.
- PETSc, `SNES` : evaluation d'une fonction non lineaire unique au candidat, comparaison des Jacobiennes et recherche lineaire. https://petsc.org/main/manual/snes/ . V49 n'est pas un solveur SNES/FGMRES global.
- *Influence of momentum interpolation methods on the accuracy and convergence of pressure-velocity coupling algorithms in OpenFOAM*, Journal of Computational and Applied Mathematics 309 (2017) 654-673, DOI https://doi.org/10.1016/j.cam.2016.03.037 : une interpolation de momentum mal controlee peut introduire une dependance aux facteurs de relaxation.

## 6. Decision technique et suite

Conserver V48 OFF et le microscope V49 **OFF** par defaut. L'oracle V49 reste read-only : il ne modifie jamais U,p,phi,k,omega, les coefficients de production, les seuils, la selection des six Schur ou le rollback.

**V50 recommandee :** (1) extraire un evaluateur F(U,p,k,omega) sur un meme flux natif coherent, y compris les mises a jour de paroi et les sources ; (2) confronter F(X) sous differents historiques et variations de directions (Jv en difference finie), (3) construire une vraie direction de descente et Armijo sur ce F ; (4) confronter les 6 Schur a empreinte A,b,x0 commune, puis BFS712 et BFS20 540, avant validation Driver-Seegmiller. Ne pas assouplir les seuils ni declarer valide un etat natif non stationnaire.

**Verdict :** progres de qualification **operateur**, pas progression demontree de la solution physique. Les six equations du BFS152 et les campagnes industrielles demeurent non qualifiees.
