# CFDX — Rejeu exact des Fp COUPLED/PCD, NACA et BFS (2026-11-08)

## Portee et statut

Ce lot implémente une capture **strictement opt-in** des matrices CSR `Fp`, seconds membres et iterations initiales de la branche PCD COUPLED; un binaire de rejeu identique pour ILU0, equilibration ILU0, AMG et Jacobi; un ordonnanceur de campagnes NACA/N9 sans timeout par defaut; et des tests de non regression. La source est l'archive `CFDX_SOURCES_COMPLETES_PCD_FP_TH_SST_2026-11-08.zip` (SHA-256 a documenter dans le manifeste de livraison). Aucun changement des seuils de convergence non lineaires, aucun changement de discretisation physique, aucune promotion automatique des cas.

### Etat de l'operateur `Fp`

La branche historique PCD utilise `-Fp^-1` et non, de maniere implicite, le commutateur complet `-Mp^-1 Fp Kp^-1`. Pour NACA, `Fp` doit etre audite avec ses conditions de pression et sa jauge; N9 Couette montre que AMG peut fonctionner sur d'autres operateurs. Le maintien d'une seule ligne ancree dans un grand Fp presque a somme de lignes nulle est un indicateur de quasi-mode nul, pas une preuve de faute unique. Aucun nouveau traitement de jauge ne doit etre active en production sans comparaison du vrai Schur et de la conservation.

### Resultats de rejeu sur le meme CSR/RHS et meme x0 (FGMRES: 100 + 100 iterations max, cible 1e-8)

| Capture | Fp lignes | ILU0: residu vrai, temps total | ILU0 equilibre: residu vrai, temps total | AMG: residu vrai, temps total |
|---|---:|---:|---:|---:|
| N64 premiere RHS | 3096 | 1.92e-9, 0.0208s | 2.53e-8, 0.0300s | 0.995, 0.0485s |
| N64 RHS difficile | 3096 | 1.247e-3, 0.1516s | 1.299e-3, 0.1107s | 0.9999, 0.0784s |
| N128 premiere RHS | 12336 | 6.16e-9, 0.1917s | 3.80e-8, 0.2193s | ~1, 0.2179s |
| N128 RHS difficile | 12336 | 1.168e-3, 0.732s | 1.160e-3, 0.715s | ~1, 0.390s |
| N256 premiere RHS | 49248 | 8.63e-9, 0.692s | 1.10e-7, 0.832s | 7.81e-5, 2.33s |
| N9 Couette premiere RHS | 192 | 3.40e-9, 0.00108s | 3.50e-9, 0.000878s | 1.78e-9, 0.0334s |

Les temps incluent setup + premiere et eventuelle seconde passe. Ils ont ete mesures sur un binaire de qualification diagnostique peu optimise; utiliser les chiffres seulement pour comparer des variantes **sur le meme echantillon**. Le rejeu calcule le residu vrai sur la **matrice originale**. `status=CONVERGED` du solveur sur la matrice equilibree ne suffit pas si le residu original depasse la tolerance. Le choix ILU0 vs AMG ne se deduit pas d'une seule RHS facile: les premiers RHS difficiles N64/N128 montrent une stagnation ILU0 autour de 1e-3 apres 200 iterations.

### Localisation du defaut AMG

- N64 `Fp` : 3095 des 3096 lignes ont une somme nulle a la precision du test, un seul ancrage; dissymetrie relative ~0.898. Factorisation LU externe: pivot minimal ~6.0e-16, norme du vecteur solution ~1.58e12, vrai residu relatif ~3.22e-4 pour un RHS difficile. Il s'agit d'un indice de quasi-singularite.
- Une application du V-cycle AMG sur la premiere RHS N64 donne `||Fp * z|| ~ 3.37e14`, contre environ 0.987 pour ILU0. D'autres familles AMG testees amplifient egalement le residu. Sur N9 Couette, le V-cycle n'amplifie pas catastrophiquement et atteint des residus internes de ~1e-14 dans le pilote de production.
- Le nouveau garde-fou signale l'amplification extreme et **refuse** une action AMG numeriquement dangereuse sans substitution silencieuse; seuil de detection d'anomalie `1e6`, distinct de toute tolerance de convergence. C'est une protection logicielle experimentale, pas la resolution du probleme scientifique.
- Le choix des conditions limites de `Fp` doit etre traite comme un sujet de methode numerique: Elman & Tuminaro (ETNA 2009) montrent qu'un mauvais choix peut deteriorer fortement les preconditionneurs a commutateur approximatif. Reference : https://etna.ricam.oeaw.ac.at/volumes/2001-2010/vol35/abstract.php?pages=257-280

### Coupled N64 / N128 / N256 / N9 sur les budgets locaux

- NACA N64 : PCD historique diverge en resolution monolithique apres 79 iterations FGMRES, vrai residu relatif ~2.585e-3; sous relaxations successives ne modifient pas l'operateur lineaire initial. Le re-test du pilote nouvellement compile a reproduit deux rejets et a ete interrompu avant la troisieme tentative par la limite de calcul locale.
- NACA N128/N256 : les binaires ont ete lances avec un budget numerique long mais les essais effectues ici ont ete interrompus apres 180 secondes chacun par un plafond d'environnement. Aucune convergence physique de N128/N256 ne peut etre declaree.
- N9 PCD : le sous-cas Couette a ete franchi et des iterations Poiseuille ont ete observees; la campagne s'est arretee a son plafond local de 300 secondes avant d'atteindre Ghia. Ne pas declarer N9 PASS.
- `scripts/run_cfdx_coupled_campaign.py` est fourni pour relancer **sans plafond de temps par defaut**; il continue a limiter le nombre d'iterations physiques a un budget affiche (1500 par defaut) et exige les verifications physiques independantes.

### Etat des essais BFS et SST/thermohydraulique

Les calculs BFS physiques ci-dessous sont des **resultats historiques consignes dans le code source**, non reexecutes avec un solveur de production recompilé dans ce lot. Ils sont visibles dans `docs/validation/14-benchmarks/09-thermohydraulics/IMPLEMENTATION_STATUS_2026-11-03.md` :

| Etude historique | Convergence | x_r/h via cisaillement simple | x_r/h via reconstruction quadratique | Qualification |
|---|---|---:|---:|---|
| BFS SST 192 cellules, Re_h = 36000, entree uniforme | oui, 145 it. | 1.76882 | 2.26478 | non |
| BFS SST 402 cellules, memes conditions | oui, 220 it. | 2.22988 | 2.81068 | non |
| BFS SST 688 cellules | non, 45 it. | non mesurable | non mesurable | non |

Le contrat experimental TH-03 fige `Re_h=37500`, profils d'entree turbulent reel, reference d'accrochage `x_r/h=6.26 +/-0.10` et comparaisons Cf/Cp/U/contraintes turbulentes. Les donnees du debogage Re_h=36000 ne sont pas une campagne equivalente. Le site NASA TMR definit une variante CFD avec Re_h ~36000; l'identite des conditions d'entree doit cependant etre verifiee explicitement. Reference : https://tmbwg.github.io/turbmodels/backstep_val.html

- Cette passe : BFS quick C++ PASS (2560 cellules, coherence patches, `solver_campaign_pending`), tests Python SST/TH03/qualifier 64/64 PASS, ciblage CTest multi-domaine 10/10 PASS (PCD, turbulence, energie, radiation P1 et S2S, VOF). Ce sont des tests de composant ou de configuration, **non** de validation physique de TH-03.
- Thermique, canal, autres SST : preparation et contrats disponibles, mais aucun cas de validation physique complet multi-grilles n'est promu avec ce lot. Cas 3D DNS canal periodique, TH-07 convection mixte, TH-08 T-junction exigent encore simulations et traces reelles.

### Contrat de captures et garde-fous

- Capture activee seulement par `CFDX_PCD_FP_CAPTURE_PREFIX` / `CFDX_PCD_FP_CAPTURE_DIFFICULT_PREFIX` ; premiere serie et serie difficile limitees a 3 fichiers chacune **pour tout le processus** et non par appel COUPLED; aucun ecrasement silencieux lors de la transition entre sous-cas N9.
- Fichier `CFDXFP01` contient matrice CSR et vecteurs **originaux** sans equilibration; le rejeu compare exactement les memes inconnues. Le fichier est atomiquement ecrit et refuse la corruption.
- Le probe AMG n'ecrase jamais la correction physique et ne fait pas de substitution automatique. Trace a chaque action seulement si `CFDX_PCD_FP_AMG_TRACE=1`, sinon journalise les echanges dangereux.
- Rapport du rejeu en JSON/CSV, hash SHA256 des captures et du binaire; si `--require-all-true-residuals` est active, un seul echec de precision fait echouer le script.

### Commandes de reproduction

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF
cmake --build build --target test_pcd_fp_replay test_pcd_fp_equilibration test_naca0012_qualification test_n9_physical_matrix
ctest --test-dir build -R 'test_pcd_fp_replay_capture|test_pcd_fp_equilibration' --output-on-failure
python scripts/pcd_fp_replay_campaign.py --capture path/to/naca_hard_n64_rhs1.cfp --executable build/test_pcd_fp_replay --output-dir out/fp-comparison
python scripts/run_cfdx_coupled_campaign.py --build build --meshes build/validation_meshes --out out/coupled-full --only naca --models pcd block_local lsc bfbt simple simplec --levels 64 128 256 --iterations 1500
python scripts/run_cfdx_coupled_campaign.py --build build --meshes build/validation_meshes --out out/n9-full --only n9
```

**Interdits de promotion :** ne jamais attribuer un PASS physique a un timeout, un resultat NACA a une valeur BFS, un calcul BFS avec entree uniforme a l'experience Driver-Seegmiller, ou un simple test de Schur a une convergence COUPLED. CI GitHub complete et GCI NACA/BFS non verifies dans cette passe.
