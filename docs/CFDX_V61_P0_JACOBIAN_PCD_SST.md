# CFDX V61 — verrou P0 : Jacobien physique, PCD, momentum–SST (13 novembre 2026)

## Périmètre de la livraison

Base : sources fusionnées CFDX V60 (rayonnement CPU HO–LO/DOM V57 et BFS/SST V59–V60).
Ce lot introduit un calcul optionnel du Jacobien explicite figé 6N sur deux pas,
un audit indépendant à trois échelles du résidu physique momentum–Rhie–Chow–SST,
et un oracle SciPy de l'opérateur de Schur exact sur les matrices réelles 4N de BFS152.
Le solveur de production est compilé, mais aucun correctif numérique de la physique
primal ni des six préconditionneurs n'est activé par défaut.

### 1. Jacobien figé V61 : Richardson sous refus des non-lissités

`FrozenPhysicalJacobianV54` conserve la version historique sans changement par défaut.
En mode `richardson_v61`, une colonne j est calculée à partir de la même fonction
physique `Fphys` aux perturbations ±h et ±h/2 (réassemblages complets), avec

`D(h) = [F(x+h e_j)-F(x-h e_j)]/(2h)`

et, seulement dans une zone suffisamment lisse,

`J_j = [4 D(h/2)-D(h)]/3`.

L'accord des pentes latérales à h/2 et entre h et h/2 doit respecter la limite
V54 de fidélité existante, sans la relever. Si la colonne est non admissible,
les évaluations sont non finies ou une branche reste ambiguë, aucun Newton
ne peut l'accepter. Cette extrapolation ne constitue pas en soi un Jacobien
semi-lisse fiable aux égalités multiples SOU/MUSCL.

La direction d'Arnoldi reste multipliée par une unique matrice figée `Amat`;
la vraie convergence linéaire est auditée indépendamment du résidu estimé GMRES.
L'acceptation finale ne peut s'appuyer que sur la véritable fonction `Fphys`.

Sur le test synthétique lisse cubique, le défaut V54 `4.04959e-6` devient
`3.21965e-15` avec Richardson ; une fonction à cassure est rejetée et aucune
direction Newton n'est acceptée. Ce ne sont PAS des résultats de BFS convergé.

### 2. Audit candidat V61 à h, h/2 et h/4

`PhysicalRemainderRatesV61` réassemble les équations momentum, continuité,
convection, transport SST et flux Rhie–Chow sur les candidats de chaque échelle.
L'audit refuse toute échelle dont le défaut d'additivité est supérieur à la
limite inchangée `1e-6`, y compris si la norme globale masque un bloc mauvais.
Le constat de défaut persistant est un indicateur, pas une preuve de cassure
mathématique : erreur de dérivation, bruit et variation de branche doivent
encore être départagés par les audits face par face V58 et V59.

Rejeu BFS152 réel : 152 cellules, 634 faces ; état initial et paramètres du
cas V60, schéma SOU, SST spatial expérimental, BFBt, six directions Krylov :

| Indicateur | Valeur |
|---|---:|
| Défaut le plus fort à h | 9.825045378947e-5 |
| Défaut le plus fort à h/2 | 1.255781091194e-4 |
| Défaut le plus fort à h/4 | 1.511173645450e-4 |
| V60 : additivité globale du reste hors convection | 2.808072846303e-7 |
| V60 : convection | 9.825045378947e-5 |
| V60 : momentum | 2.631207010124e-6 |
| V60 : continuité | 2.196891215979e-6 |
| V60 : k | 1.304044178932e-6 |
| V60 : omega | 2.804885281009e-7 |
| V60 : flux Rhie–Chow | 1.516025002379e-6 |
| V60 : rAU | 4.093069076918e-7 |
| V60 : viscosité | 5.385566324159e-7 |
| V60 : contrainte/strain | 1.053092908721e-6 |
| V52 : résidu linéaire vrai après six Krylov | 0.9579208419022 |
| V52 : défaut d'additivité Jv | 1.824151640401e-4 |
| Newton accepté | 0 |
| BFS convergé | NON |

Ces indicateurs sont préliminaires : une seule itération BFS152 est exécutée
avec un certificat volontairement négatif, pas une campagne de convergence.
La sélection `CFDX_FPHYS_JFNK_V61_RATE_AUDIT=1` nécessite l'audit V60 actif.
Le mode Richardson est séparé et exige `CFDX_FPHYS_JFNK_V54_FROZEN=1` ainsi
que l'audit du véritable résidu Krylov V57. Rien n'est actif par défaut.

### 3. PCD : identifier le défaut du préconditionnement, pas du système physique

Les captures MatrixMarket communes aux six Schur V60 ont des SHA-256 identiques
pour `A`, `b`, `x0` (608 inconnues, 7 080 coefficients non nuls, 152 pressions).
L'oracle `scripts/pcd_v61_exact_schur_oracle.py` construit en lecture seule :

`S = C - D M_u^{-1} G` (Schur exact 4N, non pas Jacobien physique 6N).

Le LU sparse de la matrice 4N initiale a un résidu vrai relatif `1.361e-15`;
le Schur pression est inversible dans ce cas, avec conditionnement spectral
`5.3954e4`. C'est une preuve algébrique limitée à cette matrice figée,
non une correction et non une solution convergée BFS.

L'audit V38 des vrais coefficients auxiliaires capturés signale :

| Variante de l'opérateur auxiliaire | Rang numérique Fp | Résidu COUPLED vrai après 2 000 itérations |
|---|---:|---:|
| PCD historique | 151/152 | 6.274552e-3 |
| Condition Dirichlet d'entrée pour Fp | 152/152 | 4.936761e-3 |
| Dirichlet + commutateur de masse lumped | 152/152 | 6.199061e-3 |

Les trois méthodes échouent à converger au seuil original. La comparaison
idéalisée au Schur exact donne, sans reproduire les itérations du solveur :

- `||S * (-Fp^+) - I||F / sqrt(N) = 864.9874` pour PCD historique ;
- `||S * (-diag(Mp)^-1 Fp Kp^+) - I||F / sqrt(N) = 28.6036` ;
- même avec le Dirichlet auxiliaire, écart ≈ 778.69 et 27.44, respectivement.

Les plus petites valeurs singulières changent selon la discrétisation et la
mise à l'échelle ; ces chiffres ne sont pas des critères d'équivalence
physique. La correction pertinente est de vérifier les BC et les nullspaces
de `Fp`/`Kp`, puis de construire un préconditionnement réellement compatible
avec `D`, `G`, `Mu` et le gauge pression sur les cas BFS. Il ne faut ni
inverser arbitrairement `Fp` ni supprimer la contrainte de continuité.

### 4. Distinction entre preuves et conditions de qualification

- `BFS152` : diagnostic réel et oracle Schur, pas de Newton global accepté.
- `PCD` : défaillance confirmée sur le 4N et non rétablie.
- `LSC`, `BFBt`, `BlockLocal`, `SIMPLE-Schur`, `SIMPLEC-Schur` : états de
  campagne V60 inchangés ; ils ne sont pas physiquement qualifiés en V61.
- `SST` : contributions `k`, `omega`, viscosité et couplages dérivés
  observés mais pas certifiés contre un oracle OpenFOAM actif sur BFS.
- `BFS712`, `BFS20540` et x_r/h de référence : NON REJOUÉS dans V61 ;
  impossible de conclure à une validation physique ou expérimentale.
- Branche radiative V57 fusionnée en V60 : conservée et optionnelle,
  mais **campagne longue radiative non relancée pour V61**.
- Pas de nouveau lancement CI complet, ni MPI/GPU, ni push GitHub.

### 5. Reproduction minimale

```bash
cmake -S . -B build-v61 -G Ninja -DCMAKE_BUILD_TYPE=Release \
  -DCFDX_BUILD_TESTS=ON -DCFDX_ENABLE_MPI=OFF \
  -DCFDX_BUILD_PYTHON=OFF -DCFDX_USE_HDF5_SERIAL=OFF \
  -DCFDX_ENABLE_RADIATION_RESEARCH_V57=ON \
  -DCFDX_DIAGNOSTIC_FAST_COMPILE=ON \
  -DCFDX_DIAGNOSTIC_FAST_PRODUCTION_COMPILE=ON \
  -DCFDX_ENABLE_PCD_V61_FORENSIC_TESTS=ON
cmake --build build-v61 --target \
  test_frozen_richardson_v61 test_physical_remainder_rates_v61 \
  test_frozen_physical_jacobian_v54 test_physical_remainder_audit_v60 \
  test_global_physical_newton_v52 test_pcd_schur \
  test_pcd_fp_bisided_basis_v13 cfdx_production_solver
ctest --test-dir build-v61 --output-on-failure -R \
 '^(test_frozen_richardson_v61|test_physical_remainder_rates_v61|test_frozen_physical_jacobian_v54|test_physical_remainder_audit_v60|test_global_physical_newton_v52|test_pcd_schur|test_pcd_fp_bisided_basis_v13|test_pcd_v61_exact_schur_oracle)$'
```

La compilation rapide ciblée à `-O0` pour certaines unités est réservée au
diagnostic, pas à la mesure de performance. La reproduction des captures
PCD et du BFS152 nécessite le HDF5 joint au paquet de preuves et les
MatrixMarket de la livraison précédente ; le fichier `preuves/README.md`
détaille les commandes et les modes qui échouent normalement.

### 6. Suite obligatoire (pas une validation acquise)

1. Sur BFS152, classifier les non-lissités par active-set SOU, Rhie–Chow,
   SST/F2 et wall-functions, avec oracles de faces et Jacobien généralisé
   figé **à l'état de base** ; reconstruire la réponse physique complète.
2. Vérifier chaque vecteur Arnoldi : fidélité, vraie linéarité, vrai résidu
   recalculé, puis la descente de `Fphys` pour le pas globalisé.
3. Corriger PCD au niveau pression `Fp,Kp,Mp`, BC et projection de
   nullspace ; comparer les six Schur avec mêmes A,b,x0 et horodatage.
4. Obtenir plusieurs contractions physiques sur BFS152 avec parois SST
   actives, puis seulement lancer BFS712 et BFS20540 et comparer
   `Cf`, `Cp`, profils, contraintes et réattachement aux références NASA.

**Verdict : V61 améliore la mesure et la protection du Newton, et précise
la racine de l'échec PCD ; le verrou P0 reste ouvert.**