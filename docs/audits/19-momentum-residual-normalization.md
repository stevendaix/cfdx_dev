# FICHE MÉTHODE — Normalisation du résidu momentum dans la convergence incompressible

**En-tête d’évidence (gouvernance §16.5)**

- Scope : critère d’acceptation non linéaire de l’équation momentum dans `solve_steady_incompressible()`.
- Statut : DRAFT — la correction candidate requiert CI et `validation-total` sur le HEAD exact.
- Registre autoritaire : N12 coupled pressure-velocity MMS / issue #905 ; statut de qualification inchangé jusqu’aux preuves.
- PR liée : #956.
- Contrat : aucune modification de discrétisation, tolérance, budget d’itérations, gate de continuité, ni oracle physique.

## 1. Références sources (obligatoire : DOI ou ISBN)

- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI: 10.1201/9781482234213 (solution itérative et contrôle de convergence).
- Y. Saad and M. H. Schultz, “GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems,” *SIAM Journal on Scientific and Statistical Computing*, 7(3), 1986, pp. 856–869 — DOI: 10.1137/0907058 (résidu algébrique et méthode de Krylov ; la référence n’impose pas le gate non linéaire décrit ici).

## 2. Formulation continue

Pour une équation de quantité de mouvement discrétisée (A_u U=b_u), le résidu algébrique est (r_u=A_uU-b_u). La convergence doit mesurer la satisfaction de cette équation ainsi que les contraintes de continuité et les variations de champs. Le rapport au résidu de la première itération décrit une progression relative au point initial, mais n’est pas invariant à un choix initial déjà exact ou proche du bruit machine.

## 3. Discrétisation Volume Fini

CFDX assemble (A_u), (b_u) et le résidu dans l’opérateur volumes finis existant. Deux normalisations sont exposées dans l’historique :

- Résidu RHS-normalisé : (eta_b=|r_u|/max(|b_u|,epsilon)).
- Résidu relatif à l’itération initiale : (eta_0=|r_u|/max(|r_u^{(0)}|,epsilon)).

- Le résidu RHS-normalisé reste toujours un gate d’acceptation. Pour les algorithmes ségrégés, le ratio à la première itération reste exigé lorsque sa référence est résolue ; si la référence est au niveau du bruit machine (référence <= 100 epsilon machine × échelle RHS), ce ratio est exclu du monitor et du gate, tout en restant disponible dans l’historique. Pour `COUPLED`, le résidu est mesuré après la résolution du système monolithique : le résidu de la première itération n’est donc pas un état initial indépendant. Le ratio `eta_0` est exclu du monitor non linéaire pour `COUPLED`, quelle que soit cette valeur de référence ; le gate RHS-normalisé et les gates physiques restent obligatoires. Aucun terme de matrice, flux de face ou source discrète n’est modifié.

## 4. Traitement des frontières

Aucune condition limite n’est modifiée par cette correction. Les conditions de pression/vitesse, la jauge de pression, les flux conservatifs et les gates de continuité gardent leur comportement existant. Un cas avec frontières mal posées ou bilan de masse incorrect doit continuer à échouer.

## 5. Pièges connus

- Si (|r_u^{(0)}|) est au niveau du bruit machine, (eta_0) peut devenir très grand malgré un résidu RHS-normalisé négligeable. Le diagnostic N12 a signalé un rapport d’environ 27,8 dans ce régime.
- Ce rapport ne doit pas être utilisé comme second gate indépendant avec la même tolérance lorsque son dénominateur est numériquement dégénéré.
- Dans `COUPLED`, le résidu de référence est capturé après la résolution monolithique ; le ratio peut alors rester proche de 1 au plancher FP64 sans signaler une stagnation physique. Il est conservé comme diagnostic mais retiré du monitor non linéaire ; les résidus physiques normalisés par le RHS, la continuité et les variations de champs continuent de décider de l’acceptation.
- La correction ne doit pas contourner le gate RHS-normalisé, les gates pression/vitesse/continuité, ni un oracle physique.
- Une convergence N12 après ce changement doit être confirmée par les métriques de résidu, les champs MMS, la continuité et les campagnes complètes sur le HEAD exact.

## 6. Ordre de précision attendu et validation exigée

Cette modification ne change pas l’ordre spatial ni la discrétisation : elle corrige uniquement la normalisation utilisée par le gate de convergence. Les tolérances et budgets restent strictement identiques.

Validation requise :
- tests unitaires et suite Release ;
- build et tests DebugSanitizers ;
- N12 coupled pressure-velocity MMS, incluant le cas initial proche de la solution ;
- N9 qualification physique et campagne ;
- `validation-total` complet, rapports et artefacts examinés sur le SHA exact ;
- aucun nouvel échec attribuable à la modification et aucun échec masqué.

## 7. Statut d’implémentation

- Fichier source : `src/cfdx/physics/steady_incompressible_solver.h`.
- Changements : le RHS-normalisé reste toujours un gate. Le ratio de première itération reste un critère du monitor pour les algorithmes ségrégés lorsque la référence est bien conditionnée ; il est exclu lorsque la référence est au bruit machine et pour l’algorithme `COUPLED`, dont le résidu de référence est post-solve.
- Implémenté : candidat dans la PR #956.
- Vérifié : CI en cours ; pas encore qualifié physiquement.
- Qualification : **PENDING** jusqu’au verdict exact-HEAD de `validation-total`.
- Décision de merge : interdite tant que les preuves ci-dessus ne sont pas terminées et examinées.
