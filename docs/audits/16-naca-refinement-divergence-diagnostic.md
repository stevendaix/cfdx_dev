# FICHE MÉTHODE 16 — Diagnostic de convergence non linéaire et divergence NACA

**En-tête d'évidence**
- Scope : protocole de diagnostic exact-HEAD des trajectoires non linéaires NACA0012, qualification N64/N128/N256 et attribution des échecs d'audit de cisaillement pariétal.
- Statut : DIAGNOSTIC (protocole uniquement ; aucune qualification numérique ni correction du solveur revendiquée).
- Issue liée : #986 ; PR de documentation #987.
- Limitation : les observations historiques ne remplacent pas une exécution sur le HEAD courant. Le protocole n'autorise aucun changement de relaxation, d'oracle, de tolérance ou de budget.

## 1. Références sources
- Dembo, R. S., Eisenstat, S. C., Steihaug, T. (1982), “Inexact Newton Methods”, SIAM Journal on Numerical Analysis, 19(2), 400–408. DOI: 10.1137/0719025.
- Ferziger, J. H., Perić, M., Street, R. L. (2020), *Computational Methods for Fluid Dynamics*, 4th ed., Springer. ISBN 978-3-662-62183-2.
- Les références définissent le contexte des méthodes non linéaires et des critères de convergence ; elles ne valident pas à elles seules les résultats CFDX.

## 2. Formulation continue
Pour une discrétisation non linéaire (F(U)=0), une trajectoire d'itération doit être évaluée à partir des résidus de l'équation discrète, de la continuité/du bilan de flux et de la variation de l'état. Un résidu relatif isolé ne suffit pas à démontrer la cohérence du correcteur de pression ni la validité physique de la solution.

La qualification de NACA0012 doit conserver les équations, conditions aux limites, discrétisation, initialisation et oracles existants. Le présent protocole observe la trajectoire ; il ne modifie pas la solution calculée.

## 3. Discrétisation Volume Fini
Exécuter la voie de production SIMPLE existante sur les maillages N64, N128 et N256, indépendamment lorsque le harnais le permet. Conserver les fichiers maillage, paramètres, révisions de code, journaux et artefacts.

À chaque itération, recueillir :
- indice d'itération et nombre de cellules ;
- résidus momentum RHS-relative et reference-relative ;
- résidu linéaire pression, statut du solveur et nombre d'itérations ;
- continuité L1, L-infinity et résidu normalisé ;
- variations L-infinity de vitesse et pression ;
- valeurs effectives de `alpha_u` et `alpha_p` ;
- vérifications finite-value et première cellule/face fautive, si disponibles.

La trace actuelle à fréquence 200 est trop espacée pour localiser un début de divergence autour de l'itération 11. Pour le diagnostic seulement, utiliser l'historique non linéaire existant ou capturer chaque itération sans modifier l'algorithme de production.

## 4. Traitement des frontières
Vérifier que la correction de pression est cohérente avec le flux de face réellement utilisé et avec les conditions de pression fixée/jauge. Évaluer indépendamment l'identité de continuité/flux avant et après correction de pression.

Contrôler les orientations des faces, les normales et les coefficients géométriques locaux. Toute incohérence de flux ou de frontière doit être investiguée avant d'attribuer l'échec au facteur de relaxation. Aucun traitement de frontière ne doit être changé dans le cadre de cette fiche.

## 5. Pièges connus
- Une fréquence de trace trop faible peut masquer l'itération d'apparition de la divergence.
- Un solveur linéaire qui annonce la convergence ne prouve pas que le système assemblé correspond à l'opérateur de flux appliqué.
- L'échec de l'audit de cisaillement pariétal peut être une conséquence d'un écoulement non qualifié ; ne le classer comme cause indépendante que si sa voie de calcul échoue encore sur une solution convergée.
- Une hausse ponctuelle du résidu n'est pas, à elle seule, une preuve de divergence.
- Aucun seuil, oracle, budget d'itérations, normalisation ou test ne doit être affaibli pour obtenir un résultat vert.

## 6. Ordre de précision attendu et validation exigée
Le protocole ne définit pas de nouvelle tolérance ni d'ordre de précision. Les critères numériques faisant autorité sont ceux des tests existants.

Pour chaque niveau de maillage, présenter séparément le verdict et les métriques. Les contrôles requis avant attribution causale sont :
1. topologie du maillage, volumes cellulaires, aires et orientations de faces ;
2. signes, valeurs finies et contributions convection/diffusion des coefficients momentum aux premières cellules en croissance ;
3. matrice de correction de pression, jauge et traitement des frontières à pression fixée ;
4. identité de flux/continuité avant et après correction ;
5. résidu pression et statut réel du solveur linéaire mis en regard de la croissance du résidu momentum ;
6. comparaison du point d'échec du wall-shear audit au verdict du solveur d'écoulement.

Les tests de classification diagnostique doivent couvrir : décroissance convergente, stagnation plate, croissance monotone, croissance exponentielle, oscillations bornées, pic transitoire suivi d'une récupération, métrique de référence nulle, NaN/Inf et priorité explicite entre divergence et budget maximal. Les seuils doivent être justifiés par les contrôles existants et les traces synthétiques, pas ajustés pour faire réussir NACA.

## 7. Statut d'implémentation CFDX
- Statut actuel : protocole de diagnostic documentaire ; implémentation et exécution A/B restent à faire.
- Cibles de validation : `test_naca0012_qualification`, `test_naca0012_wall_shear_audit`, Release, DebugSanitizers, tests NACA ciblés et `validation-total` sur le HEAD exact.
- Acceptation : N64 conserve ses résultats et gates ; le verdict N128 identifie le début mesuré de divergence sans le confondre avec une stagnation ; N256 est rapporté séparément si exécutable ; l'audit wall-shear est classé comme indépendant ou en cascade sur la base de preuves.
- Changement de production interdit dans ce lot : aucun ajustement de `alpha_u`/`alpha_p`, solveur linéaire, normalisation, tolérance, oracle, budget ou fallback.
- Étude séparée éventuelle : comparaison de facteurs `alpha_u` fixes seulement après les contrôles A/B et si les preuves isolent effectivement la relaxation comme cause ; elle exige une proposition de diff indépendante et une revue.

