# CFDX V67 — diagnostic causal des pas couplés BFS152

**Branche :** à partir de V66, sans modification de SOU/MUSCL, Rhie–Chow,
`rAU`, du SST, des conditions limites, de Schur ou des critères V13/V14.

## Question

V65 échoue à l'itération extérieure 23 du pilote `TH03_debug_only`,
152 cellules. Le pas complet diminue le mérite état-seul mais entraîne
une trop forte croissance du bloc omega ; les petits pas empêchent cette
croissance mais ne satisfont plus la réduction de mérite historique de
0,995. La séparation du facteur de pas carrier (U,p) et du transport
SST (k,omega) mérite donc d'être testée, avant de toucher à un Jacobien
ou au préconditionneur.

## Expérience (lecture seule)

`CFDX_FPHYS_V67_BLOCK_PROBE=1` exige les options expérimentales V51 et
V65, elles-mêmes conditionnées par V47/V13/V14, désactivées par défaut.
Aux reprises 0, 4, 8, 10 de l'itération 23, l'évaluateur reconstitue
la même fonction état-seul `F_phys(U,p,k,omega)` pour les neuf mélanges
indépendants `(fraction_k, fraction_omega)` : `(0,0)`, `(0.25,0.25)`,
`(0.5,0.5)`, `(0.75,0.75)`, `(1,0)`, `(1,0.25)`, `(1,0.5)`,
`(0.5,1)`, `(0,1)`. Les fractions sont appliquées aux seules inconnues
SST, avec le même carrier U,p effectivement tenté pour le retry concerné.

Chaque diagnostic recalcule la fermeture Rhie–Chow et le résidu physique
via l'oracle V51. Il emploie les dénominateurs gelés et les gardes V65
**sans les modifier** : croissance <=1.20 par bloc, ratio de mérite <=0.995.

**Avertissement fondamental :** ces états sont *contre-factuels*,
ils n'ont PAS été produits par le transport SST intérieur. Même un
mélange qui satisferait les deux gardes ne constituerait pas un état
physique accepté ou un certificat de convergence ; il faudrait encore
raccorder une mise à jour cohérente de `nu_t`, des parois et des flux,
avec un test de non-régression complet.

La sonde enregistre également `1 - 0.005 * beta` comme référence
indicative de suffisance de décroissance tenant compte du pas, mais
ce n'est PAS un véritable test Armijo du solveur CFDX : `beta` est
appliqué à U,p puis le SST est recalculé, et non à l'ensemble de la
direction à six composantes. Cette référence **ne gouverne aucune
acceptation**.

## Contrat de sécurité

- Aucun état V67 n'est engagé ni retenu dans le solveur.
- Toute violation de la positivité, de la fermeture de flux, de la forme
  des champs ou des critères de métrique déclenche un refus.
- Les traces V65 doivent rester textuellement identiques aux témoins V66.
- Les cas absents ou les faux certificats sont rejetés par
  `scripts/audit_bfs_v67_block_probe.py`.
- Les validations BFS712, BFS20 540, NASA, parois SST actives et CI
  complète restent ouvertes.

## Seconde expérience : direction globale à six champs

Les fractions `alpha=2^{-j}`, `j=0..10`, sont appliquées **à tous les
six champs** entre `X_old` et le candidat complet obtenu au retry 0 :
`X(alpha) = X_old + alpha * (X_trial_full - X_old)`. Contrairement à V47,
ce chemin est bien un segment dans l'espace des états à six composantes.
`F_phys` est recomputé de manière indépendante pour chaque alpha.
Les résultats sont comparés aux mêmes gardes stricts V65 (1.20, 0.995).
Cela n'est pas non plus une résolution SST intérieure ni une autorisation
d'acceptation de la transaction : aucune valeur reconstituée de nu_t ou
flux porteur n'est publiée. Il faut ensuite qualifier l'intégration de
la mise à jour dans la transaction, et comparer contre les 40 itérations
du COUPLED historique avant toute activation.

Les deux familles sont distinctes et publiées sous
`BFS_BLOCK_COUNTERFACTUAL_V67` et `BFS_FULL_SIX_COUNTERFACTUAL_V67`.
L'auditeur exige respectivement 36 et 11 points, rejette un point manquant,
un faux commit et la modification des seuils, et s'assure que le point
`alpha=1` reproduit bien le candidat V65 du retry 0.

## Résultats réexécutés (BFS152)

- La sonde de bloc a produit **36** contre-factuels (9 mélanges sur
  4 reprises) : **aucun** ne satisfait simultanément les gardes V13/V14.
  La plupart des mélanges partiels créent des résidus SST beaucoup plus
  grands car ils ne résolvent plus les équations de transport.
- La ligne complète `6N` a produit **11** contre-factuels : **aucun** ne
  satisfait le garde SST `<=1.20`. Même le segment `alpha=1/1024` donne
  un rapport de mérite `0.9984439`, mais une croissance maximale `6.7737`.
- Les véritables essais V65 conservent 22 commits et un rejet à
  l'itération 23. Les **149 traces V65**, les **22 commits** et les
  **22 certificats V13** sont textuellement identiques aux journaux V66.
- Une référence indicative `1-0.005*beta` laisserait seulement les
  reprises 9 et 10 passer les deux tests *indicatifs* ; elle ne remplace
  pas l'acceptation V65 car le carrier et SST ne suivent pas la même
  direction paramétrée.

## Priorité physique suivant le résultat négatif

Le SST de production résout son transport avec le flux **porteur COUPLED**,
mais `F_phys` évalue le transport avec le flux Rhie–Chow réassemblé en
état seul. Le diagnostic V66 a déjà montré un écart de flux relatif
`8.678e-4` et un écart `rAU` proche de `49.35 %` au retry=0 de
l'itération 23. Ce décalage n'est pas à lui seul un bug, puisque le
système Picard réutilise des coefficients retardés. L'expérience
prioritaire suivante devra réassembler les DEUX équations SST à état
strictement identique en ne substituant que `phi_carrier` par
`phi_Fphys`, en gardant forces, gradient, nu_t, lois de paroi et
SST sources identiques. Réfuter une hypothèse de flux exigera de
montrer que cette substitution ne représente qu'une fraction négligeable
du résidu omega observé ; toute correction restera ensuite à valider
contre le COUPLED historique sur 40 itérations et sur les autres cas.