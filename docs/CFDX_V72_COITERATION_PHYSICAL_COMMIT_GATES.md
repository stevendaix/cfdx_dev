@@ -0,0 +1,82 @@
# CFDX V72 — fermeture itérative privée COUPLED–SST et certificats de non-publication

## Situation de départ

Sur le pilote `TH03_debug_only/BFS152` de V71, la co-itération privée
`gamma=1/32`, fraction de pression `1/4` a réduit le mérite V65 à 0,92835
et le résidu brut `omega` à 7,07e-4. Sa divergence physique après
réassemblage `F_phys`, 6,25e-4 kg/s, reste supérieure au seuil de
conservation de 1e-10 kg/s. La résolution d'un sous-problème de pression
avec `rAU` gelé n'est pas un certificat de conservation du système complet.

## Hypothèse V72 — ni nouveau préconditionneur, ni nouveau Jacobien

La boucle V72 répète, **sur copies privées uniquement**, la co-itération V71 :

1. Reconstituer `F_phys(X)` sur l'état `X=(Ux,Uy,Uz,p*,k,omega)` courant.
2. Appliquer la prédiction diagonal-Picard du momentum avec **son** `rAU`.
3. Reconstruire la réponse pression–vitesse et le flux Rhie–Chow avec le
   même candidat, sans remplacer arbitrairement le flux porteur historique.
4. Résoudre le véritable transport SST à flux candidat, sur une copie de
   l'état et avec les fonctions de paroi *du pilote utilisé*.
5. Recalculer `F_phys`, incluant viscosité turbulente, flux, diagonale,
   convection, diffusion et termes sources.
6. Réévaluer indépendamment une seconde fois **exactement le même état**.
   Le rejeu doit être identique en arithmétique machine.
7. Vérifier, par sommation owner/neighbour de toutes les faces, la
   divergence réassemblée et l'identité avec le résidu de continuité signé.
8. Comparer au **même** état ancien et à ses échelles figées selon V65.
   Un cycle privé n'est retenu que si les deux gardes V65 **inchangés**
   sont satisfaits ET si la divergence **réassemblée** diminue.

À chaque cycle, un petit jeu borné de fractions est essayé depuis le
**même état privé**. Le meilleur candidat admissible est conservé dans la
mémoire privée ; la boucle s'arrête en cas d'absence de progrès, de rejet
par les gardes ou lorsque le budget de cycles est atteint.

La version initiale V72 est limitée à 8 cycles sur BFS152, avec les choix
`(gamma, fraction pression)` :
`(1/32,1/4)`, `(1/32,1/2)`, `(0,1/2)`, `(0,1/4)`, `(1/32,1)`.
Ce sont des **sondes de co-itération existantes**, pas de nouveaux schémas.
L'option `CFDX_FPHYS_V72_FIXED_POINT=1` requiert les diagnostics V51/V65,
V68/V69 et V71 ; elle reste désactivée par défaut.

## Certificats réémis sur l'état final du cycle

Pour chaque cycle privé retenu :

- Les six `signed_raw` et les six `normalized` de **l'état final**, pas les
  certificats momentum d'une itération précédente.
- Le résidu physique maximum de momentum, continuité, `k`, `omega` ;
- La divergence maximale des flux Rhie–Chow **réassemblés**, le flux net
  des frontières, l'identité owner/neighbour et l'identité de la divergence
  des faces avec `signed_raw` de la continuité ;
- La fermeture du point fixe Rhie–Chow et le rejeu bitwise de `F_phys` ;
- Le même mérite V65 avec normalisation figée sur l'ancien état ;
- L'indication explicite **non réémis** pour les certificats de solveur
  linéaire momentum et de transport SST **de production**.

Une convergence de `k` et `omega` sur une copie ne remplace pas un nouveau
certificat autoritatif de leur transport, ni un vrai nouveau solveur linéaire
momentum. Le prédicteur V71 est diagonal et ne produit pas ce certificat.
**Aucun commit n'est autorisé par V72**, même si une sonde ferme toute la masse :
une future transaction de production devra d'abord obtenir les deux
certificats natifs supplémentaires à l'état final, vérifier les résidus
physiques indépendants puis appliquer le commit atomique avec rollback.

## Critères non modifiés

- Garde V13 par bloc : croissance <= 1,20.
- Garde V14 de mérite V65 : rapport <= 0,995.
- Fermeture Rhie–Chow et conservation absolue de la sonde : 1e-10.
- Seuils physiques de convergence du cas : **inchangés**.
- Aucun choix PCD/LSC/BFBt/BlockLocal/SIMPLE/SIMPLEC/selfp modifié.
- Aucun paramètre SST ou de paroi modifié.

## Limites de qualification

La boucle reste une expérience de diagnostic sur le pilote
`TH03_debug_only/BFS152` avec traitement de paroi historique. Il ne s'agit
ni du BFS Driver–Seegmiller physiquement validé, ni d'une preuve sur
BFS712/BFS20540 ou les fonctions de paroi actives. Les certificats partiels
sont explicitement distincts d'une autorisation de publication d'état.
