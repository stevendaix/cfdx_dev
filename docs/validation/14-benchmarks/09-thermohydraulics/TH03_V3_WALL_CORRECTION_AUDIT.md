# TH-03 — Audit numérique et correctifs SST V3 (2026-11-03)

## Verdict indépendant

**NON QUALIFIÉ.** Le maillage HDF5 de débogage est volontairement trop grossier pour un SST résolu à la paroi et son profil d'entrée est uniforme. Aucune comparaison expérimentale ni GCI sur trois maillages n'est démontrée ici. La référence Driver–Seegmiller est `x_r/h = 6.26 ± 0.10`, mais les cas de débogage `Re_h=36000` ne sont pas identiques au contrat TH-03 `Re_h=37500` : ne pas utiliser une valeur expérimentale comme oracle sans correspondance de conditions d'essai.

Sources : https://tmbwg.github.io/turbmodels/backstep_val.html et https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/.

## Deux défauts avérés dans le code

1. **Flux visqueux à la paroi SST** : le prédicteur momentum, la correction de contrainte transposée et le contrôle indépendant utilisaient `mu_eff(owner) = mu + mu_t`, y compris sur une face de paroi non glissante. Or la viscosité turbulente à une paroi résolue est nulle, et l'export `Cf` utilisait déjà `mu` seul. V3 introduit une viscosité explicite par face frontière pour les opérateurs de diffusion et la reconstruction des flux, avec `mu_wall = mu` sur les parois fixes SST. Le contrôle de conservation, le résidu final et le prédicteur reçoivent le **même coefficient**. Les flux intérieurs conservent leur interpolation en viscosité effective. Le test C++ vérifie la matrice et le flux intégrés.
2. **Sous-relaxation HDF5 ignorée** : `numerics.under_relaxation_momentum` et `numerics.under_relaxation_pressure` étaient chargés mais n'étaient jamais passés à `controls.coupling`. V3 exige chacun dans `(0,1]` et affecte ces valeurs au solveur, avec leur trace imprimée. La tolérance de convergence n'est pas modifiée.

## Diagnostic des résultats réels (maillage debug de 192 cellules)

| Mesure | Ancien code | Opérateur corrigé |
|---|---:|---:|
| Itérations maximales | 250 | 250 |
| Convergence déclarée | oui, 145 itérations | **non**, 250 itérations |
| Résidu momentum final relatif (RHS) | 5.65e-12 | ~3.02e-2 |
| Résidu de continuité final | 8.73e-13 | ~8.93e-3 |
| Réattachement diagnostique premier ordre | `x_r/h=1.77` | **non déterminé** : l'état n'a pas convergé |
| Résolution de paroi estimée (avant correction) | `max(y+)≈54.63` | maillage inchangé, non qualifié |
| `max(nu_t/nu)` cellule de paroi (avant correction) | ~491.42 | nouvelle itération non qualifiée |

**Interprétation :** le précédent `PASS` numérique ne constitue pas une précision physique. Modifier le contrat visqueux est nécessaire, mais ce changement déstabilise le cas de débogage sous-résolu. Il est interdit de présenter la valeur `x_r/h` de l'état instationnaire ou non convergé comme validation. La stabilité et la fidélité physique restent des verrous séparés.

3. **Raison d'arrêt trompeuse** : quand le moniteur atteint `MAX_ITERATIONS`, le message précédent annonçait que les critères porteurs étaient satisfaits même si le dernier résidu momentum n'était pas faible. V3 distingue explicitement `carrier CONVERGED / SST non convergé` de `max itérations / transport et/ou carrier non convergés`.

## Deuxième essai de stabilité sur le même maillage

`alpha_u=0.25`, `alpha_p=0.10`, 400 itérations, **tolérances inchangées** : non convergé, résidu momentum relatif RHS `2.08e-2`, continuité normalisée `3.05e-4`, déséquilibre massique global relatif `3.75e-7` mais **max local relatif `3.82e-2`**, et `max(y+)≈79.69`. Une bonne conservation globale ne peut pas masquer un déséquilibre local. Cet état est une preuve négative conservée, pas un nouveau résultat de reattachment.

## Gates ajoutées

- Export CSV natif `wall_y_plus_estimate` et `owner_nu_t_over_nu` pour chaque face non glissante ; `wall_shear_method=first_order_one_sided_diagnostic` indique explicitement un calcul **non certifié** de `Cf`.
- `python scripts/audit_th03_wall_resolution.py --face-csv .../face_fluxes.csv --cell-csv .../cell_fields.csv --patch downstream_wall --density 1.2 --kinematic-viscosity 1.5592777777777778e-5 --output .../wall_resolution.json` bloque si `max(y+) > 1` ; tous les patchs demandés doivent être présents.
- En campagne native TH-03, toutes les parois fixes de la description HDF5 sont vérifiées ; le rapport est sauvé dans `wall_resolution.json`. Le validateur strict revérifie les valeurs contre les CSV signés SHA-256. Un `PREFLIGHT_PASS` reste **insuffisant** pour une qualification physique.
- Le générateur de maillage debug propose `--wall-stretch 3.6` pour améliorer la résolution géométrique aux parois sans changer la hauteur de marche. Il ne crée pas le profil turbulent expérimental d'entrée. La résolution `y+` doit toujours être recalculée après simulation.

## Pour qualifier TH-03

Il reste à : importer la géométrie et la couche limite mesurée du vrai cas, produire des grilles suffisamment fines près des parois et dans la zone de séparation, obtenir convergence et conservation indépendantes sur trois niveaux, vérifier `Cf` par reconstruction d'ordre 2 et comparer les profils `Cp`, `Cf`, `U/Uref` et les contraintes de Reynolds avec les données originales assorties d'incertitudes. Le SST expérimental doit être amélioré si ces contrôles échouent. Ne jamais modifier une tolérance juste pour obtenir `PASS`.

## Statut des autres benchmarks

TH-01 et TH-08 demeurent bloqués par l'absence d'une chaîne instationnaire DNS/LES/URANS de production et de statistiques/spectres. TH-07 dispose d'un chemin expérimental SST–énergie–flottabilité, sans validation sur gradients thermiques réels. Cette V3 cible les deux erreurs TH-03 observables et n'invente aucun résultat de qualification sur TH-01/07/08.
