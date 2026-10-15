# port3b — journal des expériences (2026-10-15)

Outil: diag-tools/port3.ts (twin TS du solveur head 3475826, corrigé des 3 bugs de transcription).
Résultats bruts: diag-evidence/logs/twin_campaign_2026-10-15.json (même commit).

## Expériences

1. **fix D** (suppression du computeFluxes() in-loop): convergé mais solutions parasites.
   Vérification "audit vrais flux": après convergence, recalculer phi depuis (ux,ur) convergés et
   remesurer mom/cont. Échec à ~5e-2 sur Re=100 coarse/medium => le point fixe du fix D
   ne satisfait PAS le système discret. Cd faux: 1.8227 / 4.3345 / 6.0984.
2. **jumeau fidèle C++** (recalcul in-loop conservé): Cd corrects partout
   (0.8075 / 1.1867 / 1.1914 ; Re=50 medium 1.5988), mais mom cale à 5.3e-5 (coarse @1500)
   et cycle 0.7-2.2e-5 (medium/fine) vs gate 2e-6. cont/pcorr passent.
3. **blending beta=0.2** (phi <- (1-b)*phi + b*F(u) à la place du recompute): mom dégradé
   (2.3e-4 @1500), Cd dérive. Abandonné.

## Reproduire avec port3.ts

- fix D: commenter l'appel computeFluxes() in-loop dans solveChunk (expérience D du 14/10).
- audit vrais flux: gradients(p) -> drag() -> copie phi -> computeFluxes() -> momentumResidual()/maxContinuity().
- blending: remplacer le recompute par phi[f] = (1-beta)*phi[f] + beta*F(u)[f].

## Suite décidée (Steven, 2026-10-15)

Prototyper Rhie-Chow dans le twin (flux de face avec couplage pression-vitesse),
valider Cd sur les 4 cas, puis proposer le fix dans la PR #895.
