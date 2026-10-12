<!-- PR noyau numérique : la checklist ci-dessous bloque le review (constitution technique, Règles d'Or 1–4).
     Les cases [ ] doivent être cochées dans la description. Un modérateur ne merge pas une case manquante. -->

## Fiche(s) Méthode (Règle d'Or 1)

- [ ] La ou les méthodes numériques touchées ont une fiche `docs/audits/NN-*.md’ approuvée (lien : …).
- [ ] La discrétisation implémente **exactement** la fiche (toute divergence est documentée dans la fiche, section 7).

## A. Mathématique & discrétisation

- [ ] Terme de surface et de volume cohérents et **conservatifs**.
- [ ] Couplage pression-vitesse (Rhie–Chow + correction Majumdar) correct pour le colocalisé — si applicable.
- [ ] Jauge de pression / null-space : traitement explicite et testé (tout-Neumann, RHS incompatible) — si applicable.

## B. Physique & conservation (Règle d'Or 3)

- [ ] Bilan de masse : $\\sum \\dot m_{in} - \\sum \\dot m_{out} < 10^{-10}$ (ou tolérance déclarée par le contrat de cas).
- [ ] Bilans QDM/énergie : forces/sources = variation des champs (reconstruction **indépendante**).
- [ ] Bornitude : champs dans leurs bornes physiques **sans oscillation**, détection **avant** clipping + télémetrie (floors/clipping).

## C. Numérique & performance (Règle d'Or 2)

- [ ] MMS : pente de convergence théorique reproduite en log–log sur ≥ 3 niveaux.
- [ ] Benchmark de référence : superposition graphique (analytique / Ghia / OpenFOAM / SU2) jointe à la PR.
- [ ] GCI : $r$, $p$ observé, $F_s$, erreur fine, extrapolation — si affirmation de convergence.
- [ ] Parité CPU vs MPI vs GPU sur le même maillage à 1e-12 — si applicable (Règle d'Or 4). Aucun fallback silencieux.

## D. Hygiène de campagne

- [ ] Aucune tolérance élargie, aucun test désactivé, aucune campagne raccourcie pour rendre le CI vert (`docs/development-and-validation-rules.md’ §1–§5).
- [ ] Registres mis à jour dans la même PR : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json’ **et** `QUALIFICATION_REGISTRY.md’.
- [ ] Issue de roadmap (#34 / #461 / #481 / #601…) synchronisée après merge.

## Résumé

- Issue/maturité : …
- Ce qui est implémenté / testé / qualifié / hors scope : …
- Limitations connues : …
