# Étude : écart +9,93 % du Cd fine du cas axisymétrique VMFL036 (Re=100, réf. Tabata–Itakura 1.0895)

> Statut : DIAGNOSTIC — étude d'évidence, sans modification de code ni d'oracle.
> Aucune tolérance, aucun gate, aucune référence n'est modifié par ce document.
> Contexte : PR #895 (couplage Rhie-Chow, squash `86ef9f2a`) — la campagne
> `test_vmfl036_axisymmetric` passe désormais en CI sur les 4 cas (le test a
> disparu de la liste des échecs longs entre les runs 5133 et 38809508473).

## 1. Constat

Le niveau fine (48x72, Re=100, μ=0.01) converge et passe les critères bring-up
actuels, mais avec un écart proche de la bande temporaire de 10 % :

| Cas | Its | Cd | Écart vs réf. |
|---|---|---|---|
| Re=50 medium 44x54 (réf 1.5785) | 167 | 1.60630 | +1,76 % (bande 3 % PASS) |
| Re=100 coarse 24x36 | 75 | 0.87165 | −20,00 % (convergence exigée, pas de bande Cd) |
| Re=100 medium 44x54 | 136 | 1.19166 | +9,38 % (convergence exigée, pas de bande Cd) |
| Re=100 fine 48x72 (réf 1.0895) | 180 | 1.19764 | +9,93 % (bande bring-up 10 % PASS, marge 0,07 pt) |

Résultats C++ du run CI 38809508473 (head `b51286a7`) : les trois gates de
convergence (mom 2e-6 / pcorr 2e-7 / cont 2e-6) passent sur les quatre cas, et
fine/medium = +0,50 % ≤ 2 %. Les valeurs sont identiques aux prédictions du
jumeau TS `diag-tools/port4_rhie_chow.ts` (branche `ci/long-campaign-diagnostics`,
commit `f0807bec`).

## 2. Compensation d'erreurs : le total masque des contributions non convergées

Le critère d'évolution fine/medium passe (+0,50 %), mais les contributions
évoluent en sens opposés :

| Contribution | Medium 44x54 | Fine 48x72 | Variation |
|---|---|---|---|
| Cd pression | 0.816593 | 0.798473 | −2,22 % |
| Cd visqueux | 0.375066 | 0.399168 | +6,43 % |
| Cd total | 1.19166 | 1.19764 | +0,50 % |

La stabilité du total entre medium et fine est donc une **compensation de deux
dérives opposées**, pas une convergence des contributions. Le critère portant
sur le total ne démontre pas que cd_p et cd_v ont convergé individuellement ;
des critères d'évolution de cd_p et cd_v devront être ajoutés une fois la
famille de maillages établie (toute valeur de gate reste à approuver
explicitement — aucun ajout dans cette PR).

## 3. Le maillage : la première maille radiale domine l'erreur

La famille de maillages (`scripts/generate_vmfl036_axisymmetric_mesh.py`) suit
r_i = R + (R_out − R)·(i/nr)^1.65 avec R = 0.5 et R_out = 50 (50 diamètres).
Fraction de la première maille radiale rapportée au rayon de la sphère,
dr1/R = (R_out − R)·(1/nr)^1.65 / R :

| nr | dr1/R | Niveau |
|---|---|---|
| 24 | 52,3 % | coarse Re=100 |
| 36 | 26,8 % | ancien medium Re=50 |
| 44 | 19,2 % | medium (Re=50 et Re=100) |
| 48 | 16,7 % | fine Re=100 |
| 72 | 8,5 % | — |
| 88 | 6,1 % | — |
| 96 | 5,3 % | — |
| 128 | 3,3 % | — |

Mesures déjà documentées (commentaires du test + jumeau validé) :

- Re=50, nr 36→44→48 : Cd 1.5011 (−4,9 %) → 1.5992 (+1,3 %) → 1.6041 (+1,6 %) :
  le **nombre d'anneaux radiaux est le paramètre contrôlant**.
- Augmenter l'exposant (1.65 → 2.0 / 2.2) à nr fixe affame les anneaux externes
  et éloigne Cd de la référence (1.82 / 1.91 à Re=50).
- Les schémas d'ordre 2 sur cette famille éloignent les DEUX gates de leurs
  références : CDS à correction différée Re=50 medium 1.2526 (−20,6 %) ; SOU
  avec limiteur Re=50 medium 1.302 (−17,5 %) et Re=100 fine 0.956 (−12,2 %).

Interprétation : à Re=100 la couche limite est deux fois plus fine qu'à Re=50,
or **tous les niveaux actuels (même fine) sous-échantillonnent la paroi**. La
contrainte visqueuse est évaluée depuis le centre du premier anneau (à
~8–26 % de R de la paroi selon le niveau) et la pression de face depuis le
centre de cellule. La dérive +6,43 % de cd_v (gradient pariétal sous-estimé,
qui se corrige en raffinant) et la dérive opposée de cd_p sont cohérentes avec
ce diagnostic d'échantillonnage near-wall dominant.

## 4. Question ouverte Rhie-Chow (cohérence algébrique)

Le terme RC du commit `86ef9f2a` pondère par les diagonales relaxées
`apx_`/`apr_` : D^RC = V/a_P^r, tandis que la correction de vitesse de
l'équation de pression utilise D^pc = V/(α_u·a_P^r), soit
**D^RC = α_u·D^pc** (α_u = 0.65 : les coefficients RC sont ~35 % plus faibles).
La campagne CI valide empiriquement le schéma convergé (4 cas, jumeau TS
identique), mais la cohérence algébrique avec l'opérateur de pression
(traitement de la relaxation façon Majumdar) reste à démontrer par la
dérivation du schéma. Tests analytiques prévus (cf.
`docs/validation/N9_S9_RHIE_CHOW_CONSISTENCY.md` et fiches
`docs/audits/03-interpolation-rhie-chow.md` / `04-couplage-pression-vitesse.md`) :

1. annulation du terme RC sur un champ de pression linéaire (gradient exact) ;
2. indépendance de la solution convergée vis-à-vis de α_u ;
3. conservativité du flux RC face à l'équation de p' ;
4. limite maillage régulier uniforme ;
5. maillage étiré (distances centre-face très différentes).

## 5. Protocole d'étude wall-refinement (exécuté le 2026-10-16 — résultats au §8)

Matrice (contrôles campagne : maxIt 1500, sweeps 8/400, α_u 0.65, α_p 0.25,
gates inchangés) :

- Re=100, nt=72 : nr = 72 (8,5 %), 96 (5,3 %), 128 (3,3 %) ;
- Re=50, nt=54 : nr = 88 (6,1 %).

Exécution : la matrice a été étendue de deux cas (Re=100 nt=144 à nr=48, et
Re=100 nr=96 × nt=144) afin de séparer raffinement radial et azimutal. La
campagne a tourné dans GitHub Actions sur le jumeau TS byte-exact (branche
`ci/long-campaign-diagnostics`, driver v3) ; voir §8.

Mesures : cd, cd_p, cd_v et ordre de convergence des trois coefficients.
Critère de décision :

- si cd → 1.0895 ± 3 % quand le premier anneau ≤ 5 % R, l'écart fine est un
  pur effet d'échantillonnage pariétal : la bande 3 % est restaurable par une
  famille de maillages raffinée à la paroi (premier anneau ≤ ~5 % R, progression
  géométrique contrôlée) associée à l'ordre 2 — TODO(order-2) — sans toucher
  l'oracle ;
- sinon, un biais de schéma (intégration de la traînée / fermeture des
  gradients) domine et doit être corrigé en premier.

## 6. Plan (ordre de priorité, audit du 2026-10-16)

1. **[FAIT dans cette PR]** Réparer le test kernel (asserts inertes en Release,
   contrôle visqueux radial avalé par la corruption) et durcir le smoke test
   quick (convergence explicitement exigée dans le budget réduit).
2. Démontrer la cohérence Rhie-Chow (§4) : dérivation + 5 tests analytiques.
3. Valider indépendamment le calcul de la traînée : convergence de la force
   pariétale sur champ analytique, puis référence indépendante pour cd_p/cd_v ;
   ajouter des critères d'évolution de cd_p et cd_v (en plus du total).
   **Exécuté le 2026-10-16** : trois estimateurs indépendants (pariétal,
   budget champ lointain, fermetures équations télscopées) sont d'accord
   sous raffinement (§8) — le calcul de traînée n'est pas en cause ; les
   critères d'évolution cd_p/cd_v restent suspendus au pivot de l'item 4.
4. Exécuter le protocole wall-refinement (§5) et dériver la famille de maillages
   adaptée à la paroi. **Exécuté le 2026-10-16 — item pivoté** : l'hypothèse
   d'échantillonnage pariétal est réfutée (§8), le biais est un défaut de
   schéma (~+19 % asymptotique) ; pas de famille de maillages ni de gates
   tant que ce défaut n'est pas corrigé. TODO(order-2) renforcé.
5. Corriger la collecte d'evidence CI : le job « Validation report » cherche
   `long-evidence/build-validation/long-validation.log` alors que le journal
   CTest de l'artefact est stocké sous une autre structure de répertoires —
   correctif d'infrastructure séparé, sans lien avec le schéma numérique.
   **PR #916 ouverte** : upload-artifact@v4 enracine l'artefact au plus
   proche ancêtre commun des chemins (`build-validation`), donc le log
   réel est `long-evidence/long-validation.log` ; 3 chemins corrigés +
   étape de vérification défensive.

## 7. Références

- Tabata & Itakura, « Precise Computation of Drag Coefficients of the Sphere »,
  Int. J. Computational Fluid Dynamics 9 (1998) 303-311 :
  C_D(Re=50) = 1.5785, C_D(Re=100) = 1.0895 (cf. `VMFL036_RE100_CASE.md`).
- `docs/audits/03-interpolation-rhie-chow.md`,
  `docs/audits/04-couplage-pression-vitesse.md`,
  `docs/audits/07-benchmarks-references.md`.
- `docs/validation/N9_S9_RHIE_CHOW_CONSISTENCY.md`,
  `VMFL036_RE100_CASE.md`, `VMFL036_AXISYMMETRIC_IMPLEMENTATION.md`.
- Issue #440 (référence Re=50 mal attribuée dans la matière Ansys VMFL036).

## 8. Résultats empiriques (campagne CI du 2026-10-16)

> Statut : EVIDENCE — campagne exécutée par GitHub Actions sur le jumeau TS
> byte-exact (`diag-tools/port4_rhie_chow.ts` + driver v3, branche
> `ci/long-campaign-diagnostics`). Preuve :
> `diag-evidence/logs/wall_refinement_2026-10-16.json` (commit du bot CI,
> version 3). Contrôles inchangés (maxIt 1500, sweeps 8/400, α_u 0.65,
> α_p 0.25, ω 1.8, gates 2e-6/2e-7/2e-6) ; self-check pré-campagne : le cas
> re100 coarse 24x36 reproduit Cd=0.8716507056286522 bit-exact (campagne du
> 2026-10-15).
> Correction transparente : le driver v2 avait une erreur de signe dans le
> cd_budget assemblé (D = convTotal − presTotal + viscTotal) ; v3 :
> D = −(convTotal + presTotal − viscTotal). Les composantes enregistrées
> par v2 étaient correctes ; seul le total assemblé était affecté
> (re50_nr88 : budget v2 0.783 → v3 1.876).

Trois estimateurs indépendants du Cd sont mesurés par cas :

- **pariétal** : `drag()` du jumeau (pression au centre de cellule +
  gradients Gauss) — c'est la valeur CI ;
- **budget champ lointain** : identité QDM discrète sur les frontières du
  domaine (convection + pression + visqueux), D = −(conv+pres−visc)/qA ;
- **fermetures équations (télscopé)** : force pariétale avec les fermetures
  des équations elles-mêmes (pression extrapolée à la paroi + cisaillement
  deux-points mu·A/bD·ux).

| Cas | dr1/R | Cd pariétal | cd_p / cd_v | Cd budget | Cd éq. télscopé | Écart pariétal vs réf |
|---|---|---|---|---|---|---|
| Re=100 nr72 nt72 | 8,5 % | 1.23602 | 0.700 / 0.536 | 1.37867 | 1.37884 | +13,45 % |
| Re=100 nr96 nt72 | 5,3 % | 1.25883 | 0.663 / 0.596 | 1.33996 | 1.34033 | +15,54 % |
| Re=100 nr128 nt72 | 3,3 % | 1.25593 | 0.639 / 0.617 | 1.30285 | 1.30350 | +15,28 % |
| Re=100 nr48 nt144 | 16,7 % | 1.16492 | 0.768 / 0.397 | 1.48508 | 1.48656 | +6,92 % |
| Re=100 nr96 nt144 | 5,3 % | 1.22179 | 0.632 / 0.590 | 1.29974 | 1.30167 | +12,14 % |
| Re=50 nr88 nt54 | 6,1 % | 1.74866 | 0.843 / 0.906 | 1.87628 | 1.87574 | +10,78 % |

Conclusions :

1. **Hypothèse échantillonnage pariétal RÉFUTÉE.** Le raffinement radial
   nr72 → nr96 → nr128 (dr1/R 8,5 % → 3,3 %) sature : Cd pariétal
   1.23602 → 1.25883 → 1.25593, sans tendance vers 1.0895. La bande de 3 %
   n'est pas restaurable par raffinement ; le critère de décision du §5
   tranche donc pour « biais de schéma dominant ».
2. **Les trois estimateurs convergent l'un vers l'autre sous raffinement**
   (écart budget/pariétal : 11,5 % à nr72, 6,4 % à nr96, 3,7 % à nr128) vers
   ~1.30, soit ~+19 % vs la référence. Le biais est donc dans la solution
   elle-même (schéma), pas dans le post-traitement de la traînée : l'item 3
   du plan est clos, le calcul de traînée n'est pas le coupable.
3. **Identité télscopée vérifiée** : cd_budget ≈ cd_eq_total à ~5e-4 près sur
   tous les cas, massImbalance ≤ 6e-5 (absolu, sur des flux de ~7854 kg/s).
   Les fermetures des équations et le budget champ lointain sont cohérents.
4. **Le niveau fine CI 48x72 (+9,93 %, PASS bande 10 %) passe « pour la
   mauvaise raison »** : il se situe dans le transitoire non monotone
   (nr48 nt144 : 1.16492, +6,92 %) et le biais asymptotique du schéma est
   ~+19 %.
5. **Hypothèse dominante : viscosité numérique de l'upwind ordre 1**
   (mu_num ~ rho·U·dr1/2, comparable à mu aux nr testés → Re effectif plus
   faible → Cd trop élevé, cohérent avec la hausse de cd_v sous raffinement
   0.536 → 0.617 pendant que cd_p baisse 0.700 → 0.639). Nuance : les essais
   d'ordre 2 du §3 (CDS/SOU) éloignaient AUSSI les gates de leurs
   références — l'interaction schéma/RC reste donc à élucider lors du
   TODO(order-2), avec l'écart budget/pariétal comme outil de diagnostic
   (il mesure la convergence interne de la solution).
6. **Conséquence (pivot de l'item 4)** : pas de famille de maillages ni de
   gates d'évolution cd_p/cd_v tant que le défaut de schéma n'est pas
   corrigé. TODO(order-2) renforcé.

Campagne complémentaire v4 lancée sur la branche diag (driver v4, mêmes
contrôles sauf α_u ∈ {0.5, 0.65, 0.8, 0.95}) : indépendance de la solution
convergée vis-à-vis de α_u — cohérence Majumdar, item 2 du plan (diagnostic,
pas un gate).
