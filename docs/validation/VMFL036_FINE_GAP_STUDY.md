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

## 4. Rhie–Chow axisymétrique : dérivation et qualification

### Dérivation de l'opérateur

L'équation de quantité de mouvement résolue par le solveur contient la
diagonale relaxée
`aP_relaxed = aP / alpha_u` et le terme source explicite associé à l'ancien
champ. Pour une correction de pression (p'), en gelant les termes hors
diagonale et la source, la réponse de vitesse de cette équation est

[
  mathbf{u}' = -rac{V}{a_{P,\mathrm{relaxed}}}\nabla p'.
]

La même réponse doit être utilisée dans la dérivée du flux de face qui forme
l'équation de continuité, dans la correction conservative du flux et dans la
correction de vitesse. Or, avant cette branche, `compute_fluxes()` utilisait
`V/apx_`/`V/apr_`, tandis que `face_response()` et la correction de vitesse
utilisaient `V/(alpha_u*apx_)`/`V/(alpha_u*apr_)`. Comme `apx_`/`apr_`
stockent déjà (a_P/\alpha_u), cette seconde forme annulait la relaxation une
deuxième fois et donnait un écart de facteur `alpha_u` entre les opérateurs.
Aucun terme de Majumdar explicitement dérivé ne justifiait cette différence.

La branche `fix/vmfl036-relaxed-pressure-response-audit` aligne désormais
les quatre usages sur `V/apx_` et `V/apr_`. C'est une correction algébrique
à qualifier, pas une preuve de correction physique : la campagne de sensibilité
doit confirmer la convergence et mesurer séparément (C_d), (C_{d,p}) et
(C_{d,v}).

### Preuve requise

Le test long `test_vmfl036_alpha_sensitivity` exécute le même cas Re=100,
maillage fine 48×72, mêmes budgets et mêmes gates pour
`alpha_u = 0.3, 0.5, 0.65, 0.7, 1.0`. Il conserve les résidus, le nombre
d'itérations et les trois coefficients de traînée dans
`validation-report/vmfl036-alpha-sweep.json`. Il échoue si un point ne
converge pas ou produit une traînée non finie. La dispersion est mesurée mais
ne reçoit pas de seuil arbitraire avant examen de la première preuve.

Les propriétés à qualifier restent distinctes :

1. annulation de la correction RC pour un champ de pression linéaire exact ;
2. indépendance de la solution convergée vis-à-vis de `alpha_u` ;
3. fermeture de continuité avec le flux de face autoritaire ;
4. limite de maillage régulier uniforme ;
5. maillage étiré avec distances centre-face très différentes.

**Statut :** formulation modifiée par dérivation ; résultats de la campagne
d'impact encore attendus. La référence (C_d=1.0895), les tolérances de
convergence et la bande actuelle de 10 % ne sont pas modifiées.

## 5. Protocole d'étude wall-refinement (à exécuter au jumeau port4)

Matrice (contrôles campagne : maxIt 1500, sweeps 8/400, α_u 0.65, α_p 0.25,
gates inchangés) :

- Re=100, nt=72 : nr = 72 (8,5 %), 96 (5,3 %), 128 (3,3 %) ;
- Re=50, nt=54 : nr = 88 (6,1 %).

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
4. Exécuter le protocole wall-refinement (§5) et dériver la famille de maillages
   adaptée à la paroi.
5. Corriger la collecte d'evidence CI : le job « Validation report » cherche
   `long-evidence/build-validation/long-validation.log` alors que le journal
   CTest de l'artefact est stocké sous une autre structure de répertoires —
   correctif d'infrastructure séparé, sans lien avec le schéma numérique.

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
