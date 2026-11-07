# NACA0012 — campagne forensique de recherche de la stagnation momentum

**Statut :** protocole de diagnostic uniquement. Aucun changement de schéma, coefficient, relaxation, tolérance ou critère d'acceptation. Les sorties ne sont jamais des certificats de convergence.

## 1. Question à résoudre

Sur N64, la continuation atteint 100 % de TVD à l'itération 120, puis la solution stagne et termine non convergée à l'itération 330. Valeurs observées :

- Momentum RHS relatif : 7.25783e-5
- Maximum des composantes momentum : 5.374e-3
- Résidu momentum indépendant en faces, composante y : 2.06470e-3
- Continuité normalisée : 2.09097e-7
- Résidu linéaire pression : 8.43316e-12

Le résidu linéaire de pression ne permet pas d'attribuer la stagnation à la pression, ni de conclure que le système non linéaire est convergé. Il faut identifier quelle contribution physique, quelle région et quelle composante empêchent la diminution du résidu momentum.

## 2. Instrumentation déjà disponible après le merge #1187

Ne pas réimplémenter ces mécanismes avant d'avoir exploité leurs sorties :

| Instrument | Activation / emplacement | Donnée |
|---|---|---|
| Historique physique | CFDX_NACA_PHYSICAL_HISTORY_CSV | Momentum x/y/z, continuité, résidus de faces par axe, pression linéaire, variations U/p, alpha, pseudo-CFL |
| Budget momentum cellwise | CFDX_NACA_MOMENTUM_BUDGET_CSV | Convection, viscosité, contrainte transpose, pression, force, divergence bornée, divergence masse, écart de flux collocated, résidu local |
| Budget toutes cellules | CFDX_NACA_MOMENTUM_BUDGET_ALL_CELLS=1 | Toutes les cellules au lieu des 12 pires par composante |
| Trace des champs avant assemblage | CFDX_NACA_FORENSICS=1 | Extrema de p, grad(p), U, flux, divergence, viscosité effective et non-finis |
| Trace de l'équation momentum | CFDX_NACA_FORENSICS=1 | RHS, diagonale, hors-diagonale, sommes de lignes, cellules extrêmes, forces pression/corps |
| Analyse historique | scripts/naca_convergence_diagnosis.py | Tendances et critères non satisfaits ; jamais un certificat |
| Analyse budget | scripts/naca_momentum_budget_diagnosis.py | Classement cellwise des contributions ; jamais un certificat |

Les callbacks de budget sont opt-in et observateurs : ils ne doivent pas modifier l'état, la matrice, les tolérances ou la décision d'acceptation.

## 3. Exécution reproductible

Construire le binaire de validation habituel du dépôt : test_naca0012_qualification, compilé depuis tests/validation/test_naca0012_geometry.cpp. Utiliser exactement le même fichier de maillage N64 pour tous les essais.

~~~bash
mkdir -p artifacts/naca-forensics/upwind
CASE=artifacts/naca-forensics/upwind
NACA_N64_MESH=/path/to/naca0012_n64.h5
env CFDX_NACA_PHYSICAL_HISTORY_CSV="$CASE/history.csv" \
    CFDX_NACA_MOMENTUM_BUDGET_CSV="$CASE/momentum_budget.csv" \
    CFDX_NACA_MOMENTUM_BUDGET_EVERY=10 \
    CFDX_NACA_MOMENTUM_BUDGET_ALL_CELLS=1 \
    CFDX_NACA_FORENSICS=1 \
    ./build/test_naca0012_qualification --single "$NACA_N64_MESH" \
    >"$CASE/stdout.log" 2>"$CASE/forensics.log"
~~~

Adapter le chemin du binaire et du maillage au répertoire de build réel. Créer un sous-répertoire distinct pour chaque essai. Le test peut retourner un code non nul si la convergence n'est pas atteinte : conserver les CSV et logs, c'est précisément le résultat diagnostique recherché.

Analyser les sorties sans écraser les sources :

~~~bash
python scripts/naca_convergence_diagnosis.py "$CASE/history.csv" \
  --output "$CASE/history-diagnosis.json"
python scripts/naca_momentum_budget_diagnosis.py "$CASE/momentum_budget.csv" \
  --output "$CASE/budget-diagnosis.json"
~~~

Le mode --single mesh.h5 conserve les paramètres par défaut du test (dont alpha_u=0.6, alpha_p=0.25 et 8000 itérations maximum). Ne pas comparer des historiques issus de maillages, paramètres, initialisations ou budgets différents.

## 4. Matrice A/B — une seule modification à la fois

Tous les essais utilisent N64, Re=1000, alpha=0, la même initialisation uniforme, les mêmes tolérances physiques, le même budget maximal et le même réglage de relaxation, sauf la dimension explicitement testée.

| Essai | Variables additionnelles | Question |
|---|---|---|
| A0 | aucune (UPWIND/SIMPLE) | Référence basse résolution |
| A1 | CFDX_NACA_CONVECTION_TVD=1 + CFDX_NACA_TVD_CONTINUATION=1 | La stagnation est-elle liée à la reconstruction TVD minmod ? |
| A2 | CFDX_NACA_CONVECTION_VANLEER=1 + CFDX_NACA_TVD_CONTINUATION=1 | La variante van Leer change-t-elle le défaut local ? |
| A3 | CFDX_NACA_CONVECTION_VANALBADA=1 + CFDX_NACA_TVD_CONTINUATION=1 | Même question pour van Albada |
| A4 | CFDX_NACA_ALGORITHM=SIMPLEC (avec A1) | SIMPLEC change-t-il la stagnation ? |
| A5 | CFDX_NACA_SKEW_GRADIENT=1 (avec A1) | Le gradient corrigé du skewness est-il déterminant ? |
| A6 | CFDX_NACA_PRESSURE_LINE_SEARCH=1 (avec A1) | Les corrections de pression sur-amplifient-elles le défaut de continuité ? |

Ne pas activer plusieurs variantes de convection simultanément. A4–A6 sont comparés séparément à A1. Le pseudo-temps, alpha_u=0.8 et les essais de prolongement restent hors de cette première matrice : ils ajouteraient des facteurs confondants.

Pour chaque essai, recopier la commande A0 et ajouter uniquement les variables indiquées. Chaque essai doit avoir ses propres chemins de sortie.

## 5. Lecture causale des résultats

Pour les 8 à 12 dernières itérations échantillonnées, produire par composante :

- évolution du maximum local du résidu physique et du résidu de l'opérateur borné ;
- cellules responsables, coordonnées et voisinage topologique ;
- contributions signées convection, viscosité, contrainte transpose, pression et force ;
- écart entre flux conservatif autoritaire et flux obtenu par interpolation collocated ;
- divergence de masse locale et flux de frontière ;
- diagonale, hors-diagonale, somme de ligne et RHS de l'équation assemblée ;
- extrema de U, p, grad(p), mu_eff et flux, avec compte des valeurs non finies.

Hypothèses à départager — ce ne sont pas des conclusions acquises :

1. **Limiteur / reconstruction convective :** le terme convection domine le résidu local et varie entre A1/A2/A3 à état et budget comparables.
2. **Gradient / non-orthogonalité / skewness :** les cellules fautives sont corrélées aux zones de mauvaise qualité géométrique et A5 réduit le défaut sans dégrader les autres composantes.
3. **Couplage pression-vitesse / flux collocated :** l'écart de flux collocated ou la correction pression domine ; A4/A6 améliore simultanément momentum-y et le bilan de masse.
4. **Diffusion visqueuse / contrainte transpose :** le terme visqueux domine et le défaut suit la direction d'anisotropie locale.
5. **Conditionnement momentum :** diagonale faible, ratio hors-diagonale/diagonale élevé, somme de ligne anormale ou RHS local dominant sur les mêmes cellules.
6. **Normalisation :** le résidu relatif y est élevé parce que son RHS de référence est petit ; rapporter normes absolues, relatives et échelles de référence sans remplacer le critère officiel.

Une contribution signée peut être grande et compensée par une autre. Examiner les sommes signées **et** les valeurs absolues/maxima par cellule. L'écart collocated est la différence entre flux de masse corrigé et flux interpolé directement ; ce n'est pas une mesure pure de Rhie–Chow.

## 6. Livrables avant toute correction du solveur

- Un dossier d'artefacts par essai : commande exacte, environnement, hash du binaire, hash du maillage, stdout, trace forensique, historique physique, budget cellwise et JSON d'analyse.
- Un tableau comparatif des dernières itérations acceptées, sans mélanger les normes ou leurs dénominateurs.
- Une carte des cellules fautives et leur qualité géométrique : non-orthogonalité, skewness et rapport d'aspect.
- Une conclusion causale avec au moins une observation discriminante et une contre-épreuve.
- Une proposition de correction isolée avec un test de régression reproduisant le défaut.

**Interdictions :** ne pas relâcher les tolérances, augmenter les seuils, supprimer les gates, déclarer une convergence à partir du seul résidu linéaire, ni calculer un GCI à partir de solutions non convergées. Cette campagne ne qualifie ni N64 TVD, ni N128/N256, ni les coefficients aérodynamiques.

## 7. Critère de sortie

La campagne est terminée lorsque les résultats isolent une cause dominante suffisamment étayée pour justifier un correctif de production. Si plusieurs mécanismes restent confondus, lancer une seconde matrice A/B ciblée plutôt que modifier plusieurs opérateurs simultanément.
