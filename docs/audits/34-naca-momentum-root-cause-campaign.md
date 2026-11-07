# FICHE MÉTHODE 34 — Campagne forensique de convergence NACA0012

**Statut :** protocole de diagnostic uniquement. Aucun changement de schéma, coefficient, relaxation, tolérance ou critère d'acceptation. Les sorties ne sont jamais des certificats de convergence.

## 1. Références sources

- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016, DOI 10.1007/978-3-319-16874-6.
- J. H. Ferziger, M. Perić, R. L. Street, *Computational Methods for Fluid Dynamics*, Springer, DOI 10.1007/978-3-319-99693-6.
- W. L. Oberkampf, C. J. Roy, *Verification and Validation in Scientific Computing*, Cambridge University Press, DOI 10.1017/CBO9780511760396.
- NASA NPARC, *Examining Spatial (Grid) Convergence*, https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html.

## 2. Formulation continue

Sur N64, la continuation atteint 100 % de TVD à l'itération 120, puis la solution stagne et termine non convergée à l'itération 330. Valeurs observées :

- Momentum RHS relatif : 7.25783e-5
- Maximum des composantes momentum : 5.374e-3
- Résidu momentum indépendant en faces, composante y : 2.06470e-3
- Continuité normalisée : 2.09097e-7
- Résidu linéaire pression : 8.43316e-12

Le résidu linéaire de pression ne permet pas d'attribuer la stagnation à la pression, ni de conclure que le système non linéaire est convergé. Il faut identifier quelle contribution physique, quelle région et quelle composante empêchent la diminution du résidu momentum.

## 3. Discrétisation Volume Fini

### Instrumentation disponible après le merge #1187

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

### Exécution reproductible

Construire le binaire de validation test_naca0012_qualification depuis tests/validation/test_naca0012_geometry.cpp. Utiliser exactement le même fichier de maillage N64 pour tous les essais.

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

### Matrice A/B — une seule modification à la fois

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

Ne pas activer plusieurs variantes de convection simultanément. A4–A6 sont comparés séparément à A1. Le pseudo-temps, alpha_u=0.8 et les essais de prolongement restent hors de cette première matrice : ils ajouteraient des facteurs confondants. Pour chaque essai, recopier la commande A0 et ajouter uniquement les variables indiquées. Chaque essai doit avoir ses propres chemins de sortie.

## 4. Traitement des frontières

Conserver strictement les mêmes conditions aux limites et la même référence de pression pour tous les essais. Les champs importés d'un état prolongé ne changent pas les conditions aux limites. Le flux de masse utilisé pour le budget doit rester le flux conservatif corrigé réellement consommé par le solveur ; le flux interpolé collocated n'est qu'une observable de comparaison.

Les flux de frontière, leur orientation et leur contribution au bilan de masse doivent être rapportés séparément. Une différence entre flux autoritaire et flux collocated ne doit pas être appelée « correction Rhie–Chow pure » : elle inclut aussi le partage pression-relaxation.

## 5. Pièges connus

- Un résidu linéaire pression faible n'est pas une preuve de convergence non linéaire.
- Un résidu global peut masquer une composante y dominante.
- Une contribution signée importante peut être compensée par une autre ; examiner les sommes signées **et** les maxima de valeurs absolues par cellule.
- Les diagnostics cellwise sont échantillonnés ; vérifier les itérations échantillonnées et ne pas supposer une évolution monotone entre deux échantillons.
- Le résidu relatif peut être amplifié par une échelle RHS de référence très petite ; rapporter les normes absolues, relatives et leurs dénominateurs sans remplacer le critère officiel.
- Une corrélation avec le skewness ou le rapport d'aspect ne suffit pas à démontrer la causalité ; exiger une comparaison A/B contrôlée.
- Les solutions non convergées ne peuvent alimenter un GCI qualifié.

## 6. Ordre de précision attendu et validation exigée

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

Livrables avant toute correction du solveur :

- Un dossier d'artefacts par essai : commande exacte, environnement, hash du binaire, hash du maillage, stdout, trace forensique, historique physique, budget cellwise et JSON d'analyse.
- Un tableau comparatif des dernières itérations acceptées, sans mélanger les normes ou leurs dénominateurs.
- Une carte des cellules fautives et leur qualité géométrique : non-orthogonalité, skewness et rapport d'aspect.
- Une conclusion causale avec au moins une observation discriminante et une contre-épreuve.
- Une proposition de correction isolée avec un test de régression reproduisant le défaut.

**Interdictions :** ne pas relâcher les tolérances, augmenter les seuils, supprimer les gates, déclarer une convergence à partir du seul résidu linéaire, ni calculer un GCI à partir de solutions non convergées. Cette campagne ne qualifie ni N64 TVD, ni N128/N256, ni les coefficients aérodynamiques.

## 7. Statut d'implémentation

- Le protocole est ajouté dans cette PR ; les mécanismes de collecte et d'analyse référencés sont déjà présents dans le code fusionné par #1187.
- La campagne A/B n'a pas encore été exécutée par cette PR.
- Aucune cause dominante n'est établie à ce stade.
- Aucun correctif de production, changement de tolérance, certificat de convergence ou résultat GCI n'est revendiqué.
- La campagne est terminée lorsque les résultats isolent une cause dominante suffisamment étayée pour justifier un correctif de production. Si plusieurs mécanismes restent confondus, lancer une seconde matrice A/B ciblée plutôt que modifier plusieurs opérateurs simultanément.
