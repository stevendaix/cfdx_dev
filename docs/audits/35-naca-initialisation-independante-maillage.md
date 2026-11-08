# NACA N256 - initialisation independante et bord de fuite

Date : 2026-11-07. Base : CFDX lot 4 du meme jour. Ce lot **ne modifie pas**
le solveur C++, les seuils, ni les maillages NACA historiques.

## 1. Question et objectif

Le prolongement N128 -> N256 n'est qu'une condition initiale. La validation
physique de N256 doit etre obtenue avec son propre operateur, ses bilans
conservatifs et ses criteres stricts. Le resultat final devrait etre
independant de la condition initiale si une solution stationnaire unique est
obtenue. Un test depuis un champ uniforme est requis pour verifier ce point.

## 2. Formulation continue

Pour un écoulement incompressible stationnaire, les équations gouvernantes sont :

- Continuité : \\(\\nabla\\cdot\\mathbf{u}=0\\).
- Quantité de mouvement : \\(\\rho(\\mathbf{u}\\cdot\\nabla)\\mathbf{u}=-\\nabla p+\\nabla\\cdot\\boldsymbol{\\tau}\\).

Les métriques ci-dessous mesurent l'itération non linéaire ; elles ne remplacent ni le bilan physique discret ni une étude d'indépendance au maillage.

## Tests directement reproduits (N256, Re=1000, alpha=0)

Conditions identiques : 49 248 cellules, SIMPLE, Van Albada a 100%,
alpha_u=0.60, alpha_p=0.25, 20 iterations, sans pseudo-transitoire.

| Metrique a l'iteration 20 | Uniforme | N128 interpole | Variante geometrique uniforme |
|---|---:|---:|---:|
| Momentum RHS relatif | 1.36635564e-3 | 4.12840852e-4 | 1.19223e-3 |
| Momentum transverse relatif | 2.73754058e-1 | 1.84110077e-2 | 2.86176e-1 |
| Continuite normalisee | 1.1490507e-5 | 3.3315046e-6 | 1.1481e-5 |

**Aucun des trois essais n'a converge en 20 iterations.** Le warm start
abaisse le residu initial, mais le cold start continue a decroitre; ceci ne
prouve pas que la dependance a N128 soit necessaire. Les normalisations de la
composante transverse ont des echelles RHS differentes : prudence pour
les comparaisons de ratios.

## 3. Discrétisation Volume Fini

Le raffinement tangentiel et radial doit être évalué séparément. Les rapports d'aspect extrêmes modifient le conditionnement des opérateurs discrets et peuvent amplifier les erreurs de gradient et de flux.

## Probleme de maillage localise

L'O-grid NACA cosinus serre fortement les points pres du bord de fuite :
le pas tangentiel extremal tend vers O(1/n^2), alors que le pas radial
est raffine essentiellement en O(1/n). D'ou le ratio d'aretes maximum
croissant a la premiere couche du bord de fuite.

| Grille | Max avant | Max avec TE regularisation=0.04 | p99 avant | p99 apres |
|---|---:|---:|---:|---:|
| N64 | 101.94 | 54.24 | 16.72 | 16.62 |
| N128 | 204.74 | 70.70 | 15.51 | 15.49 |
| N256 | 409.56 | 82.54 | 14.54 | 14.42 |

Nouvelle option explicite `--naca-te-regularization 0.04` : remappe les
parametres tangents de la surface et conserve exactement la loi analytique
NACA0012 et ses extremites. C'est une **nouvelle famille de maillages**, a
requalifier independamment; elle ne remplace pas la precedente. Les trois
maillages pilotes n'ont aucune cellule inversee et la porte de similarite
geometrique passe. Le facteur de croissance radial maximal reste proche de
2.5 et demande une etude separee.

Au bout de 20 iterations a froid sur le nouveau N256, le residu momentum
global est environ 13% plus bas qu'avec le maillage historique, mais
le residu transverse est legerement plus eleve. Aucun benefice de convergence
finale ni de precision des forces n'est encore prouve.

## 4. Traitement des frontières

Les conditions aux limites, l'orientation des normales et les flux de frontière doivent rester identiques entre les essais A/B. Toute comparaison de forces doit inclure les contributions de pression et visqueuses et vérifier le bilan global de masse.

## 5. Pièges connus

- Une diminution du résidu sur 20 itérations ne démontre pas la convergence stationnaire.
- Un warm start plus favorable ne démontre pas l'indépendance à l'initialisation.
- Deux familles de maillages différentes ne doivent pas être mélangées dans un même GCI.
- Un résultat de génération/similarité géométrique ne constitue pas une validation physique.

## 6. Reproduction

```bash
# Variante A/B : maillages experimentaux, sans ecraser la famille historique.
python scripts/generate_2d_validation_meshes.py naca0012 \
  --output-dir build/naca_te04 --naca-nested-radial \
  --naca-te-regularization 0.04
python scripts/naca_mesh_similarity_gate.py \
  build/naca_te04/naca0012_mesh_quality_nested_te0p04.json

# Comparer deux historiques acquis sur le meme maillage, meme schema :
python scripts/naca_initialization_ab.py \
  --cold-history cold.csv --warm-history warm.csv \
  --cold-mesh mesh_n256.h5 --warm-mesh mesh_n256.h5 \
  --report ab.json
```

L'outil A/B exige le meme contenu binaire de maillage et le meme calendrier
d'operateurs. Il marque toujours le rapport comme *diagnostic non qualifie* :
le CSV d'historique ne contient pas de preuve cryptographique de provenance.

## Priorites de correction

1. Effectuer un N256 **sans N128**, a budget complet, avec continuation
   UPWIND -> TVD Van Albada puis maintien de l'operateur complet. Comparer
   a un N256 warm-start en conservant toutes les portes physiques.
2. Si les deux convergent, comparer U(x), p(x) (reference de jauge coherente),
   Cd, Cd pression, Cd visqueux, Cl, bilans et profils Cp/Cf. Refuser une
   qualification si les solutions dependent du seed au-dela de l'incertitude.
3. Qualifier deux familles geometriques (historique et TE regularisee)
   avec *trois solutions convergentes chacune*. Ne pas melanger les deux
   familles dans un meme GCI.
4. Auditer les cellules extreme-TE, la croissance radiale, la quadrature
   des efforts et les termes Rhie-Chow/pression. Ne pas confondre un
   residu lineaire de pression faible avec l'exactitude du couplage.
5. En cas de forte sensibilite au seed apres convergence, investiguer
   solutions multiples, couplage instable et une simulation transitoire
   independante, et non simplement forcer une relaxation plus forte.

## References

- NASA grid sequencing et justification du startup grossier :
  https://www.cs.odu.edu/~mln/ltrs-pdfs/NASA-98-tm208444.pdf
- NASA convergence spatiale, donnees de GCI :
  https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html
- OpenFOAM, transferts inter-maillages :
  https://doc.openfoam.com/2312/tools/post-processing/function-objects/field/mapFields/
- OpenFOAM, initialisation alternative potentialFoam :
  https://doc.openfoam.com/2306/tools/processing/solvers/rtm/basic/potentialFoam/
- SU2, multigrille et reprise :
  https://su2code.github.io/docs_v7/Setup-Guide-RANS-External-Aero/
- PETSc, multigrille non lineaire FAS :
  https://petsc.org/release/manualpages/SNESFAS/SNESFAS/
- NACA0012 Re=1000: regime stationnaire a faible incidence, onset du
  shedding autour de 8 degres (etudes numeriques, pas oracle CFDX) :
  https://journals.sagepub.com/doi/10.1177/17568293211055656

Ce lot est une amelioration de **diagnostic et de generation du maillage**,
pas une certification N256 ni un correctif du solveur SIMPLE lui-meme.

## 7. Statut d'implémentation

**Statut :** `DIAGNOSTIC_ONLY` — résultats préliminaires, aucune qualification physique N256 et aucun GCI certifié. Les essais décrits sont des outils de diagnostic et de préparation de campagnes ; ils ne modifient pas les critères d'acceptation du solveur.
