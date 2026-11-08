# CFDX — Requalification des limites PCD Fp et pression BFS V4

**Date** : 8 novembre 2026. **Base exacte** : ZIP original SHA-256
`c3a0d40206829596dca61871c10f6079724c6d0e11bbf63a47d4e894e914e62c`.
**Statut** : correctifs logiciels et preuves partielles ; pas de convergence
COUPLED NACA, ni de qualification physique BFS. Aucun seuil physique diminué.

## 1. Opérateur auxiliaire PCD Fp

Le sous-solveur PCD agit sur l'équation auxiliaire de pression, dont les
conditions aux limites ne sont pas obligatoirement identiques à celles du
champ de pression physique. Le chemin existant assemblait Fp et Kp avec la
condition p-Neumann du farfield NACA, puis une seule ligne de jauge.

Sur la capture N64 initiale : 3 096 lignes, 15 219 nnz, pivot LU minimal
~6.006e-16 ; un RHS difficile donne une correction de norme ~1.58e12
et un résidu direct relatif ~3.2e-4.

L'option `pcd_pressure_inflow_dirichlet_experimental` assemble des
conditions **homogènes Dirichlet seulement sur les faces d'entrée**
(déterminées par `mass_flux(face) < 0`, y compris sur un farfield mixte),
pour Fp et Kp, jamais pour la matrice physique COUPLED. Les faces WALL,
SYMMETRY, EMPTY, INTERFACE et PERIODIC restent inchangées.
La valeur par défaut est `false` et le test C++ vérifie l'identité stricte
du chemin historique.

Sur le N64 modifié, les 65 lignes dont la somme est non nulle correspondent
à la jauge et aux 64 faces entrantes : pivot minimal ~0.0245037,
corrections directes modérées et résidu inférieur à 3e-15 sur trois RHS
capturés. Les applications FGMRES/ILU0 terminent fréquemment en 16–17 pas,
avec résidus vrais relatifs de l'ordre de 1e-9. Sur un essai
N64 de 90 s, **1 951 applications Fp sur 1 951** atteignent leur propre
contrat (maximum mesuré 9.981e-9, médiane 2.203e-9 et médiane 17
itérations). Ce nombre ne mesure pas les itérations non linéaires.
Le calcul extérieur NACA
COUPLED n'a **pas** terminé avant l'arrêt diagnostique de 90 s : ne pas le
présenter comme convergé. Les résultats ne portent pas encore sur N128/N256
pour cette nouvelle variante.

Activation expérimentale dans le pilote NACA :

```bash
CFDX_NACA_ALGORITHM=COUPLED CFDX_NACA_COUPLED_SCHUR=pcd \
CFDX_NACA_PCD_INFLOW_DIRICHLET=1 \
CFDX_NACA_CONVECTION_CELL_LIMITED=1 \
CFDX_NACA_FARFIELD_PRESSURE_NEUMANN=1 \
CFDX_PCD_FP_CAPTURE_PREFIX=/tmp/cfdx/naca/inflow \
  build/test_naca0012_qualification --single \
  build/validation_meshes/naca0012_laminar_n64.h5 0.6 0.25 1
```

La même option peut être fixée dans les cas de production COUPLED + Schur PCD
avec `numerics.raw_settings.pcd_pressure_inflow_dirichlet_experimental = "1"`.
Toute autre utilisation est refusée. La pression physique n'est jamais
remplacée par ces conditions auxiliaires.

## 2. BFS OpenFOAM V4 conservé

Le générateur six blocs `scripts/generate_th03_openfoam_blockmesh.py`
conserve 20 540 cellules / 82 493 faces et les huit patches. Une incohérence
`OUTLET` + `OUTLET_PRESSURE` empêchait le chargement en production ; elle est
corrigée en `PRESSURE_OUTLET` + `OUTLET_PRESSURE`.

Le fichier complet doit être nommé `*.cfdx.h5`. L'initialisation turbulence
courante du générateur utilise un ratio d'entrée `nut/nu = 0.009` et
`omega_in = k / (nu * 0.009) ≈ 7763.53 s^-1`; elle diffère du réglage
historique 181728 s^-1 et n'est pas prouvée identique au `0.orig` OpenFOAM.

### Reproduction de la vraie matrice pression

```bash
PYTHONPATH=src/cfdx/python python scripts/generate_th03_openfoam_blockmesh.py \
  --resolution full --max-iterations 3 \
  --pressure-linear-tolerance 1e-11 --pressure-gmres-restart 160 \
  --output build/th03_v4/full.cfdx.h5

build/cfdx_production_solver --mesh build/th03_v4/full.cfdx.h5 \
  --output-dir build/th03_v4/out --iterations 3 --experimental-sst \
  --benchmark-export --pressure-replay-prefix build/th03_v4/pressure \
  --pressure-replay-iterations 1

build/cfdx_linear_replay build/th03_v4/pressure_i1_c0.cfdxlin native ilu0
```

La capture de production réelle est une matrice symétrique 20 540 × 20 540
(102 034 coefficients CSR) dont 10 249 lignes sont détectées anisotropes et
le contraste maximal de couplage hors diagonale vaut 5.37e7.

| Méthode sur la capture (même b, même x0) | Vrai résidu relatif final |
|---|---:|
| FGMRES + AMG natif | 2.361583e-10 |
| FGMRES + AMG directionnel | 2.663351e-10 |
| FGMRES + ILU0 | 1.705687e-3 |
| FGMRES + Jacobi | 6.627763e-2 |
| LU SciPy directe | 3.912932e-10 |
| Critère exigé | 1.0e-11 |

La LU directe retourne une correction de pression de norme 4.292019e7,
malgré une norme du second membre 1.820455e-3. Le plancher de résidu
relatif est donc aussi observé hors CFDX : ce diagnostic doit guider la
révision du conditionnement des coefficients de bord/du correcteur, et non
un assouplissement des tolérances.

## 3. Résultats et qualification

Les trois essais sur la configuration BFS corrigée rejettent la première
itération physique à 5.220932e-10, 5.258521e-10 et 3.500920e-10. Aucune
SST convergée, distribution y+, Cf, Cp ou position de réattachement n'est
publiée pour cette grille. Le seul audit `PASS_GEOMETRY_ONLY` ne vaut pas
qualification : 0 face inversée, angle max ~9.413 degrés, et non-orthogonité
modérée n'assurent pas un opérateur de pression bien conditionné.

Les vérifications antérieures Couette, Poiseuille, Ghia, N9 et la chaîne
thermal/RANS ne sont pas automatiquement qualifiées par ces modifications.
La CI complète, le DebugSanitizers et les calculs longs restent ouverts.

## 4. Risques contrôlés

- Aucun changement de maillage, propriétés physiques ou seuil d'acceptation.
- Le nouveau comportement auxiliaire PCD n'est jamais implicite.
- La capture n'est activée que par deux options fournies ensemble ;
  les retries d'une même itération remplacent le fichier correspondant.
- Des tests algébriques et maillage/BC vérifient le défaut inchangé,
  le sens du flux et la conformité HDF5.
- Un succès Fp ne constitue pas un succès du Krylov monolithique COUPLED.

## 5. Régression CTest de l'outil de rejeu corrigée

`test_pcd_fp_replay` était enregistré en CTest sans argument, donc renvoyait
une erreur d'utilisation, même lorsque le calcul numérique était correct.
Il exécute désormais le même autotest déterministe sur un appel sans argument
ou avec `--self-test`. Le mode normal `capture.cfp [budget] [tol]` et sa
sortie numérique restent inchangés. Le CTest explicite du mode capture reste
actif. Les **6 tests CTest** ciblés ont passé, sans marquer un échec comme
attendu ni le désactiver. Les **48 tests Python** BFS/TH-03 ciblés ont passé.


## 1. Périmètre et hypothèses

Cette fiche formalise les frontières auxiliaires PCD \(F_p/K_p\) et le cas de pression BFS V4. L'opérateur auxiliaire ne remplace jamais l'opérateur physique de continuité ni le flux de masse conservatif. Les réglages de frontière expérimentaux sont désactivés par défaut.

## 2. Formulation continue

Le champ physique satisfait \(\nabla\cdot u=0\) et l'équation de quantité de mouvement incompressible. La pression est définie à une constante près en l'absence de référence Dirichlet. Dans la décomposition bloc, le complément de Schur exact est \(S=C-DA_{uu}^{-1}G\); le PCD utilise des opérateurs auxiliaires pour approcher son inverse, sans prétendre calculer \(A_{uu}^{-1}\) exactement.

## 3. Discrétisation Volume Fini

La continuité sur une cellule est le bilan des flux de masse orientés aux faces. La matrice de pression CSR est assemblée à partir des coefficients de faces et des contributions de frontières; l'opérateur de face conservatif est une représentation distincte utilisée pour vérifier les flux. Une correction algébrique high/low ne certifie pas le résidu obtenu après arrondi FP64 ni le bilan conservatif des faces. Chaque quantité doit être recalculée avec son opérateur de référence.

## 4. Traitement des frontières

Pour les opérateurs auxiliaires PCD, l'option expérimentale applique une condition homogène Dirichlet aux seules faces d'entrée détectées par \(\dot m_f<0\). Les patches WALL, SYMMETRY, EMPTY, INTERFACE et PERIODIC conservent leur traitement. Les conditions de pression physique et la jauge restent indépendantes. Toute activation doit être explicite, tracée dans les paramètres du cas et testée contre le chemin par défaut.

## 5. Pièges connus

- Une condition auxiliaire PCD n'est pas une condition physique de pression.
- La sélection des faces entrantes dépend du signe du flux orienté owner/neighbour; une convention inversée modifie le problème.
- Une jauge de pression ne doit pas être confondue avec une frontière physique.
- Un résidu CSR, un résidu d'opérateur implicite de faces et un résidu après arrondi FP64 sont trois mesures distinctes.
- Une correction high/low ne doit pas être promue comme champ physique si le champ FP64 arrondi échoue.
- Un maillage valide géométriquement peut produire un système de pression mal conditionné.

## 6. Vérifications et critères d'acceptation

Les tests doivent couvrir le sens du flux, l'ensemble des types de patches, l'absence de modification du chemin par défaut, la jauge et le null-space. Pour chaque cas, consigner séparément le résidu CSR original, le résidu face implicite, le résidu après arrondi FP64, le bilan global de masse et les erreurs locales de conservation. Le seuil algébrique est évalué sans relâchement; les preuves physiques nécessitent en plus une campagne convergée.

## 7. Statut d'implémentation

**Statut :** option Dirichlet auxiliaire PCD à l'entrée implémentée en mode expérimental opt-in; opérateurs de pression high/low et face conservatif disponibles à des fins de diagnostic. Aucun branchement expérimental n'est activé par défaut.

**Statut des preuves :** les tests de contrat et les diagnostics de pression sont partiels. Le résidu CSR après arrondi et les résultats BFS V4 ne satisfont pas encore le critère annoncé; aucune qualification physique SST n'est revendiquée.

**Références :** Kay, Loghin & Wathen (2002), DOI: 10.1137/S1064827500389549; Rhie & Chow (1983), DOI: 10.2514/3.8284; Driver & Seegmiller (1985), DOI: 10.2514/3.8890.

## 8. Contrat des frontières d'entrée dans le chemin de production

Le type de patch et l'opérateur de champ sont deux informations distinctes. Un patch topologique INLET peut porter soit une vitesse prescrite \(u=u_{in}\) (\(\texttt{FIXED}\)), soit un opérateur ouvert à gradient nul \(\partial_n u=0\) (\(\texttt{ZERO_GRADIENT}\)). Le second cas est utilisé par les canaux entraînés par force volumique; il ne doit pas être confondu avec un inlet pression/débit importé. Les sources \(\texttt{pressure-inlet}\), \(\texttt{mass-flow-inlet}\) et \(\texttt{pressure-far-field}\) restent rejetées tant que leurs opérateurs dédiés ne sont pas implémentés.

Dans les fixtures SST de couplage thermique à débit nul, le champ de vitesse est imposé à zéro et le pressure outlet fixe la jauge; le champ \(k,\omega\) non nul est alimenté par une vitesse de référence d'entrée positive distincte du vecteur de vitesse moyenne nul; cela permet de vérifier l'effet de \(Pr_t\) sur \(\lambda_{eff}\) sans transformer ce smoke test en qualification d'écoulement. Pour un futur cas à vitesse d'entrée non nulle, l'inlet FIXED devra être associé à un outlet PRESSURE_OUTLET/OUTLET_PRESSURE avec une pression de référence explicite. Les canaux entraînés par force volumique conservent leurs frontières ouvertes à gradient nul.

Un patch importé de type générique \(\texttt{UNKNOWN}\) peut recevoir un opérateur reconnu lorsque le fichier de cas déclare explicitement son type physique; les patches déjà typés (WALL, SYMMETRY, EMPTY, etc.) doivent rester compatibles et les incompatibilités sont rejetées. Cette règle évite de déduire une physique du seul nom du patch sans interdire les patches OpenFOAM génériques correctement configurés.

Le contrat de production vérifie la concordance du type topologique \(\texttt{PatchType::INLET}\) et du type physique \(\texttt{BCType::INLET}\), puis sélectionne explicitement FIXED ou ZERO_GRADIENT. Les tests unitaires couvrent les deux branches et le rejet des inlets pression/débit. Les fixtures N10 génériques déclarent explicitement leurs patches comme WALL avant d'attacher des conditions no-slip; elles ne déduisent pas une physique de l'ordre des patches VTK.

**Statut :** opérateurs FIXED et ZERO_GRADIENT pris en charge sur patch INLET; pression/débit/far-field non pris en charge et rejetés. Cette extension ne constitue pas une validation physique du cas de canal.
