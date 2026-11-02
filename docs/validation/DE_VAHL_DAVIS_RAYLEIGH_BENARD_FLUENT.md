# CFDX : convection naturelle de de Vahl Davis, Rayleigh-Benard, et benchmark Fluent

## Statut et périmètre

- **Implémenté et exécuté localement** : cavité 2-D extrudée de de Vahl Davis,
  `Ra=1e3`, `Pr=0.71`, avec le **vrai solveur de production** :
  SIMPLE/Rhie-Chow, boucle externe Boussinesq `run_production_buoyant_outer_loop`,
  solveur d'énergie et bilan thermique indépendant. Sur trois maillages,
  **PREQUALIFICATION_PASS** en Release, CTest 33.57 s.
- **Non validé / à implémenter** : Rayleigh-Benard transitoire et rouleaux au
  voisinage du seuil; balayage `Ra=1e4..1e6`; Fluent VMFL056 rayonnement DOM
  et conduction. Aucun de ces trois sujets ne doit être marqué PASS.
- L'exécutable est intentionnellement étiqueté `prequalification;long`, non
  `qualification` complète. Le statut ne doit pas être promu tant que des
  mesures fines et la convergence des profils de vitesses ne sont pas acquises.

## 1. Cavité différentiellement chauffée : de Vahl Davis (1983)

**Géométrie** : carré `L=1 m`, extrudé de `1 m` en z pour le solveur 3-D.
Pas de pénétration ni glissement sur six parois, parois z adiabatiques.
Température 305 K sur x=0, 295 K sur x=L; parois y=0/L et z=0/1 adiabatiques.
La gravité agit vers -y, `g=(0,-9.81,0)`.

**Modèle** : densité de référence `rho=1 kg/m3`, viscosité cinématique
`nu=0.01 m2/s`, capacité calorifique `cp=1 J/(kg K)`, diffusivité thermique
`alpha=nu/Pr` et conductivité `k=rho cp alpha`. Coefficient de dilatation
`beta=Ra nu alpha/(g DeltaT L^3)`, soit `beta DeltaT = 0.01436` : le modèle
Boussinesq reste dans le domaine des faibles variations de densité.
Ces propriétés non réalistes pour l'air en unités SI visent la **similitude
dimensionnelle** du cas de référence (Ra et Pr); elles ne doivent pas être
présentées comme des propriétés de l'air réel.

Definitions :

\[
Ra = \frac{g \beta \Delta T L^3}{\nu \alpha},\qquad
Pr = \frac{\nu}{\alpha},\qquad
Nu = \frac{L}{k\Delta T A}\int_{\mathrm{paroi}} q''\,dA.
\]

Le test mesure le flux de paroi par la **distance cellule-face** et calcule
indépendamment Nu chaud et Nu froid. Il vérifie : solveur couplé convergé,
valeurs finies, circulation non nulle, bilan d'énergie, accord des deux flux,
erreur Nu du maillage fin <= 5 %, écarts de vitesses échantillonnées <= 12 %,
et décroissance stricte de l'erreur absolue de Nu sur 8/12/16.
Aucun critère d'un test existant n'est modifié.

### Références de de Vahl Davis

| Ra | Nu moyen | u_max L/alpha | v_max L/alpha |
|---:|---:|---:|---:|
| 1e3 | 1.118 | 3.649 | 3.697 |
| 1e4 | 2.243 | 16.178 | 19.617 |
| 1e5 | 4.519 | 34.73 | 68.59 |
| 1e6 | 8.800 | 64.63 | 219.36 |

Source primaire : G. de Vahl Davis, *Natural convection of air in a square
cavity: a bench mark numerical solution*, IJNMF 3(3) 249-264 (1983),
https://doi.org/10.1002/fld.1650030305 .
Tableau des valeurs avec comparaison indépendante :
https://cardinal.cels.anl.gov/vv/square_cavity_natural_convection.html .

### Résultats mesurés : binaire C++ Release, 31 octobre 2026

| Maillage | Nu chaud | Nu froid | Erreur sur Nu ref. | u_max L/alpha* | v_max L/alpha* |
|---:|---:|---:|---:|---:|---:|
| 8 x 8 | 1.18373197013 | 1.18373197014 | 5.88 % | 3.24569 | 3.26107 |
| 12 x 12 | 1.14957837167 | 1.14957837168 | 2.82 % | 3.21587 | 3.22579 |
| 16 x 16 | 1.13574183264 | 1.13574183265 | 1.59 % | 3.28419 | 3.35087 |

(*) Vitesses calculées à partir des centres de cellules les plus proches de
`x=0.5` ou `y=0.5`, pas d'interpolation haute précision; les valeurs
restent environ 9-10 % sous la référence sur la grille fine. Un test plus
strict devra extraire les profils interpolés et utiliser 32/64/128 cellules
avant qualification physique complète. L'erreur Nu est décroissante, mais
**trois maillages ne suffisent pas à prouver l'ordre asymptotique**.

Tous les calculs de cette première campagne ont obtenu `converged=true`,
66 itérations de boucle externe et un bilan énergétique global relatif
inferieur a `2e-12`. La porte Nu 5 % est évaluée sur la grille fine, pas
sur le maillage 8x8 (lequel échoue volontairement à cette précision).

### Reproduction

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON
cmake --build build --target test_natural_convection_de_vahl_davis -j 2
ctest --test-dir build -V -R '^test_natural_convection_de_vahl_davis$'
# Rejeu diagnostique d'une seule grille : ne prouve PAS la convergence spatiale
build/test_natural_convection_de_vahl_davis --grid 16
```

## 2. Rayleigh-Benard : vraie campagne d'apparition des rouleaux (non executee)

Une couche horizontale de hauteur `H`, isotherme **chaude en bas** et
**froide en haut**, doit recevoir une perturbation thermique ou de vitesse
non nulle et subir une intégration temporelle. Pour deux plaques rigides
et une extension horizontale infinie (conditions periodiques), le seuil
lineaire theorique est `Ra_c~1707.76`, avec `k_c~3.117/H`.
Ce seuil ne vaut pas tel quel pour une boite carree fermee a parois laterales
no-slip : largeur, condition latérale et aspect ratio modifient le seuil.

**Contrat de test à ajouter après disponibilité du solveur transitoire
couplé et des bonnes conditions périodiques** :

1. `Pr=0.71`, domaine `Lx/H` grand et documente, resolution verticale et
   horizontale independantes, maillages raffines.
2. Initialiser `T=T_hot - DeltaT*y/H + eps*sin(kx)*sin(pi*y/H)` et `U` presque
   nul avec perturbation non nulle, pour ne pas rester artificiellement
   sur la branche conductive instable.
3. `Ra=1500` : decroissance du mode perturbe et `Nu -> 1`.
4. `Ra=2000` : croissance initiale de l'amplitude apres un regime lineaire,
   puis saturation convective; `Nu > 1`; verifier le nombre et la longueur
   d'onde des rouleaux et le bilan du flux chaud/froid.
5. Fit numerique de la croissance/décroissance pour differentes resolutions,
   temps de pas et amplitudes; confronter le seuil extrapole a la reference,
   sans assouplir les criteres. Un champ de vitesse nul a Ra=2000 **n'est
   pas une preuve** d'absence de convection : l'état de conduction reste
   une solution stationnaire formelle.

Sources :
https://doi.org/10.1103/PhysRevE.62.7987 (rouleaux en géométrie finie);
https://doi.org/10.1016/j.crme.2010.07.008 (seuil rigide/rigide 1707.76).

## 3. Référence Ansys Fluent VMFL056 (non exécutée)

Le manuel officiel Ansys Fluids Verification, cas **VMFL056**, décrit une
cavité `1 m x 1 m` couplant **conduction et rayonnement DOM** : conductivite
`k=1 W/(m K)`, coefficient d'absorption `0.228 m^-1`, coefficient de
scattering nul, paroi chaude `100 K`, trois autres parois `50 K`.
L'oracle est le profil de température adimensionnelle le long de la coupe
centrale, pas un simple statut de solveur `PASS`.

Les conditions thermiques et optiques doivent être reproduites exactement,
avec l'équivalent DOM du Fluent, puis comparées au profil publié, bilan
radiatif et convergence en maillage/directions angulaires. Les tests
existants de composants P1/S2S/DOM ne valident pas automatiquement VMFL056.

Source Ansys :
https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/fbu_vm/Hlp_VMFL056.html

## Intégration PR #1149

Le patch est construit sur la copie locale d'archive `298e4d01` incluant le
cas analytique de convection forcee. Il doit être **rebasé** sur le HEAD
courant de #1149, en ne portant que les ajouts absents : un nouveau test
C++, une ligne d'enregistrement CMake, une ligne de proprietes CTest et ce
fichier de methode. Ne pas ecraser la version deja avancee du solveur.
