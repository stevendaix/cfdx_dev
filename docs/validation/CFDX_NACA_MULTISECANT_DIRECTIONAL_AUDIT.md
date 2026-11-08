# CFDX — Accélération non linéaire TVD NACA, preuves et limites

## Périmètre

Le solveur expérimental `nonlinear_multisecant` propose un nouvel état initial au **pas SIMPLE suivant**, sans remplacer la projection pression–vitesse. Le vecteur extrapolé contient simultanément les 3 composantes de la vitesse aux cellules, la pression aux cellules et les **flux massiques conservatifs aux faces**. Les mêmes poids d'extrapolation sont appliqués aux cinq champs, puis une limite de pas scalaire commune est appliquée. Aucune interpolation indépendante des flux ne doit remplacer la reconstruction Rhie–Chow du pas accepté.

Deux modes opt-in : `anderson` (défaut de point fixe `G(X)-X`) et `ngmres` (résidu signé de l'opérateur physique TVD complet réassemblé sur `G(X)`). Ils ne sont pas des implémentations identiques aux solveurs SNES de PETSc : ils utilisent une mémoire courte, un solveur dense réduit régularisé et la projection SIMPLE propre à CFDX.

### Exécution

```bash
CFDX_NACA_CONVECTION_TVD=1 \
CFDX_NACA_MULTISECANT=anderson \
CFDX_NACA_MULTISECANT_DEPTH=4 \
CFDX_NACA_MULTISECANT_MIXING=0.7 \
CFDX_NACA_MULTISECANT_START=20 \
CFDX_NACA_PHYSICAL_HISTORY_CSV=anderson_n64.csv \
./test_naca0012_geometry --single /chemin/naca0012_laminar_n64.h5 0.6 0.25 400
```

Défaut de `CFDX_NACA_MULTISECANT` = méthode désactivée. `CFDX_NACA_MULTISECANT_STEP_CAP` contrôle la taille maximale de l'extrapolation (et non les critères physiques). Les combinaisons avec turbulence SST, régime transitoire, schémas de démarrage variables, pseudotransitoire et autres accélérateurs expérimentaux sont refusées explicitement.

### Équations et critères

- Toutes les tentatives, y compris avec prédicteur multisécant, réassemblent `A(U)` et le second membre de la convection TVD sur l'état courant, puis recalculent le prédicteur H/A, la pression et le flux Rhie–Chow.
- Un prédicteur est rejeté et le pas refait sans extrapolation s'il détériore le résidu composante ou la continuité au-delà du garde-fou conservateur. Les critères physiques de convergence finaux ne changent pas.
- L'orthogonalisation de l'historique NGMRES est basée sur les résidus signés divisés par des **échelles fixées au premier état accepté** ; il ne faut pas renormaliser par le RHS instantané de la tentative.
- Un `PASS` exige encore les résidus du momentum x/y, l'audit indépendant aux faces, la continuité, les bilans locaux et les seuils des coefficients stables. Les historiques seuls ne certifient pas une qualification.

## Jacobien directionnel TVD

`nonlinear_directional_jvp` effectue une différence finie centrale de l'opérateur **réévalué** :

\[
J_{\rm TVD}(U)v \simeq \frac{R_{\rm TVD}(U+hv)-R_{\rm TVD}(U-hv)}{2h}.
\]

`test_tvd_skew_affine_mms` compare le résultat du chemin matrice CSR à une somme **indépendante** des flux aux faces, pour minmod, van Leer et van Albada. Une direction de Newton de petit système est ensuite résolue avec une matrice directionnelle issue des reassemblages et acceptée par réduction du résidu vrai. **La pression et les flux de masse restent figés dans cette démonstration.** Elle ne qualifie PAS un Newton–Krylov couplé NACA, pour lequel une direction entièrement couplée incluant la dérivée de Rhie–Chow resterait à construire.

## Vérificateur A/B

```bash
python scripts/compare_naca_multisecant.py --manifest \
  /chemin/manifest_n64_400.json --output /chemin/comparison400.json
pytest -q tests/python/test_compare_naca_multisecant.py
```

Le manifeste fixe le SHA-256 du binaire et du maillage, les valeurs de relaxation, le nombre d'itérations et les empreintes des CSV. Le vérificateur rejette la moindre ligne manquante, la variation de relaxation, une convection non TVD finale, les NaN et les fichiers altérés. Il inscrit toujours `qualification=NOT_CERTIFIED` pour une campagne dont on ne dispose que de l'historique.

## Première campagne exécutée le 8 novembre 2026

NACA N64 laminaire, 400 itérations, TVD minmod, SIMPLE, alpha_U=0.6, alpha_p=0.25 :

| Méthode | Momentum x/RHS-x | Momentum y/RHS-y | Continuité | Certificat |
| --- | ---: | ---: | ---: | --- |
| Picard | 5.685868987e-4 | 6.140065787e-3 | 1.657767979e-7 | NON |
| Anderson, m=4, omega=0.7 | 1.034967958e-4 | 3.860996527e-5 | 1.870082971e-8 | NON |
| Anderson, m=6, omega=0.5, cap=1 | 2.352013767e-4 | 9.010452156e-5 | 5.564498070e-8 | NON |
| NGMRES, m=3, omega=0.3, cap=1 | 4.968412713e-4 | 2.766658170e-4 | 1.347705406e-7 | NON |

Un NGMRES moins borné a rencontré la protection contre la croissance sévère vers l'itération 124. Ne pas traiter cette exécution interrompue comme une campagne terminée. Aucun N64, N128 ou N256 n'est déclaré physiquement qualifié dans ce lot ; aucun calcul de GCI n'est autorisé.

### Sources méthodologiques

- PETSc SNESANDERSON : https://petsc.org/release/manualpages/SNES/SNESANDERSON/
- PETSc SNESNGMRES : https://petsc.org/main/manualpages/SNES/SNESNGMRES/
- PETSc SNES (cadre des résidus non linéaires) : https://petsc.org/main/manual/snes/
