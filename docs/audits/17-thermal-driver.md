# FICHE MÉTHODE 17 — Driver thermique autonome (A2+A3)

**En-tête d'évidence (gouvernance §16.5)**  
- Scope : driver de résolution thermique autonome pour une région maillée d'un `ThermalCase` (A1) : conduction solide stationnaire ou transitoire, diffusion fluide stagnante et advection de flux massique facial prescrit (A3), avec conditions limites thermiques typées (Dirichlet, flux imposé, convection/Robin, puissance de patch, adiabatique), conductivité k(T) traitée par Picard externe et audit énergétique indépendant. Aucune initialisation Navier–Stokes : le flux massique est nul par contrat (solide, fluide stagnant) ou fourni par un champ de faces autoritaire (PRESCRIBED_MASS_FLUX, A3).
- Statut : IMPLÉMENTÉ — A2+A3 (driver fail-closed + conditions typées + advection de flux massique prescrit + 17 cas VV enregistrés). CI exacte HEAD en cours ; aucune qualification physique revendiquée.
- Issue / PR : #1150 / (numéro de PR ajouté à l'ouverture).
- Spécification détaillée : `docs/thermal/thermal-case-data-contract.md`.
- Sources code : `src/cfdx/physics/thermal_driver.h`.
- Tests : `tests/validation/test_thermal_driver_vv.cpp`.
- Limites : une région maillée par appel ; le champ de flux massique facial est figé par appel (pas de couplage pression–vitesse) et doit être discrètement conservatif ; interfaces conjuguées reportées en A4 ; liens et nœuds 0D reportés en A5 ; conductivité tensorielle anisotrope reportée (opérateur de diffusion tensoriel) ; rho et cp figés à l'état initial déclaré de la région ; puissance de patch répartie uniformément sur l'aire intégrée.
- Convention : SI partout ; température absolue en K ; pression absolue en Pa ; aucun changement d'unité implicite.

## 1. Références sources

- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI: 10.1201/9781482234213. Référence de base pour les bilans thermiques, le traitement des conditions limites et la linéarisation par point fixe.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN: 978-3-319-16873-9.
- T. L. Bergman, A. S. Lavine, F. P. Incropera, D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*, 7e éd., Wiley, 2011 — ISBN: 978-0470501979. Condition de film convectif (Robin) et résistance série conduction/film.
- CFDX, fiches méthodes existantes : `docs/audits/01-diffusion-laplacien.md`, `docs/audits/10-thermique-cht.md`, `docs/audits/11-solveur-energie-et-convergence-lineaire.md`, `docs/audits/14-energie-fluide-production.md`, `docs/audits/15-thermal-only-case-contract.md`.

## 2. Formulation continue

Pour une région maillée du cas thermique autonome (solide immobile, fluide stagnant ou fluide à flux massique facial prescrit) :

$$
\rho c_p \frac{\partial T}{\partial t}
+ \nabla\cdot\left(\dot{m}\, c_p T\right)
= \nabla\cdot\left(k(T)\nabla T\right) + q'''
$$

- $\dot{m}$ est le flux massique facial (signe production, kg/s) : identiquement nul pour un solide ou un fluide stagnant ; pour une région PRESCRIBED_MASS_FLUX (A3), c'est le champ de faces autoritaire passé en argument.
- Stationnaire : le terme temporel est abandonné et le pas de temps du cas doit être nul (contrat A1).
- Transitoire : un appel du driver avance l'état d'un pas de temps physique (Euler implicite arrière).
- $\rho$ et $c_p$ sont évalués une seule fois à l'état initial déclaré de la région $(T_0, p_0)$ ; seule la conductivité $k(T)$ est non linéaire et traitée par point fixe.
- Conditions aux frontières de la région :
  - Dirichlet : $T = T_b$ sur le patch ;
  - flux imposé : $-k\,\partial T/\partial n = q''$ (positif vers l'intérieur) ;
  - convection (Robin) : $-k\,\partial T/\partial n = h\,(T_{amb} - T)$ avec $h > 0$ ;
  - puissance : $-k\,\partial T/\partial n = P/A_{patch}$, puissance totale $P$ répartie uniformément sur l'aire intégrée du patch ;
  - adiabatique : $\partial T/\partial n = 0$.

## 3. Discrétisation Volume Fini

L'opérateur est **inchangé** : le driver appelle l'assemblage et le solveur énergie de production (`assemble_energy_equation`, `solve_energy`, fiches 11 et 14) avec le flux massique facial résolu : identiquement nul pour un solide ou un fluide stagnant, champ de faces autoritaire pour une région PRESCRIBED_MASS_FLUX. Aucune tolérance, aucun solveur linéaire, aucun oracle existant n'est modifié.

- Conductivité par cellule : $k(T_P)$ évaluée à chaque itération externe de Picard ; la conductance de face interne reste la résistance série $A_f / (d_{Pf}/k_P + d_{Nf}/k_N)$ de l'assemblage existant.
- Robin : la condition de film est traduite en Dirichlet au point fixe discret de la face,
$$
T_b = \frac{k\,T_P + h\,d\,T_{amb}}{k + h\,d},
$$
  équivalent exact, une fois le point fixe atteint, à une conductance série film + demi-cellule $A_f\,k\,h/(k+h\,d)$ par face (même forme que le traitement CHT existant). $T_b$ est reconstruit à chaque itération externe avec la conductivité courante.
- Flux imposé et puissance : traduits en gradient normal sortant $q''/k_P$ (condition `FIXED_GRADIENT` face-résolue), intégré par l'assemblage comme flux sortant $-k A_f g$.
- Transitoire : l'état $T^n$ est la copie du champ d'entrée ; le pas $dt$ vient du cas (TRANSIENT $\Rightarrow$ $dt>0$, sinon 0). Un appel = un pas.
- Mono-passe exacte : si la loi de conductivité est CONSTANTE et qu'aucune condition CONVECTION n'est présente, le système assemblé ne dépend pas de l'état de température ; une seule passe d'énergie est exacte et la boucle externe dégénère ; le champ de flux massique prescrit ne dépend pas de $T$, la mono-passe reste applicable.
- Advection prescrite (A3) : upwind borné de production (fiche 14) avec correction $-\nabla\cdot(\dot m)\,T_P$ ; le flux est validé fail-closed avant assemblage (présent, à la taille du maillage, dimension 1, fini, conservation discrète par cellule $|\mathrm{net}_c| \le 10^{-9}\max(1, \mathrm{through}_c)$).
- Audit : après convergence, `reconstruct_energy_balance` reconstruit les watts physiques depuis l'état final avec **le même** flux massique facial, **les mêmes** conductivités par cellule et **les mêmes** conditions de faces que la passe finale d'assemblage.

## 4. Traitement des frontières

Toutes les frontières de la région sont typées par patch via `ThermalBoundaryConditions` (nom de patch $\to$ condition typée). La traduction est face-résolue (`ScalarBoundaryFaceCondition`) pour chaque face de bord, reconstruite à chaque itération externe.

Contrôles fail-closed avant tout assemblage :
- patch de type EMPTY rejeté (aucune signification thermique intégrée pour ce driver) ;
- patch sans condition thermique rejeté ; condition référençant un patch inexistant rejetée ;
- valeur Dirichlet hors plage du matériau de la région rejetée ; $h > 0$ requis pour la convection ; toutes les valeurs doivent être finies ;
- état d'entrée $T$ hors plage du matériau rejeté ; température convergée hors plage rejetée après solve ;
- région PRESCRIBED_MASS_FLUX : le champ de flux massique facial doit être présent, nommé, à la taille du maillage, de dimension 1, fini et discrètement conservatif — tout défaut est rejeté avec diagnostic explicite ;
- face de bord à flux entrant ($\dot m_f < 0$) d'une région PRESCRIBED_MASS_FLUX : condition de température exigée (FIXED_TEMPERATURE ou CONVECTION), sinon rejet fail-closed (l'upwind extrapole sinon l'état du propriétaire) ;
- interface conjuguée, lien 0D et matériau anisotrope : diagnostics explicites de report (A4 / A5 / opérateur tensoriel).

## 5. Pièges connus

- Ne jamais reconstruire un flux massique depuis un champ de vitesse : un fluide stagnant a un flux nul par contrat ; un fluide à advection prescrite exige un champ de faces autoritaire, figé par appel et discrètement conservatif — un champ déséquilibré crée de l'énergie fantôme (A3).
- Signe du flux massique facial : positif sort de la cellule propriétaire (convention production). Une face de bord à flux négatif est une entrée : sans condition de température, l'upwind y extrapole l'état du propriétaire au lieu d'advecter un état externe, d'où le rejet fail-closed.
- Le Robin est linéarisé en Dirichlet au point fixe : des itérations externes sont nécessaires et la convergence combine le solveur énergie, le changement de température externe et le bilan indépendant.
- Un problème stationnaire tout-Neumann est singulier (jauge) : le driver n'introduit aucune référence artificielle ni ancrage caché ; l'échec du solveur linéaire est rapporté (`converged=false`), jamais masqué.
- $\rho$ et $c_p$ sont figés à l'état initial déclaré : un cas à propriétés fortement variables doit utiliser des pas de temps adaptés ; aucun figeage silencieux rétroactif n'est appliqué en cours de transitoire.
- La puissance de patch est un flux uniforme sur l'aire intégrée du patch, pas une distribution spatiale.
- L'audit d'énergie doit réutiliser exactement les conductivités par cellule et les conditions de faces de la passe finale ; toute divergence mélange erreurs de stabilisation et erreurs de données.
- La valeur Dirichlet est validée contre la plage déclarée du matériau de la région, pas contre une plage universelle.
- Le gradient imposé est le gradient normal **sortant** $q''/k_P$ : une erreur de signe inverse chauffage et refroidissement.

## 6. Ordre de précision attendu et validation exigée

A2+A3 ne livrent aucun nouvel opérateur : l'ordre spatial et temporel est celui des fiches 11 et 14 (diffusion implicite orthogonale, upwind d'ordre 1, Euler arrière d'ordre 1 en temps). Cette fiche qualifie le **câblage** du driver (typage des frontières, Picard de conductivité, point fixe Robin, transitoire un-pas-un-appel, advection de flux massique prescrit), pas une nouvelle discrétisation.

Oracles exacts discrets sur maillages orthogonaux (systèmes résolus à la main), enregistrés dans `test_thermal_driver_vv` :
- conduction 1-D Dirichlet 400/300 K : profil linéaire exact 375/325 K, mono-passe ;
- génération uniforme, parois à 310 K : $T = 310 + 100/12$ ;
- flux imposé 150 W/m² : $T = 315$ K ; puissance 225 W : $T = 322{,}5$ K ;
- équilibre de convection : $T = T_{amb} = 350$ K ; point fixe Robin avec génération : $T = 325$ K exact ;
- chauffage isolé transitoire : accumulation de 100 W par pas auditée indépendamment ;
- relaxation Dirichlet transitoire : récurrence exacte $(24\,T^n + 12 \cdot 310)/36$ et forme fermée $310 - 10(2/3)^n$ ;
- diffusion fluide stagnante (k = 0,6) : $T = 310 + 50/7{,}2$ ;
- Picard de conductivité $k(T) = 1 + 0{,}001(T-300)$ : racine exacte $300 + (\sqrt{148{,}8}-12)/0{,}024$ ;
- advection d'un état constant (A3) : entrée à 320 K, flux traversant de 1 kg/s, état initial 350 K évacué exactement, $T_0 = T_1 = 320$ K, mono-passe ;
- throughflow avec génération (A3) : $\dot m = \pm 2$ kg/s, entrée 300 K, $q''' = 100$ W/m³ : $T = 300 + 100/4 = 325$ K exact ;
- lavage transitoire (A3) : $a = 10$, $b = 4$, récurrence $T^{n+1} = 300 + (5/7)(T^n - 300)$ et forme fermée $300 + 50\,(5/7)^4$ sur quatre pas.

Bilan énergétique indépendant normalisé $< 10^{-9}$ exigé sur tous les cas physiques. Quatre batteries fail-closed vérifient les diagnostics de report (A4/A5, anisotropie), le mapping de frontières (EMPTY, patch manquant, patch inconnu), les valeurs (Dirichlet hors plage, $h \le 0$, valeur non finie, état d'entrée hors plage) et le contrat du flux massique prescrit (map absente, champ introuvable, taille, valeur non finie, non-conservation discrète, entrée sans condition de température).

Exemptions documentées : MMS, benchmark graphique, GCI et parité CPU/MPI/GPU ne sont pas exigés pour A2+A3 car l'opérateur assemblé (diffusion + upwind) est inchangé et déjà couvert par les fiches 11 et 14 ; ils redeviennent exigibles si l'opérateur change. INCONCLUSIVE n'est pas PASS. CI exacte HEAD verte requise avant fusion ; aucun statut de qualification n'est promu par cette fiche.

## 7. Statut d'implémentation

- Driver : `src/cfdx/physics/thermal_driver.h`.
- Tests : `tests/validation/test_thermal_driver_vv.cpp` ; enregistrement `CMakeLists.txt` (labels `thermal;validation`).
- Réutilise tel quel : `src/cfdx/physics/thermal_case.h` (A1), `src/cfdx/physics/thermophysical_models.h`, `src/cfdx/physics/energy_solver.h`, `src/cfdx/physics/finite_volume_transport.h`.
- Aucune tolérance élargie, aucun test désactivé, aucun oracle existant modifié, aucun registre de qualification modifié.
- Reportés : A4 (interfaces conjuguées), A5 (liens et nœuds 0D), opérateur de diffusion tensoriel pour l'anisotropie.
- Prochaine étape : A4, interfaces conjuguées entre régions maillées du cas.
