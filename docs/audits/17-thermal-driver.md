# FICHE MÉTHODE 17 — Driver thermique autonome (A2+A3+A4)

**En-tête d'évidence (gouvernance §16.5)**  
- Scope : driver de résolution thermique autonome pour une région maillée d'un `ThermalCase` (A1) : conduction solide stationnaire ou transitoire, diffusion fluide stagnante et advection de flux massique facial prescrit (A3), interfaces conjuguées et résistance de contact entre régions maillées avec solveur couplé multi-régions (A4), avec conditions limites thermiques typées (Dirichlet, flux imposé, convection/Robin, puissance de patch, adiabatique), conductivité k(T) traitée par Picard externe et audit énergétique indépendant. Aucune initialisation Navier–Stokes : le flux massique est nul par contrat (solide, fluide stagnant) ou fourni par un champ de faces autoritaire (PRESCRIBED_MASS_FLUX, A3).
- Statut : IMPLÉMENTÉ — A2+A3+A4 (driver fail-closed + conditions typées + advection de flux massique prescrit + interfaces conjuguées/résistance de contact + solveur couplé multi-régions + 21 cas VV enregistrés). CI exacte HEAD en cours ; aucune qualification physique revendiquée.
- Issue / PR : #1150 / #1154.
- Spécification détaillée : `docs/thermal/thermal-case-data-contract.md`.
- Sources code : `src/cfdx/physics/thermal_driver.h`.
- Tests : `tests/validation/test_thermal_driver_vv.cpp`.
- Limites : une région maillée par appel de `solve_thermal_region` (le solveur couplé A4 itère les régions maillées du cas) ; le champ de flux massique facial est figé par appel (pas de couplage pression–vitesse) et doit être discrètement conservatif ; une face d'interface appariée ne peut pas porter de flux massique prescrit non nul ; liens et nœuds 0D reportés en A5 ; conductivité tensorielle anisotrope reportée (opérateur de diffusion tensoriel) ; rho et cp figés à l'état initial déclaré de la région ; puissance de patch répartie uniformément sur l'aire intégrée.
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

- Interface entre deux régions maillées (A4) : une face d'interface appariée relie les demi-cellules propriétaires par la résistance série exacte
$
q = A_f\,\frac{T_{P_a} - T_{P_b}}{d_a/k_a + d_b/k_b + R_c},
$
  avec $R_c$ la résistance de contact (m²·K/W) : nulle pour une interface CONJUGATE (une région solide, une région fluide), strictement positive pour CONTACT_RESISTANCE. Le flux est égal et opposé de part et d'autre de l'interface.

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
- Tolérance linéaire plafonnée : le stage linéaire du solveur énergie arrête en résidu **relatif** (échelle $\lVert b\rVert_2$), alors que l'audit indépendant clôt le bilan en **watts absolus** sur une échelle brute plancher 1 W. En stationnaire couplé équilibré, les flux de bord signés s'annulent (entrant = sortant) et l'échelle d'audit retombe à ce plancher : un critère relatif seul peut laisser un résidu cellulaire absolu de $\lVert b\rVert_2 \times \varepsilon$ que l'audit rejette. Le driver plafonne donc la tolérance linéaire du stage d'énergie à $\min(\text{tolérance}, 10^{-14})$ pour que le système assemblé soit résolu en absolu au niveau que l'audit peut certifier ; les gates de changement de température et de bilan gardent la tolérance commune configurée.
- Interface conjuguée/résistance de contact (A4) : chaque face d'interface appariée est résolue par le **même** point fixe Robin discret que la convection — côté région $a$, avec $h_f = 1/(d_b/k_b + R_c)$ (résistance série **restante**, le point fixe comptant déjà la demi-cellule propriétaire) et l'ambiance égale à la température de la cellule propriétaire appariée côté $b$. Au point fixe, la conductance effective vaut exactement $A_f/(d_a/k_a + d_b/k_b + R_c)$ : la solution couplée convergée est le système série discret exact et les flux d'interface sont égaux et opposés. L'opérateur assemblé reste inchangé (FIXED_VALUE face-résolue).
- `solve_thermal_regions_coupled` (A4) : itère en Gauss-Seidel les ambiances d'interface injectées — chaque région est résolue par `solve_thermal_region` avec les conditions d'interface reconstruites depuis l'état le plus frais du voisin ; chaque sweep repart de l'état $T^n$ sauvegardé, donc un solve couplé transitoire avance l'état d'exactement un pas physique quel que soit le nombre de sweeps ; la convergence porte sur le changement des ambiances injectées ($10^{-10}$ par défaut, 200 sweeps au plus) ; l'audit final reconstruit les flux d'interface equal/opposite avec les $h$ et ambiances réellement utilisés par la dernière passe.

## 4. Traitement des frontières

Toutes les frontières de la région sont typées par patch via `ThermalBoundaryConditions` (nom de patch $\to$ condition typée). La traduction est face-résolue (`ScalarBoundaryFaceCondition`) pour chaque face de bord, reconstruite à chaque itération externe.

Contrôles fail-closed avant tout assemblage :
- patch de type EMPTY rejeté (aucune signification thermique intégrée pour ce driver) ;
- patch sans condition thermique rejeté ; condition référençant un patch inexistant rejetée ;
- valeur Dirichlet hors plage du matériau de la région rejetée ; $h > 0$ requis pour la convection ; toutes les valeurs doivent être finies ;
- état d'entrée $T$ hors plage du matériau rejeté ; température convergée hors plage rejetée après solve ;
- région PRESCRIBED_MASS_FLUX : le champ de flux massique facial doit être présent, nommé, à la taille du maillage, de dimension 1, fini et discrètement conservatif — tout défaut est rejeté avec diagnostic explicite ;
- face de bord à flux entrant ($\dot m_f < 0$) d'une région PRESCRIBED_MASS_FLUX : condition de température exigée (FIXED_TEMPERATURE ou CONVECTION), sinon rejet fail-closed (l'upwind extrapole sinon l'état du propriétaire) ;
- interface conjuguée ou de résistance de contact touchant la région : conditions d'interface par face exigées (sinon rejet avec diagnostic « conjugate interface ») ; le patch d'interface ne peut pas porter simultanément une condition limite externe ; le vecteur de conditions doit correspondre au patch ; $h > 0$ finie et ambiance finie exigées ; une face d'interface portant un flux massique prescrit non nul est rejetée ;
- solveur couplé (A4) : entrées de régions exactes et complètes (région inconnue, région manquante, champs nuls rejetés) ; appairage fail-closed — patches trouvés et non EMPTY, aucune condition externe sur un patch d'interface, mêmes nombres de faces, appariement centroïde mutuel, aires égales, normales opposées — tout défaut est rejeté avec diagnostic explicite ;
- lien 0D et matériau anisotrope : diagnostics explicites de report (A5 / opérateur tensoriel).

## 5. Pièges connus

- Ne jamais reconstruire un flux massique depuis un champ de vitesse : un fluide stagnant a un flux nul par contrat ; un fluide à advection prescrite exige un champ de faces autoritaire, figé par appel et discrètement conservatif — un champ déséquilibré crée de l'énergie fantôme (A3).
- Signe du flux massique facial : positif sort de la cellule propriétaire (convention production). Une face de bord à flux négatif est une entrée : sans condition de température, l'upwind y extrapole l'état du propriétaire au lieu d'advecter un état externe, d'où le rejet fail-closed.
- Le Robin est linéarisé en Dirichlet au point fixe : des itérations externes sont nécessaires et la convergence combine le solveur énergie, le changement de température externe et le bilan indépendant.
- Un problème stationnaire tout-Neumann est singulier (jauge) : le driver n'introduit aucune référence artificielle ni ancrage caché ; l'échec du solveur linéaire est rapporté (`converged=false`), jamais masqué.
- $\rho$ et $c_p$ sont figés à l'état initial déclaré : un cas à propriétés fortement variables doit utiliser des pas de temps adaptés ; aucun figeage silencieux rétroactif n'est appliqué en cours de transitoire.
- La puissance de patch est un flux uniforme sur l'aire intégrée du patch, pas une distribution spatiale.
- L'audit d'énergie doit réutiliser exactement les conductivités par cellule et les conditions de faces de la passe finale ; toute divergence mélange erreurs de stabilisation et erreurs de données.
- Stationnaire équilibré : les flux de bord signés entrant = sortant s'annulent, l'échelle d'audit retombe au plancher 1 W et le gate de bilan devient un critère **absolu** en watts ; l'inexactitude linéaire relative du stage d'énergie (échelle $\lVert b\rVert_2 \sim 10^3$) doit rester $\ll 10^{-9}$ W — d'où le plafond de tolérance linéaire du driver.
- La valeur Dirichlet est validée contre la plage déclarée du matériau de la région, pas contre une plage universelle.
- Le gradient imposé est le gradient normal **sortant** $q''/k_P$ : une erreur de signe inverse chauffage et refroidissement.
- Interface A4 : le coefficient de film fourni est l'inverse de la résistance série **restante** (demi-cellule voisine + résistance de contact), le point fixe Robin comptant déjà la demi-cellule propriétaire ; fournir $1/(d_a/k_a + d_b/k_b + R_c)$ doublerait la demi-cellule locale et fausserait le flux série.
- Une face d'interface appariée ne peut pas porter de flux massique prescrit non nul : l'état advecté d'une face appariée n'est pas défini, le rejet est fail-closed.
- En transitoire couplé, chaque sweep repart de $T^n$ : itérer les sweeps n'avance jamais le temps physique — un solve couplé avance l'état d'exactement un pas.

## 6. Ordre de précision attendu et validation exigée

A2+A3+A4 ne livrent aucun nouvel opérateur : l'ordre spatial et temporel est celui des fiches 11 et 14 (diffusion implicite orthogonale, upwind d'ordre 1, Euler arrière d'ordre 1 en temps). Cette fiche qualifie le **câblage** du driver (typage des frontières, Picard de conductivité, point fixe Robin, transitoire un-pas-un-appel, advection de flux massique prescrit, couplage série des interfaces), pas une nouvelle discrétisation.

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
- lavage transitoire (A3) : $a = 10$, $b = 4$, récurrence $T^{n+1} = 300 + (5/7)(T^n - 300)$ et forme fermée $300 + 50\,(5/7)^4$ sur quatre pas ;
- deux barres conjuguées stationnaires (A4) : solide $k=2$ entre 400 K et l'interface, fluide stagnant $k=0{,}5$ jusqu'à 300 K, conductances série par unité d'aire 8, 4, 1,6, 1, 2 : $q = 40$ W exact et $(A_1, A_2, B_1, B_2) = (395, 385, 360, 320)$ K ;
- saut de résistance de contact (A4) : deux solides $k=1$, $R_c = 0{,}5$ m²·K/W entre 400 et 300 K : $q = 40$ W et $(390, 370, 330, 310)$ K ;
- pas couplé transitoire (A4) : les barres conjuguées en TRANSIENT ($dt = 1$ s, $\rho\,c_p = 10$, $V = 0{,}5$ m³ par cellule, $T^0 = 350$ K) : le système $4\times 4$ Euler arrière exact donne $A_1 = 508510/1353$, $A_2 = 486830/1353$, $B_1 = 474130/1353$, $B_2 = 456710/1353$.

Bilan énergétique indépendant normalisé $< 10^{-9}$ exigé sur tous les cas physiques, et déséquilibre equal/opposite des flux d'interface $< 10^{-9}$ exigé sur les cas couplés. Cinq batteries fail-closed vérifient les diagnostics de report (A5, anisotropie), le mapping de frontières (EMPTY, patch manquant, patch inconnu), les valeurs (Dirichlet hors plage, $h \le 0$, valeur non finie, état d'entrée hors plage), le contrat du flux massique prescrit (map absente, champ introuvable, taille, valeur non finie, non-conservation discrète, entrée sans condition de température) et le contrat des interfaces conjuguées (solve_thermal_region sans conditions d'interface, vecteur inadapté au patch, flux massique prescrit sur face d'interface ; solveur couplé : entrée de région inconnue ou manquante, champs nuls, patch inconnu, condition externe sur patch d'interface, patches non conformes en nombre de faces, aires appariées inégales).

Exemptions documentées : MMS, benchmark graphique, GCI et parité CPU/MPI/GPU ne sont pas exigés pour A2+A3+A4 car l'opérateur assemblé (diffusion + upwind) est inchangé et déjà couvert par les fiches 11 et 14 ; ils redeviennent exigibles si l'opérateur change. INCONCLUSIVE n'est pas PASS. CI exacte HEAD verte requise avant fusion ; aucun statut de qualification n'est promu par cette fiche.

## 7. Statut d'implémentation

- Driver : `src/cfdx/physics/thermal_driver.h`.
- Tests : `tests/validation/test_thermal_driver_vv.cpp` ; enregistrement `CMakeLists.txt` (labels `thermal;validation`).
- Réutilise tel quel : `src/cfdx/physics/thermal_case.h` (A1), `src/cfdx/physics/thermophysical_models.h`, `src/cfdx/physics/energy_solver.h`, `src/cfdx/physics/finite_volume_transport.h`.
- Aucune tolérance élargie, aucun test désactivé, aucun oracle existant modifié, aucun registre de qualification modifié.
- Reportés : A5 (liens et nœuds 0D), opérateur de diffusion tensoriel pour l'anisotropie.
- Prochaine étape : A5, liens et nœuds 0D.

