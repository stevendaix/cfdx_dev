# FICHE MÉTHODE — Continuation adaptative NACA / robustesse non linéaire

**En-tête d'évidence (gouvernance §16.5)**
- Scope : politique de pilotage expérimentale du solveur stationnaire (adaptation d'un indicateur CFL, activation progressive de corrections non orthogonales et d'une reconstruction d'ordre deux). Cette fiche couvre le contrôleur d'état uniquement ; elle ne décrit pas encore son intégration dans la boucle CFDX.
- Statut : DRAFT (2026-11-01) — non approuvée pour qualification de production.
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (cas) + NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json (exigences).
- Issue / PR liées : PR #1160.
- Dernier commit vérifié : cff614fc7e08532f1b19d6a3001058adc36631cb (2026-11-01).
- Limitations connues : aucun raccordement à solve_steady_incompressible ; aucun CFL pseudo-temporel réellement consommé ; les facteurs non orthogonal et second ordre sont seulement émis par le contrôleur. Aucune preuve de convergence N128/N256.

## 1. Références sources (obligatoire : DOI ou ISBN)

- C. T. Kelley et D. E. Keyes, « Convergence Analysis of Pseudo-Transient Continuation », SIAM Journal on Numerical Analysis, 35(2), 508–523, 1998 — DOI : 10.1137/S0036142996304796. Cette référence justifie l'intérêt de la continuation pseudo-transitoire pour globaliser la résolution d'états stationnaires ; elle ne valide pas les seuils heuristiques spécifiques codés ici.
- H. K. Versteeg et W. Malalasekera, *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*, 2e édition, Pearson, 2007 — ISBN : 978-0-13-127498-3. Référence de contexte pour les résidus de volumes finis et les corrections différées ; elle ne valide pas cette politique d'activation.

## 2. Formulation continue

Le problème stationnaire visé reste le système incompressible (fluide newtonien, densité constante) :

$$
\nabla \cdot \mathbf{u}=0,
\qquad
\rho(\mathbf{u}\cdot\nabla)\mathbf{u}
=-\nabla p+\nabla\cdot\boldsymbol{\tau}+\rho\mathbf{f}.
$$

Le contrôleur ne modifie pas ces EDP. Il pilote une suite d'approximations discrètes du problème complet. Soit $R_h(q)=0$ le résidu non linéaire discret complet, avec $q=(\mathbf{u},p)$. Pour deux états acceptés consécutifs, l'indicateur est :

$$
\gamma_k=\frac{\|R_h(q_k)\|_{\mathrm{phys}}}
{\|R_h(q_{k-1})\|_{\mathrm{phys}}}.
$$

La norme doit être le résidu physique complet, normalisé par des échelles déclarées et calculé après acceptation de la mise à jour. Le contrôleur augmente le CFL si $\gamma_k<\gamma_{good}$, le réduit si $\gamma_k>\gamma_{bad}$, et conserve sa valeur sinon. Les seuils sont des paramètres heuristiques à qualifier, pas des constantes universelles.

## 3. Discrétisation Volume Fini

La discrétisation spatiale est celle du solveur CFDX existant et n'est pas modifiée par cette brique. On note symboliquement le système assemblé :

$$
R_h(q;\theta)=A_h(q;\theta)q-b_h(q;\theta)=0,
$$

où $\theta=(c_{\mathrm{CFL}},\omega_{\mathrm{nonorth}},\omega_2)$ représente les contrôles émis. Le contrôleur propose les transitions suivantes :

- démarrage robuste : $\omega_{\mathrm{nonorth}}=0$, $\omega_2=0$ ;
- après $N_g$ diminutions successives suffisantes du résidu, incrément de $\omega_{\mathrm{nonorth}}$ de 0,25 jusqu'à 1 ;
- après $N_{no}$ diminutions suffisantes en correction non orthogonale complète, incrément de $\omega_2$ de 0,25 jusqu'à 1 ;
- maintien en ordre complet pendant au moins $N_f$ états acceptés avant d'autoriser le drapeau de préparation à la qualification, sous réserve du critère physique.

Le code de cette PR ne construit ni ne mélange les flux de face : il expose uniquement ces valeurs dans `Controls`. Tant que l'intégration solveur n'est pas réalisée, $\theta$ n'agit pas sur $A_h$ ni sur $b_h$. Aucun ordre de précision ne peut donc être attribué à ce contrôleur isolé.

## 4. Traitement des frontières

Aucune condition limite n'est modifiée par le contrôleur. Les conditions WALL, EMPTY, périodiques, vitesse imposée et pression imposée restent sous la responsabilité du solveur stationnaire existant. L'intégration future devra prouver que le facteur de correction non orthogonale ne supprime pas les contributions de frontière, et que le passage d'ordre ne modifie pas les flux imposés ni les bilans conservatifs.

## 5. Pièges connus

- Un CFL numérique n'est utile que s'il pilote effectivement un pseudo-pas de temps ou une stratégie de relaxation documentée ; la valeur retournée seule n'accélère pas un solveur stationnaire.
- Un facteur de correction non orthogonale ne peut pas être appliqué comme simple multiplicateur global sans vérifier les flux différés, les contributions de frontière et la conservation face par face.
- Un mélange d'ordre doit être borné et conservatif ; la borne locale ne prouve pas à elle seule la stabilité ni la précision.
- Les résidus doivent être vrais, complets et calculés sur un état accepté. Les résidus linéaires ou des composantes isolées ne sont pas des substituts.
- Une tentative rejetée doit restaurer exactement les champs et l'historique non linéaire. Le contrôleur seul ne garantit pas ce rollback.
- Le drapeau `qualification_ready` est un indicateur de préparation, pas une certification indépendante ; les oracles physiques et les bilans restent obligatoires.

## 6. Ordre de précision attendu et validation exigée

- Ordre attendu : sans objet pour le contrôleur isolé. Pour le solveur intégré, ordre spatial à mesurer séparément sur au moins trois maillages ; aucune pente ne doit être déduite des tests unitaires de cette PR.
- Oracles quantitatifs : NACA0012 Re=1000, alpha=0, comparaison de $C_l$, $C_d$, $C_p$ et $C_f$ avec les références pré-déclarées du dossier de qualification ; comparer les mêmes maillages et la même géométrie entre stratégie de référence et stratégie adaptative.
- GCI : recalcul obligatoire sur N64/N128/N256 issus d'une même famille de maillages, après convergence effective des trois calculs.
- Conservation : continuité globale et locale, déséquilibre massique, flux de paroi et efforts ; aucune tolérance existante ne doit être relâchée.
- Robustesse : enregistrer à chaque itération le résidu physique complet, le ratio, le CFL proposé et consommé, les facteurs de correction réellement appliqués, les rejets/rollbacks et le temps.
- Campagne comparative : stratégie actuelle contre stratégie adaptative, paramètres et tolérances physiques identiques ; inclure N64, N128 et N256. Une absence de convergence doit rester un échec explicite.

## 7. Statut d'implémentation

- Fichiers source concernés : `src/cfdx/physics/n256_continuation.h`.
- Tests exécutables : `tests/unit/test_n256_continuation.cpp`, `tests/python/test_compare_naca_continuation.py`.
- Écarts connus code ↔ fiche : contrôleur autonome non appelé par `solve_steady_incompressible` ; aucun résidu de production n'est encore fourni ; aucune des trois commandes de contrôle n'est consommée par l'assemblage.
- Gaps et prochaines actions : intégrer avec transaction/rollback dans la boucle de solveur ; identifier le sens physique du CFL pour le solveur stationnaire ; implémenter une modulation conservative des corrections non orthogonales et de convection ; ajouter tests d'intégration ; lancer les campagnes A/B N64/N128/N256 et recalculer GCI. Ne pas fusionner comme correctif de convergence avant ces preuves.

## 8. Journal des modifications de la fiche

- 2026-11-01 — PR #1160 — création de la fiche DRAFT accompagnant le contrôleur expérimental.
