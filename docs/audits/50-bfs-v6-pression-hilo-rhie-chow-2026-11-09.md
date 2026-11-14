# FICHE MÉTHODE — Correction de pression high/low raccordée au flux Rhie–Chow de production (BFS V6)

**En-tête d'évidence (gouvernance §16.5)**
- Scope : sous-problème de pression du solveur SIMPLE en production ; paires (high,low) physique/transport, flux Rhie–Chow compensé, raffinement globalisé. NON couvert : solveur Newton–Krylov monolithique U–p–SST, GPU, multi-CPU, restart haute précision.
- Statut : DRAFT (2026-11-09)
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (cas) + NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json (exigences)
- Issue / PR liées : #1193
- Dernier commit vérifié : 82732cbb168529edd306207bb566ad46e9b239df (2026-11-09)
- Limitations connues : chemin EXPÉRIMENTAL opt-in (`pressure_high_low_production`), désactivé par défaut ; BFS/SST complet non convergé ; x_r/h non qualifié ; restart V6 interdit (p_low non sérialisée).

## 1. Références sources (obligatoire : DOI ou ISBN)
- C. M. Rhie, W. L. Chow, « Numerical study of the turbulent flow past an airfoil with trailing edge separation », Computer Methods in Applied Mechanics and Engineering, 1983 — DOI: 10.1016/0045-7825(83)90106-9
- S. Majumdar, « Role of underrelaxation in momentum interpolation for collocated finite volume schemes », Numerical Heat Transfer, 1988 — DOI: 10.1080/10407788808913637
- S. Patankar, « Numerical Heat Transfer and Fluid Flow », Hemisphere, 1980 — ISBN 0-89116-522-3
- D. E. Knuth, « The Art of Computer Programming, Vol. 2 : Seminumerical Algorithms », 3e éd., Addison-Wesley, 1998 (somme exacte de deux flottants, algorithme TwoSum) — ISBN 0-201-89684-2
- OpenFOAM 2606, « Turbulent flow over a backward-facing step » (Driver–Seegmiller, Re_h=36000, x_r/h = 6.26 ± 0.10) — https://doc.openfoam.com/2606/examples/verification-validation/turbulent/backward-facing-step/
- PETSc SNES (acceptation d'une direction par résidu de l'opérateur physique, recherche arrière) — https://petsc.org/main/manual/snes/

## 2. Formulation continue
Équations de Navier–Stokes incompressibles stationnaires, fluide newtonien :

$$ \nabla \cdot \mathbf{u} = 0, \qquad \rho\,(\mathbf{u} \cdot \nabla)\mathbf{u} = -\nabla p + \mu\, \nabla^2 \mathbf{u} $$

Hypothèses : incompressible, newtonien, régime stationnaire, formulation colocalisée. Le point de départ du V6 n'est pas l'EDP mais le **défaut du V5** : le résidu du CSR arrondi ou des faces implicites en précision compensée était quasi nul, mais la divergence du flux Rhie–Chow complet n'était pas refermée — le flux final repartait d'une pression monolithique `double`, perdant le correcteur `low` et rompant le contrat opérateur–résidu–flux.

La méthode V6 maintient deux paires (high, low) :
1. pression **physique sous-relaxée** $p_{phys} = p_{old} + \alpha_p\, \delta p$ utilisée dans le gradient de $\mathbf{U} = H/A - (1/a_P)\, V\, \nabla p_{phys}$ ;
2. pression **complète de transport** $p_{flux} = p_{old} + \delta p$ utilisée pour construire le flux conservatif Rhie–Chow $\phi(p_{flux})$.

Le seuil de pression reste $10^{-11}$ (inchangé).

## 3. Discrétisation Volume Fini
Forme discrète : somme sur les faces $S_f \cdot \phi_f$ avec flux de masse Rhie–Chow :

$$ \phi_f = \left( \frac{H}{A} \right)_f \cdot S_f - \frac{1}{a_{P,f}}\, \nabla p_f \cdot S_f - R_f $$

Le terme de construction utilise l'opérateur de base avec $p_{cell}=0$ mais les **vraies conditions aux limites inhomogènes**, puis retranche la réponse de la paire (high, low) à **conditions homogènes** sur la correction :

$$ \phi_f(p_{high}, p_{low}) = \phi_f(p_{cell}=0, \mathrm{BC}_{reelles}) - R_f(p_{high}, p_{low}) $$

L'action $R_f$ inclut le saut de pression orthogonal et la contribution de gradient non orthogonale/différée de Rhie–Chow. La divergence signée $r_c = -\sum_f s_{cf}\, \phi_f$ est reconstruite depuis les faces en `long double`. Le Jacobien sparse V18 sert uniquement d'opérateur de direction/préconditionnement, jamais de critère d'acceptation. Le résidu CSR d'origine, le résidu des vraies faces et la continuité du flux stocké en `double` sont conservés séparément.

Raffinement : les corrections proposées par le solveur linéaire sont réacceptées seulement si la divergence Rhie–Chow **réévaluée intégralement** décroît (résidu physique) ; recherche arrière sur les pas $1, 1/2, \ldots, 1/128$ pour les directions surdimensionnées. Mise à jour `high+low` par `TwoSum` avec séparation de l'incrément sous-relaxé en deux doubles (une implémentation additionnant directement en `long double` perdait des bits de `low` à pression élevée ; défaut découvert puis éliminé par test de non-régression). L'état `low` est conservé entre itérations, snapshoté et restauré dans les quatre chemins de rejet. Le gradient de pression utilisé pour U reçoit la partie `low`.

Preuves mesurées (maillage BFS V4, 20 540 cellules, laminaire, 2 pas) : résidu face Rhie–Chow compensé $\approx 7{,}4 \times 10^{-18}$, divergence du flux stocké en double $\approx 5 \times 10^{-14}$ à $1{,}3 \times 10^{-13}$, continuité normalisée $\approx 1{-}3 \times 10^{-16}$. BFS/SST : résidu face complet après raffinement high/low $7{,}44 \times 10^{-18}$ au premier pas ; transport SST inachevé dans le budget local, `converged: false`. Ce sont des preuves de **conservation numérique du sous-problème de pression**, pas de convergence non linéaire.

## 4. Traitement des frontières
- Le flux de base est construit avec les **vraies BC inhomogènes** (y compris sorties à flux de recirculation locales) ; la réponse $R_f$ est évaluée à BC **homogènes** afin de retrancher exactement la réponse de la correction sans double-compter les contributions frontières.
- Frontières physiques inchangées : Dirichlet vitesse, sortie de pression, murs no-slip. Aucune CL modifiée par le mode V6.
- Conditions d'activation fail-closed : SIMPLE, Jacobien pression sparse V18, flux non relaxé V19, gradient GAUSS à une passe ; rejet explicite de SIMPLEC/COUPLED, skew multi-pass, multisecant, pression line-search indépendante.
- Refus explicite des plateformes sans `long double` réellement étendu (MSVC/Windows) : le mode exige une arithmétique double-double portable ; pas de fallback silencieux.

## 5. Pièges connus
- **Perte du correcteur low dans le flux final** (défaut V5) : le flux doit partir de la pression complète de transport, pas de la pression physique sous-relaxée seule.
- **Addition directe en `long double`** à pression élevée : des bits de `low` sont perdus ; la séparation TwoSum en deux doubles est obligatoire.
- **Acceptation par critère linéaire** : le Jacobien sparse V18 n'est qu'une direction ; l'acceptation exige la décroissance du résidu physique des vraies faces (principe SNES).
- **Restart** : le format DAT ne sérialise pas `p_low` ; un restart V6 haute précision est refusé en attendant un format checkpoint versionné.
- **Réduction non reproductible** : GPU/multicœur non qualifiés ; les réductions doivent rester reproductibles sans modifier le contrat physique.
- **Ne pas confondre** conservation numérique du pas de pression et convergence BFS/SST : le résidu momentum transverse reste proche de 1 au démarrage ; aucun x_r/h, Cf ou y+ n'est revendiqué. Les traitements de paroi CFDX (`wall_k_zero`, `wall_omega_sst`, facteur 6) ne sont **pas** équivalents aux lois de paroi OpenFOAM (`nutUSpaldingWallFunction`, `kqRWallFunction`, `omegaWallFunction`) ; la comparaison de réattachement sans comparabilité de modèle de paroi serait trompeuse.

## 6. Ordre de précision attendu et validation exigée
- Ordre attendu : inchangé par rapport aux fiches 03/04 (schémas non modifiés) ; le V6 n'ajoute aucune discrétisation nouvelle, il referme l'arrondi du flux existant.
- Oracles quantitatifs : aucune nouvelle tolérance ni oracle ajoutés ni élargis ; seuil de pression $10^{-11}$ inchangé.
- Validation exigée avant toute qualification (NON satisfaite à ce stade) : tous les résidus physiques (momentum, pression, continuité, SST k/omega) vérifiés positivement sur l'état accepté final ; champs k, omega, nu_t, distances de paroi et contraintes tangentielles natives enregistrés ; $u_\tau = \sqrt{|\tau_w|/\rho}$ puis $y^+$ depuis la traction physique ; $C_f = 2\tau_{w,x}/(\rho U_{ref}^2)$ avec convention d'axe et sens de normale ; zéro négatif-vers-positif unique pour $x_r/h$ en aval de la marche avec incertitude spatiale ; comparaison des profils U/Uref, Cp, budget de masse, traitement de paroi ; étude à au moins trois grilles comparables ; écart à $6{,}26 \pm 0{,}1$ reporté comme résultat de validation, pas comme seuil artificiel.
- Parallèle : parité CPU/MPI/GPU non applicable tant que le chemin reste CPU expérimental ; Règle d'Or 4 non satisfaite pour ce chemin.

## 7. Statut d'implémentation
- Fichiers source concernés : `src/cfdx/physics/steady_incompressible_solver.h` (contrôle `pressure_high_low_production`, validation fail-closed, `reconstruct_compensated_rhie_chow_pressure_step` en `std::vector<long double>` + surcharge Field conservée, `make_rhie_chow_mass_flux_high_low`, `pressure_true_face_residual`, `refine_pressure_true_face_high_low`, `gauss_gradient_pressure_high_low`, `add_pressure_high_low`, sidecar `pressure_low_state` avec snapshots/restaurations dans les quatre catch) ; `apps/production_solver/main.cpp` (lecture du raw_settings `pressure_high_low_production`, trace expérimentale) ; `CMakeLists.txt` (enregistrement `test_bfs_sst_wall_contract_v6`).
- Scripts de qualification : `scripts/qualify_bfs_sst_wall_contract_v6.py` (audit fail-closed des fonctions de paroi ; exit 2 `BLOCKED_MODEL_AND_PHYSICAL_EVIDENCE` = gate fermée, pas un échec logiciel).
- Tests exécutables : `tests/unit/test_pressure_rhie_chow_compensated.cpp` (9/9 : précision, split, vitesse, flux, signes, rejet d'une mauvaise direction, backtracking d'un pas surdimensionné, contrat), `test_pressure_conservative_face_operator` (7/7), CTest ciblés 3/3, pytest BFS/Rhie–Chow/NACA/gate SST 38 passed / 1 skipped, `tests/python/test_bfs_sst_wall_contract_v6.py` (contrat de paroi).
- Écarts connus code ↔ fiche : aucun identifié ; les six sources intégrées correspondent octet par octet au rapport V6 (fusion 3-voies sur master c4280391, blocs #1191 préservés).
- Gaps et prochaines actions : diagnostic de chaque terme de transport proche des faces de la paroi supérieure (résidu momentum transverse ~1) ; système globalisé de mise à jour U–p (puis SST) avec normes physiques et états transactionnels ; format checkpoint versionné incluant p_low ; qualification GPU/multi-CPU ; qualification SST bas-Re ou intégration de traitements de paroi équivalents OpenFOAM sur une famille de maillages adaptée.

## 8. Journal des modifications de la fiche
- 2026-11-09 — PR #1193 — création (DRAFT) à partir du rapport intégral V6 (raccordement high/low au flux Rhie–Chow de production et gate SST).
