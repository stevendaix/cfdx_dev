# FICHE MÉTHODE 10 — Thermique, CHT multi-régions et approximation de Boussinesq

**En-tête d'évidence**
- Scope : équation d'énergie en température pour fluide et solide, couplage conjugué fluide-solide (CHT) multi-régions, approximation de Boussinesq pour la flottabilité. Non couvert : rayonnement (fiche 11), changement de phase, milieux poreux.
- Statut : APPROUVÉE (2026-10-12, PR #861) — contrat mathématique approuvé ; le code existant doit être audité contre cette fiche avant toute promotion (Règle d'Or 1).
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
- Issues liées : #34 (roadmap, Phase 11) ; docs/validation/THERMOPHYSICAL_MODELS.md
- Dernier commit vérifié : f5a728f205cb45ad6186c9fab92554ae4af06bdd (master, 2026-10-12 — présence des sources et tests confirmée)
- Limitations connues : bilans d'énergie d'interface non exécutés depuis l'environnement d'audit ; statut déclaratif jusqu'à campagne sur HEAD exact.

## 1. Références sources
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI : 10.1201/9781482234213. (couplage conjugué, moyenne harmonique, linéarisation)
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN : 978-3-319-16873-9.
- J. H. Ferziger, M. Peric, *Computational Methods for Fluid Dynamics*, 3e éd., Springer, 2002 — DOI : 10.1007/978-3-642-56026-2.

## 2. Formulation continue
Équation d'énergie du fluide (température, propriétés à $\rho c_p$ constantes par région) :
$$\frac{\partial(\rho c_p T)}{\partial t}+\nabla\cdot(\rho c_p\,\mathbf{u}T)=\nabla\cdot(k\nabla T)+S_T$$
Région solide (conduction pure, $\mathbf{u}=0$) :
$$\frac{\partial(\rho c_p T)}{\partial t}=\nabla\cdot(k\nabla T)+S_T$$
Couplage conjugué à l'interface fluide-solide $\Gamma_{fs}$ : continuité de la température et du flux :
$$T_f=T_s,\qquad q_{\Gamma}=-k_f\frac{\partial T_f}{\partial n}=-k_s\frac{\partial T_s}{\partial n}$$
Flottabilité (approximation de Boussinesq, valide pour $\beta\Delta T\ll 1$) : la masse volumique reste constante partout **sauf** dans le terme gravitaire :
$$\rho=\rho_0\left[1-\beta(T-T_0)\right]\;\text{dans}\;\rho\mathbf{g}\;\text{uniquement}$$

## 3. Discrétisation Volume Fini
Forme discrète sur la cellule $P$ :
$$a_P T_P=\sum_{nb}a_{nb}T_{nb}+b,\qquad \sum_f \dot m_f c_p (\rho)_f \ldots$$
- **Flux convectif** : porte $\rho c_p$ — leçon d'audit #419 **B2** : quand $c_p$ varie (ou par prudence systématique), le produit $\rho c_p$ doit être interpolé et appliqué au flux convectif ; jamais de flux convectif en $\rho$ seul additionné d'un $c_p$ séparé après coup. Le bilan d'énergie doit être vérifié exact sur ce point.
- **Flux diffusif** : conductivité de face $k_f$ par **moyenne harmonique** des valeurs de cellules (Patankar) pour les coefficients discontinus à l'interface fluide-solide ; coefficients assemblés dans la matrice $A$.
- **Interface CHT conservatif** : le flux d'interface $q_{\Gamma}$ est calculé **une seule fois** et débité/credité symétriquement aux deux régions ($-q_{\Gamma}$ d'un côté, $+q_{\Gamma}$ de l'autre) — conservation machine de l'énergie d'interface ; jamais deux flux calculés indépendamment.
- **Terme source** $S_T$ : linéarisé en $S_T=S^{exp}+S^{lin}T_P$ quand il dépend de T (source radiatif, voir fiche 11) ; signe vérifié par bilan.

## 4. Traitement des frontières
- Dirichlet $T=T_w$ : modifie $a_P$ et $b$ explicitement (valeur limite) ; Neumann $-k\partial T/\partial n=q_w$ : flux imposé entrant dans $b$ ; Robin $-k\partial T/\partial n=h(T-T_\infty)$ : coefficients répartis entre $a_P$ et $b$.
- Interface CHT : Robin croisé conservatif entre deux régions (élimination du flux d'interface par substitution) — ou équivalent matriciel exact documenté.
- Symétrie/adiabatique : flux nul. EMPTY (2D) : aucun flux.
- **Gauge de température** : un cas tout-Neumann (adiabatique complet, seul régime permanent) a un null-space de T (solution définie à une constante près) — un point de référence doit être fixé, comme pour la pression (fiche 04) ; le test doit l'établir explicitement.

## 5. Pièges connus
- **B2** : $\rho c_p$ dans la convection d'énergie — pas $\rho$ seul.
- Bilan d'énergie d'interface : l'énergie qui quitte le fluide doit entrer dans le solide à la tolérance machine — sinon le couplage fuit.
- Boussinesq : seulement pour $\Delta T$ petits ; hors validité, utiliser la flottabilité pleine (variable $\rho$) et documenter le choix.
- $c_p$ variable : la formulation température peut perdre la conservativité — la formulation enthalpie est l'alternative documentée.
- Moyenne harmonique obligatoire pour $k$ à l'interface (l'arithmétique sur-lisse les discontinuités et casse le flux exact).
- Régime permanent tout-Neumann : null-space de T — point de référence explicite.
- Rayonnement couplé : le source radiatif doit être linéarisé en Newton pour rester implicite (fiche 11, leçon B3) — jamais de source purement explicite instable.

## 6. Ordre de précision attendu et validation exigée
- Ordre attendu : $O(h^2)$ pour la diffusion (faces orthogonales), $O(h)$ à $O(h^2)$ selon le schéma convectif déclaré ; temps $O(\Delta t)$ implicite Euler / $O(\Delta t^2)$ selon le schéma. MMS sur au moins 3 maillages (pente log-log).
- GCI sur les QoI thermiques ($Nu$, $T_{max}$, flux d'interface).
- Oracles quantitatifs : mur composite 1D analytique (résistance en série), plaque en conduction transitoire (solution de Fourier), bilans d'énergie machine à l'interface fluide-solide et en volume ; comparaison OpenFOAM sur même maillage (fiche 08).
- Conservation : énergie totale (fluide + solide) à la tolérance machine ; bornitude T physique (pas d'oscillations non physiques) avec télémétrie.
- Parallèle : parité CPU/MPI à 1e-12 sur le même maillage (fiche 12).

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : src/cfdx/physics/thermal.h, energy_solver.h, thermophysical_models.h, cht_solver.h, boussinesq.h.
- Tests exécutables : tests/validation/test_cht_validation.cpp, test_thermal_vv.cpp, test_thermal_challenge_vv.cpp, test_thermal_radiation_model_matrix.cpp.
- Écarts connus code <-> fiche : à auditer ligne à ligne lors de la Phase 3 (leçon B2 = premier point de contrôle).
- Gaps et prochaines actions : bilan d'énergie d'interface multi-région comme oracle machine ; campagne GCI ; comparaison OpenFOAM même maillage ; promotion BLOCKED derrière Phase 1.
- Étape 4 du plan #885 (couplage inverse énergie→fluide, flottabilité Boussinesq) : `src/cfdx/physics/production_buoyancy_stage.h` (test `tests/validation/test_production_buoyancy_outer_loop_vv.cpp`) — le contrat Boussinesq de la fiche (§2 : ρ = ρ₀[1−β(T−T₀)] dans le seul terme gravitaire) est implémenté comme source de quantité de mouvement per-cellule a_i = −β(T_i−T_ref)·g_vec (m/s², la part constante ρ₀·g_vec absorbée dans la pression modifiée), injectée par le canal cellwise du porteur (fiche 04 §7) ; déclarations fail-closed au chargement du cas (exactement 1 matériau, β strictement positif via `extra_properties["thermal_expansion_coefficient"]`, `gravity_vector` strictement 3 composantes finies non toutes nulles, `ref_temperature` finie > 0, pas de `body_force` uniforme déclaré) ; boucle externe flow-énergie (sous-relaxation fixe ω = 0,25 de la température vue par la source, budget 200 itérations fail-closed, une résolution complète du porteur par itération externe, warm-start) dont le verdict est le point fixe couplé exact — stationnarité T/U/p de la réponse NON relaxée ≤ 1e-10 normalisée par l'échelle brute, critères SÉPARÉS, INCONCLUSIVE n'est jamais PASS — avec l'étage énergie (fiche 14) et, si déclaré, le point fixe énergie-rayonnement S2S de l'étape 3 (fiche 11) comme opérateur de réponse interne ; oracles VV : équivalence bit-exacte source uniforme ↔ `body_force` uniforme, identité T = T_ref (source nul), signes (paroi chaude monte), replay indépendant du point fixe, catalogue fail-closed complet.
- Convergence-at-entry (étape 4, correctif de la PR #1016) : lorsque le porteur redémarré à chaud est classé `STAGNATED` alors que son résidu est déjà au plancher machine — résidus de pression/continuité/flux corrigé ≤ 1e-14, résidu d'équation momentum RHS-normalisé ≤ 1e-12 — la boucle flottabilité accepte le porteur d'entrée (restauration U/p, conservation du porteur précédemment accepté) ; comptabilité via `flow_converged_at_entry` sur l'itération et le résultat ; aucune tolérance/gate/oracle/budget de la fiche modifiée.

## 8. Journal des modifications de la fiche
- 2026-10-12 — création (Phase 3 du Plan Directeur, PR empilée sur governance/constitution-technique).
- 2026-10-12 — approbation du contrat (PR d'approbation #861) ; l'audit du code contre cette fiche reste requis avant toute promotion (Règle d'Or 1).
- 2026-10-21 — étape 4 du plan #885 : stage de production flottabilité `production_buoyancy_stage.h` + routage `buoyant` dans `cfdx_production_solver` — le contrat Boussinesq de §2 (source a = −β(T−T_ref)·g_vec per-cellule, part constante absorbée dans la pression modifiée, valide βΔT ≪ 1) est implémenté via le canal cellwise du porteur (fiche 04 §7) dans une boucle externe flow-énergie (ω = 0,25 sur la température vue par la source, budget 200 fail-closed, critères bloquants séparés sur la réponse non relaxée ≤ 1e-10 : stationnarité T/U/p normalisée ; verdict = point fixe couplé exact, INCONCLUSIVE n'est jamais PASS) ; oracles : équivalence bit-exacte avec le `body_force` uniforme équivalent, identité T = T_ref, signes, replay du point fixe, catalogue fail-closed (test `test_production_buoyancy_outer_loop_vv`) ; aucune tolérance/gate/oracle/budget de la fiche modifiée ; la flottabilité pleine ρ variable reste hors validité documentée (§5) et exigerait une nouvelle fiche.
- 2026-10-21 — correctif de la PR #1016 (échecs A/C/D de `test_production_buoyancy_outer_loop_vv`) : convergence-at-entry ajoutée à la boucle externe flottabilité — acceptation du porteur warm classé `STAGNATED` au plancher machine (résidus pression/continuité/flux corrigé ≤ 1e-14, momentum RHS-normalisé ≤ 1e-12), restauration U/p d'entrée, conservation du porteur précédemment accepté, comptabilité `flow_converged_at_entry` ; cas A à la 3e itération (itérations 1-2 sans at-entry, itération 3 at-entry, porteur inchangé vs itération 2) ; échecs C/D relèvent de la calibration de la gate locale de l'étage énergie (fiche 14 §7/§8) ; aucune tolérance/gate/oracle/budget de la fiche modifiée.
