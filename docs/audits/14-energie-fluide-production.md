# FICHE MÉTHODE 14 — Énergie fluide en écoulement gelé : étage de production couplé fluide→énergie

**En-tête d'évidence (gouvernance §16.5)**
- Scope : transport d'énergie stationnaire `div(ṁ cp T) − div(k grad T) = q'''` sur le flux de masse **autoritaire** du solveur incompressible convergé (couplage unidirectionnel « frozen flow », `energy_model = "conjugate"` dans `cfdx_production_solver`). Couvert : puissance volumique (W/m^3), Dirichlet entrée/paroi, flux de paroi imposé (`heat_flux`), adiabatique, audit de conservation indépendant. Non couvert : flottabilité (`buoyant` rejeté fail-closed), couplage inverse énergie→fluide, CHT multi-régions fluide/solide (fiche 10), transitoire, itération externe (déléguée aux étages turbulence/S2S de #885), correction non orthogonale (fiche 06).
- Statut : APPROUVÉE (2026-10-13, PR #888) — contrat rédigé sur l'implémentation du driver T2 ; l'exécution réelle de la campagne T2 sur un HEAD exact mergé reste à la charge du propriétaire avant toute promotion (gouvernance #118).
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (cas `T2-FLUID-ENERGY`) + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
- Issues / PR liées : #880 (T1), #885 (étage 2), #118
- Limitations connues : couplage unidirectionnel gelé (exact sans flottabilité, voir §2) ; stationnaire ; une seule propriété fluide par cas (exactement 1 matériau) ; `0` = « non déclaré » pour `temperature`/`wall_temperature` ; ≥ 1 patch Dirichlet requis (sinon rejet fail-closed) ; ordre 1 en convection (l'oracle est discret, pas continu).

## 1. Références sources
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI: 10.1201/9781482234213. (traitement upwind de la convection scalaire, half-cell Dirichlet, convection gelée)
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN: 978-3-319-16873-9. (assemblage transport scalaire convectif, correction bornée `div(φ)T`)
- J. H. Ferziger, M. Peric, *Computational Methods for Fluid Dynamics*, 3e éd., Springer, 2002 — DOI: 10.1007/978-3-642-56026-2. (conservation discrète, audits indépendants de la matrice)

## 2. Formulation continue
Régime stationnaire, fluide incompressible à propriétés constantes, vitesse $\mathbf{u}$ **gelée** (champ porteur convergé, indépendant de $T$) :
$$\rho c_p\,\nabla\cdot(\mathbf{u}T) - \nabla\cdot(k\nabla T) = q''' \quad \text{dans } \Omega,\qquad k>0,\ c_p>0$$
Le flux d'enthalpie autoritaire par face est $\dot m_f c_p$ avec $\dot m_f = \rho\,\mathbf{u}\cdot\mathbf{S}_f$ le flux de masse conservatif du solveur incompressible. Hypothèses : convection forcée **sans flottabilité** (pas de terme $\rho\mathbf{g}(T-T_0)$) ; $q'''$ volumique uniforme fini (nul autorisé) ; propriétés constantes (pas de $k(T)$). Sous ces hypothèses le problème est **linéaire** et le couplage unidirectionnel gelé est le point fixe exact du problème couplé complet : itérer jusqu'à convergence externe ne changerait rien (le champ porteur ne dépend pas de $T$). Toute itération externe est donc différée aux étages où elle porte physiquement (flottabilité, turbulence, S2S), et non simulée ici artificiellement.

## 3. Discrétisation Volume Fini
Schéma volumes finis two-point sur cellules $P$, convection upwind ordre 1 sur le flux de masse autoritaire :
$$a_P T_P = \sum_{nb} a_{nb} T_{nb} + b,\qquad F_{h,f} = \dot m_f c_p \ \text{(flux d'enthalpie, W/K)}$$
- **Faces internes** : $D_f = k A_f/d_f$ avec $d_f=\|\mathbf{C}_n-\mathbf{C}_o\|$ ; $a_{on}=D_f+\max(F_{h,f},0)$, $a_{no}=D_f+\max(-F_{h,f},0)$ (upwind conservatif) ; débit symétrique `$-a_{no}$`/`$-a_{on}$` en hors-diagonale (conservation machine entre cellules adjacentes).
- **Convection bornée** : en régime stationnaire le terme $-\mathrm{div}(F_h)\,T_P$ est retranché de la diagonale ; il s'annule exactement quand le flux porteur est à divergence nulle discrète et préserve les constantes sinon.
- **Terme source** : $b_P \mathrel{+}= q''' V_P$ (explicite, indépendant de $T$).
- **Linéarité / point fixe** : l'équation assemblée est linéaire pour flux, source et $k$ fixés → relaxation $=1$ et un solveur linéaire pleinement convergé **sont** le point fixe exact ; toute relaxation $<1$ briserait cette identité (l'état « convergé » ne résoudrait plus le système assemblé).
- **Ordre** : ordre 1 en convection, ordre 2 en diffusion ; la solution est exacte pour l'oracle **discret** reconstruit à la main (le gate de vérification est l'oracle discret, jamais une solution continue à la précision machine).

## 4. Traitement des frontières
- **Dirichlet** : `WALL_THERMAL` avec `wall_temperature` ≠ 0, ou `INLET`/`FIXED` avec `temperature` ≠ 0 → conductance half-cell $G_b = k_o A_f/d_b$ ($d_b$ distance normale centre→face) ; $a_P \mathrel{+}= G_b + \max(F_{h,f},0)$, $b \mathrel{+}= (G_b + \max(-F_{h,f},0))\,T_b$. `0` signifie « non déclaré » (convention étage 1 / fiche 13), pas un Dirichlet 0 K.
- **Flux imposé** (`raw_params["heat_flux"]` $q''$ W/m^2, positif = entrant) → gradient fixé $g = q''/k_o$ : $a_P \mathrel{+}= F_{h,f}$, $b \mathrel{+}= k_o A_f g$ (le terme $-F_{h,f} g d_b$ s'annule pour un flux de paroi nul).
- **Adiabatique / sortie** (`ZERO_GRADIENT`) : convection seule, $a_P \mathrel{+}= F_{h,f}$ (bornée).
- **EMPTY** (2D/axisym.) : aucune contribution, patch ignoré.
- **Null-space (jauge de température)** : un système tout-Neumann a un null-space constant → rejet fail-closed (≥ 1 patch Dirichlet exigé), jamais résolu silencieusement.

## 5. Pièges connus
- `energy_model = "buoyant"` (ou toute valeur inconnue) est **rejeté** : l'étage ne doit jamais dégrader silencieusement en isotherme — le refus porte la sémantique, pas un avertissement. Exception : `energy_model = "compressible"` (produit par la conversion SU2 des cas Euler / Navier-Stokes compressibles, cf. `case_schema.h` / `adapters/su2.py`) est accepté comme déclaration valide mais ne déclenche **jamais** l'étage (seul `conjugate` le fait) — le solveur de production exécute alors l'écoulement seul.
- Le flux consommé doit être le flux de masse **autoritaire conservatif** (`IncompressibleSolveResult.authoritative_mass_flux`), jamais un flux reconstruit depuis la vitesse cellulaire : la conservativité discrète du porteur est ce qui rend le bilan d'énergie indépendant exact à la machine ; un flux reconstruit la casse.
- Flux non fini, température initiale non finie, champ de mauvaise taille → rejet fail-closed avant l'assemblage.
- Tout-Neumann → singulier → rejet (miroir fiche 13) ; 0 = « non déclaré » pour `wall_temperature`/`temperature`.
- Relaxation ≠ 1 : voir §3 — interdit par construction dans cet étage (contrôles gelés à 1.0).
- Normalisation du bilan d'énergie par le **flux brut** (fiche 13), pas le flux net.
- Ordre 1 : la solution continue n'est approchée qu'à $O(h)$ en convection ; ne jamais gater sur un oracle continu à la précision machine — l'oracle exact est le système discret reconstruit indépendamment (tridiagonal, Thomas).
- L'étage ne s'exécute que sur un écoulement porteur **convergé** ; un résultat d'énergie porté par un écoulement divergé n'est pas une évidence (exit code global = fluide **et** énergie).

## 6. Ordre de précision attendu et validation exigée
- Ordre attendu : $O(h)$ convection (upwind ordre 1), $O(h^2)$ diffusion ; exact pour l'oracle discret à raffinement fixé.
- Vérification implémentée (`test_production_energy_stage_vv`, labels qualification;validation;phase13) : plug flow (parois glissantes, entrée uniforme) avec gates |u−u_in| < 1e-9 et flux de masse autoritaire uniforme et signé vers l'extérieur ; 3 cas physiques (puissance volumique + entrée Dirichlet ; parois chauffées Dirichlet ; flux de paroi imposé) comparés à un oracle tridiagonal **reconstruit à la main** (formules de §3–§4, Thomas) à < 1e-9 ; fermeture énergétique indépendante < 1e-12 ; identités de flux par patch (parois adiabatiques nulles ; flux imposé = q''·aire) ; principe du maximum discret (bornes Dirichlet) ; montée bulk $q'''V/(\dot m c_p)$ en sanity physique ; ~15 rejets fail-closed (modèle, matériaux, propriétés, volumetric_power, heat_flux, tout-Neumann, flux non fini, tailles).
- Campagne T2 (à exécuter par le propriétaire sur un HEAD exact mergé) : MMS 3 niveaux sur cas `.cfdx.h5` réels via `cfdx_production_solver` (`energy_model=conjugate`), GCI sur les QoI (flux par patch, $T_{max}$, montée bulk), exécution réelle + accord quantitatif archivés — aucun PASS sans cette évidence (gouvernance #118).
- Conservation : bilan d'énergie reconstruit **indépendamment** de la matrice assemblée (`reconstruct_energy_balance`, puis ×$c_p$) ; tolérance machine sur le flux brut.
- Parallèle : N/A à ce stade (solveur séquentiel ; fiche 12 quand le chemin sera parallélisé).

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : src/cfdx/physics/production_energy_stage.h (étage), apps/production_solver/main.cpp (app de production `cfdx_production_solver` : validation fail-closed du modèle, étage après un écoulement convergé, artéfacts `energy_convergence.json` / `energy_summary.json` / `energy_result_0.vtu`, exit code global = fluide ∧ énergie).
- Tests exécutables : tests/validation/test_production_energy_stage_vv.cpp (`test_production_energy_stage_vv`).
- Conventions utilisateur : `energy_model = "conjugate"` déclenche l'étage ; exactement 1 matériau fluide (`density`, `thermal_conductivity`, `specific_heat` > 0 finis) ; `extra_properties["volumetric_power"]` (W/m^3, défaut 0) ; BC température : `wall_temperature` (WALL_THERMAL) ou `temperature` (INLET/FIXED) ≠ 0 → Dirichlet, sinon `raw_params["heat_flux"]` (W/m^2, positif entrant) → gradient fixé, sinon adiabatique ; 0 = « non déclaré » ; ≥ 1 Dirichlet requis ; température initiale = `initial_condition.temperature` sinon `ref_temperature`.
- Écarts connus code ↔ fiche : aucun identifié à la rédaction (la fiche décrit l'implémentation de la PR #888) ; audit ligne à ligne à refaire à chaque évolution (Règle d'Or 1).
- Gaps et prochaines actions : campagne T2 sur HEAD exact (MMS 3 niveaux + GCI + exécution réelle archivée) ; itération externe et couplage inverse différés aux étages turbulence/S2S de #885 ; CHT fluide/solide multi-régions (fiche 10) ; transitoire et non orthogonale hors scope présent.

## 8. Journal des modifications de la fiche
- 2026-10-13 — création (PR #888, T2 #885 étage 2) : contrat rédigé sur l'implémentation, approbation par fusion de la PR par le propriétaire ; campagne d'exécution T2 reste à produire avant toute promotion.
