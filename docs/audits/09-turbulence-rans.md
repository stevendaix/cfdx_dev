# FICHE MÉTHODE 09 — Turbulence RANS (k-epsilon, k-omega, SST, Spalart-Allmaras) et lois de paroi

**En-tête d'évidence**
- Scope : modèles RANS à viscosité tourbillonnaire (k-epsilon standard, k-omega, SST de Menter), modèle à une équation de Spalart-Allmaras, lois de paroi (couche logarithmique), distance de paroi. Non couvert : LES/hybrides, RSM, transition, correction de compressibilité (à documenter séparément si activée).
- Statut : contrat mathématique approuvé (PR #861). Le runtime SST transactionnel est présent dans le solveur incompressible de production sur le `master` courant ; cela ne ferme pas la qualification physique RANS, qui reste bloquée tant que les cas de référence de la section 6 ne sont pas validés.
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (cas TUR-CHANNEL, TUR-FLATPLATE, TUR-NACA0012) + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
- Issues liées : #34 (roadmap, Phase 10) ; docs/validation/TURBULENCE_QUALIFICATION_MATRIX.md
- Dernier commit vérifié : f5a728f205cb45ad6186c9fab92554ae4af06bdd (master, 2026-10-12 — présence des sources et tests confirmée)
- Limitations connues : les campagnes TUR-* n'ont pas pu être exécutées depuis l'environnement d'audit ; les statuts de test restent déclaratifs jusqu'à exécution sur un HEAD exact.

## 1. Références sources
- W. P. Jones, B. E. Launder, « The prediction of laminarization with a two-equation model of turbulence », *Int. J. Heat Mass Transfer* 15:301-314, 1972. (k-epsilon original ; DOI non confirmé — citation bibliographique complète exigée)
- B. E. Launder, D. B. Spalding, « The numerical computation of turbulent flows », *Comput. Methods Appl. Mech. Engrg.* 3(2):269-289, 1974 — DOI : 10.1016/0045-7825(74)90029-2. (constantes k-epsilon standard et lois de paroi)
- D. C. Wilcox, « Reassessment of the scale-determining equation for advanced turbulence models », *AIAA Journal* 26(11):1299-1310, 1988 — DOI : 10.2514/3.10041. (k-omega)
- F. R. Menter, « Two-equation eddy-viscosity turbulence models for engineering applications », *AIAA Journal* 32(8):1598-1605, 1994 — DOI : 10.2514/3.12149. (SST)
- P. R. Spalart, S. R. Allmaras, « A One-Equation Turbulence Model for Aerodynamic Flows », *AIAA Paper* 92-0439, 30th Aerospace Sciences Meeting, 1992 — DOI : 10.2514/6.1992-439.
- D. C. Wilcox, *Turbulence Modeling for CFD*, 3e éd., DCW Industries, 2006. (référence de synthèse)

## 2. Formulation continue
Moyenne de Reynolds de Navier-Stokes incompressible, hypothèse de Boussinesq pour le tenseur de Reynolds :

$$-\langle u_i' u_j' \rangle = \nu_t\left(\frac{\partial U_i}{\partial x_j}+\frac{\partial U_j}{\partial x_i}\right)-\frac{2}{3}k\delta_{ij}, \qquad \nu_t = C_\mu\, k^2/\varepsilon \;\;(k\text{-}\varepsilon)$$

Équations transportées (k-epsilon) :
$$\frac{\partial k}{\partial t}+\nabla\cdot(\mathbf{u}k)=P_k-\varepsilon+\nabla\cdot\left[\left(\nu+\frac{\nu_t}{\sigma_k}\right)\nabla k\right]$$
$$\frac{\partial \varepsilon}{\partial t}+\nabla\cdot(\mathbf{u}\varepsilon)=C_{\varepsilon 1}\frac{\varepsilon}{k}P_k-C_{\varepsilon 2}\frac{\varepsilon^2}{k}+\nabla\cdot\left[\left(\nu+\frac{\nu_t}{\sigma_\varepsilon}\right)\nabla \varepsilon\right]$$
avec $P_k=2\nu_t\,\mathbf{S}:\mathbf{S}$, $\mathbf{S}=\tfrac{1}{2}(\nabla\mathbf{u}+\nabla\mathbf{u}^T)$ ; constantes Launder-Spalding : $C_\mu=0.09$, $C_{\varepsilon 1}=1.44$, $C_{\varepsilon 2}=1.92$, $\sigma_k=1.0$, $\sigma_\varepsilon=1.3$.

k-omega (Wilcox 1988) : $\nu_t=k/\omega$, même équation k avec $\sigma_k=0.5$, et
$$\frac{\partial \omega}{\partial t}+\nabla\cdot(\mathbf{u}\omega)=\alpha\frac{\omega}{k}P_k-\beta\omega^2+\nabla\cdot\left[\left(\nu+\frac{\nu_t}{\sigma_\omega}\right)\nabla \omega\right]$$

SST (Menter 1994) : mélange de k-omega près de la paroi et k-epsilon transformé à distance via $F_1$ ; limiteur :
$\nu_t=\frac{a_1 k}{\max(a_1\omega,\; S\,F_2)},\qquad a_1=0.31$

Pour SST-2003, la production utilisée dans les transports de `k` et `omega` doit être cohérente :
$\widetilde P_k=\min\left(\rho\nu_t S^2,\;10\beta^*\rho k\omega\right),\qquad P_\omega=\begin{cases}\gamma\widetilde P_k/\nu_t,&\nu_t>0,\\0,&\nu_t=0,\end{cases}$
avec $\gamma=F_1\gamma_1+(1-F_1)\gamma_2$. La branche $\nu_t=0$ évite une division indéfinie ; les tests algébriques ne remplacent pas la qualification sur des références DNS/expérimentales.

Spalart-Allmaras (variable de travail $\tilde{\nu}$, $\nu_t=\tilde{\nu} f_{v1}$, $f_{v1}=\chi^3/(\chi^3+C_{v1}^3)$, $\chi=\tilde{\nu}/\nu$) :
$$\frac{\partial \tilde{\nu}}{\partial t}+\nabla\cdot(\mathbf{u}\tilde{\nu})=c_{b1}\tilde{S}\tilde{\nu}+\nabla\cdot\left[(\nu+\tilde{\nu})\nabla\tilde{\nu}\right]+c_{b2}|\nabla\tilde{\nu}|^2-c_{w1}f_w\left(\frac{\tilde{\nu}}{d}\right)^2$$

## 3. Discrétisation Volume Fini
Forme discrète pour chaque scalaire turbulent $\phi\in\{k,\varepsilon,\omega,\tilde{\nu}\}$ sur la cellule $P$ :
$$a_P\phi_P=\sum_{nb}a_{nb}\phi_{nb}+b,\qquad a_P=\sum_{nb}a_{nb}+\frac{\rho V_P}{\Delta t}\;(+\,S^{lin}_{P})$$
- **Convection** : somme sur les faces $\sum_f \dot m_f \phi_f$ avec schémas bornés (UPWIND au minimum ; convection bornée obligatoire pour k, epsilon, omega, nu-tilde — pas d'overshoot non physique, Règle d'Or 3).
- **Diffusion turbulente** : les coefficients de face utilisent $\nu_{eff}=\nu+\nu_t$ (et $\sigma$ du modèle) et sont **assemblés dans la matrice $A$** — leçon d'audit #419 B1 : jamais de diffusion turbulente reportée au second membre $b$, sous peine de perdre la conservativité et la stabilité.
- **Production** $P_k=2\nu_t\mathbf{S}:\mathbf{S}$ : linéarisée quand possible (traiter $\nu_t$ explicite de l'itération externe), **bornée** (limiteur realizable ou Kato-Launder — choix documenté par modèle) ; jamais de source non bornée entrant dans la matrice.
- **Destruction** : linéarisation de Patankar implicite ($-C\phi_P^{k-1}\phi_P$) pour la robustesse ; le signe de chaque contribution doit être vérifié par le test de bilan.
- **Bornitude** : $k\ge 0$, $\varepsilon\ge 0$, $\omega\ge 0$, $\tilde{\nu}\ge 0$ vérifiés **avant tout clipping**, avec télémétrie de clipping (Règle d'Or 3).

## 4. Traitement des frontières
- **Paroi résolue** ($y^+\lesssim 1$) : $k=0$ (Dirichlet) ; $\varepsilon_P=2\nu k_P/y_P^2$ à la cellule adjacente ; $\omega_{wall}=6\nu/(\beta_1 y_1^2)$ (Wilcox, $\beta_1$ du modèle) ; SA : $\tilde{\nu}=0$ à la paroi.
- **Loi de paroi** ($y^+\gtrsim 30$) : $U^+=\frac{1}{\kappa}\ln y^+ + B$, $\kappa=0.41$, $B\simeq 5.2$ (paroi lisse) ; production amortie et dissipation imposées dans la cellule de paroi ($k_P$ du $u_\tau$ local, $\varepsilon_P=u_\tau^3/(\kappa y_P)$). **Leçon B7** : les mêmes constantes $\kappa$, $E$, $B$ doivent être utilisées de façon cohérente entre la loi de vitesse, la production amortie et le calcul de $\tau_w=\rho u_\tau^2$ — aucun mélange de jeux de constantes.
- **Entrée** : $k=\tfrac{3}{2}(I\,|\mathbf{u}|)^2$, $\varepsilon=C_\mu^{3/4}k^{3/2}/\ell$, $\omega=\sqrt{k}/(C_\mu^{1/4}\ell)$ — intensité $I$ et longueur $\ell$ déclarées par cas.
- **Sortie** : gradient nul sur tous les scalaires turbulents. **Freestream** : sensibilité documentée du k-omega (motivation même du SST de Menter) — les valeurs freestream doivent être enregistrées dans le rapport de cas.
- **SA** : termes de trip désactivés par défaut, activation explicite documentée par cas.

## 5. Pièges connus
- **B1** — diffusion turbulente assemblée dans la matrice, jamais au second membre.
- **B4** — SA : la fonction $f_w=g\left[(1+g_{w3}^6)/(1+g^6)\right]^{1/6}$, $g=r+c_{w2}(r^6-\dots)$ doit être implémentée exactement (formule continue), pas par morceaux discontinus : une discontinuité tue la convergence près de la paroi.
- **B7** — cohérence $\kappa$/$E$ entre loi de paroi, production amortie et frottement pariétal.
- Anomalie du point d'arrêt : surproduction k-epsilon ; limiteur (realizable, Durbin, ou Kato-Launder) — choix documenté par modèle dans le registre.
- Sensibilité du k-omega au freestream : déclarer et geler les valeurs ; en comparer l'effet (étude de sensibilité) pour toute qualification.
- $\omega$ pariétal très grand : raideur numérique ; formule de la cellule adjacente, jamais de division par $y_1\to 0$ non gardée.
- Régime tampon ($5<y^+<30$) : ni paroi résolue ni log-law valide — politique explicite par cas (rejeter ou déclarer), pas de zone grise silencieuse.
- Négativité post-solve de k/epsilon/omega : correction bornée **télémetrée** ou rejet fort ; jamais de clipping silencieux (Règle d'Or 3).

## 6. Ordre de précision attendu et validation exigée
- Ordre attendu : $O(h)$ UPWIND, $O(h^2)$ centré/correction différée pour k, epsilon, omega, nu-tilde ; temps : ordre du schéma temporel déclaré. À valider par MMS sur au moins 3 niveaux de maillage (pente log-log).
- GCI (Roache, DOI : 10.1115/1.2910291 ; Celik, DOI : 10.1115/1.2960953) sur les QoI : $u_\tau$, $C_f$, profils $u^+$.
- Oracles quantitatifs : TUR-CHANNEL (profil $u^+$ vs log-law et données DNS de canal ; écart pré-déclaré), TUR-FLATPLATE ($C_f(x)$ vs corrélations), TUR-NACA0012 (Cp/Cf vs SU2, matrice N14). Tolérances pré-déclarées dans docs/validation/TURBULENCE_QUALIFICATION_MATRIX.md.
- Bornitude : k, epsilon, omega >= 0 sur tous les cas qualifiés, avant clipping, avec télémétrie.
- Parallèle : parité CPU/MPI/GPU à 1e-12 sur le même maillage (fiche 12).

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : src/cfdx/physics/turbulence.h, turbulence_models.h, turbulence_solver.h, turbulence_transport.h, sst_solver.h, spalart_allmaras.h, wall_functions.h, wall_distance.h, wall_distance_turbulence.h, apps/production_solver/main.cpp.
- Tests exécutables : tests/validation/test_channel_qualification.cpp, tests/validation/test_phase10_turbulence_hardening.cpp et les autres tests test_phase10_* ; cas TUR-CHANNEL, TUR-FLATPLATE, TUR-NACA0012 (docs/validation/TUR-*_QUALIFICATION_CASE.md).
- État actuel : les API de transport k-epsilon, SST et Spalart–Allmaras existent ; le runtime SST transactionnel est raccordé à la boucle externe de `solve_steady_incompressible` via `IncompressibleSolverControls::sst_runtime`, de manière explicite et opt-in. Cela ne signifie pas que les sélecteurs RANS génériques ni les autres modèles sont activés/qualifiés.
- Correction de performance dans `sst_solver.h` : `cell_gradient_dot(k, omega, mesh)` est calculé une fois par itération externe, au lieu d'être recalculé dans chaque cellule (suppression du travail O(N²) redondant sans changement de formule ni tolérance).
- Écarts connus code <-> fiche : la qualification physique SST, la validation des résidus sur l'état accepté après bornage, les campagnes DNS/couche limite/NACA, et l'intégration de production des modèles RANS autres que SST restent ouverts. Les statuts `COMPONENT_READY` du registre C++/JSON ne valent pas `SOLVER_READY`.
- Gaps et prochaines actions : résidus non linéaires recalculés indépendamment sur l'état SST accepté ; MMS transport indépendantes sur au moins trois maillages ; TUR-CHANNEL contre DNS Moser–Kim–Mansour, puis TUR-FLATPLATE et TUR-NACA0012 ; qualification MPI ; GCI et conservation indépendante sur un HEAD exact. Promotion physique toujours BLOCKED selon la gouvernance #118.

## 8. Journal des modifications de la fiche
- 2026-10-12 — création (Phase 3 du Plan Directeur, PR empilée sur governance/constitution-technique).
- 2026-10-12 — approbation du contrat (PR d'approbation #861) ; l'audit du code contre cette fiche reste requis avant toute promotion (Règle d'Or 1).


## 9. Contrat d'itération couplée SST en production (PR #1028)
- `advance_sst_production_iterate_transactionally` orchestre un pas de transport SST sans activer à lui seul le modèle RANS dans le solveur de production.
- Avant le transport, la distance à la paroi et le taux de déformation sont recalculés depuis le maillage et l'itéré courant de vitesse, avec reconstruction tenant compte des conditions aux limites de vitesse fixes.
- Le contrat d'acceptation est atomique pour les champs porteurs `U`, `p`, `phi` et l'état SST (`k`, `omega`, `nu_t`, distance à la paroi et taux de déformation). Un rejet du transport restaure l'état antérieur exact ; une exception restaure l'état avant propagation.
- Le commit n'est autorisé qu'après convergence du transport et validation de l'état SST. Le helper transactionnel est appelé par le runtime SST explicite de la boucle externe de `solve_steady_incompressible` sur le `master` courant ; cela ne constitue pas une qualification physique.
- Le sélecteur RANS générique et les modèles autres que SST doivent rester fail-closed jusqu'à leur intégration complète. Pour SST, le runtime opt-in existe, mais les tests MMS multi-maillages et la qualification DNS/couche limite/NACA de la section 6 restent obligatoires.

## 10. Journal des modifications de la fiche
- 2026-10-22 — contrat d'itération SST transactionnelle documenté pour la PR #1028 ; la vérification du code et les tests exact-head restent obligatoires.


## 11. Outer-iteration rollback regression (PR pending)
- A regression now exercises a nested transaction: the inner SST/carrier transaction commits, then a later outer-iteration rejection rolls back `U`, `p`, authoritative `phi`, and all SST fields to the original outer snapshot.
- This verifies the rollback primitive needed by the production outer loop; it does **not** itself wire the transaction into `solve_steady_incompressible`, refresh `mu_eff` per iteration, or qualify SST.


## 12. Production outer-loop SST runtime (PR pending)
- `IncompressibleSolverControls::sst_runtime` is an explicit opt-in runtime seam. When absent, the legacy fixed-viscosity path is unchanged.
- When present, the caller must provide a validated `SSTProductionState` and SST `TurbulenceTransportControls`; initialization from declared case data or DAT restart remains explicit through the existing initialization/restart helpers. Carrier and SST density/kinematic viscosity must agree exactly. No implicit SST initial values are fabricated.
- Before the first momentum assembly, the solver refreshes wall distance, strain rate and `nu_t` from the actual initialized/restarted carrier field. The outer iteration then advances SST using the pressure-corrected authoritative face mass flux and the explicitly supplied omega face values when configured. A failed transport step raises the existing retryable-iteration path; the enclosing retry restores the SST snapshot together with carrier rollback.
- After the carrier/SST transaction commits, the refreshed `nu_t` is converted to `mu_eff = rho * (nu + nu_t)` for the next momentum assembly. SST convergence is an additional hard gate based on the transport solver's convergence status, finite equation residuals, zero non-finite bound events, and dimensionless relative `nu_t` change no larger than the configured carrier relative tolerance. This prevents declaring the carrier solution converged before momentum has been assembled with the accepted SST closure. Raw dimensional SST residuals are not mixed into the normalized carrier convergence metric; carrier, continuity, and QoI gates are not relaxed.
- This is integration scaffolding, not permission to enable the general RANS selector or a claim of physical qualification. Restart/state lifecycle and production regression coverage must be demonstrated by exact-head CI before merge.
