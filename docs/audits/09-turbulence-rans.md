# FICHE MÉTHODE 09 — Turbulence RANS (k-epsilon, k-omega, SST, Spalart-Allmaras) et lois de paroi

**En-tête d'évidence**
- Scope : modèles RANS à viscosité tourbillonnaire (k-epsilon standard, k-omega, SST de Menter), modèle à une équation de Spalart-Allmaras, lois de paroi (couche logarithmique), distance de paroi. Non couvert : LES/hybrides, RSM, transition, correction de compressibilité (à documenter séparément si activée).
- Statut : APPROUVÉE (2026-10-12, PR #861) — contrat mathématique approuvé ; le code Phase 10 existant doit être audité contre cette fiche avant toute promotion (Règle d'Or 1).
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
$$\nu_t=\frac{a_1 k}{\max(a_1\omega,\; S\,F_2)},\qquad a_1=0.31$$

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
- Audit production (PR #941, base `120070dedebecd70e9f78f30691e9e5172a71925`) : les API `solve_kepsilon_transport`, `solve_sst_transport`/`solve_sst_transport_dynamic_blending` et `solve_spalart_allmaras_transport` existent, mais la recherche des appelants et l'inspection de `apps/production_solver/main.cpp` n'ont trouvé aucun branchement de ces transports dans la boucle de quantité de mouvement de production. Les cas qui déclarent un modèle autre que laminaire sont donc rejetés explicitement par le garde-fou ajouté dans cette PR ; aucun RANS de production n'est revendiqué.
- Correction de performance dans `sst_solver.h` : `cell_gradient_dot(k, omega, mesh)` est calculé une fois par itération externe, au lieu d'être recalculé dans chaque cellule (suppression du travail O(N²) redondant sans changement de formule ni tolérance).
- Écarts connus code <-> fiche : le couplage de la viscosité turbulente à la matrice de quantité de mouvement, l'initialisation des variables transportées et la convergence externe restent à implémenter et qualifier. Les statuts `COMPONENT_READY` du registre C++/JSON ne valent pas `SOLVER_READY`.
- Gaps et prochaines actions : MMS transport indépendantes sur au moins trois maillages ; intégration SST de bout en bout ; TUR-CHANNEL contre DNS Moser–Kim–Mansour, puis TUR-FLATPLATE et TUR-NACA0012 ; qualification MPI des modèles RANS ; GCI et conservation indépendante sur un HEAD exact. Promotion physique toujours BLOCKED selon la gouvernance #118.

## 8. Journal des modifications de la fiche
- 2026-10-12 — création (Phase 3 du Plan Directeur, PR empilée sur governance/constitution-technique).
- 2026-10-12 — approbation du contrat (PR d'approbation #861) ; l'audit du code contre cette fiche reste requis avant toute promotion (Règle d'Or 1).
