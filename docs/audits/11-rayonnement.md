# FICHE MÉTHODE 11 — Rayonnement thermique (Rosseland, P1, S2S, DOM)

**En-tête d'évidence**
- Scope : modèles de rayonnement gris couplés à l'énergie : Rosseland (diffusion), P1 (approximation différentielle), surfaces-to-surfaces (S2S par facteurs de forme), méthode des ordonnées discrètes (DOM). Non couvert : gaz non gris (WSGG), rayonnement particulaire (MC), semi-transparents.
- Statut : DRAFT (2026-10-12) — contrat mathématique à faire approuver dans la PR ; le code Phase 12 existant doit être audité contre cette fiche avant toute promotion (Règle d'Or 1).
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
- Issues liées : #34 (roadmap, Phase 12) ; docs/validation/THERMAL_RADIATION_REGRESSION_MATRIX.md
- Dernier commit vérifié : f5a728f205cb45ad6186c9fab92554ae4af06bdd (master, 2026-10-12 — présence des sources et tests confirmée)
- Limitations connues : bilans radiatifs non exécutés depuis l'environnement d'audit ; statut déclaratif jusqu'à campagne sur HEAD exact.

## 1. Références sources
- M. F. Modest, *Radiative Heat Transfer*, 3e éd., Academic Press, 2013 — ISBN : 978-0-12-386944-9. (référence de référence pour tous les modèles ci-dessous)
- W. A. Fiveland, « Discrete Ordinates Solutions of the Radiative Transport Equation for Enclosures with Participating Media », *ASME J. Heat Transfer* 106(4):699-706, 1984. (DOM ; DOI non confirmé — citation bibliographique complète)
- S. Rosseland, « Theoretical Astrophysics: Atomic Theory and the Analysis of Line Atmospheres », Oxford University Press, 1936. (approximation de Rosseland ; pas de DOI — ouvrage historique)
- R. Siegel, J. R. Howell, *Thermal Radiation Heat Transfer* (référence complémentaire, à éditer avec ISBN à la première utilisation).

## 2. Formulation continue
Équation de transport radiatif (RTE), milieu gris, absorbant-émettant, non diffusant :
$$\frac{dI(\mathbf{r},\mathbf{s})}{ds}=\kappa I_b(T)-\kappa I(\mathbf{r},\mathbf{s})$$
Quantités dérivées : irradiation $G=\int_{4\pi}I\,d\Omega$, flux radiatif net :
$$q_r=\int_{4\pi}I\,\mathbf{s}\,d\Omega=(1-\pi\text{ selon convention})\;\text{— la convention } G, \; \mathbf{q}_r \text{ doit être gelée par implémentation (leçon B13)}$$

**Rosseland** (milieu optiquement épais, $\tau\gg 1$) : le rayonnement devient un terme de diffusion effectif :
$$\mathbf{q}_r=-k_R\nabla T,\qquad k_R=\frac{16\sigma n^2 T^3}{3\beta_R}$$

**P1** (approximation différentielle) : $I\approx\frac{1}{4\pi}(G+3\,\mathbf{s}\cdot\mathbf{q}_r)$ :
$$\nabla\cdot\left(\frac{1}{3\kappa}\nabla G\right)=\kappa G-4\kappa\sigma T^4$$
**Leçon B13** : le facteur $\pi$ (irradiation $G=4\sigma T^4$ du corps noir intégral vs $I_b=\sigma T^4/\pi$) est le piège classique de la mise en œuvre P1 — la convention intégrale vs intensive doit être cohérente de la source à la frontière.

**S2S** : échanges radiatifs entre surfaces grises diffuses par facteurs de forme $F_{ij}$ :
$$q_i=\frac{\varepsilon_i A_i}{1-\varepsilon_i}\left(E_{bi}-J_i\right)+A_i\sum_j (J_i-J_j)F_{ij},\qquad \sum_j F_{ij}=1$$
**DOM (Fiveland)** : RTE discrétisée en $N$ ordonnées $\mathbf{s}_k$ avec quadrature $w_k$ :
$$\mu_k\frac{\partial I_k}{\partial x}+\eta_k\frac{\partial I_k}{\partial y}+\xi_k\frac{\partial I_k}{\partial z}=\kappa(I_b-I_k)$$

Terme source couplé à l'énergie :
$$S_r=\kappa\left(G-4\sigma T^4\right)$$

## 3. Discrétisation Volume Fini
- **Rosseland** : le terme $-k_R\nabla T$ entre la matrice de conduction/convection comme conductivité effective totale $k_{eff}=k+k_R$ — **leçon B18** : la conductivité **moléculaire k ne doit jamais être perdue** (aucun remplacement, uniquement l'addition) ; $k_R$ est assemblé dans la matrice, itéré de Picard externe sur T.
- **P1** : équation de diffusion elliptique pour $G$ discrétisée comme le Laplacien (fiche 01) avec $\Gamma=1/(3\kappa)$ par face ; maillages non orthogonaux : correction de non-orthogonalité (fiche 01 §3).
- **S2S** : facteurs de forme par ray tracing ou contours (leçon **B5** : réciprocité $A_iF_{ij}=A_jF_{ji}$ et fermeture $\sum_jF_{ij}=1$ vérifiées à la tolérance machine par test dédié — un facteur de forme erroné viole la conservation d'énergie radiative).
- **DOM** : balayage par direction, schéma UPWIND spatial borné ; symétrie des quadratures conservée (aucune direction lâchée) ; verrouillage (ray effect) documenté par nombre d'ordres.
- **Source couplé** $S_r=\kappa(G-4\sigma T^4)$ — **leçon B3** : linéarisation de Newton pour l'implicite :
$$\frac{\partial S_r}{\partial T}=\kappa\left(\frac{\partial G}{\partial T}-16\sigma T^3\right)$$
le terme $-16\kappa\sigma T^3$ entre au coefficient diagonal $a_P$ (stabilise), $\partial G/\partial T$ selon le modèle (0 pour S2S, implicite pour P1/Rosseland).

## 4. Traitement des frontières
- Surfaces grises diffuses : condition Robin radiative sur $G$ (P1,Marshak) ou contribution directe au système S2S ; émissivités $\varepsilon_i$ déclarées par patch.
- Frontières noires ($\varepsilon=1$) : $J=E_b$ ; réfléchissantes ($\varepsilon=0$) : symétrie radiative.
- Ouvertures (vers l'extérieur à $T_\infty$) : flux net $\sigma(T_w^4-T_\infty^4)$ linéarisé de la même façon (B3).
- **Gauge** : le rayonnement S2S n'a pas de null-space (bilan fermé), P1 a le null-space du Laplacien si toutes les frontières sont adiabatiques radiativement — point de référence à la P1 (fiche 04) si le cas se présente.

## 5. Pièges connus
- **B13** : facteur pi du P1 (convention $G$ intégrale vs $I_b$ intensive) — cohérence de bout en bout exigée.
- **B3** : signe du source radiatif dans l'équation d'énergie ($S_r$ chauffe < 0 ou refroidit > 0 selon $G$ vs $4\sigma T^4$) + linéarisation de Newton obligatoire pour la stabilité du couplage implicite.
- **B18** : Rosseland — conserver $k$ moléculaire en plus de $k_R$ ; le $c_p$ du bilan d'énergie reste celui de la fiche 10 (aucune confusion $k_R$/$c_p$).
- **B5** : S2S — réciprocité et fermeture des facteurs de forme à la tolérance machine ; petit nombre de surfaces : conditionnement du système (facteurs de forme) à surveiller.
- Domaine de validité : Rosseland uniquement $\tau\gg 1$ ; P1 pour géométries sans obstruction forte ; DOM/S2S sinon — le modèle utilisé doit être déclaré avec sa justification dans le rapport de cas (pas de Rosseland silencieux sur cas optiquement mince).
- Couplage instable si le source radiatif est purement explicite (oscillation T) — toujours linéariser (B3).
- Anisotropie : hors scope ici (milieu non diffusant) ; toute extension devra avoir sa fiche.

## 6. Ordre de précision attendu et validation exigée
- Ordre : P1/Rosseland = ordre du Laplacien ($O(h^2)$ orthogonales) ; DOM : $O(h)$ en espace par direction, erreur angulaire dominée par la quadrature (à documenter) ; S2S : exact par paires, l'erreur vient des facteurs de forme (tolérance machine exigée sur réciprocité/fermeture).
- MMS pour le P1 (éq. de diffusion sur $G$) sur au moins 3 maillages ; GCI sur les QoI ($q_r$ pariétal, $T_{max}$, fraction radiative du flux total).
- Oracles quantitatifs : **cavité chauffée** (benchmark classique : parois à températures imposées, solution S2S exacte en géométrie simple), bilan d'énergie radiatif fermé à la tolérance machine (ce qui entre = ce qui sort + ce qui est absorbé), comparaison croisée entre modèles P1/S2S/DOM sur le même cas (doivent converger vers la même physique dans leur domaine de validité) ; matrice docs/validation/THERMAL_RADIATION_REGRESSION_MATRIX.md.
- Conservation : bilan radiatif fermé machine ; bornitude de T avec télémétrie.
- Parallèle : parité CPU/MPI à 1e-12 sur le même maillage (fiche 12).

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : src/cfdx/physics/radiation.h, radiation_models.h, radiation_advanced.h, radiation_s2s.h, radiation_solver.h, radiation_streaming.h.
- Tests exécutables : tests/validation/test_phase12_radiation_vv.cpp, test_s2s_radiation_vv.cpp, test_radiation_vv.cpp, test_thermal_radiation_regression.cpp, test_thermal_radiation_model_matrix.cpp.
- Écarts connus code <-> fiche : à auditer ligne à ligne lors de la Phase 3 — B3, B5, B13 et B18 sont les quatre premiers points de contrôle.
- Gaps et prochaines actions : bilan radiatif fermé comme oracle machine systématique ; cavité chauffée benchmark avec superposition graphique dans le rapport ; comparaison croisée P1/S2S/DOM documentée ; promotion BLOCKED derrière Phase 1.

## 8. Journal des modifications de la fiche
- 2026-10-12 — création (Phase 3 du Plan Directeur, PR empilée sur governance/constitution-technique).
