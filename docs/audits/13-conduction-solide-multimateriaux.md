# FICHE MÉTHODE 13 — Conduction solide multi-matériaux : zones thermiques, puissance volumique, résistances de contact

**En-tête d'évidence (gouvernance §16.5)**
- Scope : conduction solide stationnaire `div(k grad T) + q''' = 0` avec conductivité par zone (un matériau = une zone), puissance volumique par zone (W/m^3) et résistance de contact thermique (m^2.K/W) aux faces internes séparant deux zones différentes. Couvert : maillages orthogonaux, two-point, Dirichlet / gradient fixé / adiabatique. Non couvert : convection (fiche 10 pour l'énergie fluide et le CHT), rayonnement (fiche 11), transitoire, correction non orthogonale (fiche 06), couplage fluide-solide (étage 2 de #885).
- Statut : APPROUVÉE (2026-10-13, PR #887) — contrat rédigé sur l'implémentation du driver T1 ; l'exécution réelle de la campagne T1 sur un HEAD exact mergé reste à la charge du propriétaire avant toute promotion (gouvernance #118).
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (cas `T1-SOLID-COND`) + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
- Issues / PR liées : #880 (T1), #885 (étage 1), #118, PR #887
- Dernier commit vérifié : branche `feature/thermal-solid-production-path` (présence des sources et tests confirmée dans la PR ; à re-vérifier sur le SHA mergé)
- Limitations connues : two-point sans correction non orthogonale ; `heat_flux` (gradient fixé) requiert un patch mono-zone ; R_c déclarée par zone (les deux côtés d'une interface doivent déclarer la même valeur, sinon rejet) ; stationnaire.

## 1. Références sources
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI: 10.1201/9781482234213. (résistance série / moyenne harmonique aux faces, traitement half-cell de Dirichlet)
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN: 978-3-319-16873-9. (assemblage diffusif multi-matériaux)
- J. H. Ferziger, M. Peric, *Computational Methods for Fluid Dynamics*, 3e éd., Springer, 2002 — DOI: 10.1007/978-3-642-56026-2. (conservation discrète, ordre)

## 2. Formulation continue
Régime stationnaire, solide immobile ($\mathbf{u}=0$), propriétés constantes par zone $\Omega_z$ :
$$\nabla\cdot(k_z\nabla T) + q_z''' = 0 \quad \text{dans } \Omega_z,\qquad k_z>0$$
Aux faces internes entre deux zones $a,b$ sans résistance de contact : continuité de la température et du flux. Avec résistance de contact $R_c$ (m^2.K/W) :
$$q_f = \frac{T_a - T_b}{R_c},\qquad q_f=-k_a\frac{\partial T}{\partial n}\Big|_a=-k_b\frac{\partial T}{\partial n}\Big|_b$$
Hypothèses : résistance de contact en série (saut de température localisé à la face, flux continu), $R_c\ge 0$ finie, puissance volumique finie (nulle autorisée).

## 3. Discrétisation Volume Fini
Schéma volumes finis two-point sur cellules $P$ :
$$a_P T_P = \sum_{nb} a_{nb} T_{nb} + b,\qquad a_P=\sum_{nb}a_{nb}+\sum_{f\in\partial P_{Dirichlet}}G_{b,f}$$
- **Faces internes** : conductance en résistance série exacte
$$G_f=\frac{A_f}{\left(d_f/(2k_o)+R_c+d_f/(2k_n)\right)}$$
avec $d_f=\|\mathbf{C}_n-\mathbf{C}_o\|$ (distance entre centres), $A_f$ aire de face. $R_c=0$ redonne la moyenne harmonique (Patankar) ; $R_c>0$ s'applique uniquement aux faces séparant deux zones **différentes** déclarées (une face intra-zone porte $R_c=0$).
- **Terme source** : $b_P \mathrel{+}= q_z''' V_P$ (explicite, indépendant de $T$ : pas de linéarisation requise).
- **Conservativité** : le flux de face est débité symétriquement ($-G_fT_o$, $-G_fT_n$ en hors-diagonale) : conservation machine de l'énergie entre cellules adjacentes.
- **Ordre** : exact pour les profils linéaires par morceaux aux centres (maillage uniforme) ; sinon $O(h^2)$ sur maillage orthogonal.

## 4. Traitement des frontières
- **Dirichlet** ($T=T_w$ sur un patch WALL_THERMAL avec `wall_temperature`) : conductance half-cell $G_b=k_o A_f/d_{b}$ avec $d_b$ distance normale centre→face ; contribution $a_P \mathrel{+}= G_b$, $b \mathrel{+}= G_b T_w$.
- **Gradient fixé / flux imposé** (`heat_flux` $q''$, positif = entrant) : $b \mathrel{+}= k_o A_f g$ avec $g=q''/k_o$ ; requiert un patch **mono-zone** (sinon $k_o$ ambigu — rejet).
- **Adiabatique / symétrie** : contribution nulle.
- **EMPTY** (2D/axisym.) : aucune contribution, patch ignoré.
- **Null-space (jauge de température)** : un système tout-Neumann a un null-space constant ; rejet fail-closed à l'assemblage (diagonale nulle → exception), jamais résolu silencieusement ; l'app exige ≥ 1 Dirichlet.

## 5. Pièges connus
- Une parabole (source volumique uniforme) n'est **pas** exacte aux centres avec Dirichlet aux extrémités : erreur $O(h^2)$ liée au traitement half-cell — ne pas gater sur une fausse « solution exacte », mais sur la convergence d'ordre 2 et des identités discrètes exactes.
- Identité de flux cumulé : le flux à travers un plan interne n'est égal à la puissance générée en amont **que si** la frontière amont est adiabatique ; avec Dirichlet amont, le flux de paroi participe.
- Normalisation du bilan d'énergie : utiliser le **flux brut** ($\sum_{patch}|q_{patch}|+|source|$), jamais le flux net — un cas équilibré (flux entrant = flux sortant, sans source) a un flux net ~0 au niveau du résidu linéaire et une normalisation par le flux net transformerait le bruit d'arrondi en « déséquilibre » $O(1)$.
- $R_c$ déclarée par les deux côtés : incohérence → rejet (fail-closed), pas de moyenne silencieuse.
- Ids de zones : cellule référençant un id non déclaré, ou ids dupliqués → rejet.
- Moyenne harmonique (pas arithmétique) obligatoire à l'interface entre conductivités différentes — l'arithmétique sur-lisse la discontinuité et casse le flux exact.
- Bornitude : l'opérateur assemblé est une M-matrice sur maillage orthogonal (principe du maximum discret garanti) ; aucune correction non orthogonale dans ce driver.

## 6. Ordre de précision attendu et validation exigée
- Ordre attendu : $O(h^2)$ espace sur maillages orthogonaux (two-point), exact pour profils linéaires par morceaux et réseau de résistances.
- Vérification implémentée (`test_thermal_zones_vv`, labels qualification;validation;phase13) : mur composite 3 zones + $R_c$ exacts vs réseau de résistances (1e-9) ; série de raffinement 3 niveaux (8/16/32) sur la parabole source avec ordre observé log-log borné (1.5, 2.5) ; identité de flux cumulé au plan médian (paroi amont adiabatique) à la précision machine ; fermeture énergétique $<$ 1e-12 normalisée par le flux brut ; principe du maximum discret (bornes Dirichlet) sur les cas sans source ; 6 rejets fail-closed.
- Campagne T1 (à exécuter par le propriétaire sur un HEAD exact mergé) : MMS ≥ 3 niveaux avec pente log-log sur cas .cfdx.h5 réels via `cfdx_thermal_solver`, GCI sur les QoI (flux par patch, $T_{max}$), exécution réelle + accord quantitatif archivés — aucun PASS sans cette évidence (gouvernance #118).
- Conservation : bilan d'énergie reconstruit **indépendamment** de la matrice assemblée (`audit_thermal_boundary_flux`) ; tolérance machine sur le flux brut.
- Parallèle : N/A à ce stade (solveur séquentiel ; parité CPU/MPI/GPU à la fiche 12 quand le chemin sera parallélisé).

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : src/cfdx/physics/thermal_zones.h (noyau), apps/thermal_solver/main.cpp (app de production `cfdx_thermal_solver`).
- Tests exécutables : tests/validation/test_thermal_zones_vv.cpp (`test_thermal_zones_vv`).
- Conventions utilisateur : 1 matériau = 1 zone ; `extra_properties` par matériau : `zone_id` (requis si > 1 matériau), `volumetric_power` (W/m^3), `contact_resistance` (m^2.K/W) ; carte cellules→zones par champ scalaire `material_id` (ou `zone_id`) dans le .cfdx.h5 ; BC : `wall_temperature` → Dirichlet, `heat_flux` (W/m^2, positif entrant) → gradient fixé, sinon adiabatique.
- Écarts connus code ↔ fiche : aucun identifié à la rédaction (la fiche décrit l'implémentation du PR #887) ; audit ligne à ligne à refaire à chaque évolution (Règle d'Or 1).
- Gaps et prochaines actions : campagne T1 sur HEAD exact (MMS 3 niveaux + GCI + exécution réelle archivée) ; étage 2 de #885 (couplage fluide + énergie dans le production solver, échange CHT fluide/solide) ; extension transitoire et non orthogonale hors scope présent.

## 8. Journal des modifications de la fiche
- 2026-10-13 — création (PR #887, driver T1 #880 / #885 étage 1) : contrat rédigé sur l'implémentation, approbation par fusion de la PR par le propriétaire ; campagne d'exécution T1 reste à produire avant toute promotion.
