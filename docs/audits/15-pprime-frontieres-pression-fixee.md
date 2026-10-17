# FICHE MÉTHODE — Frontières à pression fixée dans l'équation de correction de pression (SIMPLE/SIMPLEC)

> Fiche 15. Créée à l'occasion du correctif du bug #894 (terme de flux de pression explicite au second membre de p′ aux frontières à pression imposée).

**En-tête d'évidence (gouvernance §16.5)**
- Scope : assemblage de l'équation de correction de pression p′ aux frontières à pression imposée (`FIXED_VALUE`) du solveur incompressible stationnaire — SIMPLE et SIMPLEC ; 2D plan, 3D, axisymétrique. Non couvert : problèmes tout-Neumann et jauge de pression (cas canal, cf. #889), chemin COUPLED+PCD (Ghia 32×32), tolérances de qualification.
- Statut : UPDATED (2026-10-15, correctif #904 ; création DRAFT 2026-10-14)
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (cas) + NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json (exigences)
- Issue / PR liées : #894 (branche de diagnostic — ne jamais merger), #904 (facteur de volume dAU — branche fix/870-planar-mms-flux-closure), #870, #893, #889
- Dernier commit vérifié : 6dff867b (2026-10-13)
- Limitations connues : validation des correctifs #894/#904 en attente du run validation-total au HEAD mergé.

## 1. Références sources (obligatoire : DOI ou ISBN)
- S. V. Patankar, « Numerical Heat Transfer and Fluid Flow », Hemisphere Publishing, 1980 — DOI: 10.1201/9781482234213 (chap. 6 : algorithme SIMPLE, traitement des frontières de pression).
- C. M. Rhie, W. L. Chow, « Numerical Study of the Turbulent Flow Past an Airfoil with Trailing Edge Separation », AIAA Journal 21(11), 1983, pp. 1525–1532 — DOI: 10.2514/3.8297 (flux de face cohérent, anti-checkerboarding).
- F. Moukalled, L. Mangani, M. Darwish, « The Finite Volume Method in Computational Fluid Dynamics », Springer, 2016 — ISBN 978-3-319-16873-9 (chap. 15 : SIMPLE/SIMPLEC, équation de pression et conditions aux limites).

## 2. Formulation continue
Écoulement incompressible newtonien stationnaire :
$$ \nabla \cdot \mathbf{u} = 0, \qquad \rho(\mathbf{u} \cdot \nabla)\mathbf{u} = -\nabla p + \mu \nabla^2 \mathbf{u} $$
Hypothèses : fluide newtonien, densité constante. La contrainte de continuité fait de la pression le multiplicateur de Lagrange de la divergence : toute frontière à pression imposée doit pouvoir évacuer le débit nécessaire à la satisfaction de la contrainte de masse.

## 3. Discrétisation Volume Fini
Contrainte de continuité discrétisée sur la cellule $P$ :
$$ \sum_f \phi_f = 0, \qquad \phi_f = \phi^{HbyA}_f - \rho \, dAU_f \, |S_f| / d_f \, (p_F - p_P) $$
- $\phi^{HbyA}_f$ : flux prédicteur interpolé (Rhie–Chow, cf. fiche 03).
- $dAU_f = V_f \, r_f$ : réponse de vitesse directionnelle de la face — $r_f = 1/A_P$ pour SIMPLE ($r_{AU}$), dénominateur cohérent pour SIMPLEC ($r_{AtU}$), multiplié par le volume ($V_f$ interpolé sur les faces internes ; $V_P$ sur les faces de bord). Le même opérateur dAU = V·rAU doit apparaître dans la matrice p′, dans le flux de bord autoritaire (`make_rhie_chow_mass_flux`) et dans le coefficient des faces internes (`directional_face_coefficient`).
- L'équation de correction $p'$ hérite de cette structure : les contributions des faces internes sont assemblées dans la matrice, celles des faces de bord dans la diagonale et le second membre.

## 4. Traitement des frontières
Frontière à pression imposée (`FIXED_VALUE`, ex. sortie de pression) : la correction de pression y est nulle, $p'_b = 0$. La contribution de la face de bord à l'équation de $p'$ est un demi-Dirichlet :
- matrice : $a_{PP} \mathrel{+}= \rho \, dAU_b \, (\mathbf{S}_b \cdot \mathbf{e}) / d_b$, avec $dAU_b = V_P \, r_{AU,P}$ (resp. $r_{AtU}$ en SIMPLEC), $\mathbf{e}$ la direction cellule→centre de face — le facteur de volume manquait dans l'implémentation d'origine (**bug #904**, corrigé le 2026-10-15) ;
- second membre : $b_P \mathrel{+}= \rho \, dAU_b \, (\mathbf{S}_b \cdot \mathbf{e}) / d_b \, (p_{bc} - p_P)$, avec exactement le même coefficient que la matrice — **terme explicite qui manquait (bug #894)**.
Le flux de bord cohérent s'écrit $\phi_b = \phi^{HbyA}_b - \rho \, dAU_b \, (\mathbf{S}_b \cdot \mathbf{e}) / d_b \, (p_{bc} - p_P)$ — l'opérateur exact reconstruit par `make_rhie_chow_mass_flux` : sans le terme au second membre, l'équation discrétisée force $\phi_b \to 0$ et le déséquilibre de masse au bord ne peut plus être corrigé. Le second membre $b$ étant assemblé ($b(c) = -\text{continuité}[c]$) avant la boucle des faces à pression fixée, le terme doit être ajouté à $b$ directement.
Frontières WALL / EMPTY / périodiques : flux nul ou miroir, traitées dans les fiches 03 et 04 ; non modifiées ici.

## 5. Pièges connus (Gotchas)
- **Demi-Dirichlet sans source (bug #894)** : résidu momentum au zéro machine mais plancher de continuité figé et `flux_velocity_mismatch_linf` strictement constant au fil des itérations — signature d'une incohérence matrice/second membre au bord, pas d'un défaut de convergence du solveur linéaire.
- **Demi-Dirichlet sans facteur de volume (bug #904)** : un opérateur de bord $r_f = r_{AU}$ sans $V_P$ met la matrice p′ en désaccord avec le flux autoritaire reconstruit par `make_rhie_chow_mass_flux` (qui applique dAU = V·rAU en interne). La continuité gèle sur un plancher non nul avec momentum au zéro machine (naca0012 : 1.2828e-5 figé sur 8000 it. ; bfs : 5.55785e-5 figé sur 5000 it.). Preuve FP64 (réplique #872, duct type-NACA 16×8) : opérateur défectueux figé à 5.95283444e-4 bit-identique des itérations 100 à 1200 ; opérateur corrigé convergé à l'itération 85 (continuité 4.89e-11 ≤ gate 1e-7). Le défaut était invisible sur les cas à pression constante/linéaire (Couette/Poiseuille : $p \approx p_{bc}$ annule le terme de bord), c'est pourquoi les tests rapides restaient verts.
- $b$ est assemblé avant la boucle des faces à pression fixée : le terme explicite doit modifier $b$, pas la variable de continuité déjà consommée.
- $dAU_b = V_P \, r_f$ doit être le même opérateur directionnel (V·rAU/SIMPLE, V·rAtU/SIMPLEC) pour la matrice, la source et le flux autoritaire, sinon l'opérateur de bord diffère de l'opérateur interne.
- Ne pas confondre avec la jauge de pression : un problème tout-Neumann sans ancre est singulier (cas canal, #889) — problème distinct, non couvert par ce terme.

## 6. Ordre de précision attendu et validation exigée
- Ordre attendu : $O(h^2)$ en espace pour le gradient (Gauss) ; le demi-Dirichlet de bord est une extrapolation d'ordre 1 en $d_f$, bornée par la physique de la CL. Inchangé par le correctif (il rétablit la cohérence, pas l'ordre).
- Oracles quantitatifs : qualification BFS et NACA0012 (labels `qualification;validation;long`), tolérances pré-déclarées, inchangées par le correctif.
- Conservation : le flux de bord évacue désormais le débit prédicteur corrigé ; bilan de masse global à vérifier dans le run validation-total (tolérance du contrat de cas).
- Parallèle : assemblage séquentiel inchangé ; pas de nouvelle surface de divergence CPU/MPI/GPU.

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : src/cfdx/physics/steady_incompressible_solver.h (assemblage p′, boucle des faces à pression fixée).
- Tests exécutables : test_bfs_qualification, test_naca0012_qualification, test_channel_qualification (labels qualification/validation/long ; le masquage de leur code de retour de test, retiré par la PR selon la politique #893) ; ghia76_cavity (régression).
- Écarts connus code ↔ fiche : le terme explicite de la section 4 est ajouté par la PR liée (fix/pprime-fixed-pressure-rhs) ; le facteur de volume dAU = V·rAU (sections 3–4) est corrigé par #904 sur fix/870-planar-mms-flux-closure ; à confirmer par CI verte puis run validation-total au HEAD mergé.
- Gaps et prochaines actions : cas canal mal posé (tout-Neumann + ancre sur cellule d'entrée — audit séparé, #889) ; blocage COUPLED+PCD Ghia 32×32 (audit séparé) ; cas axisymétriques (#895).

## 8. Journal des modifications de la fiche
- 2026-10-14 — PR fix/pprime-fixed-pressure-rhs (diagnostic #894) — création.
- 2026-10-15 — Issue #904 — correction du facteur de volume manquant (dAU = V·rAU) dans le demi-Dirichlet de bord : §3/§4 réécrits avec $dAU_b = V_P \, r_f$ et $(\mathbf{S}_b \cdot \mathbf{e})/d_b$ ; gotcha §5 ajouté ; preuve FP64 (plancher figé 5.95283444e-4 bit-identique it. 100→1200 vs convergence it. 85, continuité 4.89e-11).