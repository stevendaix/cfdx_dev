# FICHE MÉTHODE 24 — Diffusion non orthogonale en transport scalaire

**En-tête d'évidence**
- Scope : flux diffusif dans l'assembleur de transport scalaire et reconstruction des flux de transport ; correction tangente non orthogonale avec conductance normale cohérente avec la matrice. Non couvert : qualification complète sur polyèdres généraux, ordre observé multi-maillages, tous les appels Poisson/thermique/turbulence.
- Statut : PARTIEL — correction partagée implémentée ; campagne CI sur HEAD exact en attente.
- Registre autoritaire : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` et `docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json`.
- Issues liées : F06 / #461 (N3 diffusion et correction non orthogonale).
- Limitations connues : le premier oracle porte sur un champ affine sur maillage hexaédrique incliné et non uniforme ; la validation polyédrique et GCI restent requises.

## 1. Références sources
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI : 10.1201/9781482234213. Décomposition flux/matrice et correction différée.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN : 978-3-319-16873-9. Diffusion non orthogonale, interpolation de face et correction explicite.
- J. H. Ferziger, M. Perić, *Computational Methods for Fluid Dynamics*, 3e éd., Springer, 2002 — DOI : 10.1007/978-3-642-56026-2. Cohérence entre flux de face et discrétisation du volume fini.

## 2. Formulation continue
Pour un scalaire \(\psi\) avec coefficient diffusif \(\Gamma\) :
\[
\nabla\cdot(\mathbf{u}\psi)-\nabla\cdot(\Gamma\nabla\psi)=S_\psi.
\]
Le flux diffusif orienté selon le vecteur de surface \(\mathbf S_f\) est
\[
F_{diff,f}=\Gamma_f\nabla\psi_f\cdot\mathbf S_f.
\]
Une décomposition orthogonale pure ne représente pas le flux exact lorsque le vecteur entre centres \(\mathbf d=\mathbf C_N-\mathbf C_P\) n'est pas aligné avec \(\mathbf S_f\). La composante tangente au terme implicite doit être conservée ; l'omettre rend le flux de diffusion incohérent avec le gradient physique.

## 3. Discrétisation Volume Fini
Le terme implicite utilise la même conductance normale que l'assembleur :
\[
D_f=\frac{|\mathbf S_f|}{d_P/\Gamma_P+d_N/\Gamma_N},
\]
avec \(D_f=\Gamma |\mathbf S_f|/(d_P+d_N)\) pour un coefficient constant. Pour la décomposition sur-relaxée associée à cette conductance :
\[
\alpha_f=\frac{|\mathbf S_f|^2}{\mathbf S_f\cdot\mathbf d},\qquad
\mathbf S_{t,f}=\mathbf S_f-\alpha_f\mathbf d.
\]
Le flux corrigé est
\[
F_{diff,f}=D_f(\psi_N-\psi_P)
+\Gamma_f\nabla\psi_f\cdot\mathbf S_{t,f},\qquad
\Gamma_f=D_f/\alpha_f.
\]
- Le flux tangentiel est explicite et différé, évalué avec le champ scalaire courant et un gradient aux moindres carrés, linéaire-exact lorsque le stencil est de rang plein. Green–Gauss avec interpolation deux-points n'est pas utilisé pour cette correction, car il est biaisé lorsque les centres de cellules et les centroïdes de faces ne sont pas alignés.
- La correction est assemblée de façon conservative : ajouter le flux orienté au second membre du propriétaire et le soustraire au voisin. Les contributions d'une face interne s'annulent globalement.
- `compute_nonorthogonal_diffusion_face_fluxes` fournit le flux de correction par face ; `assemble_scalar_equation` et `reconstruct_scalar_transport_flux` réutilisent le même calcul et la même conductance normale. Le gradient de correction utilise `compute_gradient_least_squares`.
- Si une face diffusive non orthogonale exige un champ scalaire courant et que celui-ci n'est pas fourni, l'assembleur échoue explicitement au lieu de retomber silencieusement sur la diffusion orthogonale.

## 4. Traitement des frontières
- EMPTY : aucune contribution de flux ni de correction.
- Valeur imposée : le gradient utilise la valeur de face résolue, avec priorité aux conditions par face puis aux valeurs de face puis à la condition de patch.
- Gradient imposé : la valeur de face est extrapolée à partir du gradient normal et de la distance normale ; la correction tangente emploie ce même état de frontière.
- Gradient nul : le flux diffusif de frontière reste nul ; aucune correction tangente indépendante n'est ajoutée.
- Les conditions et valeurs de frontière sont résolues par `resolve_scalar_boundary_condition`, partagé entre reconstruction du gradient et assemblage.

## 5. Pièges connus
- Utiliser la distance euclidienne entre centres dans la conductance à la place des distances normales rend la matrice incohérente avec le flux corrigé.
- Ajouter la correction tangente avec le mauvais signe rompt l'exactitude affine ; pour \(-\nabla\cdot F_{diff}=S\), la divergence de correction est déplacée au second membre avec le signe positif.
- Une correction calculée pour l'assemblage mais omise dans la reconstruction produit des bilans de flux incompatibles.
- Une face quasi singulière avec \(\mathbf S_f\cdot\mathbf d\approx0\) rend la décomposition instable et doit être rejetée.
- Le gradient LS est linéaire-exact lorsque le stencil est de rang plein ; les maillages dégénérés ou les stencils insuffisants restent à couvrir par des tests explicites.
- Le test hexaédrique incliné ne constitue pas une preuve d'ordre sur un maillage polyédrique arbitraire.

## 6. Ordre de précision attendu et validation exigée
- Exactitude affine : \(\nabla^2\psi=0\) pour un champ linéaire, y compris sur un maillage incliné non uniforme, dans les cellules intérieures utilisées comme oracle.
- Parité opérateur/flux : le résidu de l'assembleur et la divergence du flux reconstruit doivent correspondre cellule par cellule ; les corrections internes doivent s'annuler globalement.
- Régression N3 : le résidu corrigé doit être inférieur à \(10^{-10}\), tandis que le résidu deux-points non corrigé doit rester supérieur à \(10^{-3}\) sur le cas choisi, pour prouver que la correction est active.
- Qualification restante : au moins trois maillages pour l'ordre observé, un cas polyédrique incliné, et des contrôles de parité pour les appels Poisson, énergie, turbulence, reconstruction et correction différée.
- Aucun élargissement de tolérance et aucun repli silencieux vers une méthode orthogonale ne sont autorisés.

## 7. Statut d'implémentation CFDX
- Fichiers source : `src/cfdx/physics/finite_volume_transport.h`.
- Entrées partagées : `resolve_scalar_boundary_condition`, `compute_scalar_gradient_with_boundary`, `compute_nonorthogonal_diffusion_face_fluxes` ; la correction utilise le gradient LS.
- Gradient : `compute_nonorthogonal_diffusion_face_fluxes` utilise le gradient aux moindres carrés pour préserver l'exactitude affine sur maillages déformés.
- Assembleur : `assemble_scalar_equation` ajoute le flux tangentiel au second membre de façon conservative.
- Périmètre : la correction est activée par défaut pour le transport scalaire ; les wrappers de quantité de mouvement la désactivent explicitement afin de préserver la qualification COUPLED existante. Son activation sur la quantité de mouvement doit faire l'objet d'une qualification dédiée.
- Reconstruction : `reconstruct_scalar_transport_flux` utilise la même correction et la même conductance normale.
- Test exécutable : `tests/validation/test_nonorthogonal_laplacian_campaign.cpp` — champ affine sur maillage hexaédrique incliné et non uniforme, résidu assemblé, divergence du flux reconstruit, rejet si le champ courant manque.
- Registre : `docs/validation/CFDX_QUALIFICATION_REGISTRY.json` (cas `SCALAR-TRANSPORT-NONORTH-AFFINE`, statut DIAGNOSTIC), `NUMERICAL_METHOD_MATURITY_AUDIT.json/.md` et `NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json`.
- Gaps : CI exacte en attente ; MMS polyédrique, GCI, et audit exhaustif de tous les appels Poisson/reconstruction/correction différée requis avant promotion.


## Journal de correction post-merge
- 2026-10-26 — audit après F06 / PR #1126 : le chemin N9 PCD construisait des opérateurs utilisés uniquement sous forme de matrice, mais l'assembleur tentait aussi d'évaluer la correction non orthogonale explicite sans champ courant. Le correctif désactive cette correction différée uniquement dans les deux assemblages PCD matrix-only ; la conductance normale implicite et la diffusion de production restent inchangées. L'audit énergie transmet explicitement la température courante. Les assertions de résidu avant/après relaxation sont documentées dans la fiche 11. Statut : validation CI du HEAD exact en cours ; aucune promotion de qualification n'est revendiquée.
