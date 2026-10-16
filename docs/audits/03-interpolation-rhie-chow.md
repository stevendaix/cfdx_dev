# FICHE MÉTHODE 03 — Interpolation de Rhie–Chow (maillage colocalisé)

**En-tête d'évidence**
- Scope : calcul du flux de masse aux faces pour un arrangement colocalisé (vitesse et pression aux centres de cellules). Couvre le terme de stabilisation, la sensibilité à la relaxation (Majumdar) et aux forces volumiques.
- Statut : APPROVED (2026-10-12)
- Issues liées : #601 (N9, S9 — `docs/validation/N9_S9_RHIE_CHOW_CONSISTENCY.md’), #34 ; PRs #821, #851
- Dernier commit vérifié : 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : audit code ↔ fiche (section 7) à exécuter en Phase 1 (N9).

## 1. Références sources
- C. M. Rhie, W. L. Chow, « Numerical study of the turbulent flow past an airfoil with trailing-edge separation », *AIAA Journal* 21(11):1525–1532, 1983. (interpolation de quantité de mouvement)
- S. Majumdar, « Role of underrelaxation in momentum interpolation for calculation of flow with nonstaggered grids », *Numerical Heat Transfer, Part B* 13:125–132, 1988 — DOI : 10.1080/10407798808551376.
- Moukalled, Mangani, Darwish (2016), ISBN 978-3-319-16873-9, ch. 15 ; Ferziger & Perić (2002), DOI 10.1007/978-3-642-56026-2, § 7.5.

## 2. Formulation continue
Le flux de masse aux faces découle de l'équation de quantité de mouvement ; en colocalisé, interpoler simplement $\\mathbf{u}_f = \\overline{\\mathbf{u}}$ (moyenne des centres) produit un champ de pression découplé (checkerboarding).

## 3. Discrétisation Volume Fini (forme exacte attendue dans le code)
$$ \\dot m_f = \\rho_f\\,\\mathbf{S}_f\\cdot\\overline{\\mathbf{u}}_f \\;-\\; \\\\overline{\\left(\\frac{V}{a_P}\\right)}_f\\,(\\nabla p)_f\\cdot\\mathbf{S}_f \\;+\\; \\\\overline{\\left(\\frac{V}{a_P}\\right)}_f\\,\\overline{(\\nabla p)}_f\\cdot\\mathbf{S}_f \\;+\\; \\\\text{(termes de correction Majumdar)} $$
Contrats :
1. La moyenne $\\overline{(V/a_P)}_f$ est pondérée par les coefficients **au moment de l'assemblage**, et les écarts de relaxation sont corrigés par le facteur de Majumdar $\\alpha/(1+(1-\\alpha)(V/a_P))$ — sans lui, la solution dépend de $\\alpha$ (leçon #419 : indépendance à α vérifiée à 3,95e–9).
2. La différence $(\\nabla p)_f - \\\\overline{(\\nabla p)}_f$ s'annule pour une pression linéaire : **la stabilisation ne doit pas ajouter d'erreur sur les champs exactement linéaires** (test obligatoire : Rhie–Chow nul sur $p$ linéaire).
3. Corps force (gravité/Boussinesq) : traiter les termes de force exactement comme la pression dans la reconstruction, sinon artefacts de colonne hydrostatique.

## 4. Traitement des frontières
- WALL : $\\dot m_f = 0$ exact (aucune reconstruction).
- Entrée/sortie : flux imposé ou reconstruit à partir de la valeur frontière ; la contribution Rhie–Chow aux faces frontières doit être identique à celle utilisée dans l'équation de continuité.
- EMPTY/périodique : reconstruction au niveau de la face virtuelle avec même formule.

## 5. Pièges connus
- **Checkerboarding** si le terme de stabilisation est omis ; **sur/sous-correction** si $a_P$ relaxé est utilisé à la place du $a_P$ non relaxé dans le $rAU$ (leçon A-B2 de l'audit #419 : facteur de sur-amplification $a_P/(a_P-\\Sigma|a_{nb}|)$).
- Dépendance à la relaxation sans correction Majumdar.
- Flux de continuité ≠ flux de conservation : N9.0 exige que le même flux serve aux deux, ou que l'écart soit documenté (PR #851, test d'antisymétrie des flux Poiseuille).

## 6. Ordre de précision attendu et validation exigée
- $O(h^2)$ sur vitesse si la reconstruction est consistante ; le terme de stabilisation est $O(h^2)$ par construction (différence de gradients).
- Validation exigée : Poiseuille (solution exacte + antisymétrie des flux) ; indépendance à α ; indépendance au pas de temps (transitoire) ; N9.7 négatifs : rejet d'une épingle de pression, détection de null-space faux, etc. (PR #851) ; bilan de masse $< 10^{-10}$.

## 7. Statut d'implémentation CFDX
- Sources : `src/cfdx/physics/steady_incompressible_solver.h’, `finite_volume_transport.h’ (HbyA, rAU/rAtU), `src/cfdx/core/numerics/…’.
- Tests : N9 S9 (`docs/validation/N9_S9_RHIE_CHOW_CONSISTENCY.md’), `test_n9_negative_diagnostics’ (#851) ; campagne fix4b/fix4c.
- Écarts code ↔ fiche : **à auditer en Phase 1** (N9) : reconstruction HbyA, gradient de pression aux faces, correction Majumdar, cohérence flux continuité/conservation.
- Gaps : campagne physique complète (Couette, Poiseuille, Ghia Re=100, skewed) sur le HEAD exact avec vrais résidus et chemin `.cfdx.h5’ — critère de fermeture #601.

## 8. Journal
- 2026-10-12 — création ; intègre les leçons #419 (A-B2/A-B4), #821, #851.

- 2026-10-16 — audit 3D du couplage pression–vitesse : le coefficient implicite doit être identique dans le flux Rhie–Chow et la matrice de correction. Avec e = d/|d| et dAU_f le tenseur diagonal de réponse volumique interpolé, D_f = rho * (dAU_f : (e⊗e)) * (Sf·e)/|d|. L'utilisation de |Sf|/|d| et d'une projection sur la normale de face n'est équivalente que sur une face orthogonale. L'écart a été trouvé dans les chemins SIMPLE/SIMPLEC, COUPLED, et aux frontières à pression fixée, y compris dans le terme de flux SIMPLEC. PR #923 mutualise ce coefficient, vérifie sa limite orthogonale et teste directement le flux Rhie–Chow de production sur une face skewée avec réponse directionnelle anisotrope ; validation CI et campagnes physiques encore en attente. Cette correction aligne l'opérateur implicite mais ne qualifie pas, à elle seule, la précision globale sur maillages fortement non orthogonaux.
