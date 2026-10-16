# FICHE MÉTHODE 03 — Interpolation de Rhie–Chow (maillage colocalisé)

**En-tête d'évidence**
- Scope : calcul du flux de masse aux faces pour un arrangement colocalisé (vitesse et pression aux centres de cellules). Couvre le terme de stabilisation, la sensibilité à la relaxation (Majumdar) et aux forces volumiques.
- Statut : AUDIT EN COURS — PR `fix/rhie-chow-operator-consistency-audit`; merge interdit avant validation numérique complète
- Issues liées : #601 (N9, S9 — `docs/validation/N9_S9_RHIE_CHOW_CONSISTENCY.md’), #34 ; PRs #821, #851
- Dernier commit vérifié : 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : les modifications d'opérateur sont en revue; les campagnes physiques et le balayage d'indépendance à α doivent être exécutés sur le SHA final de la PR.

## 1. Références sources
- C. M. Rhie, W. L. Chow, « Numerical study of the turbulent flow past an airfoil with trailing-edge separation », *AIAA Journal* 21(11):1525–1532, 1983. (interpolation de quantité de mouvement)
- S. Majumdar, « Role of underrelaxation in momentum interpolation for calculation of flow with nonstaggered grids », *Numerical Heat Transfer, Part B* 13:125–132, 1988 — DOI : 10.1080/10407798808551376.
- Moukalled, Mangani, Darwish (2016), ISBN 978-3-319-16873-9, ch. 15 ; Ferziger & Perić (2002), DOI 10.1007/978-3-642-56026-2, § 7.5.

## 2. Formulation continue
Le flux de masse aux faces découle de l'équation de quantité de mouvement ; en colocalisé, interpoler simplement $\\mathbf{u}_f = \\overline{\\mathbf{u}}$ (moyenne des centres) produit un champ de pression découplé (checkerboarding).

## 3. Discrétisation Volume Fini (forme exacte attendue dans le code)
$$ \\dot m_f = \\rho_f\\,\\mathbf{S}_f\\cdot\\overline{\\mathbf{u}}_f \\;-\\; \\\\overline{\\left(\\frac{V}{a_P}\\right)}_f\\,(\\nabla p)_f\\cdot\\mathbf{S}_f \\;+\\; \\\\overline{\\left(\\frac{V}{a_P}\\right)}_f\\,\\overline{(\\nabla p)}_f\\cdot\\mathbf{S}_f \\;+\\; \\\\text{(termes de correction Majumdar)} $$
Contrats :
1. La réponse tensorielle et son interpolation utilisent les diagonales **au moment de l'assemblage**. La décomposition implicite + différée doit reconstruire le flux tensoriel interpolé, avec un coefficient implicite $((\\mathbf D_f\\mathbf e)\\cdot\\mathbf S_f)/d_{PN}$ identique dans la matrice et le flux.
2. Sur une face orthogonale et un champ de pression linéaire exact, le résidu RC s'annule. Sur une face oblique, le reste différé n'est pas forcé à zéro : il doit fermer exactement la différence entre le flux tensoriel interpolé et la contribution implicite.
3. La correction de Majumdar reste un helper testable, mais **n'est pas activée en production** par cette PR tant que son interaction avec ce split RC et l'indépendance à $\\alpha_u$ ne sont pas démontrées.
4. Corps force (gravité/Boussinesq) : traiter les termes de force exactement comme la pression dans la reconstruction, sinon artefacts de colonne hydrostatique.

## 4. Traitement des frontières
- WALL : $\\dot m_f = 0$ exact (aucune reconstruction).
- Entrée/sortie : flux imposé ou reconstruit à partir de la valeur frontière ; la contribution Rhie–Chow aux faces frontières doit être identique à celle utilisée dans l'équation de continuité.
- EMPTY/périodique : reconstruction au niveau de la face virtuelle avec même formule.

## 5. Pièges connus
- **Checkerboarding** si le terme de stabilisation est omis ; incohérence si la diagonale relaxée/non relaxée est mélangée sans dérivation.
- Correction de Majumdar ajoutée sans preuve : elle peut déplacer le point fixe si le split RC inclut déjà un résidu explicite. La dérivation et le balayage `alpha_u` sont bloquants avant activation en production.
- Flux de continuité ≠ flux de conservation : N9.0 exige que le même flux serve aux deux, ou que l'écart soit documenté (PR #851, test d'antisymétrie des flux Poiseuille).

## 6. Ordre de précision attendu et validation exigée
- $O(h^2)$ sur vitesse si la reconstruction est consistante ; le terme de stabilisation est $O(h^2)$ par construction (différence de gradients).
- Validation exigée : Poiseuille (solution exacte + antisymétrie des flux) ; indépendance à α ; indépendance au pas de temps (transitoire) ; N9.7 négatifs : rejet d'une épingle de pression, détection de null-space faux, etc. (PR #851) ; bilan de masse $< 10^{-10}$.

## 7. Statut d'implémentation CFDX

### Correctifs proposés dans la PR d'audit Rhie–Chow
- L'opérateur interne utilise le même coefficient tensoriel `(dAU e)·Sf / d` dans le flux, la matrice de pression ségrégée et l'approximation de Schur couplée. La réponse différée est `(D ∇p)_f·Sf - ((D e)·Sf)(p_N-p_P)/d`; elle est appariée au même saut de pression que le coefficient implicite.
- Le flux conservatif reconstruit depuis la vitesse cellulaire utilise l'interpolation produit `(D∇p)_f`, puis retranche la réponse implicite `((D e)·Sf)(p_N-p_P)/d`. Cela évite le remplacement incohérent `interp(D∇p) -> interp(D) interp(∇p)` lorsque la réponse varie entre cellules.
- Les coefficients aux frontières à pression fixée sont alignés sur la réponse centre-face et incluent `dAU = V/A_P`.
- La correction de Majumdar reste disponible dans le helper testé, mais n'est pas appliquée par la boucle de production : son interaction avec le split tensoriel doit d'abord être qualifiée par un balayage `alpha_u`.
- Le helper PISO rejette désormais un dénominateur effectif non positif au lieu de le remplacer silencieusement par une petite valeur.
- Le helper historique `assemble_pressure_correction` utilise désormais la convention `dAU = V/A_P` et la projection centre-centre, au lieu d'une formule supposant implicitement volume unitaire et face orthogonale.

### Statut de validation — non clos
- Tests ajoutés : découpage tensoriel sur face oblique/réponse anisotrope, correction de Majumdar, dénominateur PISO invalide, coefficient du helper historique.
- À exécuter sur le SHA final : build C++/tests unitaires, test N9 S9, Poiseuille piloté en pression, Couette, Ghia Re=100, cavité, cas skewed, ainsi que `validation-total`.
- À mesurer impérativement : balayage `alpha_u = 0.3, 0.5, 0.7, 1.0` avec comparaison des champs pression/vitesse (pression alignée à la jauge), flux massiques de face et résidus de quantité de mouvement/continuité recalculés indépendamment.
- Le merge reste bloqué tant que les contrôles ci-dessus ne sont pas verts sur le HEAD exact et que l'indépendance à `alpha_u` n'est pas démontrée numériquement.

## 8. Journal
- 2026-10-12 — création ; intègre les leçons #419 (A-B2/A-B4), #821, #851.
