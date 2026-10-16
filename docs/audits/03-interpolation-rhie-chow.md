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
1. La moyenne $\\overline{(V/a_P)}_f$ est pondérée par les coefficients **au moment de l'assemblage**. Avec la relaxation implicite de la quantité de mouvement, le flux de face reçoit la correction de Majumdar $\\phi_f \\mathrel{+}= (1-\\alpha_u)(\\phi_f^{old}-\\rho\\,\\overline{\\mathbf{u}}_f^{old}\\cdot\\mathbf{S}_f)$, calculée à partir du flux conservatif et de la vitesse précédents. Elle est nulle pour $\\alpha_u=1$. L'indépendance de la solution convergée vis-à-vis de $\\alpha_u$ doit être démontrée sur les cas physiques, pas déduite de la formule seule.
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

### Correctifs proposés dans la PR d'audit Rhie–Chow
- L'opérateur interne utilise le même coefficient tensoriel `(dAU e)·Sf / d` dans le flux, la matrice de pression ségrégée et l'approximation de Schur couplée. La réponse différée est `(D ∇p)·Sf - ((D e)·Sf)(∇p·e)`, et non une projection scalaire qui n'est équivalente que pour une réponse isotrope.
- Le flux conservatif reconstruit depuis la vitesse cellulaire utilise le résidu Rhie–Chow `((p_N-p_P)/d - ∇p_f·e)`; ce terme s'annule pour un champ de pression linéaire exact.
- Les coefficients aux frontières à pression fixée sont alignés sur la réponse centre-face et incluent `dAU = V/A_P`.
- La correction de Majumdar est réintroduite via le flux conservatif précédent, l'interpolation de la vitesse précédente et `alpha_u`; elle est nulle lorsque `alpha_u = 1`.
- Le helper PISO rejette désormais un dénominateur effectif non positif au lieu de le remplacer silencieusement par une petite valeur.
- Le helper historique `assemble_pressure_correction` utilise désormais la convention `dAU = V/A_P` et la projection centre-centre, au lieu d'une formule supposant implicitement volume unitaire et face orthogonale.

### Statut de validation — non clos
- Tests ajoutés : découpage tensoriel sur face oblique/réponse anisotrope, correction de Majumdar, dénominateur PISO invalide, coefficient du helper historique.
- À exécuter sur le SHA final : build C++/tests unitaires, test N9 S9, Poiseuille piloté en pression, Couette, Ghia Re=100, cavité, cas skewed, ainsi que `validation-total`.
- À mesurer impérativement : balayage `alpha_u = 0.3, 0.5, 0.7, 1.0` avec comparaison des champs pression/vitesse (pression alignée à la jauge), flux massiques de face et résidus de quantité de mouvement/continuité recalculés indépendamment.
- Le merge reste bloqué tant que les contrôles ci-dessus ne sont pas verts sur le HEAD exact et que l'indépendance à `alpha_u` n'est pas démontrée numériquement.

## 8. Journal
- 2026-10-12 — création ; intègre les leçons #419 (A-B2/A-B4), #821, #851.
