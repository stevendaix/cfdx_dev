# FICHE MÉTHODE 01 — Diffusion / Laplacien (Volume Fini)

**En-tête d'évidence**
- Scope : terme de diffusion dans les équations de transport (quantité de mouvement, énergie, scalaire, variables turbulentes). Couvre l'opérateur, l'orthogonalité du maillage, les corrections de non-orthogonalité et de skewness. Non couvert : diffusion turbulente (traitée comme $\\mu + \\mu_t$ dans le même opérateur — voir la fiche turbulence à venir).
- Statut : APPROVED (2026-10-12)
- Issues liées : #461 (N3), #34 ; PRs #845, #852 (campagnes polyédriques exact-HEAD)
- Dernier commit vérifié : 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : fiche rédigée a posteriori pour le code existant ; l'audit code ↔ fiche (section 7) reste à exécuter en Phase 1.

## 1. Références sources
- J. H. Ferziger, M. Perić, *Computational Methods for Fluid Dynamics*, 3e éd., Springer, 2002 — DOI de l'ouvrage : 10.1007/978-3-642-56026-2 (ISBN 978-3-540-42607-3 / ebook 978-3-642-56026-2).
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN 978-3-319-16873-9 (ch. diffusion et faces non orthogonales).
- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI (rééd.) : 10.1201/9781482234213.

## 2. Formulation continue
Transport d'un scalaire $\\phi$ avec diffusivité $\\Gamma$ :
$$ \\frac{\\partial(\\rho\\phi)}{\\partial t} + \\nabla\\cdot(\\rho\\,\\mathbf{u}\\,\\phi) = \\nabla\\cdot(\\Gamma \\nabla\\phi) + S_\\phi $$
Hypothèses : $\\Gamma$ peut être variable (par ex. $\\mu_\\text{eff}$) ; le terme de diffusion est intégré sur les faces et n'exige pas la continuité de $\\nabla\\phi$.

## 3. Discrétisation Volume Fini
Intégration sur la cellule $P$, somme sur les faces $f$ :
$$ \\sum_f \\Gamma_f\\, |S_f|\\,(\\nabla\\phi)_f\\cdot\\mathbf{n}_f \\;\\approx\\; \\sum_f \\Gamma_f\\,D_f\\,(\\phi_N-\\phi_P) + \\sum_f \\Gamma_f\\,[\\text{corrections}] $$
avec $D_f = |S_f \\cdot \\mathbf{e}_{PN}| / |\\mathbf{d}_{PN}|$ la part orthogonale. Les corrections de non-orthogonalité (décomposition $\\mathbf{S}_f = (\\mathbf{S}_f\\cdot\\mathbf{e}_{PN})\\mathbf{e}_{PN} + \\mathbf{T}_f$) et de skewness sont **explicites**, évaluées par itération différée, jamais implicites : la matrice reste un M-matrice de type Laplacien (diagonale ≥ somme des extra-diagonales en valeur absolue, décroissance pour $\\phi=\\text{cste}$).
Terme source : linéarisation $S = S_c - S_p\\,\\phi_P$ pour préserver la bornitude/discretisation.

## 4. Traitement des frontières
- **Dirichlet** ($\\phi=\\phi_b$) : flux exact $\\Gamma_f |S_f| (\\phi_b-\\phi_P)/d_{Pf}$ — coefficient diagonal renforcé, pas de ligne de cellule fantôme sans coefficient.
- **Neumann** (flux imposé $q$) : contribution au second membre $b_P\\, +=\\, q\\,|S_f|$ seulement.
- **Robin** ($a\\phi + b\\,\\partial\\phi/\\partial n = c$) : élimination du flux → coefficient diagonal et second membre explicites.
- Symétrie : flux nul. EMPTY (2D/axisym.) : flux nul, aucune contribution matricielle.

## 5. Pièges connus
- Non-orthogonalité ≥ ~70–75° sans correction → perte d'ordre ; les corrections doivent être **comptées dans le budget de convergence** (itérations explicites).
- Skewness : l'intersection $\\mathbf{e}_{PN}$–face doit être corrigée séparément de la non-orthogonalité ; confondre les deux détruit l'ordre 2.
- Terme source sous forme purement explicite sans linéarisation → découplage/oscillations (leçon de l'audit « B1 » : la diffusion turbulente de k/ε doit transiter par la matrice, pas par le second membre).
- Maillage à fort rapport d'aspect : perte de diagonalité → préconditionnement à adapter (voir fiche 05).

## 6. Ordre de précision attendu et validation exigée
- $O(h^2)$ sur maillage orthogonal ; $O(h)$ garanti sur maillages généraux tant que les corrections ne sont pas convergées (le contrat de #845/#852 impose la **linéarité exacte à chaque niveau** + décroissance stricte de L2 pondéré par les volumes).
- Validation exigée : MMS diffusion ≥ 3 niveaux ; Poisson polyédrique avec vrai résidu recalculé ≤ 1e-8 ; bilan de flux global indépendant à la tolérance machine ; parité CPU/MPI (même maillage, 1e-12).
- Familles de maillage obligatoires : orthogonal, non orthogonal, skewed, fort rapport d'aspect (N10/#639).

## 7. Statut d'implémentation CFDX
- Sources : `src/cfdx/physics/finite_volume_transport.h’ (assemblage `assemble_scalar_equation’), `src/cfdx/core/…’ opérateurs.
- Tests : campagne Poisson polyédrique (PR #852), N3 (`docs/validation/N3_DIFFUSION_QUALIFICATION.md’).
- Écarts code ↔ fiche : **à auditer en Phase 1** (vérifier que la part non orthogonale et la part skew sont bien séparées et explicites).
- Gaps : familles skew/high-aspect-ratio incomplètes (cf. #34 N10) ; qualification MPI limitée.

## 8. Journal
- 2026-10-12 — création (PR gouvernance, Plan Directeur 2026-10-12).
