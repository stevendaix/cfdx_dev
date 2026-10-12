# FICHE MÉTHODE 07 — Oracles de validation : analytique, Ghia, MMS, GCI

**En-tête d'évidence**
- Scope : oracles quantitatifs utilisés par toutes les autres fiches : solutions analytiques (Couette, Poiseuille), benchmark publié (Ghia Re=100), MMS, et estimation d'incertitude de discrétisation (GCI).
- Statut : APPROVED (2026-10-12)
- Issues liées : #34, #461, #601, #639 ; PRs #851 (négatifs), #852 (Poisson polyédrique)
- Dernier commit vérifiéà: 00f5ff8aedeafa9ab802ee302496f3eae217937b (master, 2026-10-12)
- Limitations : les données Ghia tabulées doivent être gélées dans `docs/validation/reference_codes/’ avec leur provenance.

## 1. Références sources
- U. Ghia, K. N. Ghia, C. T. Shin, « High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method », *J. Comput. Phys.* 48(3):387–411, 1982 — DOI : 10.1016/0021-9991(82)90058-4. (cavité entraînée, tables de u/v sur les lignes médianes)
- P. J. Roache, « Perspective: A method for uniform reporting of grid refinement studies », *ASME J. Fluids Eng.* 116(3):405–413, 1994 — DOI : 10.1115/1.2910291. (GCI)
- I. B. Celik, U. Ghia, P. J. Roache, C. J. Freitas, H. Coleman, « Procedure for estimation and reporting of uncertainty due to discretization in CFD applications », *ASME J. Fluids Eng.* 130(7):078001, 2008 — DOI : 10.1115/1.2960953.
- P. J. Roache, *Verification and Validation in Computational Science and Engineering*, Hermosa, 1998 — ISBN 978-0-913478-08-3. (MMS)
- ASME V&V 20, *Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer*, 2009.

## 2. Oracles analytiques (contrats exacts)
- **Couette plan** : $u(y) = U_w\\,y/H$, profil exact linéaire ; gradient de pression nul → tout bug de couplage P-V s'y masque (« Couette ne voit pas le gradient de pression ») → Couette **jamais suffisant seul** pour le couplage P-V.
- **Poiseuille plan** : $u(y) = -(dp/dx)\\,(y(H-y))/(2\\mu)$, débit exact ; cas de référence **minimal obligatoire** de toute modification du couplage (le bug A-B1 #419 a été détecté par lui) ; test d'antisymétrie des flux (N9.7/#851).
- **Canal/entrée développée** : même contrat avec sortie de pression.

## 3. Benchmark publié — Ghia Re=100
Grille de comparaison gélée : profils $u$ sur la ligne médiane verticale et $v$ sur l'horizontale, aux abscisses tabulées de la publication, maillage de comparaison spécifié ; métriques L2/L∞ **pré-déclarées** ; graphique de superposition joint à la PR (script PyVista/Matplotlib, cf. Partie 5 du Plan Directeur).

## 4. MMS (Method of Manufactured Solutions)
Choisir des champs $(\\mathbf{u}, p)$ analytiques lisses, injecter le terme source exact, mesurer les erreurs L2 et L∞ sur ≥ 3 niveaux de maillage ; pente log–log = ordre observé, comparé à l'ordre théorique de la fiche correspondante. Le champ MMS, la famille de maillage, le rapport de raffinement, la norme et la procédure de fit sont **gelés** avant l'exécution (gouvernance §5).

## 5. GCI (Grid Convergence Index)
$$ GCI_\\text{fine} = F_s\\,\\frac{|\\epsilon|}{r^p-1}, \\qquad \\epsilon = \\frac{f_2-f_1}{f_1} $$
avec $r$ le rapport de raffinement, $p$ l'ordre observé, $F_s = 1.25$ (3 grilles, Celik 2008). Toute affirmation « convergence vérifiée » doit citer : $r$, $p$ observé, $F_s$, erreur fine, valeur extrapolée — et vérifier que les hypothèses (range asymptotique, $p$ cohérent) tiennent. GCI = vérification de solution, **pas** une métrique de validation (gouvernance §5).

## 6. Tests négatifs obligatoires
Toute la gamme N9.7 (PR #851) : singularité tout-Neumann, faux null-space, RHS incompatible, coefficients non positifs, antisymétrie de flux, épingle de pression — la machinerie doit **rejeter**, jamais rattraper. Un oracle n'est qualifié que si ses tests négatifs sont verts.

## 7. Statut d'implémentation CFDX
- Data : `docs/validation/reference_codes/’ (provenance Ghia/OpenFOAM/SU2), `EXTENDED_BENCHMARK_MATRIX.md’, `M1_M4_ANALYTICAL_BENCHMARKS.md’.
- Scripts : `scripts/run_validation.py’, `analyze_benchmark.py’, `compare_validation_results.py’, `validation_report.py’.
- Gaps : superposition graphique automatique dans le rapport CI à généraliser ; comparaison OpenFOAM/SU2 sur même maillage exact à déployer (fiche 08).

## 8. Journal
- 2026-10-12 — création.
